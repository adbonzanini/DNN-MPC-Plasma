%
%   This script uses results from main_dnn.m to run a controller which
%   substitutes the OCP with a DNN. Relevant variables from main_dnn.m are
%   saved in the DNN_training.mat file.
%

% Clear workspace
clear all

%% User-defined inputs

% Experiments or simulations
runExperiments = 0;

% Number of simulations/time-steps
Nsim = 40;

% Define reference
Sdes = [zeros(nx, 20), [4.9;2].*ones(nx, Nsim+1-20)];



%% Load relevant inputs for DNN training
load('Supporting-Data-Files/DNN_training.mat')



%% Project into maximal robust control invariant set

% Bounds on w (not needed here)
w_upper = [0.8; 2.0]';
w_lower = [-0.8; -2.0]';
W = Polyhedron('lb',w_lower','ub',w_upper');

% Calculate robust control invariant set
sys = ULTISystem('A',A,'B',B,'E',eye(nx));
sys.x.min = x_min';
sys.x.max = x_max';
sys.u.min = u_min';
sys.u.max = u_max';
sys.d.min = w_lower';
sys.d.max = w_upper';
Cinf = sys.invariantSet('maxIterations',50);

% Define problem to project into Cinf
Cinf_next = Cinf - W;
Cinf_next.computeVRep();
xcurrent = sdpvar(nx,1);
uexplicit = sdpvar(nu,1);
uproject = sdpvar(nu,1);
constraints = [];
objective = 0;
xnext = A*xcurrent + B*uproject;
constraints = [constraints, Cinf_next.A*xnext <= Cinf_next.b];
constraints = [constraints, U.A*uproject <= U.b];
objective = objective + (uexplicit - uproject)'*(uexplicit - uproject);

% Add constraints on the explicit variable to bound the size of the mp map
constraints = [constraints, -2*(u_max-u_min)' + u_min' <= uexplicit <= u_max' + 2*(u_max-u_min)'];
constraints = [constraints, Cinf.A*xcurrent <= Cinf.b];

% Create optimizer object
ops = sdpsettings('verbose',0);
explicit_controller = optimizer(constraints,objective,ops,[xcurrent;uexplicit],[uproject]);

% Calculate the explicit solution using yalmip
[mptsol,diagn,Z,Valuefcn,Optimizer] = solvemp(constraints,objective ,ops,[xcurrent;uexplicit],[uproject]);



%% Perform simulations/experiments to check results

% If the option to run experiments is turned on, send the initial point to
% the plasma jet
if runExperiments==1
    % Set up tcp/ip socket client
    host = 'pi3.dyn.berkeley.edu';
    port = 2223;
    t = tcpip(host, port, 'NetworkRole', 'client');
    fopen(t)
    print('Connection established...')
    
    % Send initial point as a comma-delimited string
    u_opt=[1.5;1.5];
    Usend = sprintf('%.1f,', u_opt(:));
    Usend = Usend(1:end-1);
    disp('Sending initial point...') 
    fwrite(t, Usend)
    disp(['[', Usend, ']'])
end

% calculate offset gain
Hd = C*inv(eye(nx)-(A+B*K));
lambdaf = 0.5;

% initialize
Xsim = zeros(nx,Nsim+1);
Ysim = zeros(ny,Nsim+1);
Usim = zeros(nu,Nsim);
Wsim = zeros(nx,Nsim);
What = zeros(nx,Nsim);

% initial states
Xsim(:,1) = [1;0];
Ysim(:,1) = C*Xsim(:,1);

% Initialize vectors for plotting
Yplot = Ysim;
Uplot = Usim;

% reset random seed
rng(200, 'twister')

% run loop over time
for k = 1:Nsim
    % evaluate the explicit controller
    xscaled = ([Xsim(:,k);Sdes(:,k)-Hd*What(:,k)] - xscale_min')./(xscale_max-xscale_min)';
    tscaled = net(xscaled)';
    uexp = (tscale_min+(tscale_max-tscale_min).*tscaled)';

    % specify to use the projection or just the DNN
    useProj = 1;
    if useProj == 1
        assign(xcurrent, Xsim(:,k));
        assign(uexplicit, uexp);
        value(Optimizer)
        Usim(:,k) = value(Optimizer);        
    else
        Usim(:,k) = uexp;
    end
    % this calls the original offset-free mpc
%     [sol,errorcode] = controller{[Xsim(:,k);Sdes(:,k)-Hd*What(:,k)]};
%     Usim(:,k) = double(sol(1:nu));
    
    if runExperiments == 0
        
        % get most recent disturbance
        Wsim(:,k) = [0.5;0.5];
        
        % provide values to plant
        Xsim(:,k+1) = A*Xsim(:,k) + B*Usim(:,k) + Wsim(:,k);
        Ysim(:,k+1) = C*Xsim(:,k+1);
        
    else
        % Send optimal input to the set-up
        Usend = Usim(:,k)+[3.0;3.0];
        Usend=[Usend(:)', Sdes(1,k)+38.0,(Sdes(2,k)+9.7)*10.0];
        Usend = sprintf('%.1f,', Usend(:));
        Usend = Usend(1:end-1);
        disp('Sending inputs...')        
        fwrite(t, Usend)
        disp(['[', Usend, ']'])  
        
        % Receive measurement
        disp ('Receiver started');
        t=tcpip(host, port,'NetworkRole','server');

        % Wait for connection
        disp('Waiting for connection');
        fopen(t);
        disp('Connection OK');
        y_meas=fread(t,t.BytesAvailable);
        try
            % Convert string to array
            y_m = str2num(y_meas);
        catch
            % Assign previous measurements if no measurement is obtained
            warning('FAILURE IN OBTAINING MEASUREMENT')
            disp('Assigning previous measurement')
            y_m = Yplot(:,k-1)';
        end
        
        % Temperature
        Yplot(1,k+1) = y_m(1);
        Ysim(1,k+1) = y_m(1)-38.0;

        % Intensity
        Yplot(2,k+1) = y_m(2);
        Ysim(2,k+1) = y_m(2)*0.1-9.7;
        
        % State feedback
        Xsim(:,k) = Ysim(:,k);

        % Optimal inputs
        Uplot(:,k) = Usim(:,k)+[3.0;3.0];
        
        % Measure plant-model mismatch (state feedback case)
        Wsim = Xsim(:,k+1)-A*Xsim(:,k)-B*Usim(:,k);
    end

    % estimate disturbance
    What(:,k+1) = lambdaf*What(:,k) + (1-lambdaf)*Wsim(:,k);
end



%% Plot results

% Plot phase plot
figure; hold on
X.plot('wire',1,'edgecolor','r','linestyle','--','linewidth',2)
plot(Xsim(1,:),Xsim(2,:),'-ok','MarkerFaceColor','k','linewidth',1.5,'MarkerSize',6)
set(gcf,'color','w');
set(gca,'FontSize',16)
axis([-10, 10, -20, 25])

% Plot output
figure; hold on;
time = 0:Nsim;
stairs(time,Sdes(1,:),'k')
plot(time,Ysim(1,:),'r')
plot([0,Nsim],[x_max(1),x_max(1)],'--r')
set(gcf,'color','w');
set(gca,'FontSize',12)

