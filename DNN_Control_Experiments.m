%
%   This script uses results from main_dnn.m to run a controller which
%   substitutes the OCP with a DNN. Relevant variables from main_dnn.m are
%   saved in the DNN_training.mat file.
%

% Two files before 14:48
% Two files before 15:06
% 5-->    <16:10; 16:11, 
% 5.5 --> <16:17; 16:19; || 16:23; 16:27 || 16:32; 16:36 ||

% <13:34; 13:36; 13:43 [2;-2]
% <14:10;<14:25; <14:30; 14:34 [2;3]

% Re-visit model --> collect data
% Reduce sampling time? --> is this going to affect data collection?

% Clear workspace
clear all

%% Load relevant inputs for DNN training
load('Supporting-Data-Files/DNN_training.mat');
model_ID=load('Supporting-Data-Files/MIMOmodelGlass.mat');
Tss = model_ID.steadyStates(1);
Iss = model_ID.steadyStates(2);
Pss = model_ID.steadyStates(3);
Qss = model_ID.steadyStates(4);

%% User-defined inputs

% Experiments or simulations
runExperiments = 0;

% Switch for projection to a safe set
useProj = 0;

% Number of simulations/time-steps
Nsim = 160;

% Define reference
tChange = 80;
% Sdes = [[0;0].*ones(nx, tChange), [5.5;-2].*ones(nx, Nsim+1-tChange)];
% Sdes = [[-1;0].*ones(nx, tChange), [2;2.8].*ones(nx, Nsim+1-tChange)];
Sdes = [1;0].*ones(nx, Nsim+1);



%% Project into maximal robust control invariant set

% Bounds on w (not needed here)
w_upper = [1.1; 0]'; % try robustifying one at a time if too conservative and you know where you are going to operate
w_lower = [0; 0]';
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
    fopen(t);
    disp('Connection established...')
    pause(6);
    
    % Send initial point as a comma-delimited string
    u_opt=[2.5;1.5];
    Usend = sprintf('%6.1f, %6.1f ', [u_opt(1), u_opt(2)]);
    disp('Sending initial point...') 
    pause(0.2)
    fwrite(t, Usend)
    disp(['Sent inputs (Q, P) = ', '[', Usend, ']'])

end

% calculate offset gain
Hd = C*inv(eye(nx)-(A+B*K));
lambdaf = 0.95;

% initialize
Xsim = zeros(nx,Nsim+1);
Ysim = zeros(ny,Nsim+1);
Usim = zeros(nu,Nsim);
Wsim = zeros(nx,Nsim);
What = zeros(nx,Nsim);

% initial states
Xsim(:,1) = [0;0];
Ysim(:,1) = C*Xsim(:,1);

% Initialize vectors for plotting
Yplot = Ysim;
Uplot = Usim;

% reset random seed
rng(200, 'twister')

figure(1)
% run loop over time
for k = 1:Nsim
    % evaluate the explicit controller
    xscaled = ([Xsim(:,k);Sdes(:,k)-Hd*What(:,k)] - xscale_min')./(xscale_max-xscale_min)';
    tscaled = net(xscaled)';
    uexp = (tscale_min+(tscale_max-tscale_min).*tscaled)';
    


    % specify to use the projection or just the DNN
    if useProj == 1
        assign(xcurrent, Xsim(:,k));
        assign(uexplicit, uexp);
        value(Optimizer)
        Usim(:,k) = value(Optimizer);        
    else
        Usim(:,k) = uexp;
    end

    if isnan(Usim(:,k))==[1;1]
        warning('NaN in Usim. Assigning previous value...');

        Usim(:,k) = [0;-1];
        What(:,k) = Wsim(:,k-1);
        Ysim(:,k) = Ysim(:,k-1);
        Xsim(:,k) = Ysim(:,k);
        
    end
    % this calls the original offset-free mpc
%     [sol,errorcode] = controller{[Xsim(:,k);Sdes(:,k)-Hd*What(:,k)]};
%     Usim(:,k) = double(sol(1:nu));
    
    if runExperiments == 0
        
        % get most recent disturbance
        Wsim(:,k) = [0;0];
        
        % provide values to plant
        Xsim(:,k+1) = A*Xsim(:,k) + B*Usim(:,k) + Wsim(:,k);
        Ysim(:,k+1) = C*Xsim(:,k+1);
        
    else
        % Send optimal input to the set-up
        Usend = [Usim(1,k)+Qss;Usim(2,k)+Pss];
        Usend=[Usend(:)', Sdes(1,k)+Tss,(Sdes(2,k)+Iss)*10.0];
        Usend = sprintf('%.1f,', Usend(:));
        Usend = Usend(1:end-1);
        disp('Sending inputs...') 
        pause(0.3)
        fwrite(t, Usend)
        disp(['Sent inputs (q,P,Tss, Iss) = ','[', Usend, ']'])
        
        pause(0.9)
        
        % Receive measurement
        disp('Receive Measurement...')
        measurements = fread(t, [1, t.BytesAvailable]);
        try
            y_meas = char(measurements(end-12:end));
            % Convert string to array
            y_m = str2num(y_meas);
            
        catch
            % Assign previous measurements if no measurement is obtained
            disp('Assigning previous measurement')
            
            if Yplot(1,k)==0
                y_m = [Tss, Iss*10];
            else
                y_m = Yplot(:,k)'; 
            end
        end
        
        disp(y_m)
        
        % Temperature
        Yplot(1,k+1) = y_m(1);
        Ysim(1,k+1) = y_m(1)-Tss;

        % Intensity
        Yplot(2,k+1) = y_m(2);
        Ysim(2,k+1) = y_m(2)*0.1-Iss;
        
        % State feedback
        Xsim(:,k+1) = Ysim(:,k+1);

        % Optimal inputs
        Uplot(:,k) = Usim(:,k)+[Qss;Pss];
        
        % Measure plant-model mismatch (state feedback case)
        Wsim(:,k) = Xsim(:,k+1)-A*Xsim(:,k)-B*Usim(:,k);
        if isnan(Wsim(:,k))==[1;1]
            Wsim(:,k) = Wsim(:,k-1);
        end
        
        % Live plotting
        subplot(2,1,1)
        hold on
        plot(Yplot(1,1:k), 'r')
        plot(Sdes(1,1:k)+Tss, 'k-')
        plot([0,k],[x_max(1)+Tss,x_max(1)+Tss], 'k--')
        ylim([25, 50])
        xlabel('Time Step')
        ylabel('Temperature/oC')
        
        subplot(2,1,2)
        hold on
        plot(Yplot(2,1:k), 'r')
        plot((Sdes(2,1:k)+Iss)*10, 'k-')
        plot([0,k],([x_max(2),x_max(2)]+Iss)*10, 'k--')
        ylim([0, 255])
    end


    % estimate disturbance
    What(:,k+1) = lambdaf*What(:,k) + (1-lambdaf)*Wsim(:,k);
end
%%

% Close pcp ip connection
if runExperiments==1
    restorePlasma
%     fclose(t)
end

%%
% Plot results

% Plot phase plot
figure; hold on
X.plot('wire',1,'edgecolor','r','linestyle','--','linewidth',2)
plot(Xsim(1,:),Xsim(2,:),'-ok','MarkerFaceColor','k','linewidth',1.5,'MarkerSize',6)
ylim([-5, 10])
set(gcf,'color','w');
set(gca,'FontSize',16)
axis([-10, 10, -20, 25])

% Plot output
figure; hold on;
time = 0:Nsim;
subplot(2,1,1)
hold on
stairs(time,Sdes(1,:),'k')
plot(time,Ysim(1,:),'r')
plot([0,Nsim],[x_max(1),x_max(1)],'--r')
ylim([-8, 8])
set(gcf,'color','w');
set(gca,'FontSize',12)

subplot(2,1,2)
hold on
stairs(time,Sdes(2,:),'k')
plot(time,Ysim(2,:),'r')
plot([0,Nsim],[x_max(2),x_max(2)],'--r')
set(gcf,'color','w');
set(gca,'FontSize',12)

