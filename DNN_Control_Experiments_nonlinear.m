close all
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

% <15:08; <15:10; <15:33 
% <15:15; <15:19 <15:23; <15:28

% <13:39
% <13:37; 13:41

% 13:56; 14:01; 14:06
% 13:55; 13:57; 14:03

%16:10
%16:09 16:12

% 14:39; <14:41; 14:43 || 14:52; 14:53; 14:54

% Re-visit model --> collect data
% Reduce sampling time? --> is this going to affect data collection?

% Clear workspace
clear all

%% Load relevant inputs for DNN training
load('Supporting-Data-Files/DNN_training.mat');
model_ID=load('Supporting-Data-Files/MIMOmodelGlass.mat');
Tss = model_ID.steadyStates(1);
Iss = round(model_ID.steadyStates(2),1);
Pss = round(model_ID.steadyStates(3),1);
Qss = round(model_ID.steadyStates(4),1);

%User-defined inputs

% Experiments or simulations
runExperiments = 1;

% Switch for projection to a safe set
useProj = 0;

% Number of simulations/time-steps
Nsim = 40;

% Define reference
% tChange = 10;
% Sdes = [[0;0].*ones(nx, tChange), [5.5;-2].*ones(nx, Nsim+1-tChange)];
% Sdes = [[-1;0].*ones(nx, tChange), [2;2.8].*ones(nx, Nsim+1-tChange)];
Sdes = 1.5*ones(1, Nsim+1);
KcemThreshold = 35;
Kcem = 0.5;

%% Project into maximal robust control invariant set
if useProj==1
% Bounds on w
w_upper = [2.5; 0]'; %2.5 (0.8sim) try robustifying one at a time if too conservative and you know where you are going to operate
w_lower = -[0; 0]';
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
constraints = [constraints, -100*(u_max-u_min)' + u_min' <= uexplicit <= u_max' + 100*(u_max-u_min)'];
constraints = [constraints, Cinf.A*xcurrent <= Cinf.b];

% Create optimizer object
ops = sdpsettings('verbose',0);
explicit_controller = optimizer(constraints,objective,ops,[xcurrent;uexplicit],[uproject]);

% Calculate the explicit solution using yalmip
[mptsol,diagn,Z,Valuefcn,Optimizer] = solvemp(constraints,objective ,ops,[xcurrent;uexplicit],[uproject]);
end


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
Ysim = zeros(nx,Nsim+1);
Usim = zeros(nu,Nsim);
Wsim = zeros(nx,Nsim);
What = zeros(nx,Nsim);
CEMcurr = zeros(ny, Nsim+1);

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
%     xscaled = ([Xsim(:,k);Sdes(:,k)-Hd*What(:,k)] - xscale_min')./(xscale_max-xscale_min)';
    xscaled = ([Xsim(:,k);CEMcurr(k)] - xscale_min')./(xscale_max-xscale_min)';
    
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

    if any(isnan(Usim(:,k)))==1
        disp('-----------------------------------------');
        disp('NaN in Usim. Assigning nominal value...');
        disp('-----------------------------------------');

        Usim(:,k) = [4;-1];
        if k==1
            What(:,k) = [0;0];
            Ysim(:,k) = [0;0];
        else
            What(:,k) = Wsim(:,k-1);
            Ysim(:,k) = Ysim(:,k-1);
        end
        Xsim(:,k) = Ysim(:,k);
        
    end
    % this calls the original offset-free mpc
%     [sol,errorcode] = controller{[Xsim(:,k);Sdes(:,k)-Hd*What(:,k)]};
%     Usim(:,k) = double(sol(1:nu));
    
    if runExperiments == 0
        
        % get most recent disturbance
        Wsim(:,k) = [0;0];
        
        % provide values to plant
        Areal =1.2*A;
        Breal = 1.3*B;
        Xsim(:,k+1) = Areal*Xsim(:,k) + Breal*Usim(:,k) + Wsim(:,k);
        Ysim(:,k+1) = C*Xsim(:,k+1);
        
        if Ysim(1,k)+Tss<= KcemThreshold
            CEMcurr(k+1) = CEMcurr(k);
        else
            CEMcurr(k+1) = CEMcurr(k)+Kcem.^(6-Ysim(1,k+1));
        end
        
    else
        % Send optimal input to the set-up
        Usend = [Usim(1,k)+Qss;Usim(2,k)+Pss];
        Usend=[Usend(:)', Sdes(1,k)+Tss,0];
        Usend = sprintf('%.1f,', Usend(:));
        Usend = Usend(1:end-1);
        disp('Sending inputs...')
        
        pause(0.3)
        
        fwrite(t, Usend)
        disp(['Sent inputs (q,P,Tss, Iss) = ','[', Usend, ']'])
        
        pause(1.0)
        
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
        if Yplot(1,k+1)<= KcemThreshold
            CEMcurr(k+1) = CEMcurr(k);
        else
            CEMcurr(k+1) = CEMcurr(k)+Kcem.^(43-Yplot(1,k+1));
        end

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
        plot(CEMcurr(1,1:k), 'r')
        plot(Sdes(1,1:k), 'k-')
        ylim([0, 1.2*max([Sdes(1), max(CEMcurr(1,1:k))])])
        xlabel('Time Step')
        ylabel('CEM/min')
        subplot(2,1,2)
        hold on
        plot(Yplot(1,1:k), 'r')
        plot([0, k], [x_max(1)+Tss, x_max(1)+Tss], 'k--')
        ylim([KcemThreshold-1, x_max(1)+1+Tss])
        
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
hold on
stairs(time,Sdes(1,:),'k')
plot(time,CEMcurr,'r')
ylim([0, 2*max(Sdes(1,:))])
xlabel('Time Step')
ylabel('CEM/min')
set(gcf,'color','w');
set(gca,'FontSize',12)

figure;
subplot(2,1,1)
hold on
plot(Ysim(1,:)+Tss, 'r')
plot([0,Nsim],[x_max(1),x_max(1)]+Tss,'--k')
plot([0,Nsim],[x_min(1),x_min(1)]+Tss,'--k')
set(gcf,'color','w')
set(gca,'FontSize',12)
xlabel('Time Step')
ylabel('Temperature/ ^{\circ}C')

subplot(2,1,2)
hold on
plot(10*Ysim(2,:)+Iss,'r')
plot([0,Nsim],10*([x_max(2),x_max(2)]+Iss),'--k')
set(gcf,'color','w')
set(gca,'FontSize',12)
xlabel('Time Step')
ylabel('Intensity/ a.u')

