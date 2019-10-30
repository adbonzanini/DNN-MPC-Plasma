
%
%   This script uses results from main_dnn.m to run a controller which
%   substitutes the OCP with a DNN. Relevant variables from main_dnn.m are
%   saved in the DNN_training.mat file.
%


% Clear workspace
clear all

%% Load relevant inputs for DNN training
load('Supporting-Data-Files/DNN_training.mat');
model_ID=load('Supporting-Data-Files/MIMOmodelGlass.mat');
Tss = model_ID.steadyStates(1);
Iss = round(model_ID.steadyStates(2),1);
Pss = round(model_ID.steadyStates(3),1);
Qss = round(model_ID.steadyStates(4),1);

%% User-defined inputs

for scenario = [1,2,3]
rng(300)

if scenario==1
    useProj=0;
    approximate=0;
elseif scenario==2
    useProj=0;
    approximate=1;
elseif scenario==3
    useProj=1;
    approximate=1;
end

% Number of simulations/time-steps
Nsim = 50;

% Sampling time
Tsampling=1.3;

% Define reference
Sdes = 1.5*ones(1, Nsim+1);
KcemThreshold = 35;
Kcem = 0.5;

% Stylistic choices depending on case
if useProj==0 && approximate==0
    color='r';
elseif useProj==1
    color = 'g';
elseif useProj==0 && approximate==1
    color = 'b';
end

lambdaf=0.95;

%% Project into maximal robust control invariant set
if useProj==1
    % Bounds on w
    w_upper = [1.75; 0]'; %2.5 (0.8sim) try robustifying one at a time if too conservative and you know where you are going to operate
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


if approximate==0
    [solver, args, f] = getNMPCSolver(Tsampling, N, x_init, x_min, x_max, u_init, u_min, u_max, CEM_init, CEM_min, CEM_max, CEMsp);
end


%% Perform simulations to check results

% initialize
Xsim = zeros(nx,Nsim+1);
Ysim = zeros(nx,Nsim+1);
Usim = zeros(nu,Nsim);
% Wsim = zeros(nx,Nsim);
Wsim = normrnd(0, 0.75, nx, Nsim);
What = zeros(nx,Nsim);
CEMcurr = zeros(ny, Nsim+1);
Trun = zeros(Nsim,1);

% initial states
Xsim(:,1) = [0;0];
Ysim(:,1) = C*Xsim(:,1);

% Initialize vectors for plotting
Yplot = Ysim;
Uplot = Usim;

% reset random seed
rng(200, 'twister')

% run loop over time
for k = 1:Nsim
    tStart = tic;
    if approximate == 1
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
            try
                Usim(:,k) = Usim(:,k-1);
                What(:,k) = Wsim(:,k-1);
                Ysim(:,k) = Ysim(:,k-1);
            catch
                Usim = [0;0];
                What(:,k) = [0;0];
                Ysim(:,k) = [0;0];
            end
            Xsim(:,k) = Ysim(:,k);
        end
    else
        data_rand = [Xsim(:,k)', CEMcurr(k)];
        [U_mpc, Feas, V_opt] = solveSampleNMPC(solver, args, data_rand);
        Usim(:,k) = U_mpc;
    end
    % this calls the original offset-free mpc
%     [sol,errorcode] = controller{[Xsim(:,k);Sdes(:,k)-Hd*What(:,k)]};
%     Usim(:,k) = double(sol(1:nu));
        
    % get most recent disturbance
%     Wsim(:,k) = [0;0];

    % provide values to plant
    Areal =1.1*A;
    Breal = 0.9*B;
    Xsim(:,k+1) = Areal*Xsim(:,k) + Breal*Usim(:,k) + Wsim(:,k);
    Ysim(:,k+1) = C*Xsim(:,k+1);

    if Ysim(1,k)+Tss<= KcemThreshold
        CEMcurr(k+1) = CEMcurr(k);
    else
        CEMcurr(k+1) = CEMcurr(k)+Kcem.^(6-Ysim(1,k+1));
    end
        


    % estimate disturbance
    What(:,k+1) = lambdaf*What(:,k) + (1-lambdaf)*Wsim(:,k);
    Trun(k) = toc(tStart);
end

%% Plot results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ALL TRAJECTORIES ON SEPARATE GRAPHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Plot phase plot
figure; hold on
X.plot('wire',1,'edgecolor','r','linestyle','--','linewidth',2)
plot(Xsim(1,:),Xsim(2,:),'-ok','MarkerFaceColor','k','linewidth',1.5,'MarkerSize',6)
ylim([-5, 10])
set(gcf,'color','w');
set(gca,'FontSize',16)
axis([-10, 10, -20, 25])
%}
time = 0:Tsampling:Nsim*Tsampling;

% Plot CEM
figure(3)
hold on
stairs(time,Sdes(1,:),'k')
hCEM{scenario} = plot(time,CEMcurr,color, 'Linewidth', 2);
ylim([0, 2*max(Sdes(1,:))])
xlabel('Time Step')
ylabel('CEM/min')
set(gcf,'color','w');
set(gca,'FontSize',12)

% Plot states
figure(4);
subplot(2,1,1)
hold on
hT{scenario} = plot(time, Ysim(1,:)+Tss, color, 'Linewidth', 2);
plot([time(1),time(end)],[x_max(1),x_max(1)]+Tss,'--k')
plot([time(1),time(end)],[x_min(1),x_min(1)]+Tss,'--k')
set(gcf,'color','w')
set(gca,'FontSize',12)
xlabel('Time/ s')
ylabel('Temperature/ ^{\circ}C')
set(gca,'FontSize',12)

subplot(2,1,2)
hold on
hI{scenario} = plot(time, 10*(Ysim(2,:)+Iss),color, 'Linewidth', 2);
plot([time(1),time(end)],10*([x_max(2),x_max(2)]+Iss),'--k')
plot([time(1),time(end)],10*([x_min(2),x_min(2)]+Iss),'--k')
set(gcf,'color','w')
xlabel('Time/ s')
ylabel('Intensity/ a.u')
set(gca,'FontSize',12)

% Plot inputs
figure(5);
subplot(2,1,1)
hold on
hQ{scenario} = plot(time(1:end-1), Usim(1,:)+Qss, color, 'Linewidth', 2);
plot([time(1),time(end)],[10, 10],'--k')
plot([time(1),time(end)], [0.8, 0.8],'--k')
set(gcf,'color','w')
xlabel('Time/ s')
ylabel('He Flowrate/ slm')
set(gca,'FontSize',12)

subplot(2,1,2)
hold on
hP{scenario} = plot(time(1:end-1), Usim(2,:)+Pss,color, 'Linewidth', 2);
plot([time(1),time(end)], [5, 5],' --k')
plot([time(1),time(end)],[0.5, 0.5], '--k')
set(gcf,'color','w')
xlabel('Time/ s')
ylabel('Applied Power/ W')
set(gca,'FontSize',12)

if approximate==1
    disp('Approximate NMPC ')
else
    disp('Full NMPC')
end
disp(['Average Run Time = ', num2str(mean(Trun)), ' seconds'])

end



figure(5)
subplot(2,1,1)
legend([hQ{1}, hQ{2}, hQ{3}], 'Full NMPC', 'Approximate NMPC', 'Approximate NMPC with Projection')
set(gca,'FontSize',12)
figure(4)
legend([hT{1}, hT{2}, hT{3}], 'Full NMPC', 'Approximate NMPC', 'Approximate NMPC with Projection')
set(gca,'FontSize',12)
figure(3)
legend([hCEM{1}, hCEM{2}, hCEM{3}], 'Full NMPC', 'Approximate NMPC', 'Approximate NMPC with Projection')
set(gca,'FontSize',12)


