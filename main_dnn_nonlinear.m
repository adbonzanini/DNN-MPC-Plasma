
% Clear variables
yalmip('clear')
clear

% Start uqlab
uqlab

% Fix random seed
rng(100, 'twister')



%% Setup tracking mpc problem
Ts=1.3;

% Cost matrices
% Tuning for contstraint violation, but with some oscillations:
Q = diag([10, 0.2]);
R = diag([20, 10]);
% Tuning without oscillations
% Q = diag([5, 0.001]);
% R = diag([40, 20]);
PN = 1*Q;

% Load identified model matrices for 
% x(k+1) = Ax(k) + Bu(k)
% y(k) = Cx(k) + Du(k)
% Reduced variables around xss = [37.4; 9.7] and uss = [3;3]
directory = pwd;
modelGlass = load([directory, '/Supporting-Data-Files/MIMOmodelGlass']);
A = modelGlass.A; B=modelGlass.B; C=modelGlass.C; D=0;
Tss = modelGlass.steadyStates(1);Iss = modelGlass.steadyStates(2);
qss = modelGlass.steadyStates(3); Pss = modelGlass.steadyStates(4);


% Sizes
nx = size(A,2);
nu = size(B,2);
ny = 1;

% Horizon length
N = 5;

% Infinite horizon LQR
[K,Pinf,Eig] = dlqr(A,B,Q,R);
K = -K;                         % Matlab vs. literature conventions

% Observer parameters
Hs = (C+D*K)*inv(eye(nx)+(A+B*K));


% State constraints
Xcon = [eye(2,3);-eye(2,3)];
XlimReal = [41;250;33;0];
Xcon(:,3) = [XlimReal(1)-Tss; XlimReal(2)*0.1-Iss;-(XlimReal(3)-Tss);-(XlimReal(4)*0.1-Iss)];

% Input constraints
Ucon = [1, 0, 0;-1, 0, 0;0, 1, 0;0, -1, 0];
Ucon(:,3) = [7;1.2;2;2.5];

% Terminal constraints
Xf = Xcon;

% Create constraints in MPT 
X = Polyhedron('A',Xcon(:,1:2),'b',Xcon(:,3));
U = Polyhedron('A',Ucon(:,1:2),'b',Ucon(:,3));

% Marginals for states, inputs, and outputs
x_min = -[Xcon(3,3), Xcon(4,3)];
x_max = [Xcon(1,3), Xcon(2,3)];
x_init = [0,0];
u_min = -[Ucon(2,3), Ucon(4,3)];
u_max = [Ucon(1,3), Ucon(3,3)];
u_init = [0,0];
% Current CEM bounds
CEM_min = 0;
CEM_max = 2;
CEM_init = 0;
CEMsp = CEM_max-0.5;

%%
% Setup the mpc problem
[solver, args, f] = getNMPCSolver(Ts, N, x_init, x_min, x_max, u_init, u_min, u_max, CEM_init, CEM_min, CEM_max, CEMsp);

%%

%{
% Setup the mpc problem
if N == 1
    u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
    u = {u};
else
    u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
end
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
xinit = sdpvar(nx,1);
xs = sdpvar(nx,1);
us = sdpvar(nu,1);
ys = sdpvar(ny,1);
constraints = [];
constraints = [constraints, xs == A*xs + B*us];
constraints = [constraints, x{1} == xinit];
constraints = [constraints, X.A*xs <= X.b];
constraints = [constraints, U.A*us <= U.b];
objective = 0;

for k = 1:N
    
    %Define Kcem
    if value(x{k})<=2
        Kcem=0;
    else
        Kcem = 0.25;
    end
    
%    objective = objective + (Kcem^(43)/log(Kcem))*(Kcem^(-x{k}(1)-37)-Kcem^(-x{k+1}(1)-37));
%    objective = objective + Kcem^(43-(x{k}(1)+37));
%     objective = objective + Kcem^6-Kcem^6*log(Kcem)*x{k}(1) + 0.5*Kcem^6*(log(Kcem))^2*x{k}(1).^2;
%     objective = objective + sum((Q*(x{k}-xs)).*(x{k}-xs)) + sum((R*(u{k}-us)).*(u{k}-us));
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}]; % dynamic equality constraints
    constraints = [constraints, U.A*u{k} <= U.b]; % input constraints
    constraints = [constraints, X.A*x{k+1} <= X.b]; % state constraints
end
objective = objective+Kcem^6-Kcem^6*log(Kcem)*ys(1) + 0.5*Kcem^6*(log(Kcem))^2*ys(1).^2;
% objective = objective + sum((PN*(x{N+1}-xs)).*(x{N+1}-xs)); % terminal cost
% objective = objective + 1000*(ys - C*xs)'*(ys - C*xs);
%}

% Calculate feasible region of mpc problem
DoA = X;
for i = 1:N
    pre = inv(A) * (DoA + (-B*U));
    DoA1 = pre & X;
    if DoA1 == DoA
        i
        break    
    else
        DoA = DoA1;
    end
end

%{
% Create optimizer object
ops = sdpsettings('verbose',0);
controller = optimizer(constraints,objective,ops,[xinit;ys],[u{1}]);
%}


%% Create input objects


for i = 1:nx
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [x_min(i), x_max(i)];    
end

% Create state "input" object
myInput_X = uq_createInput(Input);

% Marginals for the reference
for i = 1:ny
    Input.Marginals(nx+i).Type = 'Uniform';
    Input.Marginals(nx+i).Parameters = [CEM_min(i), CEM_max(i)];
end

% Create total "input" object (both states and reference)
myInput_P = uq_createInput(Input);



%% Sample within the domain of attraction

% Specify number of samples
Nsamp = 2500;

% Sample the state/reference space
Psamp = uq_getSample(myInput_P, 10*Nsamp, 'MC');

% Check which points are inside DoA
index = DoA.contains(Psamp(:,1:nx)');
data_rand = Psamp(index,:);
data_rand = data_rand(1:Nsamp,:);

%{
% Solve tube mpc problem over samples
U_mpc = zeros(Nsamp,nu);
Feas = zeros(Nsamp,1);
for i = 1:Nsamp
    xcurr = data_rand(i,1:nx)';
    ytcurr = data_rand(i,nx+1:nx+ny)';
    
    [sol,errorcode] = controller{[xcurr;ytcurr]};
    uopt = double(sol(1:nu));
    
    if errorcode == 1
        fprintf('QP infeasible!\n')
        Feas(i) = 1;
    else
        if mod(i,10)==0
            fprintf('%g\n', i)
        end
    end

    U_mpc(i,:) = uopt';
end
%}

[U_mpc, Feas, V_opt] = solveSampleNMPC(solver, args, data_rand);

target_rand = [U_mpc];

% Combine data
data = [data_rand];
target = [target_rand];

% Scale variables
xscale_min = [x_min, CEM_min];
xscale_max = [x_max, CEM_max];
x = (data - repmat(xscale_min,[size(data,1),1]))./(repmat(xscale_max-xscale_min,[size(data,1),1]));
x = x';
tscale_min = u_min;
tscale_max = u_max;
t = (target - repmat(tscale_min,[size(data,1),1]))./(repmat(tscale_max-tscale_min,[size(data,1),1]));
t = t';

% List of nodes and layers
Nlayers_list = 5;%3:1:10;
Nnodes_list = 5; %5:1:20;
mse_tol = 1e-4;

% Fit deep neural network for each hyperparameter
net_list = cell(length(Nlayers_list),length(Nnodes_list));
mse_list = zeros(length(Nlayers_list), length(Nnodes_list));
Memory_dnn_kb = zeros(length(Nlayers_list), length(Nnodes_list));
index = 0;
for i = 1:length(Nlayers_list)
    for j = 1:length(Nnodes_list)
        fprintf('\n Training %d of %d ...', [length(Nnodes_list)*(i-1)+j, length(Nlayers_list)*length(Nnodes_list)])
        index = index + 1;
        net = feedforwardnet(Nnodes_list(j)*ones(1,Nlayers_list(i)), 'trainlm');
        for l = 1:Nlayers_list(i)
            net.layers{1}.transferFcn = 'poslin';
        end
        [net,tr] = train(net, x, t);
        net_list{i,j} = net;
        mse_list(i,j) = tr.best_perf;
        
        % Calculate memory
        M = Nnodes_list(j);
        L = Nlayers_list(i);
        ninput = nx+ny;
        noutput = nu;
        Mdnn = (ninput+1)*M + (L-1)*(M+1)*M + (M+1)*noutput;
        Memory_dnn_kb(i,j) = Mdnn*8/1e3;   
        fprintf('Done!')
    end
end


% Save variables used later so that we can run the DNN-controller in a separate file
save('Supporting-Data-Files/DNN_training.mat','net','xscale_min','xscale_max', 'tscale_min', 'tscale_max', ...
    'x_min', 'x_max', 'u_min', 'u_max', 'A', 'B', 'C', 'nx', 'nu', 'ny', 'X', 'U', 'Q', 'R', 'PN', 'K', 'mse_list')

% save(['Supporting-Data-Files/MSE_Ns_', num2str(Nsamp), '.mat'], 'Nlayers_list', 'Nnodes_list', 'mse_list')



%% Simulations/experiments performed in a different script to avoid repeating the training every time
%{
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




%% Perform simulations to check results



% calculate offset gain
Hd = C*inv(eye(nx)-(A+B*K));
lambdaf = 0.5;

% number of simulations
Nsim = 40;

% initialize
Xsim = zeros(nx,Nsim+1);
Ysim = zeros(ny,Nsim+1);
Usim = zeros(nu,Nsim);
Wsim = zeros(nx,Nsim);
What = zeros(nx,Nsim);
Sdes = zeros(ny,Nsim+1);

% define reference
Sdes = [zeros(nx, 20), [4.9;2].*ones(nx, Nsim+1-20)];

% initial states
Xsim(:,1) = [1;0];
Ysim(:,1) = C*Xsim(:,1);

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
    if useProj
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

    % get most recent disturbance
%     Wsim(:,k) = [(w_upper-w_lower).*rand(1,nx)+w_lower]';
    Wsim(:,k) = [0.5;0.5];
    
    % provide values to plant
    Xsim(:,k+1) = A*Xsim(:,k) + B*Usim(:,k) + Wsim(:,k);
    Ysim(:,k+1) = C*Xsim(:,k+1);
    
    % estimate disturbance
    What(:,k+1) = lambdaf*What(:,k) + (1-lambdaf)*Wsim(:,k);
end

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

%}
