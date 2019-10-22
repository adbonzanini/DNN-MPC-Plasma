%
%
%   This script implements an offset-free MPC algorithm based on the paper 
%   by Alvarado et al. (2007). It uses a model for an atmospheric pressure
%   plasma jet derived on glass. The purpose is to be able to solve the OCP
%   for various initial states to eventually substitute this by a deep 
%   neural network.
%
%
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Activate offset-free action
OFswitch =1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cost matrices
Q = diag([15,1]);
R = 0.01*diag([0.5, 0.2]);
PN = Q;
nx = size(Q,1);
nu = size(R,1);

% Bounds on w (not needed here)
% w_upper = [1.8530; 2.2558];
% w_lower = [-1.3929; -3.4765]

% Load identified model matrices for 
% x(k+1) = Ax(k) + Bu(k)
% y(k) = Cx(k) + Du(k)
% Reduced variables around xss = [37.4; 9.7] and uss = [3;3]
directory = pwd;
modelGlass = load([directory, '/Supporting-Data-Files/MIMOmodelGlass']);
A = modelGlass.A; B=modelGlass.B; C=modelGlass.C; D=0;

% Initial point(s)
xi = [1;0];
yi = C*xi;

% Prediction/simulation horizon and number of starts
Np = 4;
N = 40;
Nstart = size(yi,2);

% Setpoint trajectory
ysp = [zeros(nx, 20), [4;2].*ones(nx, N+1-20)];

% Add offset in "real" plant to test offset-free action
offset = 0.2*ones(2,1);

% Infinite horizon LQR
[K,Pinf,Eig] = dlqr(A,B,Q,R);
K = -K;                         % Matlab vs. literature conventions

% Observer parameters
Hs = (C+D*K)*inv(eye(nx)+(A+B*K));
lf = 0.5;

% Load constraints (loading robust constraints just as a placeholder. 
% We will only use the first, i.e non-tightened, constraints of these)
Constraints = load([directory,'/supporting-data-files/robustConstraintsWorstCase']);
% State constraints
Xcon = Constraints.Xtight{1};
% Input constraints
Ucon = Constraints.Ucon;
% Terminal constraints
Xf = Constraints.Xtight{1};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE/PRE-ALLOCATE MPC MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uOptSim = zeros(nu, N);
fopt = zeros(1, N);
xk = zeros(nx, N+1);
yTraj = zeros(nx, N+1);
yTraj(:,1) = yi;
what = zeros(nx, N+1);
exitflag = zeros(1, N+1);

% Start timer
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPC LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:N
    
    % Update initial condition
    xk(:,k) = xi;
    
    % Update setpoint
    ysp_k = ysp(:,k)-Hs*what(:,k);

    % Define variables
    x = sdpvar(nx, Np+1);
    u = sdpvar(nu, Np);
    
    % Re-initialize cost and constraints
    J = 0;
    Cons = [x(:,1)-xk(:,k) == 0];
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OCP LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:Np
        
        [xNext, Lstage] = dynamics(x(:,j), u(:,j),  ysp_k, A, B, C, Q, R);
        
        % Update cost
        J = J + Lstage;
        
        % Update constraints
        Cons =[Cons; 
               Xcon(:, 1:2)*x(:,j)<=Xcon(:,3);
               Ucon(:,1:2)*u(:,j)<=Ucon(:,3);
               xNext-x(:,j+1)==0
               Xcon(:, 1:2)*x(:,j+1)<=Xcon(:,3)];
    end
    % Terminal cost and constraints
    yN = C*x(:,Np+1);
    J = J + (yN-ysp_k)'*PN*(yN-ysp_k);
    Cons = [Cons; Xf(:, 1:2)*x(:,Np+1)<=Xf(:,3)]; 
    
    % Create nlp
    options = sdpsettings('solver','quadprog', 'verbose', 0);
    sol = optimize(Cons, J, options);
    
    exitflag(k) = sol.problem;
    
    % Assign solution to appropriate vectors
    xopt = value(x);
    uopt = value(u);
    
    
    % Keep only the first control input
    uOptSim(:,k) = uopt(:,1);
    
    % Save the optimal cost
    fopt(k) = value(J);
    
    % One step ahead simulation to find new initial point
    xi = A*xk(:,k) + B*uopt(:,1)+[0.2;0.2];
    yk = C*xi;
    
    if OFswitch ==1
        % Estimate uncertainty
        what(:,k+1) = lf*what(k)+(1-lf)*(xi-A*xk(:,k)-B*uopt(:,1));
        color = 'r';
    else
        % Test what happens when no uncertainty is estimated by setting what=0
        what(:,k+1) = zeros(2,1);
        color = 'b';
    end
    
    % Store real trajectory
    yTraj(:,k+1) = yk;
    
end

% End timer
tMPC = toc;
disp(['Total MPC time = ', num2str(tMPC), 's'])
fprintf('\n')

% Check whether any step was infeasible
for k = 1:N+1
    if exitflag(k) ~= 0
        warning('Problem at step %2.f', (k))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ub = Xcon(1:end/2,1:2)\Xcon(1:end/2,3);
lb = Xcon(end/2+1:end,1:2)\Xcon(end/2+1:end,3);


figure(1)
% Tempearature trajectory
subplot(2,1,1)
plot(yTraj(1,:), color)
hold on
plot([1, N+1], [ub(1), ub(1)], 'k--')
plot([1, N+1], [lb(1), lb(1)], 'k--')
stairs(ysp(1,:), 'k')
xlim([1, N+1])
ylim([lb(1)+4, 1.1*ub(1)])
xlabel('Simulation Step')
ylabel('Temperature Deviation')
% Intensity trajectory
subplot(2,1,2)
plot(yTraj(2,:), color)
hold on
plot([1, N+1], [ub(2), ub(2)], 'k--')
plot([1, N+1], [lb(2), lb(2)], 'k--')
stairs(ysp(2,:), 'k')
xlim([1, N+1])
ylim([lb(2)+10, 1.1*ub(2)])
xlabel('Simulation Step')
ylabel('Intensity Deviation')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS (LINEAR MODEL AND COST)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xNext, Lstage] = dynamics(x, u, ysp, A, B, C, Q, R)

xNext = A*x + B*u;

% State feedback => x=y
xsp = ysp;                      
y = C*x;
usp = B\(-(A-eye(2,2))*xsp);    % Control input steady-state

Lstage = (y-ysp)'*Q*(y-ysp) + (u-usp)'*R*(u-usp);
end





