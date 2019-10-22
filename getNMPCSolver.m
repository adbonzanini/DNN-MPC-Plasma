function [solver, args, f] = getNMPCSolver(Ts, N, x_init, x_min, x_max, u_init, u_min, u_max, CEM_init, CEM_min, CEM_max)

% Import casadi
import casadi.*

% Load identified model matrices for 
% x(k+1) = Ax(k) + Bu(k)
% y(k) = Cx(k) + Du(k)
% Reduced variables around xss = [37.4; 9.7] and uss = [2;2]
directory = pwd;
modelGlass = load([directory, '/Supporting-Data-Files/MIMOmodelGlass']);
A = modelGlass.A; B=modelGlass.B; C=modelGlass.C; D=0;

% Sizes
nx = size(A,2);
nu = size(B,2);

x = SX.sym('x', nx);
u = SX.sym('u', nu);

%Model Parameters
Kcem = 0.25;


% Model equations
xdot = A*x+B*u;

% Continuous time dynamics
f = Function('f', {x, u}, {xdot});

% Output function
CEM = Kcem^(43-x(1)-37);
h = Function('h', {x}, {CEM});
ny = length(CEM);

% Stage cost function (deviation variables)
R = 0;
CEMtarget = 1;
l = Kcem.^(43-x(1)-37) + u'*R*u;
% l = 1*x(1)^2 + 1/100*x(2)^2 + 1/10*u(1)^2;
L = Function('L', {x,u}, {l});

% Terminal cost function (deviation variables)
m = 0.0;
M = Function('M', {x}, {m});

% Offset cost function (deviation from target)
vO = 0*CEM;
VO = Function('VO', {x}, {vO});

% Start with an empty NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g = {};
lbg = [];
ubg = [];

% "Lift" the initial conditions
Xk = MX.sym('X0', nx);
w = [w, {Xk}];
lbw = [lbw, x_init];
ubw = [ubw, x_init];
w0 = [w0, x_init];

% Add NLP variables to represent the steady state
CEMcurr = MX.sym('CEMcurr', ny);
w = [w, {CEMcurr}];
lbw = [lbw, CEM_min];
ubw = [ubw, CEM_max];
w0 = [w0, CEM_init];

%{
Xss = MX.sym('Xss', nx);
w = [w, {Xss}];
lbw = [lbw, x_min];
ubw = [ubw, x_max];
w0 = [w0, x_init];
Uss = MX.sym('Uss', nu);
w = [w, {Uss}];
lbw = [lbw, u_min];
ubw = [ubw, u_max];
w0 = [w0, u_init];
Yss = MX.sym('Yss', ny);
w = [w, {Yss}];
lbw = [lbw, y_min];
ubw = [ubw, y_max];
w0 = [w0, CEM_init];

% Add steady state equations
g = [g, {Xss-f(Xss,Uss)}];
g = [g, {Yss-h(Xss,Uss)}];
lbg = [lbg, zeros(1,nx+ny)];
ubg = [ubg, zeros(1,nx+ny)];
%}

% Formulate the NLP
for k = 0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], nu);
    w = {w{:}, Uk};
    lbw = [lbw, u_min];
    ubw = [ubw, u_max];
    w0 = [w0, u_init];

    % Integrate till the end of the interval
    Xk_end = f(Xk, Uk);
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx);
    w = [w, {Xk}];
    lbw = [lbw, x_min];
    ubw = [ubw, x_max];
    w0 = [w0, x_init];

    % Add equality constraint
    g = [g, {Xk_end-Xk}];
    lbg = [lbg, zeros(1,nx)];
    ubg = [ubg, zeros(1,nx)];
    
    % Add contribution of stage cost to objective
    J = J + L(Xk,Uk);
end

% Add terminal and offset cost to objective
% J = J + M(Xk-Xss) + VO(Yss-CEMtarget);
J = (J+CEMcurr-CEMtarget).^2;


% Add terminal equality constraint
%{
g = [g, {Xk-Xss}];
lbg = [lbg, zeros(1,nx)];
ubg = [ubg, zeros(1,nx)];
%}
% Set solver options
options = struct;
options.ipopt.tol = 1e-8;
options.ipopt.print_level = 0;
options.print_time = 0;

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob, options);

% Store arguments
args.w0 = w0;
args.lbw = lbw;
args.ubw = ubw;
args.lbg = lbg;
args.ubg = ubg;
% args.offset_mpc = nx+ny+nx+nu+ny;
args.offset_mpc = nx+ny;
args.offset_x0 = 0;
args.offset_ytar = nx;
args.nx = nx;
args.ny = ny;
args.nu = nu;
args.warm_start = 1;

end