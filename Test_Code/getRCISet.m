
clear

% Load identified model matrices for 
% x(k+1) = Ax(k) + Bu(k)
% y(k) = Cx(k) + Du(k)
directory = pwd;
modelGlass = load([directory, '/Supporting-Data-Files/MIMOmodelGlass']);
A = modelGlass.A; B=modelGlass.B; C=modelGlass.C; D=0;

% Load constraints (loading robust constraints just as a placeholder. 
% We will only use the first, i.e non-tightened, constraints of these)
Constraints = load([directory,'/supporting-data-files/robustConstraintsWorstCase']);
% State constraints
Xcon = Constraints.Xtight{1};
% Input constraints
Ucon = Constraints.Ucon;
% Terminal constraints
Xf = Constraints.Xtight{1};

% Bounds on w
w_upper = [0.8; 2.0];
w_lower = [-0.8; -2.0];

% Create constraints in MPT 
X = Polyhedron('A',Xcon(:,1:2),'b',Xcon(:,3));
U = Polyhedron('A',Ucon(:,1:2),'b',Ucon(:,3));
W = Polyhedron('lb',w_lower,'ub',w_upper);

% Compute robust control invariant set
sys = ULTISystem('A',A,'B',B,'E',eye(2));
sys.x.with('setConstraint');
sys.x.setConstraint = X;
sys.u.with('setConstraint');
sys.u.setConstraint = U;
sys.d.min = -Inf*ones(2,1);
sys.d.max = Inf*ones(2,1);
sys.d.with('setConstraint');
sys.d.setConstraint = W;
Cinf = sys.invariantSet;
