
%% Import required packages

% Clear variables
uqlab
clearvars

% Import casadi
import casadi.*



%% Setup nmpc problem and input objects

% NMPC problem parameters
Ts = 1.3;    % sampling time
N = 10;       % number of control intervals

% State constraints
x_min = [-6, 0];
x_max = [6, 15];
x_init = [0, 0];

% Input constraints
u_min = [-1.2, -2.5];
u_max = [7, 2];
u_init = [0, 0];

% Output constraints
y_min = x_min;
y_max = x_max;
y_init = x_init;

% Get NMPC solver
[solver, args, f] = getNMPCSolver(Ts, N, x_init, x_min, x_max, u_init, u_min, u_max, y_init, y_min, y_max);
args.warm_start = 0;
nx = args.nx;
ny = args.ny;
nu = args.nu;

% Marginals for states
for i = 1:nx
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [x_min(i), x_max(i)];    
end

% Create state "input" object
myInput_X = uq_createInput(Input);

% Marginals for the reference
for i = 1:ny
    Input.Marginals(nx+i).Type = 'Uniform';
    Input.Marginals(nx+i).Parameters = [y_min(i), y_max(i)];
end

% Create total "input" object (both states and reference)
myInput_P = uq_createInput(Input);



%% Build svm classifier for domain of attraction iteratively

% User specified parameters
N_init = 250;
N_add = 250;
missclass_tol = 0.01;
prob_tol = 0.005;
relVol_tol = 0.05;
ytarget = 3;
N_max = 2500;

% Get initial samples 
X0_train = uq_getSample(myInput_X, N_init, 'Sobol');
Ytarget_train = repmat(ytarget, [N_init,1]);
CEMtarget_train = 1*ones(N_init,1);

% Solve nmpc problem over samples
[~, Feas_train, ~] = solveSampleNMPC(solver, args, X0_train, CEMtarget_train);

if isempty(find(Feas_train~=1))
    Feas_train(end)=0;
else
end
% Train svm to represent doa
svm_doa = fitcsvm(X0_train,Feas_train,'KernelFunction','rbf','Standardize',true);

% Fit posterior
[svm_doa,~] = fitPosterior(svm_doa); 

% Estimate uncertainty in decision region
X0_posterior = uq_getSample(myInput_X, 10000, 'MC');
[~,prob_posterior] = predict(svm_doa,X0_posterior);
index = find(abs(prob_posterior(:,1)-0.5) < 0.5-missclass_tol);
relVol_missclass = length(index)/size(X0_posterior,1)

% Iterate until tolerance is met
X0_list = X0_train;
N_iter = 1;
while relVol_missclass > relVol_tol && size(X0_train,1) < N_max
    % Get new samples
    X0_new = uq_enrichSobol(X0_list, N_add, myInput_X);
    X0_list = [X0_list ; X0_new];
    
    % Find samples near decision boundary
    [~,prob_new] = predict(svm_doa, X0_new);
    index = find(abs(prob_new(:,1)-0.5) < 0.5-prob_tol);
    X0_add = X0_new(index,:);
    Ytarget_add = repmat(ytarget, [length(index),1]);
    
    % Solve nmpc problem at new samples
    [~, Feas_add, ~] = solveSampleNMPC(solver, args, X0_add, Ytarget_add);
    
    % Add to training set
    X0_train = [X0_train ; X0_add];
    Feas_train = [Feas_train ; Feas_add];
    
    % Retrain svm with new dataset
    svm_doa = fitcsvm(X0_train,Feas_train,'KernelFunction','rbf','Standardize',true);
    
    % Fit posterior
    [svm_doa,~] = fitPosterior(svm_doa);

    % Estimate uncertainty in decision region
    [~,prob_posterior] = predict(svm_doa,X0_posterior);
    index = find(abs(prob_posterior(:,1)-0.5) < 0.5-missclass_tol);
    relVol_missclass = length(index)/size(X0_posterior,1)
    
    % Add to counter
    N_iter = N_iter + 1;
end

X = X0_train;
Y = Feas_train;

[x1Grid,x2Grid] = meshgrid(0:0.01:1,280:1:370);
[~,PosteriorRegion] = predict(svm_doa,[x1Grid(:),x2Grid(:)]);
figure;
contourf(x1Grid,x2Grid,...
        reshape(PosteriorRegion(:,2),size(x1Grid,1),size(x1Grid,2)),20);
h = colorbar;
h.Label.String = 'P({\it{versicolor}})';
h.YLabel.FontSize = 16;
caxis([0 1]);
colormap jet;

hold on
index_class1 = find(Y==-1);
index_class2 = find(Y==1);
sorted_data = [X(index_class1,:) ; X(index_class2,:)];
sorted_class = [-1*ones(length(index_class1),1) ; 1*ones(length(index_class2),1)];
h(1:2) = gscatter(sorted_data(:,1),sorted_data(:,2),sorted_class,'rb','.',10);
index = find(abs(prob_posterior(:,1)-0.5) < 0.5-missclass_tol);
scatter(X0_posterior(index,1), X0_posterior(index,2),'gx')



%% Fit explicit controller using samples in the svm set

% User specified parameters
N_init = 3000;
N_add = 3000;

% Sample the state/reference space
P_full = uq_getSample(myInput_P, N_init, 'Sobol');

% Find the samples inside the doa
[~,prob_posterior] = predict(svm_doa,P_full(:,1:nx));
index = find(prob_posterior(:,2) > 1-missclass_tol);
P_mpc = P_full(index,:);

% Solve nmpc problem over samples
[U_mpc, Feas_mpc, ~] = solveSampleNMPC(solver, args, P_mpc(:,1:nx), P_mpc(:,nx+1:end));

% Define new names
x = P_mpc';
t = U_mpc';
x(2,:) = (x(2,:)-280)/(370-280);
x(3,:) = (x(3,:)-280)/(370-280);
t(1,:) = (t(1,:)-280)/(370-280);

% Fit dnn
NL = 5;
Nnodes = 10;
net = feedforwardnet(Nnodes*ones(1,NL), 'trainlm');
net = train(net, x, t);
y = net(x);

% Run closed-loop simulation
Nsim = 1000;
args.warm_start = 1;
X = zeros(nx,Nsim+1);
Xexp = zeros(nx,Nsim+1);
U = zeros(nu,Nsim);
Uexp = zeros(nu,Nsim);
Yt = zeros(ny,Nsim);
Feas_cl = zeros(1,Nsim);
X(:,1) = [0.95 ; 281];
Xexp(:,1) = X(:,1);
for i = 1:Nsim
    % Print start statement
    fprintf('Step %g of %g...', i, Nsim)
    tic
    
    % Define target
    if i < 500
        Yt(:,i) = 350;
    else
        Yt(:,i) = 370;
    end
    
    % Get feedback
    [U(:,i), Feas_cl(i), ~, args, ~] = getFeedback(solver, args, X(:,i), Yt(:,i));
    
    % Use the explicit controller
    xscale = Xexp(:,i);
    xscale(2) = (xscale(2)-280)/(370-280);
    yscale = (Yt(:,i)-280)/(370-280);
    uscale = net([xscale ; yscale]);
    Uexp(:,i) = 280 + (370-280)*uscale;
    
    % Print end statement
    fprintf('took %g seconds\n', toc)
    
    % Give optimal input to system
    X(:,i+1) = full(f(X(:,i), U(:,i)));
    Xexp(:,i+1) = full(f(Xexp(:,i), Uexp(:,i)));
end

% Plots
figure
plot(X(1,:), X(2,:), '-k')
hold on
plot(Xexp(1,:), Xexp(2,:), '--k')

figure
stairs(1:Nsim,U(1,:), '-b')
hold on
stairs(1:Nsim,Uexp(1,:), '-r')


% % Sample wrt states only (domain of attraction is independent of ytarget)
% N = 1000;
% X0 = uq_getSample(myInput_X, N, 'Sobol');
% Ytarget = 0.5*ones(N,ny);
% 
% % Solve nmpc problem over samples
% [U_mpc, Feas, V_opt] = solveSampleNMPC(solver, args, X0, Ytarget);
% 
% % Train SVM classifier
% X = X0;
% Y = Feas;
% SVMModel = fitcsvm(X,Y,'KernelFunction','rbf','Standardize',true);
% 
% % Fit posterior
% [SVMModel,ScoreParameters] = fitPosterior(SVMModel); 
% 
% % Cross-validation
% CVModel = crossval(SVMModel);
% misclass1 = kfoldLoss(CVModel)
% 
% % Plot classes
% index_class1 = find(Y==-1);
% index_class2 = find(Y==1);
% sorted_data = [X(index_class1,:) ; X(index_class2,:)];
% sorted_class = [-1*ones(length(index_class1),1) ; 1*ones(length(index_class2),1)];
% figure;
% h(1:2) = gscatter(sorted_data(:,1),sorted_data(:,2),sorted_class,'rb','.',10);
% hold on
% 
% % Run closed-loop simulation
% Nsim = 100;
% args.warm_start = 1;
% X = zeros(nx,Nsim+1);
% U = zeros(nu,Nsim);
% Yt = zeros(ny,Nsim);
% Feas_cl = zeros(1,Nsim);
% X(:,1) = [0.95 ; 281];
% for i = 1:Nsim
%     % Print start statement
%     fprintf('Step %g of %g...', i, Nsim)
%     tic
%     
%     % Define target
%     if i < 500
%         Yt(:,i) = 350;
%     else
%         Yt(:,i) = 370;
%     end
%     
%     % Get feedback
%     [U(:,i), Feas_cl(i), ~, args] = getFeedback(solver, args, X(:,i), Yt(:,i));
%     
%     % Print end statement
%     fprintf('took %g seconds\n', toc)
%     
%     % Give optimal input to system
%     X(:,i+1) = full(f(X(:,i), U(:,i)));    
% end
% 
% % Plots
% plot(X(1,:), X(2,:), '-k')
