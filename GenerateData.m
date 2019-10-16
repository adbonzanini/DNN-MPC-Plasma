clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE PARAMETERS AND DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steady-states
% y1ss = 34.5; %T
% y2ss = 75;   %I
% u1ss = 1.5;  %q
% u2ss = 2.0;  %P

% Pick a value for the lag
L = 1;
% Choose if you want to overwrite the datasets
overwriteData =1;

% Pick the standard deviation of the noise
x1sd = 0.1; %0.5
x2sd = 0.1; %0.6

for Case=[1,2]
rng(Case+1)
%1 for test data 2 for training data
TestOrTrain = Case;
% Time between inputs
DTinput = 200;
% How many data points to skip in the dataset (to reduce the number of points)
Nskip = 25;
% Pick constant offset (if any)
offset = [0;0];

%Load identified model  (A, B, C matrices). This is the "real"/plant model
load('Supporting Data Files/MIMOmodelMetal')
Areal = A;
Breal = B;

%Modified model (used for MPC predictions)
load('Supporting Data Files/MIMOmodelGlass')
A = A;
B = B;

% Check that the system is stabilizable
fprintf('\n')
lambda = eig(Areal);
if abs(lambda) >= 1
    disp(['eig(A) >= 1 => Expect Full Rank = ', num2str(size(A,1))])
elseif abs(lambda)<1
    disp('eig(A) < 1 => Deficient Rank OK')
end
cMat = [Areal-lambda*eye(1), Breal];
disp(['rank = ', num2str(rank(cMat))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE TRAINING OR TEST INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create appropriate file (test vs. train)
if TestOrTrain == 1
    filename = 'TestDataSimulation.csv';
    % Create inputs

%     u1 = [0, -0.5, 0.5, 1, 2, 3];
%     u2 = [0, 1, -0.5];
    %aa=load('/users/adbonzanini/Box Sync/Berkeley/Research/Plasma Code/Collect GP Data/inputTrajectories20180503');
    u1 = [0, 1, -1];
    u2 = [0, 1, -1, 2];
    u = combvec(u1, u2);
    % Shuffle while keeping first and last columns
    shuffleVec = eye(size(u,2)-2);
    u(:,2:end-1) = u(:,2:end-1)*shuffleVec(randperm(size(shuffleVec,1)),:);
   
    % Number of data points
    Ndata = size(u,2)*DTinput;
    udata = imresize(u, [2 Ndata], 'nearest');
    % Skip some entries
    udata = udata(:,end:-1:1);
    
    % Initial state
    % syms xx1 xx2
    % ss = solve([xx1;xx2]-A*[xx1;xx2]-B*[3;-0.5], [xx1;xx2])
    % double(ss.xx1)
    % double(ss.xx2)
    ydata(:,1) = [-0.0886;-0.1555]; %By solving x-Ax-Bu==0
    
%     Nskip=1;
    
else
    filename = 'TrainingDataSimulation.csv';
    % Create inputs
    
    %{
    u1 = [0, 0.5, 1, -0.5, 2, 3];
    u2 = [0, 1, -0.5];
    u = combvec(u1, u2);
    % Shuffle while keeping first and last columns
    shuffleVec = eye(size(u,2)-2);
    u(:,2:end-1) = u(:,2:end-1)*shuffleVec(randperm(size(shuffleVec,1)),:);
    %}
    
    u1 = [0, 2.5, -1, 1];
    u2 = [0, -1, 1, 1.5, 2.5];
    u = combvec(u2, u1);
    u = u(end:-1:1,:);
    % Shuffle while keeping first and last columns
    shuffleVec = eye(size(u,2)-2);
    u(:,2:end-1) = u(:,2:end-1)*shuffleVec(randperm(size(shuffleVec,1)),:);
    
    % Number of data points
    Ndata = size(u,2)*DTinput;
    udata = imresize(u, [2 Ndata], 'nearest');
    
    % Skip some entries
    udata = udata(:,1:1:end);
    
    % Initial state
%     ydata(:,1) = [0;0];
    ydata(:,1) = [0;0]; %By solving x-Ax-Bu==0
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLANT MODEL TO GENERATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate Noise
rng(1)   
w1 = normrnd(0, x1sd, 1, Ndata);
w2 = normrnd(0, x2sd, 1, Ndata);
w = [w1;w2];

%Run the model
xModel = zeros(2, Ndata);
xModel(:,1)= C\ydata(:,1);
yModel = zeros(2, Ndata);
yModel(:,1) = C*xModel(:,1);
yPseudo = zeros(2,Ndata);
yPseudo(:,1) = yModel(:,1) + w(:,1);
yReal = zeros(2, Ndata);
yReal(:,1) = C*xModel(:,1)+offset;

for j=1:Ndata-1
    xModel(:,j+1)=A*xModel(:,j)+B*udata(:,j);
    yModel(:,j+1)=C*xModel(:,j+1);
    yPseudo(:,j+1) = yModel(:,j+1) + w(:,j+1);
    yReal(:,j+1) = C*(Areal*xModel(:,j)+Breal*udata(:,j))+w(:,j+1)+offset;
end


%Select larger timesteps between data
% yPseudo = yPseudo(:,1:Nskip:end);
% yModel = yModel(:,1:Nskip:end);
% yReal = yReal(:,1:Nskip:end);
% udata = udata(:,1:Nskip:end);
%  
% Ndata = size(yPseudo,2);


yPseudo1 = [yPseudo(1,1:Nskip:end-1); yPseudo(1,2:Nskip:end)];
yPseudo1 = yPseudo1(:)';
yPseudo2 = [yPseudo(2,1:Nskip:end-1); yPseudo(2,2:Nskip:end)];
yPseudo2 = yPseudo2(:)';
yPseudo = [yPseudo1;yPseudo2];
  
yModel1 = [yModel(1,1:Nskip:end-1);yModel(1,2:Nskip:end)];
yModel1 = yModel1(:)';
yModel2 = [yModel(2,1:Nskip:end-1);yModel(2,2:Nskip:end)];
yModel2 = yModel2(:)';
yModel = [yModel1;yModel2];
 
yReal1 = [yReal(1,1:Nskip:end-1);yReal(1,2:Nskip:end)];
yReal1 = yReal1(:)';
yReal2 = [yReal(2,1:Nskip:end-1);yReal(2,2:Nskip:end)];
yReal2 = yReal2(:)';
yReal = [yReal1;yReal2];
 
udata1 = [udata(1,1:Nskip:end-1);udata(1,2:Nskip:end)];
udata1 = udata1(:)';
udata2 = [udata(2,1:Nskip:end-1);udata(2,2:Nskip:end)];
udata2 = udata2(:)';
udata = [udata1;udata2];
 
% Update the number of data points
Ndata = size(yPseudo,2)/2;



% Plot model identification fitting
 figure(Case)
 subplot(2, 2, 1)
 plot(yPseudo(1,:), 'b')
 hold on
 plot(yReal(1,:), 'r')
 plot(yModel(1,:), 'g')
 legend('MPC data', 'Plant Model', 'MPC Model')
 ylabel('T-T_{ss}/ \circ C')
 xlabel('Time step') 
 title(filename)
 
 
 subplot(2, 2, 2)
 plot(yPseudo(2,:), 'b')
 hold on
 plot(yReal(2,:), 'r')
 plot(yModel(2,:), 'g')
 legend('MPC data', 'Plant Model', 'MPC Model')
 ylabel('I-I_{ss}/ \circ C')
 xlabel('Time step')
 
 subplot(2, 2, 3)
 plot(udata(1,:))
 ylim([min(udata(1,:))-0.5, max(udata(1,:))+0.5])
 ylabel('q-q_{ss}/slm')
 xlabel('Time step')
 
 subplot(2, 2, 4)
 plot(udata(2,:))
 ylim([min(udata(2,:))-0.5, max(udata(2,:))+0.5])
 ylabel('P-P_{ss}/W')
 xlabel('Time step')
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STORE THE DATA IN AN APPROPRIATE FORM FOR ITERATIVE TRAINING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define deviation from the model and the data (to be captured by ML)
% x(k+1) = Ax(k) + Bu(k) + g(x,u)
yDev = (yReal-yModel)';
yDev = yDev(1:2:end, :);

% Pre-allocation for speed
X = zeros(Ndata-L, 4*L);
count=1;

% Store the data along with L previous data points
for j=1+L:Ndata
    %First L columns are y and the remaining are u
    X(count, :) = [yDev(j-L:j-1, 1)', yDev(j-L:j-1, 2)', udata(1,j-L:j-1), udata(2, j-L:j-1)];
    count = count+1;
end

rng(1)
X = [normrnd(0, x1sd, 10,1), normrnd(0, x2sd, 10,1), 0.3*ones(10,1), 0.4*ones(10,1);X];
Y = [normrnd(0, x1sd, 10,1), normrnd(0, x2sd, 10,1); yDev(L+1:end, 1), yDev(L+1:end, 2)];

time = L+1:1:Ndata+10;

Ndata = size(X,1)+1;
time = L+1:1:Ndata;

%Combine the data into a single matrix
exportData = [-999*ones(1,1+4*L+2);time', X, Y];

fprintf('\n')
disp(['Number of data points = ', num2str(Ndata)])
fprintf('\n')

filename = ['/users/adbonzanini/Box Sync/Berkeley/Research/Gaussian Process/MIMO GP State Feedback Substrates/Supporting Data Files/' filename];
% Save datasets
if overwriteData==1
    csvwrite(filename, exportData)
    warning('Datasets overwritten!')
    fprintf('\n')
end

end



