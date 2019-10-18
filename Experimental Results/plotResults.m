%
%   This script loads the csv file with the measurements collected
%   from the setup and plots the state and input trajectories.
%

%%
clear all

directory = pwd;

%% Load files
fnameProj = ['PI_Server_Out_2019-10-17_161504.889852-explicit2';
             'PI_Server_Out_2019-10-17_162329.360513-explicit2';
             'PI_Server_Out_2019-10-17_163158.799811-explicit2'];

fnameNoProj= ['PI_Server_Out_2019-10-17_161915.729781-explicit2';
              'PI_Server_Out_2019-10-17_162743.399614-explicit2';
              'PI_Server_Out_2019-10-17_163613.639705-explicit2'];

% Initialize
dataWithProjection = cell(3,1);
dataNoProjection = cell(3,1);
N = zeros(3,1);
for j = 1:3
    dataWithProjection{j} = csvread(fnameProj(j,:),1,0);
    dataNoProjection{j} = csvread(fnameNoProj(j,:),1,0);
    N(j) = min(size(dataWithProjection{j},1), size(dataNoProjection{j},1));
end
N = min(N);


%% Column legend for reference
%{
(1) time,(2) Tset,(3) Ts,(4) Ts2,(5) Ts3, (6) P, (7) Imax, (8) Ip2p, 
(9) O777, (10) O845, (11) N391, (12) He706, (13) sum_int, 
(14, 15, 16, 17) *U_m --> (V, freq, q, dsep), (18) q_o, (19) D_c, (20) x_pos, 
(21) y_pos, (22) T_emb, (23) Pset, (24) P_emb, (25) Prms, 
(26) Rdel, (27) Is, (28, 29) sig --> (1 and 2), (30) subs_type, (31) Trot, 
(32) tm_el
%}

% Eliminate small differences in length of datasets
for j=1:3
    dataWithProjection{j} = dataWithProjection{j}(1:N,:);
    dataNoProjection{j} = dataNoProjection{j}(1:N,:);
end

%% Extract relevant data
varP = cell(4,1);
varNP = cell(4,1);
count = 1;
for j =[3, 27, 16, 23]
    varP{count} = mean([dataWithProjection{1}(:,j), dataWithProjection{2}(:,j), dataWithProjection{3}(:,j)],2);
    varNP{count} = mean([dataNoProjection{1}(:,j), dataNoProjection{2}(:,j), dataNoProjection{3}(:,j)],2);
    count = count+1;
end
% Temperature
Tp = varP{1};
Tnp = varNP{1};

% Intensity
Ip = varP{2};
Inp = varNP{2};

% He Flow
qp = varP{3};
qnp = varNP{3};

% Power
Pp = varP{4};
Pnp = varNP{4};

%% Constraints
cd ../Supporting-Data-Files
load('DNN_training.mat');
model_ID=load('MIMOmodelGlass.mat');
cd(directory)


%% Desired Reference
tChange = 80;
Sdes = [zeros(nx, tChange), [5.5;2].*ones(nx, N-tChange)];

%% Other parameters
Tsampling = 1.3;
tPlot = 1:Tsampling:N*Tsampling;

%% Plot states
figure(1)
subplot(2,1,1)
hold on
h1 = plot(tPlot, Tp, 'b', 'Linewidth', 2);
h2 = plot(tPlot, Tnp, 'r', 'Linewidth', 2);
h3 = stairs(tPlot, Sdes(1,:)+model_ID.steadyStates(1), 'k', 'Linewidth', 1);
h4 = plot([tPlot(1), tPlot(end)], [x_max(1), x_max(1)]+model_ID.steadyStates(1), 'k--');
h5 = plot([tPlot(1), tPlot(end)], [x_min(1), x_min(1)]+model_ID.steadyStates(1), 'k--');
legend([h1, h2, h3, h4], 'With Projection', 'Without Projection', 'Setpoint', 'Constraints', 'Location', 'southeast')
ylim([35, 45])
xlabel('Time/s')
ylabel('Temperature/ ^{\circ}C')
set(gca,'FontSize',15)

subplot(2,1,2)
hold on
h1 = plot(tPlot, Ip, 'b', 'Linewidth', 2);
h2 = plot(tPlot, Inp, 'r', 'Linewidth', 2);
h3 = stairs(tPlot, 10*(Sdes(2,:)+model_ID.steadyStates(2)), 'k', 'Linewidth', 1);
h4 = plot([tPlot(1), tPlot(end)], 10*([x_max(2), x_max(2)]+model_ID.steadyStates(2)), 'k--');
h5 = plot([tPlot(1), tPlot(end)], 10*([x_min(2), x_min(2)]+model_ID.steadyStates(2)), 'k--');
legend([h1, h2, h3, h4], 'With Projection', 'Without Projection', 'Setpoint', 'Constraints', 'Location', 'northwest')
ylim([0, 200])
xlabel('Time/s')
ylabel('Intensity/ a.u.')
set(gca,'FontSize',15)



%% Plot inputs
figure(2)
subplot(2,1,1)
hold on
h1 = plot(tPlot, qp, 'b', 'Linewidth', 2);
h2 = plot(tPlot, qnp, 'r', 'Linewidth', 2);
h4 = plot([tPlot(1), tPlot(end)], [u_max(1), u_max(1)]+model_ID.steadyStates(3), 'k--');
h5 = plot([tPlot(1), tPlot(end)], [u_min(1), u_min(1)]+model_ID.steadyStates(3), 'k--');
legend([h1, h2, h4], 'With Projection', 'Without Projection', 'Constraints', 'Location', 'northwest')
xlabel('Time/s')
ylabel('He Flowrate/ slm')
set(gca,'FontSize',15)


subplot(2,1,2)
hold on
h1 = plot(tPlot, Pp, 'b', 'Linewidth', 2);
h2 = plot(tPlot, Pnp, 'r', 'Linewidth', 2);
h4 = plot([tPlot(1), tPlot(end)], [u_max(2), u_max(2)]+model_ID.steadyStates(4), 'k--');
h5 = plot([tPlot(1), tPlot(end)], [u_min(2), u_min(2)]+model_ID.steadyStates(4), 'k--');
legend([h1, h2, h4], 'With Projection', 'Without Projection', 'Constraints', 'Location', 'southeast')
xlabel('Time/s')
ylabel('Applied Power/ W')
set(gca,'FontSize',15)


