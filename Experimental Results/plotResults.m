%
%   This script loads the csv file with the measurements collected
%   from the setup and plots the state and input trajectories.
%

%%

directory = pwd;

% extract file names
files = dir(directory);

%% Constraints
cd ../Supporting-Data-Files
load('DNN_training.mat');
model_ID=load('MIMOmodelGlass.mat');
cd(directory)
steadyStates = round(model_ID.steadyStates, 1);
Tss = steadyStates(1); Iss = steadyStates(2); qss = steadyStates(3); Pss = steadyStates(4); 
% u_max = [10,11]-[qss,Pss];
% u_min = [0.5, 1]-[qss,Pss];

%% Load files
% Keep only PI_Server[...] files
idx=[];
for j=1:length(files)
    try
        if files(j).name(1:5)=='PI_Se'
            idx=[idx;j];
        end
    catch
    end
end
files=files(idx, :);
Nfiles = length(files);


data = cell(Nfiles,1);
T = cell(Nfiles,1);I = cell(Nfiles,1);q = cell(Nfiles,1);P = cell(Nfiles,1);
N = zeros(Nfiles,1);
count=1;
CEMall = cell(1,3);

for j = jvals
    data{j} = csvread(files(j).name,1,0);
    data{j} = data{j}(1:35,:);
    N = size(data{j}, 1);
    
    
    % Column legend for reference
    %{
    (1) time,(2) Tset,(3) Ts,(4) Ts2,(5) Ts3, (6) P, (7) Imax, (8) Ip2p, 
    (9) O777, (10) O845, (11) N391, (12) He706, (13) sum_int, 
    (14, 15, 16, 17) *U_m --> (V, freq, q, dsep), (18) q_o, (19) D_c, (20) x_pos, 
    (21) y_pos, (22) T_emb, (23) Pset, (24) P_emb, (25) Prms, 
    (26) Rdel, (27) Is, (28, 29) sig --> (1 and 2), (30) subs_type, (31) Trot, 
    (32) tm_el
    %}
    varIdx = [3, 27, 16, 23]; %[T, I, q, P]
    variables = data{j}(:,varIdx);
    T{j} = variables(:,1); I{j} = variables(:,2); q{j} = variables(:,3); P{j} = variables(:,4);
    
    
    %% Other parameters
    Tsampling = 1.3;
    tPlot = 1:Tsampling:N*Tsampling;
    %% Desired CEM Reference
    Sdes = 1.5*ones(1,N);


    % Calculate CEM
    CEM = zeros(1, N);
    for k=1:N-1
        if T{j}(k)<35
            CEM(k+1) = CEM(k);
        else
            CEM(k+1) = CEM(k)+0.5.^(43-T{j}(k));
        end
    end

    CEMall{count} = CEM;


    if plotWithinScript==1
    %% Plot states
    figure(1)
    hold on
    h1 = plot(tPlot, CEM, 'Linewidth', 2);
    h2 = plot(tPlot, Sdes, 'k', 'Linewidth', 1);
    xlabel('Time/ s')
    ylabel('CEM/ min')
    xlim([tPlot(1), tPlot(end)])
    legend([h2], 'CEM setpoint', 'Location', 'southeast')
    set(gca,'FontSize',15)
    ylim([0, 4.5])
    
    figure(2)
    subplot(2,1,1)
    hold on
    h1 = plot(tPlot, T{j}, 'Linewidth', 2);
    h2 = plot([tPlot(1), tPlot(end)], [x_max(1), x_max(1)]+Tss, 'k--');
    h3 = plot([tPlot(1), tPlot(end)], [x_min(1), x_min(1)]+Tss, 'k--');
    xlabel('Time/s')
    ylabel('T/ ^{\circ}C')
    xlim([tPlot(1), tPlot(end)])
    legend([h2], 'Constraints', 'Location', 'southwest')
    set(gca,'FontSize',15)
    
    figure(2)
    subplot(2,1,2)
    hold on
    h1 = plot(tPlot, I{j}, 'Linewidth', 2);
    h2 = plot([tPlot(1), tPlot(end)], 10*([x_max(2), x_max(2)]+Iss), 'k--');
    h3 = plot([tPlot(1), tPlot(end)], 10*([x_min(2), x_min(2)]+Iss), 'k--');
    xlabel('Time/s')
    ylabel('Intensity/ a.u.')
    xlim([tPlot(1), tPlot(end)])
    legend([h2], 'Constraints', 'Location', 'northwest')
    set(gca,'FontSize',15)
    
    figure(3)
    subplot(2,1,1)
    hold on
    h1 = plot(tPlot, q{j}, 'Linewidth', 2);
%     h2 = plot([tPlot(1), tPlot(end)], [u_max(1), u_max(1)]+qss, 'k--');
%     h3 = plot([tPlot(1), tPlot(end)], [u_min(1), u_min(1)]+qss, 'k--');
    xlabel('Time/s')
    ylabel('He Flowrate/ slm')
    xlim([tPlot(1), tPlot(end)])
%     legend([h2], 'Constraints', 'Location', 'southeast')
    set(gca,'FontSize',15)
    
    figure(3)
    subplot(2,1,2)
    hold on
    h1 = plot(tPlot, P{j}, 'Linewidth', 2);
%     h2 = plot([tPlot(1), tPlot(end)], [u_max(2), u_max(2)]+Pss, 'k--');
%     h3 = plot([tPlot(1), tPlot(end)], [u_min(2), u_min(2)]+Pss, 'k--');
    xlabel('Time/s')
    ylabel('Applied Power/ W')
    xlim([tPlot(1), tPlot(end)])
%     legend([h2], 'Constraints', 'Location', 'northwest')
    set(gca,'FontSize',15)
    end
    
    count=count+1;
end


