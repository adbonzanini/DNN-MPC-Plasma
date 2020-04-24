%
%   This script needs the plotResults.m script to run. Make sure they
%   are in the same folder!
%

% clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT AVERAGE TRAJECTORIES ON THE SAME GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontSz = 15;
cases = [1,2];
plotWithinScript=0;
color = {'r', [0 0.4470 0.7410]};
hC = cell(2,1);
hT = cell(2,1);
hQ = cell(2,1);

figure(4)
hold on
xlabel('Time/s')
ylabel('Average CEM/min')
set(gcf,'color','w');




for ii=cases
    
    if ii==1
        jvals = [1,2,3];
    else
        jvals =[4,5,6];
    end
    
    plotResults;
    
    tPlot = [0, tPlot(1:end-1)]; %shift to start from zero
    
    % Take averages
    CEMavg = mean([CEMall{1}; CEMall{2}; CEMall{3}]);
    CEMmax = max([CEMall{1}; CEMall{2}; CEMall{3}]);
    CEMmin = min([CEMall{1}; CEMall{2}; CEMall{3}]);
    
    Tavg = mean([T{jvals(1)}, T{jvals(2)}, T{jvals(3)}],2)';
    Tmax = max([T{jvals(1)}, T{jvals(2)}, T{jvals(3)}],[], 2)';
    Tmin = min([T{jvals(1)}, T{jvals(2)}, T{jvals(3)}],[], 2)';
    
    Iavg = mean([I{jvals(1)}, I{jvals(2)}, I{jvals(3)}],2)';
    Imax = max([I{jvals(1)}, I{jvals(2)}, I{jvals(3)}],[], 2)';
    Imin = min([I{jvals(1)}, I{jvals(2)}, I{jvals(3)}],[], 2)';
    
    qavg = mean([q{jvals(1)}, q{jvals(2)}, q{jvals(3)}],2)';
    qmax = max([q{jvals(1)}, q{jvals(2)}, q{jvals(3)}],[], 2)';
    qmin = min([q{jvals(1)}, q{jvals(2)}, q{jvals(3)}],[], 2)';
    
    Pavg = mean([P{jvals(1)}, P{jvals(2)}, P{jvals(3)}],2)';
    Pmax = max([P{jvals(1)}, P{jvals(2)}, P{jvals(3)}],[], 2)';
    Pmin = min([P{jvals(1)}, P{jvals(2)}, P{jvals(3)}],[], 2)';
    

    % Determine when the CEM setpoint is reached (would switch off the plasma
    % in practice)
    idx = find(CEMavg>=1.5);

%    CEMavg(idx) = CEMavg(idx(1));
%     Tavg(idx) = 35;
%     Iavg(idx) = Iavg(idx(1));
%     qavg(idx) = qavg(idx(1));
%     Pavg(idx) = Pavg(idx(1));
    idx = idx(1);

    % If the CEM exceeds the setpoint by a lot, interpolate
    if CEM(idx)>1.6
        CEMavg(idx) = (CEMavg(idx)+CEMavg(idx-1))/2;
        Tavg(idx) = (Tavg(idx)+Tavg(idx-1))/2;
        Iavg(idx) = (Iavg(idx)+Iavg(idx-1))/2;
        qavg(idx) = (qavg(idx)+qavg(idx-1))/2;
        Pavg(idx) = (Pavg(idx)+Pavg(idx-1))/2;
        tPlot(idx) = (tPlot(idx)+tPlot(idx-1))/2;
    end
    

    % Plot CEM
    figure(4)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [CEMmin(1:idx), fliplr(CEMmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
%     plot(tPlot(1:idx),CEMavg(1:idx), color{ii}, 'LineWidth', 2)
    hC{ii} = plot(tPlot(1:idx),CEMavg(1:idx), '-o', 'color', color{ii}, 'LineWidth', 2);
    stairs(tPlot(1:idx),Sdes(1,1:idx),'k', 'LineWidth', 2)
    ylim([0, 2])
    set(gca,'FontSize',fontSz)

    
    
    % Plot states
    figure(5)
    subplot(2,1,1)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Tmin(1:idx), fliplr(Tmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    hT{ii} = plot(tPlot(1:idx),Tavg(1:idx), '-o', 'color', color{ii}, 'LineWidth', 2);
%     plot(tPlot,Tavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [x_max(1), x_max(1)]+Tss, 'k--')
    plot([tPlot(1), tPlot(idx)], [x_min(1), x_min(1)]+Tss, 'k--')
    hold on
    xlabel('Time/s')
    ylabel('Average T/ ^{\circ}C')
    set(gcf,'color','w');
    set(gca,'FontSize',fontSz)

    
    subplot(2,1,2)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Imin(1:idx), fliplr(Imax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),Iavg(1:idx), '-o', 'color', color{ii}, 'LineWidth', 2)
%     plot(tPlot, Iavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], 10*([x_max(2), x_max(2)]+Iss), 'k--')
    plot([tPlot(1), tPlot(idx)], 10*([x_min(2), x_min(2)]+Iss), 'k--')
    ylim(10*[x_min(2)+Iss-2, x_max(2)+Iss+2])
    xlabel('Time/s')
    ylabel('Average I/ a.u.')
    set(gcf,'color','w');
    set(gca,'FontSize',fontSz)

    
    
    % Plot inputs
    figure(6)
    subplot(2,1,1)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [qmin(1:idx), fliplr(qmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    hQ{ii} = plot(tPlot(1:idx),qavg(1:idx), '-o', 'color', color{ii}, 'LineWidth', 2);
%     plot(tPlot, qavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [10, 10], 'k--')
    plot([tPlot(1), tPlot(idx)], [0.5, 0.5], 'k--')
    xlabel('Time/s')
    ylabel('Average q/ slm')
    set(gcf,'color','w');
    set(gca,'FontSize',fontSz)

    
    subplot(2,1,2)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Pmin(1:idx), fliplr(Pmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),Pavg(1:idx), '-o', 'color', color{ii}, 'LineWidth', 2)
%     plot(tPlot,Pavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [5, 5], 'k--')
    plot([tPlot(1), tPlot(idx)], [1,1], 'k--')
    xlabel('Time/s')
    ylabel('Average Appled Power/ W')
    set(gcf,'color','w');
    set(gca,'FontSize',fontSz)

end

legend([hC{1}, hC{2}], 'PNN-based NMPC', 'DNN-based NMPC', 'Location', 'southeast')
legend([hT{1}, hT{2}], 'PNN-based NMPC', 'DNN-based NMPC', 'Location', 'southeast')
legend([hQ{1}, hQ{2}], 'PNN-based NMPC', 'DNN-based NMPC', 'Location', 'southeast')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MSE SCORES AND MEMORY FOOTPRINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../Supporting-Data-Files/MSE_Ns_5000.mat')
% Rows: Layers; Columns: Nodes

NlayersUnique = unique(Nlayers_list);
Nl = length(NlayersUnique);
count=1;
idx = cell(Nl,1);
mse_avg = zeros(Nl, length(Nnodes_list));
mse_sd = mse_avg;
memory_avg = mse_avg;
memory_sd = mse_avg;

for j = NlayersUnique
	idx{count}= find(Nlayers_list==j);
	mse_avg(count,:) = mean(mse_list(idx{count}, :), 1);
    mse_sd(count,:) = std(mse_list(idx{count}, :), 0, 1);
    
    memory_avg(count,:) = mean(Memory_dnn_kb(idx{count}, :), 1);
    memory_sd(count,:) = std(Memory_dnn_kb(idx{count}, :), 0, 1);
    
    count = count+1;
end

mse_min = mse_avg - 2*mse_sd;
mse_max = mse_avg + 2*mse_sd;

memory_min = memory_avg - 2*memory_sd;
memory_max = memory_avg + 2*memory_sd;

legendVec = [];
h = cell(Nl, 1);
color = ['r';'b';'m';'k'];

figure(7)
hold on
%{
for j = 1:Nl
    fx = [Nnodes_list, fliplr(Nnodes_list)];
    fy = [mse_min, fliplr(mse_max)];
    alpha(0.5)
    fill(fx, fy, [1,1,1]*0.9, 'LineStyle', ':')

end
%}
subplot(2,1,1)
hold on
for j = 1:Nl
    hold on
    h{j} = semilogy(Nnodes_list, mse_avg(j,:), '-o', 'color', color(j,:), 'Linewidth', 2);
end
legend([h{1}, h{2}, h{3}, h{4}], '1 Layer', '3 Layers', '5 Layers', '7 Layers')
xlabel('Nodes')
ylabel('MSE')
set(gca,'FontSize',fontSz)



figure(7)
subplot(2,1,2)
hold on
for j = 1:Nl
    h{j} = plot(Nnodes_list, memory_avg(j,:),  '-o', 'color', color(j,:), 'Linewidth', 2);
end
legend([h{1}, h{2}, h{3}, h{4}], '1 Layer', '3 Layers', '5 Layers', '7 Layers', 'Location', 'northwest')
xlabel('Nodes')
ylabel('Memory/kb')
set(gca,'FontSize',fontSz)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT AVERAGE COMPUTATION TIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Nhorizon = [5;10;20;40;60;80;100];
Tapprox = [0.013;0.014;0.013;0.014;0.014;0.014;0.013];
Tfull = [0.024;0.024;0.027;0.039;0.055;0.082;0.093];
TapproxProj = round([0.021356;0.021907;0.022512;0.021304;0.020411;0.022059;0.021807], 3);
figure(9)
hold on
plot(Nhorizon, Tfull*1000, '-o', 'color', 'r', 'Linewidth', 2)
plot(Nhorizon, Tapprox*1000, '-o', 'color', [0 0.4470 0.7410], 'Linewidth', 2)
plot(Nhorizon, TapproxProj*1000, '-o', 'color',[0.4940 0.1840 0.5560], 'Linewidth', 2)
legend('NMPC', 'DNN-based NMPC', 'PNN-based NMPC', 'Location', 'northwest')
xlabel('Prediction Horizon')
ylabel('Average computation time/ ms')
set(gca,'FontSize',fontSz)
%}

%%
Nplot = [5, 10, 20, 40 ,60, 80, 100];
tNMPC = [0.05, 0.057, 0.067, 0.087, 0.114, 0.145, 0.176];
tDNN = [52, 34, 42, 34, 33, 35, 49]*10^(-5);
tPNN = [24, 20, 24, 20, 20, 20, 21]*10^(-4);
colorvec = [0 0.4470 0.7410; 0.4940 0.1840 0.5560];

hf1 = figure(9);
hold on
plot(Nplot, tNMPC*1000,'r-o', 'Linewidth', 2)
plot(Nplot, tDNN*1000, '-o', 'color', colorvec(1,:), 'Linewidth', 2)
plot(Nplot, tPNN*1000, '-o', 'color', colorvec(2,:), 'Linewidth', 2)
xlabel('Prediction Horizon')
ylabel('Average computation time/ ms')
legend('Full NMPC', 'DNN-based NMPC', 'PNN-based NMPC', 'Location', 'northwest')
box on
set(gca,'FontSize',15)
axes('parent',hf1,'position',[0.35 0.2 0.5 0.2]);
hold on
plot(Nplot, tDNN*1000, '-o', 'color', colorvec(1,:), 'Linewidth', 2)
plot(Nplot, tPNN*1000, '-o', 'color', colorvec(2,:), 'Linewidth', 2)
box on
set(gca,'FontSize',fontSz), 




%% Change figures
figure(4)
title('(a)')
xlim([0, 27])
ylim([0, 1.6])
set(gca,'FontSize',fontSz)
box on

figure(5)
subplot(2,1,1)
title('(b)')
set(gca,'FontSize',fontSz)
xlim([0, 27])
ylim([32.5, 43])
box on
subplot(2,1,2)
xlim([0, 27])
box on

figure(6)
subplot(2,1,1)
title('(c)')
set(gca,'FontSize',fontSz)
xlim([0, 27])
box on
subplot(2,1,2)
xlim([0, 27])
subplot(2,1,2)
xlim([0, 27])
box on

figure(7)
xlim([2,12])
subplot(2,1,1)
title('(a)')
set(gca,'FontSize',fontSz)
box on
subplot(2,1,2)
title('(b)')
set(gca,'FontSize',fontSz)
box on
% figure(8)
% xlim([2,12])
% title('(b)')
% set(gca,'FontSize',fontSz)
box on
figure(9)
box on




