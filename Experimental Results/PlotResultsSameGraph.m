%
%   This script needs the plotResults.m script to run. Make sure they
%   are in the same folder!
%

% clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT AVERAGE TRAJECTORIES ON THE SAME GRAPH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cases = [1,2];
plotWithinScript=0;
color = {'r', 'b'};

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

    CEMavg(idx) = CEMavg(idx(1));
%     Tavg(idx) = 35;
%     Iavg(idx) = Iavg(idx(1));
%     qavg(idx) = qavg(idx(1));
%     Pavg(idx) = Pavg(idx(1));
    idx = idx(1);

    

    % Plot CEM
    figure(4)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [CEMmin(1:idx), fliplr(CEMmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
%     plot(tPlot(1:idx),CEMavg(1:idx), color{ii}, 'LineWidth', 2)
    plot(tPlot,CEMavg, color{ii}, 'LineWidth', 2)
    stairs(tPlot(1:idx),Sdes(1,1:idx),'k')
    ylim([0, 2*max(Sdes(1,:))])
    set(gca,'FontSize',12)

    
    
    % Plot states
    figure(5)
    subplot(2,1,1)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Tmin(1:idx), fliplr(Tmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),Tavg(1:idx), color{ii}, 'LineWidth', 2)
%     plot(tPlot,Tavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [x_max(1), x_max(1)]+Tss, 'k--')
    plot([tPlot(1), tPlot(idx)], [x_min(1), x_min(1)]+Tss, 'k--')
    hold on
    xlabel('Time/s')
    ylabel('Average T/ ^{\circ}C')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

    
    subplot(2,1,2)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Imin(1:idx), fliplr(Imax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),Iavg(1:idx), color{ii}, 'LineWidth', 2)
%     plot(tPlot, Iavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], 10*([x_max(2), x_max(2)]+Iss), 'k--')
    plot([tPlot(1), tPlot(idx)], 10*([x_min(2), x_min(2)]+Iss), 'k--')
    ylim(10*[x_min(2)+Iss-2, x_max(2)+Iss+2])
    xlabel('Time/s')
    ylabel('Average T/ ^{\circ}C')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

    
    
    % Plot inputs
    figure(6)
    subplot(2,1,1)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [qmin(1:idx), fliplr(qmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),qavg(1:idx), color{ii}, 'LineWidth', 2)
%     plot(tPlot, qavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [10, 10], 'k--')
    plot([tPlot(1), tPlot(idx)], [0.5, 0.5], 'k--')
    xlabel('Time/s')
    ylabel('Average q/ slm')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

    
    subplot(2,1,2)
    hold on
%     fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
%     fy = [Pmin(1:idx), fliplr(Pmax(1:idx))];
%     fill(fx, fy, [1,1,1]*0.9)
%     alpha(0.5)
    plot(tPlot(1:idx),Pavg(1:idx), color{ii}, 'LineWidth', 2)
%     plot(tPlot,Pavg, color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [5, 5], 'k--')
    plot([tPlot(1), tPlot(idx)], [1,1], 'k--')
    xlabel('Time/s')
    ylabel('Average Appled Power/ W')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MSE SCORES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../Supporting-Data-Files/MSE_Ns_2500.mat')
% Rows: Layers; Columns: Nodes

NlayersUnique = unique(Nlayers_list);
count=1;
idx = cell(3,1);
mse_avg = zeros(3, length(Nnodes_list));
mse_sd = mse_avg;
for j = NlayersUnique
	idx{count}= find(Nlayers_list==j);
	mse_avg(count,:) = mean(mse_list(idx{count}, :), 1);
    mse_sd(count,:) = std(mse_list(idx{count}, :), 0, 1);
    count = count+1;
end

mse_min = mse_avg - 3*mse_sd;
mse_max = mse_avg + 3*mse_sd;

legendVec = [];
h = cell(length(NlayersUnique), 1);
color = ['r';'b';'c'];

figure(7)
hold on
for j = 1:length(NlayersUnique)
    fx = [Nnodes_list, fliplr(Nnodes_list)];
    fy = [mse_min, fliplr(mse_max)];
    alpha(0.5)
    fill(fx, fy, [1,1,1]*0.9, 'LineStyle', ':')

end

for j = 1:length(NlayersUnique)
    h{j} = plot(Nnodes_list, mse_avg(j,:), color(j,:), 'LineWidth', 2);
end
legend([h{1}, h{2}, h{3}], '1 Layer', '2 Layers', '3 Layers')
xlabel('Nodes')
ylabel('MSE')
set(gca,'FontSize',12)



figure(8)
hold on
for j = 1:length(NlayersUnique)
    plot(Nnodes_list, mse_max(j,:), [color(j,:),':'], 'LineWidth', 1);
    plot(Nnodes_list, mse_min(j,:), [color(j,:),':'], 'LineWidth', 1);

end

for j = 1:length(NlayersUnique)
    h{j} = plot(Nnodes_list, mse_avg(j,:), color(j,:), 'LineWidth', 2);
end
legend([h{1}, h{2}, h{3}], '1 Layer', '2 Layers', '3 Layers')
xlabel('Nodes')
ylabel('MSE')
set(gca,'FontSize',12)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT AVERAGE COMPUTATION TIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nhorizon = [5;10;20;40;60;80;100];
Tapprox = [0.013;0.014;0.013;0.014;0.014;0.014;0.013];
Tfull = [0.024;0.024;0.027;0.039;0.055;0.082;0.093];


figure(8)
hold on
plot(Nhorizon, Tapprox, 'Linewidth', 2)
plot(Nhorizon, Tfull, 'Linewidth', 2)
legend('Approximate NMPC', 'Full NMPC', 'Location', 'northwest')
xlabel('Prediction Horizon')
ylabel('Average computation time/ s')
set(gca,'FontSize',15)

