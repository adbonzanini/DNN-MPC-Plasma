%
%   This script needs the plotResults.m script to run. Make sure they
%   are in the same folder!
%

clear all

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
    idx = idx(1);


    % Plot CEM
    figure(4)
    hold on
    fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
    fy = [CEMmin(1:idx), fliplr(CEMmax(1:idx))];
    fill(fx, fy, [1,1,1]*0.9)
    alpha(0.5)
    plot(tPlot(1:idx),CEMavg(1:idx), color{ii}, 'LineWidth', 2)
    stairs(tPlot(1:idx),Sdes(1,1:idx),'k')
    ylim([0, 2*max(Sdes(1,:))])
    set(gca,'FontSize',12)

    
    
    % Plot states
    figure(5)
    subplot(2,1,1)
    hold on
    fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
    fy = [Tmin(1:idx), fliplr(Tmax(1:idx))];
    fill(fx, fy, [1,1,1]*0.9)
    alpha(0.5)
    plot(tPlot(1:idx),Tavg(1:idx), color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [x_max(1), x_max(1)]+Tss, 'k--')
    plot([tPlot(1), tPlot(idx)], [x_min(1), x_min(1)]+Tss, 'k--')
    hold on
    xlabel('Time/s')
    ylabel('Average T/ ^{\circ}C')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

    
    subplot(2,1,2)
    hold on
    fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
    fy = [Imin(1:idx), fliplr(Imax(1:idx))];
    fill(fx, fy, [1,1,1]*0.9)
    alpha(0.5)
    plot(tPlot(1:idx),Iavg(1:idx), color{ii}, 'LineWidth', 2)
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
    fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
    fy = [qmin(1:idx), fliplr(qmax(1:idx))];
    fill(fx, fy, [1,1,1]*0.9)
    alpha(0.5)
    plot(tPlot(1:idx),qavg(1:idx), color{ii}, 'LineWidth', 2)
    plot([tPlot(1), tPlot(idx)], [10, 10], 'k--')
    plot([tPlot(1), tPlot(idx)], [0.5, 0.5], 'k--')
    xlabel('Time/s')
    ylabel('Average q/ slm')
    set(gcf,'color','w');
    set(gca,'FontSize',12)

    
    subplot(2,1,2)
    hold on
    fx = [tPlot(1:idx), fliplr(tPlot(1:idx))];
    fy = [Pmin(1:idx), fliplr(Pmax(1:idx))];
    fill(fx, fy, [1,1,1]*0.9)
    alpha(0.5)
    plot(tPlot(1:idx),Pavg(1:idx), color{ii}, 'LineWidth', 2)
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

mse_list = round(mse_list,3);
figure(7)
hold on
contourf(Nnodes_list, Nlayers_list, mse_list,4, 'ShowText','on')
xlabel('Number of Nodes')
ylabel('Number of Layers')
set(gca,'FontSize',12)
xlim([5,8])
% ylim([3,5])
