%   Run this script after running DNN_Control_Experiments.m so
%   that the plasma is restored to its previous inputs
u_opt=[1.5;1.5;];
Usend = sprintf('%6.1f, %6.1f ', [u_opt(1), u_opt(2)]);
disp('Sending initial point...') 
fwrite(t, Usend)
disp(['Sent inputs (Q, P) = ', '[', Usend, ']'])
fclose(t)
close all