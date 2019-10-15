%Test pcp ip

t = tcpip(host, port, 'NetworkRole', 'client');
UU = [1.5;1.5].*ones(2, 30);
UU(:, 1:10) = [1.5;3.0].*ones(2, 10);
UU(:, 10:15) = [2;1.5].*ones(2,6);

y_meas=cell(10, 1);
fopen(t);
pause(0.1);
disp('Connection established...')
for i=1:20
disp('Start Loop')
% Send initial point as a comma-delimited string
u_opt=UU(:,i);
Usend = sprintf('%6.1f, %6.1f ', [u_opt(1), u_opt(2)]);
% Usend = sprintf('%.1f,', u_opt(:));
% Usend = Usend(1:end-1);
disp('Sending initial point...') 
fwrite(t, Usend)
disp(['[', Usend, ']'])

pause(3)
% Receive measurement
disp('Receive Measurement...')
measurements = fread(t, [1, t.BytesAvailable]);
try
    y_meas{i} = char(measurements(end-12:end));
catch
    ymeas{i} = [];
end
y_meas{i}
disp('Done!\n')
end


fclose(t);