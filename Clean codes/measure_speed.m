clear all; clc; close all
[table, ~]=importfile(['Z:\Private Folders\Souheil Serroud\Test Data\With compliant bearing\8 springs\continous 1000.xlsx']);

set(groot,'defaultAxesTickLabelInterpreter','latex','DefaultAxesFontSize',12);  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
wn = 7;
t=table.data(1:end,5);
dt=1/2000;

% t0=min(t);
% 
% 
% dt = 1/2000;
% if t0==0
%     p0=1;
% else
%     p0=t0/dt;
% end
% 
% pf=p0+max(t)/dt;
% t=t(p0:pf)-(t0-dt);
% N=0;
% laser1=(table.data(p0+N:pf+N,2)*2.735978);
% prox = (table.data(p0:pf,4));

t0 = 2.5;
time_window=max(t)-t0;


if t0==0
    p0=1;
else
    p0=t0/dt;
end

pf=p0+time_window/dt;

t=t(p0:pf)-(t0-dt);

prox=table.data(p0:pf,4);



% Parameters
time_interval = 0.5; % Change this value to the interval you desire, in seconds

% Number of points to check before and after
ptc_before = 3; % change according to your needs
ptc_after = 1;
threshold = 0.03866;
% Initialize variables
counts = 0;
n = length(prox);
indices = [];

% Initialize speed storage
rot_speeds = [];
rot_speed_times = [];

% Initialize time interval tracking
last_time_interval = 0;

for j=ptc_before+1:n-ptc_after
    if prox(j) <= threshold && all(prox(j-ptc_before:j-1) == prox(j)) && any(prox(j+1:j+ptc_after) ~= prox(j))
        counts = counts + 1;
        indices = [indices j];
    end

    % Check if the current time exceeds the next time interval
    if t(j) - last_time_interval >= time_interval
        % Compute the rotational speed in RPM
        rpm = (counts / 4) * 60 / time_interval;
        % Store speed and corresponding time
        rot_speeds = [rot_speeds rpm];
        rot_speed_times = [rot_speed_times t(j)];
        % Reset counts for the next interval
        counts = 0;
        % Update the last time interval
        last_time_interval = t(j);
    end
end

% Plot the rotational speed over time
figure('Position', [100, 0, 1000, 200]);
plot(rot_speed_times, rot_speeds/60);
ylabel('   $\Omega$', 'FontSize', 16);
xlabel('$t$ (s)', 'FontSize', 16)
xlim([0 60]); 
text(-0.1, 0.5,'b)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
hold off
