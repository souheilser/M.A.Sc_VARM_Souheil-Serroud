clear all;clc;close all

% This script plots timeseries and FFTs of constant speed tests. The fata
% files should be named by the speed in RPM at which they were recorded so that a
% int. vector "n" can be used to call each excel file. See examples below.

% n=[180:60:600, 630:30:840, 900:60:1200]; %shaft only basic 60s Z:\Public Folders\Souheil Serroud\Final tests\Shaft only 3\excel
% n=[900:60:1140];
% n=900;
% n=[180:60:480 870 900 960]; %Shaft only repeat 3 Z:\Public Folders\Souheil Serroud\Final tests\Shaft only - repeat 3\excel
% n=[180:60:600, 630:30:840, 900:60:1140]; %misalignement both Z:\Public Folders\Souheil Serroud\Final tests\misal both
% n= [300 360 630 750:30:840 900:60:1140]; %misalignement 1.8cm 180s 

% n=[1 2:0.5:6.5]; %shaft and rotor 120 sec Z:\Public Folders\Souheil Serroud\Final Final\Shaft rotor\Aligned\Continuous speed\excel

n = 1200;   % Here a single data file is read.

%Setup properties
E = 110300000000; % Pa Youngs Modulus
b2b = 803; %mm Distance between the two bearings
d=6.35; %shaft diameter in mm
I = pi*(d)^4/64; %mm^4 Moment of inertia
% wn = (23.6+26.4)/2*2*pi; %rad/s
% fn = (23.6+26.4)/2; %hz
fn = 25; %Experimental natural frequency used for non-dimensional values
wn = fn*2*pi; %rad/s
M = 0.281; % Kg Total mass ;
disp_height = 60; %mm Height of the lasers from the lower bearing
dt=1/2000;  %Data acquisition frequency

rms_max=[];
fft_max1=[];
fft_max2=[];
fft_max3=[];

zeta = 0.01;    %Estimated damping coefficient
Ks = 192*E*(I*0.001^4)/(b2b*0.001)^3;   %Shaft equivalent stiffness
Ks2 = wn^2*M;   %System equivalent stiffness from the measured nat. freq.
b = zeta*2*sqrt(M*Ks2);

%y=2.735978112175103x+14.000000000000002
% set(gca,'TickLabelInterpreter','latex');

%Latex style axis functions
set(groot,'defaultAxesTickLabelInterpreter','latex','DefaultAxesFontSize',12);  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

bigpeakfft1 = [];
bigpeakfft11 = [];
bigpeakfft111 = [];
max_r = [];
max_x = [];
max_y = [];

for i=1:length(n)
    disp(i);disp(length(n))
    t2=[];
    [table, table2]=importfile(['Z:\Public Folders\Souheil Serroud\Final tests\Shaft only 3\excel\' num2str(n(i)) '.xlsx']);

%     n(i) = n(i) *60; Uncomment if the file names are in Hz
    t=table.data(1:end,4);
%     t2=table2.data(1:end,9);
    
    data_pt=max(size(t));
    
    time_window=max(t);
    t0=min(t);
    
    if t0==0
        p0=1;
    else
        p0=t0/dt;
    end

    pf=p0+time_window/dt;
    
    t=t(p0:pf)-(t0-dt);
    
    accel1=(table.data(p0:pf,1))*9.81/((d*0.001)*(wn^2));
    N=0;
    laser1=(table.data(p0+N:pf+N,2)*2.735978)/d; 
    laser2=(table.data(p0+N:pf+N,3)*2.735978)/d;
    
    prox = (table.data(p0+N:pf+N,4));


    %FFT__to_find_series_period_______________________________________________________
    Fs=1/dt;
    Fs2=10;
    T=1/Fs;
    
    L=length(t);
    L2=length(t2);
    
    Y=fft(laser1);
    Y2=fft(laser2);
    Y3=fft(accel1);
        
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    P22 = abs(Y2/L);
    P11 = P22(1:L/2+1);
    P11(2:end-1) = 2*P11(2:end-1);

    P222 = abs(Y3/L);
    P111 = P222(1:L/2+1);
    P111(2:end-1) = 2*P111(2:end-1);
     
    f = Fs*(0:(L/2))/L;
    f2=Fs2*(0:(L2/2))/L2;
    
% Get rid of low frequencies Peaks
%     indices_1=[0 1 find(ceil(f2)==1)];
%     
%     P1=[zeros(indices_1(end),1); P1(indices_1(end)+1:end)];
%     P11=[zeros(indices_1(end),1); P11(indices_1(end)+1:end)];
%     P111=[zeros(indices_1(end),1); P111(indices_1(end)+1:end)];
%     
    %_____Measure speed_________
    speeds = [];
    threshold = 0.03866;
    ptc_before = 3;
    
    counts = 0;
    last_time = t(1);
    
    for j=ptc_before+1:length(prox)-1
        if prox(j)<threshold && min(prox(j-ptc_before:j-1) == prox(j))==1 && prox(j+1)>threshold
            counts = counts+1;
            
            % Calculate speed every 5 seconds
            if t(j) - last_time >= 5
                speeds = [speeds, counts/4/(t(j)-last_time)];
                counts = 0;
                last_time = t(j);
            end
        end
    end
    
    % Calculate speed for the remaining time period if it's less than 5 seconds
    if t(end) - last_time > 0
        speeds = [speeds, counts/4/(t(end)-last_time)];
    end
    time_array = 0:5:(length(speeds)-1)*5;

    % Bearings reaction force estimation from disp data
    radial_disp = (sqrt((laser1*d*0.001).^2+(laser2*d*0.001).^2)); %m
    
    max_r = [max_r; max(radial_disp)];
    max_x = [max_x; max(laser1)];
    max_y = [max_y; max(laser2)];
    
    
    %_______PLotting figures_______________________________________________________________________
    fig(i)=figure('Position', [100, 0, 1200, 900]);
    
    %nextile_________________________________________________________

%     plot(t,accel1,'.k')
%     hold on
%     if n2(i)==450
%     
%         %title(['\fontsize{16} Time series at ' num2str(n(i)) ' RPM (' num2str(n(i)/60) ' Hz)'])
%     
%     end
    subplot(6,1,1);
    plot(t*wn,accel1,'-k')
%   ylim([-0.1 0.1]);
    xlim([min(t)*wn max(t)*wn]);
    xlabel('$\tau$', 'FontSize', 16);
    %title(['\fontsize{11} Time series at ' num2str(n(i)) ' RPM (' num2str(n(i)/60) ' Hz)'])
    ylabel('$\alpha$','FontSize',16);
    text(-0.1, 0.5,'a)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    %nextile_________________________________________________________
    subplot(6,1,2);
    plot(t*wn,laser1,'k')
    xlim([min(t)*wn max(t)*wn]);
%     ylim([-0.07 0.07]);
    ylim([min(laser1)*1.2 max(laser1)*1.2]);
    ylabel('$\delta_1$','FontSize',16)
    xlabel('$\tau$', 'FontSize', 16);
    text(-0.1, 0.5,'b)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    %nextile_________________________________________________________
    subplot(6,1,3);
    plot(t*wn,laser2,'k')
    xlim([min(t)*wn max(t)*wn]);
%     ylim([-0.07 0.07]);
    ylim([min(laser2)*1.2 max(laser2)*1.2]);
    ylabel('$\delta_2$','FontSize',16)
    xlabel('$\tau$', 'FontSize', 16);
    text(-0.1, 0.5,'c)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    %___Speed every 5 seconds____
%     subplot(6,1,4);
%     plot(time_array, speeds);
%     ylabel('\f_e','FontSize',16)
%     xlabel(sprintf('  \\tau'), 'FontSize', 16);
%     text(-0.1, 0.5,'d)','Interpreter','latex','Units', 'Normalized', 'VerticalAlignment', 'Top');
%     
    hold off
    
   % Orbit______________________________________
    subplot(6,1,[4 5 6]);
    x = laser1';
    y = laser2';
    z = t'*wn;
    patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp'); 
    c=colorbar;colormap parula;
    c.Label.String = '\fontsize{14} \tau';
    hold on
    plot3(0,0,9999,'ok')
    plot3(mean(x),mean(y),9999,'+r')
    
    xlabel('$\delta_1$','FontSize',16);
    ylabel('$\delta_2$','FontSize',16);
    xlim([min([min(x) min(y)])*1.2 max([max(x) max(y)])*1.2])
    axis('square');
    text(-0.35, 0.5,'d)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    hold off
   % FFT plots________________________________________________________________
   
    fig(i+length(n))=figure('Position', [100, 0, 1200, 900]);
    subplot(2,1,1);
    semilogy(abs(f/(n(i)/60)),P1,'r','LineWidth',1.5)
    hold on
    semilogy(abs(f/(n(i)/60)),P11,'b')
    %title(['\fontsize{11}Displacement FFT']) %' num2str(n(i)) ' RPM (' num2str(n(i)/60) ' Hz)'])
    %xlabel('frequency (Hz)')
    ylabel('$A$')
    xlim([0 20]);
%     xlim([0 max(abs(f/(n(i)/60)))]);
%     ylim([10^-6 0.4]);
    ylim([10^-5 max(findpeaks(P1))*1.2]);
    xlabel('   $\omega/\omega_e$', 'FontSize', 16);
    legend('Laser 1', 'Laser 2');
    text(-0.1, 0.5,'d)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    subplot(2,1,2);
    semilogy(abs(f/(n(i)/60)),P111,'k');
    %title(['\fontsize{10}Acceleration FFT']); %@ ' num2str(n(i)) ' RPM (' num2str(n(i)/60) ' Hz)'])
    xlabel('   $\omega/\omega_e$', 'FontSize', 16);
    ylabel('$\Delta$')
    xlim([0 20]);
%     xlim([0 max(abs(f/(n(i)/60)))]);    
%     ylim([0 1.5]);
%     ylim([0 max(findpeaks(P111(1:length(P111)/10)))*1.2]);
    text(-0.1, 0.5,'e)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    
    
    rms_max = [rms_max [rms(laser1); rms(laser2); rms(accel1)]];
    
    % Find FFT peak at excitation frequency
    % norm_f = f/fn;
    
    norm_f = f/(n(i)/60);

%     idx = find(norm_f >= 0.95 & norm_f <= 1);
%     % For laser 1
%     peak_window = P1(min(idx):max(idx));
%     [fft_peaks, ind] = findpeaks(peak_window);
%     [peak_w, loc] = max(fft_peaks); 
%     fft_max1 = [ fft_max1; peak_w norm_f(idx(ind(loc)))*(n(i)/60)];
%     % For laser 2
%     peak_window = P11(min(idx):max(idx));
%     [fft_peaks, ind] = findpeaks(peak_window);
%     [peak_w, loc] = max(fft_peaks); 
%     fft_max2 = [ fft_max2; peak_w norm_f(idx(ind(loc)))*(n(i)/60)];
%     % For aceel 1
%     peak_window = P111(min(idx):max(idx));
%     [fft_peaks, ind] = findpeaks(peak_window);
%     [peak_w, loc] = max(fft_peaks); 
%     fft_max3 = [ fft_max3; peak_w norm_f(idx(ind(loc)))*(n(i)/60)];
%     
%     %Find max peak in FFT
%     [bigpeak1, bigpeakpos1] = max(P1);
%     [bigpeak11, bigpeakpos11] = max(P11);
%     [bigpeak111, bigpeakpos111] = max(P111);
%     
%     bigpeakfft1 = [bigpeakfft1; bigpeak1, bigpeakpos1];
%     bigpeakfft11 = [bigpeakfft11; bigpeak11, bigpeakpos11];
%     bigpeakfft111 = [bigpeakfft111; bigpeak111, bigpeakpos111];
    
    
    %nextile_______ORBIT__________________________________________
    
%     subplot(6,1,[4 5 6]);
%     x = laser1';
%     y = laser2';
%     z = t'*n(i)/60;
%     
%     ppo = round(1/(norm_f(idx(ind(loc)))*n(i)/60*dt));
%     avg_orbitx = [];
%     avg_orbity = [];
%     for j=1:ppo:length(t)-ppo
%         
%         avg_orbitx = [avg_orbitx x(j:j+ppo)'];
%         avg_orbity = [avg_orbity y(j:j+ppo)'];
%         
%     end
%     orbitx = mean(avg_orbitx,2);
%     orbity = mean(avg_orbity,2);
%     
%     
%     patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp'); 
%     c=colorbar;colormap parula;
%     c.Label.String = '\tau';
%     hold on
%     plot(0,0,'ok')
%     plot(mean(x),mean(y),'.r')
%     plot3(orbitx,orbity,10^6*ones(length(orbitx),1),'k');
    %plot(mean_x,mean_y,'r')
%     ylim([-0.07 0.07]);
    %title(["\fontsize{10}Orbital path"])
%     xlim([min([min(x) min(y)]) max([max(x) max(y)])])
%     xlabel('\delta_1','FontSize',16)
%     ylabel('\delta_2','FontSize',16)
%     axis('equal')

end

hold off
% figure(i+length(n)+1)
% plot(abs(n)/60/fn, rms_max(1,1:end)/d,'.-r', abs(n)/60/fn, rms_max(2,1:end)/d,'o-b');
% %title('\fontsize{10}Maximum displacement RMS for multiple excitation frequencies');
% xlabel('\f_e/\f_n')
% ylabel('\delta RMS')
% legend('Laser 1', 'Laser 2')

% Fb = sqrt((rms_max(1,1:end)*d*0.001).^2+(rms_max(2,1:end)*d*0.001).^2)*48*E*I/(3*L-4*0.16);

fig(i+length(n)+1)=figure('Position', [100, 0, 1000, 500]);
% plot(abs(n)/60/fn, rms_max(3,1:end),'.-r')
% xlabel('\f_e/\f_n')
% ylabel('\alpha RMS')

ax1 = axes; % create the first set of axes
hold(ax1, 'on'); % hold the axes to allow multiple plots
xlabel('$\Omega$'); % label x-axis
xlim([0 0.85]);
ylabel(ax1, '$\alpha$ RMS'); % label left y-axis
p1 = plot(abs(n)/60/fn, rms_max(3,1:end),'+-k');

ax2 = axes; % create the second set of axes
hold(ax2, 'on'); % hold the axes to allow multiple plots
p2 = plot(abs(n)/60/fn, rms_max(1,1:end),'.-r', abs(n)/60/fn, rms_max(2,1:end),'o-b');
% Position the second y-axis on the right
set(ax2, 'YAxisLocation', 'right', 'Color', 'none');
ylabel(ax2, '$\delta$ RMS'); % label right y-axis
set(ax2, 'XTick', [], 'Box', 'off');
linkaxes([ax1, ax2], 'x');
ax1.YColor = 'k'; % Set the color of the y-axis and y-label of the first curve
ax2.YColor = 'k'; % Set the color of the y-axis and y-label of the second curve
legend([p1, p2(1,:), p2(2,:)], {'Acceleration', 'Laser 1','Laser 2'},'Location','northwest'); 

% figure(i+length(n)+3)
% hold on
% plot(fft_max1(:,2)/fn, fft_max1(:,1),'.-b')
% plot(fft_max2(:,2)/fn, fft_max2(:,1),'o-r')
% %title('\fontsize{10}Amplitude of peak displacement at excitation frequency from FFT data for multiple speeds');
% xlabel('\f_e/\f_n')
% ylabel('\Delta')
% legend('Laser 1', 'Laser 2')
hold off

% fig(i+length(n)+2)=figure('Position', [100, 0, 1000, 500]);
% 
% % plot(fft_max3(:,2)/fn, fft_max3(:,1),'.-b')
% % xlabel('\f_e/\f_n')
% % ylabel('A')
% 
% ax1 = axes; % create the first set of axes
% hold(ax1, 'on'); % hold the axes to allow multiple plots
% xlabel('$\Omega$'); % label x-axis
% ylabel(ax1, '$A$'); % label left y-axis
% p1 = plot(fft_max3(:,2)/fn, fft_max3(:,1),'+-k');
% % p1 = plot(bigpeakfft111(:,2),bigpeakfft111(:,1),'+-k');
% ax2 = axes; % create the second set of axes
% hold(ax2, 'on'); % hold the axes to allow multiple plots
% p2 = plot(fft_max1(:,2)/fn, fft_max1(:,1),'.-b',fft_max2(:,2)/fn, fft_max2(:,1),'o-r');
% % p2 = plot(bigpeakfft1(:,2),bigpeakfft1(:,1),'.-b',bigpeakfft11(:,2),bigpeakfft11(:,1),'o-r');
% % Position the second y-axis on the right
% set(ax2, 'YAxisLocation', 'right', 'Color', 'none');
% ylabel(ax2, '$\Delta$'); % label right y-axis
% set(ax2, 'XTick', [], 'Box', 'off');
% linkaxes([ax1, ax2], 'x');
% ax1.YColor = 'k'; % Set the color of the y-axis and y-label of the first curve
% ax2.YColor = 'k'; % Set the color of the y-axis and y-label of the second curve
% legend([p1, p2(1,:), p2(2,:)], {'Acceleration', 'Laser 1','Laser 2'},'Location','northwest'); 


speeds_adim = abs(n)/60/fn;
we = abs(n')/60*2*pi;
fe = we/(2*pi);
% F = max_r*48*E*I/(disp_height^2*(3*b2b-4*disp_height));
% me = F ./ ((n'/60*2*pi).^2);

max_mid_x = (48/192)*((max_x*d*b2b^3)/(disp_height^2*(3*b2b-4*disp_height))); %mm
max_mid_y = (48/192)*((max_y*d*b2b^3)/(disp_height^2*(3*b2b-4*disp_height))); %mm
max_mid_r = (48/192)*((max_r*b2b^3)/(disp_height^2*(3*b2b-4*disp_height))); %m


% Fx = max_mid_x*0.001*m.*(wn^2-((abs(n')/60*2*pi)).^2)*0.5;
% Fy = max_mid_y*0.001*m.*(wn^2-((abs(n')/60*2*pi)).^2)*0.5;

% Fx = sqrt(m^2*(we.^2-wn^2).^2+b^2*we.^2).*max_mid_x*0.001; %N
% Fy = sqrt(m^2*(we.^2-wn^2).^2+b^2*we.^2).*max_mid_y*0.001; %N

Fr = max_mid_r*192*E*(I*10^-12)/(b2b*0.001)^3;%N

% Fr = sqrt(m^2*(we.^2-wn^2).^2+b^2*we.^2).*max_mid_r; %N
% Fr3 = max_mid_r.*m.*sqrt((2*zeta*we./wn).^2+(1 - we.^2/wn^2).^2);

me = Fr./((n'/60*2*pi).^2.*max_mid_r);
m = (192*E*(I*10^-12))./((b2b*0.001)^3*we.^2);
m2 = max_mid_r*M.*(wn^2-we.^2)./(2*we.^2);

fig(i+length(n)+3)=figure('Position', [100, 0, 1000, 500]);
plot(speeds_adim,Fr/2)
hold on
plot(speeds_adim,127*ones(1,length(n)),'-r');
xlabel('$\Omega$'); % label x-axis
ylabel('$F_B$ (N)'); % label left y-axis
legend('Bearing load', 'Minimal load');

fig(i+length(n)+4)=figure('Position', [100, 0, 1000, 500]);
plot(abs(n)/60/fn,me)
hold on
plot(abs(n)/60/fn,mean(me)*ones(1,length(n)),'-r');
xlabel('$\Omega$'); % label x-axis
ylabel('$me$ (kg m)'); % label left y-axis
legend('Unbalance', 'Average');
hold off

% figure(i+length(n)+10)
% plot(Fr)
% 
% figure(i+length(n)+20)
% plot(me)
plot(fe,max_r)
