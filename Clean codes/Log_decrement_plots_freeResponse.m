clear all;clc;close all


dt=0.0005;
rms_max=[];
fft_max1=[];
fft_max2=[];
fft_max3=[];

%y=2.735978112175103x+14.000000000000002

for i=1:1
    t2=[];
    [table, table2]=importfile(strcat( 'Z:\Private Folders\Souheil Serroud\Test Data\With compliant bearing\8 springs\free response_flick', '.xlsx' ));
    t=table.data(1:end,5);
%     t2 = table2.data(1:end,7);
    
    data_pt=max(size(t));
    dt=1/2000;
    
%     CHOSE TIME WINDOW TO ANALYSE
%__________________________________________________
%     shaft rotor impact.xlsx
    time_window= 5;
    t0 = 3.3;
%__________________________________________________   
%     shaft rotor free response.xlsl
%     time_window= 2.5;
%     t0 = 27.5;
%__________________________________________________   
%     shaft only impactlast.xlsl
%     time_window= 2;
%     t0 = 24.7+0.05+0.288;
% %__________________________________________________   
%     time_window= max(t);
%     t0 =min(t);
%     
    if t0==0
        p0=1;
    else
        p0=t0/dt;
    end

    pf=p0+time_window/dt;
    
    t=t(p0:pf)-(t0-dt);
    
    accel1=table.data(p0:pf,1);
    N=0;
    laser1=table.data(p0+N:pf+N,2)*2.735978; 
    laser2=table.data(p0:pf,3)*2.735978;
    
    accel1_2=table2.data(1:end,1);
    step_position = table2.data(1:end,end);
    laser1_2=table2.data(1:end,2)*2.735978; 
    laser2_2=table2.data(1:end,3)*2.735978;
    
    accel1 = accel1 - mean(accel1);
    radial_disp = sqrt(laser1.^2+laser2.^2);
    %FFT__to_find_series_period_______________________________________________________
    Fs=1/dt;
    Fs2=10;
    T=1/Fs;
    
    L=length(t);
    L2=length(t2);
    
    Y=fft(laser1);
    Y2=fft(laser2);
    Y3=fft(accel1);
    Y4=fft(radial_disp);
    
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    P22 = abs(Y2/L);
    P11 = P22(1:L/2+1);
    P11(2:end-1) = 2*P11(2:end-1);

    P222 = abs(Y3/L);
    P111 = P222(1:L/2+1);
    P111(2:end-1) = 2*P111(2:end-1);
    
    P2222 = abs(Y4/L);
    P1111 = P2222(1:L/2+1);
    P1111(2:end-1) = 2*P1111(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    f2=Fs2*(0:(L2/2))/L2;
    % Get rid of low frequencies Peaks
    indices_1=[0 1 find(ceil(f2)==1)];
    
    P1=[zeros(indices_1(end),1); P1(indices_1(end)+1:end)];
    P11=[zeros(indices_1(end),1); P11(indices_1(end)+1:end)];
    P111=[zeros(indices_1(end),1); P111(indices_1(end)+1:end)];
    P1111=[zeros(indices_1(end),1); P1111(indices_1(end)+1:end)];
%     % Estimate damping ratio
%     % Accelerometer
%     wn_accel1 = 5.87*2*pi;
%     low_freq=5.5;
%     high_freq=6.5;
%     filtered_accel1 = bandpass_filter(accel1, low_freq, high_freq,Fs, 1);
%     
%     damping_ratio = estimate_damping_ratio(t, filtered_accel1)
%     
%     curve = exp(-damping_ratio*wn_accel1*t);
%     % Laser 1
%     wn_laser1 = 5.87*2*pi;
%     low_freq=5.5;
%     high_freq=6.5;
%     filtered_laser1 = bandpass_filter(laser1, low_freq, high_freq,Fs, 1);
%     
%     damping_ratio_l1 = estimate_damping_ratio(t, filtered_laser1)
%     
%     curve_l1 = exp(-damping_ratio_l1*wn_laser1*t);
%     % Laser 2
%     wn_laser2 = 5.87*2*pi;
%     low_freq=5.5;
%     high_freq=6.5;
%     filtered_laser2 = bandpass_filter(laser2, low_freq, high_freq,Fs, 1);
%     
%     damping_ratio_l2 = estimate_damping_ratio(t, filtered_laser2)
%     
%     curve_l2 = exp(-damping_ratio_l2*wn_laser2*t);
    %______________________________________________________________________________
%     fig(i)=figure('Position', [10, 1, 1800, 950]);
% 
%     subplot(6,2,[1 3]);
%     hold on
%     plot(t,accel1,'-k')
% %     plot(t,filtered_accel1,'r')
% %     plot(t,curve,'b')
%     xlim([min(t) max(t)]);
%     title(['Impact test response'])
%     ylabel('Acceleration (g)','FontSize',11)
    
    %nextile_________________________________________________________
%     subplot(6,2,3);
%     hold on
%     plot(t,laser1,'k')
% %     plot(t,filtered_laser1,'r')
% %     plot(t,curve_l1,'b')
%     xlim([min(t) max(t)]);
%     ylim([min(laser1)*1.2 max(laser1)*1.2]);
%     ylabel('Displacement (mm)','FontSize',11)

    %nextile_________________________________________________________
    
    % Compute the spline interpolant
    [yupper,ylower] = envelope(laser1);
    windowSize = 100000; % Size of the moving average window
    yupper_smooth = smooth(yupper, windowSize);
    ylower_smooth = smooth(ylower, windowSize);
    
    figure(1)
    figure('Position', [100, 0, 900, 300]);
    ax1 = axes; % create the first set of axes
    hold(ax1, 'on'); % hold the axes to allow multiple plots
    xlabel('Time (s)'); % label x-axis
    ylabel(ax1, 'Displacement (mm)','Interpreter','latex'); % label left y-axis
    text(-0.15, 0.5,'b)','Interpreter','latex','FontSize',16,'Units', 'Normalized', 'VerticalAlignment', 'Top');
    ylim([-0.3 0.3]);
    p1 = plot(t,laser1,'b','LineWidth',0.5);
%     plot(t,yupper,t,ylower)
    
%     line([0.37 0.37], ylim, 'Color', 'r');
%     line([1.65 1.65], ylim, 'Color', 'r');
%     line([2.37 2.37], ylim, 'Color', 'r');
%     
%     ax2 = axes; % create the second set of axes
%     hold(ax2, 'on'); % hold the axes to allow multiple plots
%     p2 = plot(t,laser2,'b');
%     % Position the second y-axis on the right
%     set(ax2, 'YAxisLocation', 'right', 'Color', 'none');
%     ylabel(ax2, 'Displacement (mm)'); % label right y-axis
%     set(ax2, 'XTick', [], 'Box', 'off');
%     linkaxes([ax1, ax2], 'x');
%     ax1.YColor = 'k'; % Set the color of the y-axis and y-label of the first curve
%     ax2.YColor = 'k'; % Set the color of the y-axis and y-label of the second curve
    
%_______BEATS PERIODS LINES___________________________
%     line([0.45 0.45],[0.3 -0.3])
%     line([0.45+1/(7.19-6.29) 0.45+1/(7.19-6.29)],[0.3 -0.3])
%     line([0.45+2/(7.19-6.29) 0.45+2/(7.19-6.29)],[0.3 -0.3])
%     line([0.45+3/(7.19-6.29) 0.45+3/(7.19-6.29)],[0.24 -0.3])
    
    legend([p1, p2], {'Acceleration', 'Displacement'}); 
    xlim([0 time_window])
    
    
    %nextile_________________________________________________________
%     subplot(6,2,[7 9 11]);
%     x = laser1';
%     y = laser2';
%     z = t';
%     patch([x nan],[y nan],[z nan],[z nan], 'edgecolor', 'interp'); 
%     c=colorbar;colormap parula;
%     c.Label.String = 'Time (s)';
%     hold on
%     plot(0,0,'.k')
%     plot(mean(x),mean(y),'.r')
%     %plot(mean_x,mean_y,'r')
%     title(["\fontsize{10}Orbital path"])
%     xlabel('x (mm)','FontSize',14)
%     ylabel('y (mm)','FontSize',14)
%     axis('equal')
%     
   % FFT plots________________________________________________________________
    hold off
    figure(2)
    figure('Position', [100, 0, 900, 600]);
    ax1 = axes; % create the first set of axes
    hold(ax1, 'on'); % hold the axes to allow multiple plots
    xlabel('Frequency'); % label x-axis
    ylabel(ax1, 'Acceleration amplitude (-)'); % label left y-axis
    p1 = plot(abs(f),P11,'b','LineWidth',1);
    
    ax2 = axes; % create the second set of axes
    hold(ax2, 'on'); % hold the axes to allow multiple plots
    p2 = plot(abs(f),P11,'b');
    % Position the second y-axis on the right
    set(ax2, 'YAxisLocation', 'right', 'Color', 'none');
    ylabel(ax2, 'Displacement amplitude (-)'); % label right y-axis
    set(ax2, 'XTick', [], 'Box', 'off');
    linkaxes([ax1, ax2], 'x');
    ax1.YColor = 'k'; % Set the color of the y-axis and y-label of the first curve
    ax2.YColor = 'k'; % Set the color of the y-axis and y-label of the second curve
    xlim(ax1, [0 100]); % set x-axis limits for the first axis
    xlim(ax2, [0 100]); % set x-axis limits for the second axis
    legend([p1, p2], {'Acceleration', 'Displacement'});    
    
    
%%%%%%% Log decrement technique see https://josephcslater.github.io/MechanicalVibrationLab/labnotes/Matlab1.html
    
%     figure(i+2)
%     plot(t,accel1)
%     title(['Impact test response acceleration'])
%     ylabel('Acceleration \alpha','FontSize',11)
%     xlabel('Time (s)')
    
    
    
%     accel1 = accel1*9.81;
%     acc_avg = mean(accel1);
%     acce1l = accel1 - acc_avg;
%     
%     laser2 = laser2*9.81;
%     lase2_avg = mean(laser2);
%     laser2 = laser2 - lase2_avg;
    
    m=0.281%+1.45;
%     x1 = 0.79561;
%     t1 = 0.1115;
%     x2 = 0.18413;
%     t2 = 1.092;
%     n = 23;
    x1 = -0.23195;
    t1 = 0.6405;
%     x2 = 0.055334;
%     t2 = 0.981;

    x2 = -0.089353;
    t2 = 2.459; 
    
    n = 34;
    time = t2-t1;
    
    d = (1/n)*log(x1/x2) % delta
    z= d/sqrt(4*pi^2+d^2) % zeta
    Td = (time/n)%damped time period
    wd = 1/Td % damped natural frequency
    wne = wd/sqrt(1-z^2) % natural frequency in rad/sec
    k2=(wne^2)*m % one way to use experimental zeta and nat freq to calculate C
    c_cr = 2*sqrt(k2*m); % critical
    c=c_cr*z % damping constantwd
    figure(3)
    plot(abs(f),P11,'r')
    xlim([0 50]);
    
    plot(t,laser2)
    period = 1.3333;
%     line([0.397 0.397], ylim, 'Color', 'r');
%     line([0.397+period 0.397+period], ylim, 'Color', 'r');
%     line([0.397+period*2 0.397+period*2], ylim, 'Color', 'r');

end
