clear all;clc;close all

n=[180:60:600, 630:30:840, 900:60:1200]; %shaft only basic 60s
% n=[900:60:1140];
% n=900;
% n=[180:60:480 870 900 960]; %Shaft only repeat 3
n=[180:60:600, 630:30:840, 900:60:1140]; %misalignement repeat 1.8cm

% n2=[1000 2000 3000 4000 5000 6000 7000 8000 9000 9500 10000 10500 11000 12000];
dt=0.0005;
n = 960;
d=6.35;
rms_max=[];
fft_max1=[];
fft_max2=[];
fft_max3=[];

d=6.35; %diameter in mm
wn = (23.6+26.4)/2*2*pi;
fn = (23.6+26.4)/2;
%y=2.735978112175103x+14.000000000000002

for i=1:length(n)
    t2=[];
    [table, table2]=importfile(['Z:\Public Folders\Souheil Serroud\Final tests\Shaft only - repeat 3\excel\' num2str(n(i)) '.xlsx']);
    %shaft impact
    %     [table, table2]=importfile(strcat( 'Z:\Private Folders\Souheil Serroud\Nat frequencies_structure\SHAFT\','Impact_accel_on_rotor2', '.xlsx' ))%importfile(['Z:\Private Folders\Souheil Serroud\Test Data\Latest_05052023_newshaft1220mm_rotor\' num2str(n(i)) '.xlsx']);
    %ramping speed
    %     table=importfile(['\\triton.meca.polymtl.ca\usagers\sosera\Documents\TEST_2xIG_accel - Copie\0_450.xlsx']);
%     we=n(i)/60*2*pi;
    t=table.data(1:end,4);
    t2=table2.data(1:end,7);
    
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
    
    N=0;
        
%     accel1=(table.data(p0:pf,1))*9.81/((d*0.001)*(wn^2));
%     laser1=(table.data(p0+N:pf+N,2)*2.735978)/d; 
%     laser2=(table.data(p0:pf,3)*2.735978)/d;
%     
    accel1=(table.data(p0:pf,1))*9.81;
    laser1=(table.data(p0+N:pf+N,2)*2.735978); 
    laser2=(table.data(p0:pf,3)*2.735978);
    
    laser1 = laser1 - mean(laser1);
    
    accel1_2=table2.data(1:end,1);
    step_position = table2.data(1:end,end);
    laser1_2=table2.data(1:end,2)*2.735978; 
    laser2_2=table2.data(1:end,3)*2.735978;
    
    % Time derivative for Poincare plot
    supportlength = 3;
    modelorder = 2;
    dydt = movingslope(laser1,40,modelorder,dt);

    figure(i)
    plot(t,laser1)
    
    figure(i+length(n))
    %plot(t,dydt)
    hold on
    plot(t,dydt)
    legend('dydt2','dydt3')
        %FFT__to_find_series_period_______________________________________________________
    Fs=1/dt;
    Fs2=10;
    T=1/Fs;
    
    L=length(t);
    L2=length(t2);
    
    Y=fft(laser1);
    Y2=fft(laser2);
    Y3=fft(accel1);
    Y4=fft(accel1_2);
    
    
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    P22 = abs(Y2/L);
    P11 = P22(1:L/2+1);
    P11(2:end-1) = 2*P11(2:end-1);

    P222 = abs(Y3/L);
    P111 = P222(1:L/2+1);
    P111(2:end-1) = 2*P111(2:end-1);
    
    P2222 = abs(Y4/L2);
    P1111 = P2222(1:L2/2+1);
    P1111(2:end-1) = 2*P1111(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    f2=Fs2*(0:(L2/2))/L2;
% Get rid of low frequencies Peaks
    indices_1=[0 1 find(ceil(f2)==1)];
    
    P1=[zeros(indices_1(end),1); P1(indices_1(end)+1:end)];
    P11=[zeros(indices_1(end),1); P11(indices_1(end)+1:end)];
    P111=[zeros(indices_1(end),1); P111(indices_1(end)+1:end)];
    P1111=[zeros(indices_1(end),1); P1111(indices_1(end)+1:end)];

    
    % Find FFT peak at excitation frequency
    norm_f = f/(n(i)/60);
    
    idx = find(norm_f >= 0.95 & norm_f <= 1.05);
    % For laser 1
    peak_window = P1(min(idx):max(idx));
    [fft_peaks, ind] = findpeaks(peak_window);
    [peak_w, loc] = max(fft_peaks); 
    fft_max1 = [ fft_max1; peak_w norm_f(idx(ind(loc)))*n(i)/60 ];
    % For laser 2
    peak_window = P11(min(idx):max(idx));
    [fft_peaks, ind] = findpeaks(peak_window);
    [peak_w, loc] = max(fft_peaks); 
    fft_max2 = [ fft_max2; peak_w norm_f(idx(ind(loc)))*n(i)/60 ];
    % For aceel 1
    peak_window = P111(min(idx):max(idx));
    [fft_peaks, ind] = findpeaks(peak_window);
    [peak_w, loc] = max(fft_peaks); 
    fft_max3 = [ fft_max3; peak_w norm_f(idx(ind(loc)))*n(i)/60 ];
    
    ppo = round(1/(norm_f(idx(ind(loc)))*n(i)/60*dt));
    pcmap = [];
    for j = 1:ppo:length(t)-ppo
        
        [pc1, k] = max(abs(laser1(j:j+ppo)));
        pc2 = dydt(k);
        pcmap = [pcmap; pc1 pc2];
        
        
    end
    figure(i+2*length(t))
    plot(pcmap(:,1),pcmap(:,2))
end