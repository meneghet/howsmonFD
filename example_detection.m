clc
clear
close all force
restoredefaultpath;
addpath(genpath('src'))

%% Load dataset

% datasets_folder = 'simulated_data';
% fname = 'sim2018_PID_pump_faults.csv';

datasets_folder = 'clinical_data';
fname = 'h_data.csv';

% load data from csv
loaded_data = readtable(fullfile(datasets_folder, fname));
for n = 1:100
    T = loaded_data(loaded_data.ID == n,:);
    TT{n,1} = table2timetable(T);
end
clear('loaded_data');

%% select patient
N = 1;
pat = TT{N};

%% extract GFM, IFM, glucose slope
LW = 24; % h
SW = 1; % h
Ts = minutes(pat.time(2)-pat.time(1));

CGM  = pat.CGM;
insulin = pat.basal/60 + pat.bolus/60;

K = 0.02;
[GFM, IFM, CGMslope] = howsmon_features(CGM, insulin, K, LW, SW, Ts);

%% apply threshold values
GFM_th = 100;
IFM_th = 0.4;
slope_th = 0.3;

fd = (GFM > GFM_th) + (IFM > IFM_th) + (CGMslope > slope_th);
fd(fd < 3) = 0;
fd = boolean(fd);

%% plot results
fault_start = find(diff(pat.fault_pump) > 0);
fault_end = find(diff(pat.fault_pump) < 0);
N_fault = 2; % select fault to visualize

figure('Color','w')
ax(1) = subplot(5,1,1);
hold on
c = plot(pat.time,pat.CGM,'k','DisplayName','CGM (mg/dL)');
f1 = xline(pat.time(fault_start(N_fault)),'r','DisplayName','fault start','LineWidth',2);
% f2 = xline(pat.time(fault_end(N_fault)),'g','DisplayName','fault end','LineWidth',2);
ylabel('Glucose [mg/dL]')
% ======== GFM ========
ax(2) = subplot(5,1,2);
plot(pat.time,GFM,'k')
ylabel('GFM')
hold on
yline(GFM_th,'r-');
% ======== IFM ========
ax(3) = subplot(5,1,3);
plot(pat.time,IFM,'k')
ylabel('IFM')
yline(IFM_th,'r-');
% ======== slope of CGM in short window ========
ax(4) = subplot(5,1,4);
plot(pat.time,CGMslope,'k')
ylabel('CGM slope')
yline(slope_th,'r-');
% ==== FD ====
subplot(5,1,5)
ax(5) = subplot(5,1,5);
plot(pat.time,fd,'b');
ylabel('Alarm')

linkaxes(ax,'x');
xlim([pat.time(fault_start(N_fault))-days(2),pat.time(fault_start(N_fault))+days(0.3)])

%% save
saveas(gcf, 'fig_example_detection.png')

%%
rmpath(genpath('src'))


