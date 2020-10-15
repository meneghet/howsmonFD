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
% separate subjects
subj_id_list = unique(loaded_data.ID);
for n = subj_id_list(:)'
    T = loaded_data(loaded_data.ID == n,:);
    TT{n,1} = table2timetable(T);
end
clear('loaded_data');

%% select patient
N = 1;
pat = TT{N};

%% Get CGM slope
SW = 1; % Short Window duratino in hours
ts = minutes(pat.time(2)-pat.time(1));
CGMslope = moving_slope(pat.CGM, SW, ts);

%% Plot
fault_start = find(diff(pat.fault_pump) > 0);
fault_end = find(diff(pat.fault_pump) < 0);
N_fault = 2; % select fault to visualize

figure('Color','w')
t = pat.time;

ax(1) = subplot(2,1,1);
hold on
c = plot(pat.time, pat.CGM, 'k','DisplayName','CGM (mg/dL)');
L(1) = xline(t(fault_end(N_fault)),'r','DisplayName','set replacement','LineWidth',1.5);
ylabel('Glucose [mg/dL]')

ax(2) = subplot(2,1,2);
hold on
L2(1) = plot(t, CGMslope,'r','DisplayName','CGM_{SW}');
yline(0,'k');
legend(L2, 'Location', 'southeast')
ylabel('Glucose slope [mg/dL/min]')
linkaxes(ax,'x');

days_before = 1;
xlim([t(fault_start(N_fault))-days(days_before), t(fault_start(N_fault))+days(0.3)])


%%
rmpath(genpath('src'))