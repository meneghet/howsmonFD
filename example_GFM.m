clc
clear
% close all force
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
N_patient = 1;
pat = TT{N_patient};
ts = minutes(pat.time(2) - pat.time(1)); % sampling time

%% GFM
LW = 24; % Long Window durations in hours
SW = 1; % Short Window durations in hours
[GFM, CGM_SW, CGM_LW] = gfm(pat.CGM, LW, SW, ts);

% in the article they divide by 60 (by deduction)
GFM = GFM / 60;

%% plot

fault_start = find(diff(pat.fault_pump) > 0);
fault_end = find(diff(pat.fault_pump) < 0);
N_fault = 2; % select the fault to visualize

f = figure('Color','w','Units','normalized','Position',[0.1 0.1 0.3 0.7]);
t = pat.time;

ax(1) = subplot(3,1,1);
L = [];
hold on
c = plot(t,pat.CGM,'k','DisplayName','CGM (mg/dL)','LineWidth',1.5);
L(1) = xline(t(fault_end(N_fault)),'r','DisplayName','set replacement','LineWidth',2);
% L(2) = xline(t(fault_start(N_fault)),'r','DisplayName','fault start','LineWidth',2);
ylim([0 400])
legend(L, 'Location', 'northwest')
ylabel('Glucose [mg/dL]')

ax(2) = subplot(3,1,2);
hold on
L2(1) = plot(t, CGM_SW,'r','DisplayName','CGM_{SW}','LineWidth',1.5);
L2(2) = plot(t, CGM_LW,'b','DisplayName','CGM_{LW}','LineWidth',1.5);
legend(L2, 'Location', 'northwest')
ylabel('Average Glucose [mg/dL]')

ax(3) = subplot(3,1,3);
plot(t, GFM,'k','LineWidth',1.5)
ylabel('GFM [mg/dL/min]')

% figure settings
linkaxes(ax,'x');
days_before = 2;
xlim([t(fault_start(N_fault))-days(days_before),t(fault_start(N_fault))+days(0.3)])

%% save
saveas(f, 'fig_example_GFM.png')

%%
rmpath(genpath('src'))