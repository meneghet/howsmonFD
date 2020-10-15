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
ts = minutes(pat.time(2) - pat.time(1));

%% IFM
LW = 24; % Long Window durations in hours
SW = 1; % Short Window durations in hours

K = 0.02;
insulin = pat.basal/60 + pat.bolus/60;
[PIE, IOB] = twoCompSysD(insulin, K, ts);
[IFM, PIE_LW, PIE_SW] = ifm(PIE, LW, SW, ts);

%% Plot
fault_start = find(diff(pat.fault_pump) > 0);
fault_end = find(diff(pat.fault_pump) < 0);
N_fault = 2; % select fault to visualize

f = figure('Color','w','Units','normalized','Position',[0.1 0.1 0.3 0.7]);
t = pat.time;

ax(1) = subplot(3,1,1);
hold on
c = plot(t,pat.CGM,'k','DisplayName','CGM (mg/dL)','LineWidth',1.5);
L(1) = xline(t(fault_end(N_fault)),'r','DisplayName','set replacement','LineWidth',1.5);
% L(2) = xline(t(fault_start(N_fault)),'r','DisplayName','fault start','LineWidth',2);
ylim([0 400])
legend(L, 'Location', 'northwest')
ylabel('Glucose [mg/dL]')

ax(2) = subplot(3,1,2);
hold on
L1(1) = plot(t,pat.basal,'k','DisplayName','basal (U/h)','LineWidth',1.5);
bolus = pat.bolus;
bolus(bolus <= 0) = nan;
L1(2) = stem(t,bolus/60*ts,'b','DisplayName','bolus (U)','LineWidth',1.5,'MarkerFace','b','MarkerSize',3);
L1(3) = plot(t,PIE,'r-.','DisplayName','PIE','LineWidth',1.5);
legend(L1, 'Location', 'northwest')
ylabel('Insulin')

ax(3) = subplot(3,1,3);
plot(t,IFM,'k','LineWidth',1.5)
ylabel('IFM')

% figure settings
linkaxes(ax,'x');
days_before = 1;
xlim([t(fault_start(N_fault))-days(days_before),t(fault_start(N_fault))+days(0.3)])

%% save
saveas(f, 'fig_example_IFM.png')

%%
rmpath(genpath('src'))