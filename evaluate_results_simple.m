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

%% settings
LW = 24; % h
SW = 1; % h
K = 0.02;

GFM_th = 100;
IFM_th = 0.4;
slope_th = 0.3;

%% make alarms
for N = 1:length(TT)
    subj = TT{N};
    CGM  = subj.CGM;
    insulin = subj.basal/60 + subj.bolus/60;
    Ts = minutes(subj.time(2)-subj.time(1));
    % extract features
    [GFM, IFM, CGMslope] = howsmon_features(CGM, insulin, K, LW, SW, Ts);
    % apply threshold
    fd = (GFM > GFM_th) + (IFM > IFM_th) + (CGMslope > slope_th);
    fd(fd < 3) = 0;
    alarms{N,1} = alarm_engineering(fd,Ts);
end

%% make experiment and evaluate
TPmax = 6; % hours
time = TT{1}.time;

% T1D data and label
for N = 1:length(TT)
    subj = TT{N};
    
    t1d_data{N,1}.CGM = subj.CGM;
    t1d_data{N,1}.meal = subj.meal;
     
    L = rising_edge_trigger(subj.fault_pump);
    label{N,1} = L * TPmax;
end
experiment = FDexperiment('pump',label,alarms,time,t1d_data,GFM_th,IFM_th,slope_th);

%% evaluate
R = experiment.evaluate().report();

disp(R)


%%
rmpath(genpath('src'))





