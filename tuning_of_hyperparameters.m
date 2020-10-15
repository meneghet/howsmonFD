clc
clear
close all force
restoredefaultpath;
script_start = tic;
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
subj_id_list = subj_id_list(:);

TT = {};
for n = subj_id_list(:)'
    T = loaded_data(loaded_data.ID == n,:);
    TT{n,1} = table2timetable(T);
end
clear('loaded_data');

%% settings
LW = 24; % h
SW = 1; % h
K = 0.02;
TPmax = 6; % hours
time = TT{1}.time;

%% extract features
for N = 1:length(TT)
    subj = TT{N};
    
    CGM  = subj.CGM;
    insulin = subj.basal/60 + subj.bolus/60;
    ts = minutes(subj.time(2)-subj.time(1));
    
    [GFM{N}, IFM{N}, CGMslope{N}] = howsmon_features(CGM, insulin, K, LW, SW, ts);
end

%% grid search to find the optimal values for the thresholds

T = table(); % this table will store the results obtained with every threshold
k = 1; % a counter

% grid search over these values
GFM_th_selection = 50:25:400;
IFM_th_selection = 0:0.05:1;
slope_th_selection = 0:0.1:2;


for GFM_th = GFM_th_selection
    for IFM_th = IFM_th_selection
        for slope_th = slope_th_selection
            s = fprintf('GFM: %g, IFM: %g, der: %g \n',GFM_th,IFM_th,slope_th);
            
            % for every subject
            for N = 1:length(TT)
                subj = TT{N};
                ts = minutes(subj.time(2) - subj.time(1));
                % save subject data
                t1d_data{N,1}.CGM = subj.CGM;
                t1d_data{N,1}.meal = subj.meal;
                L = rising_edge_trigger(subj.fault_pump);
                label{N,1} = L * TPmax;
                % apply threshold to generate alarm
                y = (GFM{N} > GFM_th) + (IFM{N} > IFM_th) + (CGMslope{N} > slope_th);
                y(y < 3) = 0; % all thresholds must be exceeded
                alarms{N,1} = alarm_engineering(y,ts);
            end
            % make FD experiment object to perform evaluation
            experiment = FDexperiment('pump',label,alarms,time,t1d_data,GFM_th,IFM_th,slope_th);
            % evaluate and report
            T(k,:) = experiment.evaluate().report();
            k = k + 1;
        end
    end
end
disp('done!')

%%
close all
for par = ["GFM_th","IFM_th","slope_th"]
    pROC(T.SE, T.FPday_avg, T.(char(par)),'msize',5,'labels',{'SE','FP/day',par},'connect',0);
    set(gcf,'Units','normalized','Position',[0.2 0.1 0.5 0.7])
end

%%
x = T.FPday_avg;
y = T.SE;

T.distance_from_optimal = (x.^2 + (1-y).^2);
[B, I] = sort(T.distance_from_optimal);
fprintf('Optimal parameters \n GFM: %g \n IFM: %g \n slope: %g \n',T.GFM_th(I(1)),T.IFM_th(I(1)),T.slope_th(I(1)))

figure('Color','w')
sz = 3;
col = [0.0 0.0 0.0];
scatter(x, y, sz,'MarkerEdgeColor',col, 'MarkerFaceColor',col, 'LineWidth',1.5)
hold on
scatter(x(I(1)), y(I(1)), 50, 'o', 'MarkerEdgeColor','r', 'LineWidth',1.5)
xlim([-0.05 1.25])
ylim([-0.05 1.05])
yticks(0:0.2:1)
xticks(0:0.2:1.2)
ylabel('Sensitivity')
xlabel('FP/day')
grid on

%%
T2 = sortrows(T,'distance_from_optimal');
T2(1:5,:)

%% save
save_dir = fullfile('tuning_results', datasets_folder);
if ~isfolder(save_dir)
    mkdir(save_dir)
end
save(fullfile(save_dir, 'tuning_results'), 'T2')

%%
fprintf('Elapsed time: %g minutes \n', script_start/60);
rmpath(genpath('src'))





