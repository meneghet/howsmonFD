function [KPI, detection_times, single_cm, FP_pos, classification, FP_source] = ...
    evaluate_experiment(experiment, day_portion, removeCalibrations)
%
% TODO: description
%

% silence alarms if requested
all_alarms = experiment.get_alarms('day_portion',day_portion,'removeCalibration',removeCalibrations);
% time is common for all subjects
time = datenum(experiment.time);
% sampling time
Ts = round((time(2)-time(1))*24*60);

% perform evaluation for each subject ---------------------
patients_list = 1:length(all_alarms);
for N = patients_list
    
    % input data for current patient
    labels = experiment.labels{N};
    alarms = all_alarms{N};
    
    % evaluation on subject ------------------------------------------
    [cm, det_time, clf, fp_pos, fp_source] = evaluate_subject(alarms, labels, Ts);
    
    % confusion matrix
    single_cm(N,1) = cm;
    % detection time
    detection_times{N,1} = det_time;
    % FP positions on alarms vector
    FP_pos{N,1} = fp_pos;
    FP_source{N,1} = fp_source;
    % classification of events
    classification{N,1} = clf;
end

single_cm = struct2table(single_cm);

% complicated way to make a table
total_cm = array2table(sum(table2array(single_cm),1));
total_cm.Properties.VariableNames = single_cm.Properties.VariableNames;
% get Key Performance Indices
n_days = time(end)-time(1);
KPI = get_KPI(total_cm, single_cm, n_days, FP_source);

end


%% evaluate each subject
function [cm, det_time, classification, FP_pos, FP_source] = evaluate_subject(alarms, label, Ts)
% Perform evaluation for each subject

% read data
if size(label,2) > 1
    label = sum(label,2);
end
label = label(:);
alarms = alarms(:);

% performance evaluation  
[cm, det_time, FP_pos, classification, FP_source] = evaluate_alarms(label, alarms, Ts);

end


%% evaluate alarms
function [cm, detection_time, FP_pos, classification, FP_source] = evaluate_alarms(label, alarms, Ts)
%
% Evaluate alarms using event-based evaluation
%
% TP: alarm near a fault
% FP: alarm away from a fault
% FN: fault without alarm
% 
% TN: no way to count that
%
%
% INPUTS
% label: Mx1 vector, ground truth of fault portions
% alarms: Mx1 vector, portions that were labeled by FD algorithm
% 
% RETURNS
% cm: confusion matrix
% detection_time: list of detection times
% FP_pos: position of FPs in alarms vector
% classification: classification of detected events
%


% times of fault event, alarm events and meals
fault_time = find(label);
FD_pos = find(alarms);

% no fault are present -> set time of fault to inf
if isempty(fault_time)
    fault_time = NaN;
else
    n_faults = length(fault_time);
end

% initialization
TP = 0;
FP = 0;
detection_time = nan * ones(size(fault_time));
classification = nan * ones(size(fault_time));

FP_pos = [];
FP_source = [];

% for each alarm raised, check if it was raised correctly
for ind_alarm = 1:length(FD_pos)
    
    % ------- determine if it is TP -------------------------------------- 
    % determine distance from previous fault
    fault_relative_distance = FD_pos(ind_alarm) - fault_time;
    time_from_previous_faults = fault_relative_distance;
    time_from_previous_faults(time_from_previous_faults <= 0) = NaN;
    [D, previous_fault] = min(time_from_previous_faults);
    % was the alarm raised after a fault?
    isAfterFault = ~isnan(D);
    % was the alarm raised in time?
    if isAfterFault
        TPmax = label(fault_time(previous_fault))*60/Ts;
        foundInTime = D < TPmax;
    else
        TPmax = 0;
        foundInTime = 0;
    end
    % was the alarm raised after a minimum reasonable time?
    TPmin = floor(TPmax/60/Ts);
    isPossible = (D > TPmin);
    % is TP?
    isTP = (isAfterFault && foundInTime && isPossible);
    % -------------------------------------------------------------------
    
    if isTP % alarm raised correctly -------------------
        
        % increase TP counter
        TP = TP+1;
        % record detection time
        if isnan(detection_time(previous_fault))
            detection_time(previous_fault) = D*Ts;
        end
        % record classification based on alarm
        if isnan(classification(previous_fault))
            classification(previous_fault) = alarms(FD_pos(ind_alarm));
        end
        
    else % alarm raised incorrectly --------------------
        
        % we suspend the evaluation after a fault for a time equal to TPmax*2
        Tsuspend = TPmax * 2;
        % if the alarm was raised too late but within Tsuspend -> don't register FP
        isForgivable = (D < TPmax + Tsuspend);
        
        if isAfterFault && isForgivable
            % don't count FP
        else
            % increase FP count by 1
            FP = FP+1;
            % save FP position
            FP_pos = [FP_pos FD_pos(ind_alarm)];
            FP_source = [FP_source alarms(FD_pos(ind_alarm))];
        end        
    end
end
% -----------------------------------------------------------------

% TP and FN
if isempty(find(label,1))
    TP = 0;
    FN = 0;
else
    TP = min(n_faults,TP);
    FN = n_faults - TP;
end

% TN count (empirical calculation)
TPtime = 12; % totally random
total_time = length(alarms);
total_observations = floor(total_time/TPtime);
TN = total_observations - TP - FN - FP;

% make confusion matrix
cm = struct('TP',TP,'FP',FP,'FN',FN,'TN',TN);

end


%%
function KPI = get_KPI(cm, single_cm, n_days, FP_source)
% Get Key Performance Index

TP = cm.TP;
TN = cm.TN;
FP = cm.FP;
FN = cm.FN;

% Key Performance Indices
KPI.sensitivity = TP / (TP + FN);
KPI.precision = TP / (TP + FP);
KPI.F1score = 2 * TP / (2 * TP + FP + FN);

% nonsense
KPI.specificity = TN / (TN + FP);
KPI.accuracy = (TP + TN) / (TP + FN + TN + FP);

% FPs
KPI.FPday_med = median(single_cm.FP / n_days);
KPI.FPday_avg = mean(single_cm.FP/n_days);

% FPs by alarm that generated them
FP1 = sum([FP_source{:}] == 1);
FP2 = sum([FP_source{:}] == 2);
N = height(single_cm); % number of subjects
KPI.FPday_avg_1 = FP1/N/n_days;
KPI.FPday_avg_2 = FP2/N/n_days;

% round everything to 3rd figure -------
fnames = string(fieldnames(KPI));
for s = fnames(:)'
    KPI.(s) = round(KPI.(s),3);
end

% Confusion matrix
KPI.cm = cm;

end




