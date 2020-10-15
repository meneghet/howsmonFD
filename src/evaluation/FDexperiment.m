classdef FDexperiment
    properties %(Access = private)
        
        % type of fault simulated
        fault_type
        
        labels
        alarms
        
        % simulation time
        time
        
        % T1D data
        t1d_data
        
        % FD settings
        GFM_th
        IFM_th
        slope_th
        
        % results
        KPI
        detection_times
        single_cm
        FP_pos
        FP_source
        clf
    end
    
    methods
        
        
        function obj = FDexperiment(fault_type,labels,alarms,time,t1d_data,GFM_th,IFM_th,slope_th)
            
            obj.fault_type = fault_type;
            
            obj.labels = labels;
            obj.alarms = alarms;
            
            obj.time = time;
            
            obj.t1d_data = t1d_data;
            
            obj.GFM_th = GFM_th;
            obj.IFM_th = IFM_th;
            obj.slope_th = slope_th;
            
            obj.KPI = struct();
            obj.detection_times = [];
            obj.single_cm = table();
            obj.FP_pos = [];
            obj.FP_source = [];
            obj.clf = [];
        end
        
        
        function obj = evaluate(obj, varargin)
            % Evaluate the experiment outcome
            
            p = inputParser;
            addParameter(p,'calibrations',1,@isnumeric)
            addParameter(p,'day_portion','all_day',@ischar)
            addParameter(p,'snooze_alarms',0,@isnumeric)
            parse(p,varargin{:});
            calibrations = p.Results.calibrations;
            day_portion = p.Results.day_portion;
            snooze_alarms = p.Results.snooze_alarms;
            
            % apply alarm snooze, if requested
            obj.alarms = obj.get_alarms('snooze_alarms', snooze_alarms);
            
            [obj.KPI, obj.detection_times, obj.single_cm, obj.FP_pos, obj.clf, obj.FP_source] = ...
                evaluate_experiment(obj, day_portion, calibrations);
        end
        
        
        function T = report(obj)
            % Report results in table
            
            if isempty(fieldnames(obj.KPI))
                disp('did not perform evaluation of this experiment')
                T = [];
            else
                T = table(obj.KPI.sensitivity,...
                    obj.KPI.precision,...
                    obj.KPI.FPday_avg,...
                    obj.KPI.FPday_avg_1,...
                    obj.KPI.FPday_avg_2,...
                    obj.KPI.FPday_med,...
                    round(nanmean(cell2mat(obj.detection_times)),2),...
                    obj.GFM_th, obj.IFM_th, obj.slope_th,...
                    'VariableNames',...
                    {'SE','PR','FPday_avg','FPday_avg1','FPday_avg2','FPday_med','Tdet','GFM_th','IFM_th','slope_th'});
            end
        end
        
        
        function FP_stats = FP_stats_report(obj, varargin)
            % Make a report on some FP statistics
            
            p = inputParser;
            addParameter(p,'plot_flag',0,@isnumeric)
            parse(p,varargin{:});
            plot_flag = p.Results.plot_flag;
            
            % initalize output
            FP_stats = struct();
            
            % distance from fault time unit
            tu = (obj.time(2)-obj.time(1))*24; % [hours]
            
            if isempty(fieldnames(obj.KPI))
                disp('did not perform evaluation of this experiment')
            else
                % for each patient
                for p = 1:length(obj.FP_pos)
                    FP = obj.FP_pos{p}; % every FP position
                    
                    % FP time of the day
                    t = obj.time(FP);
                    FP_stats(p,1).ToD = mod(t,1)*24;
                    
                    % get position of meals
                    meals = obj.t1d_data{p}.meal;
                    meal_pos = find(meals);
                    
                    for k = 1:length(FP) % for each FP
                        time_from_meals = FP(k) - meal_pos;
                        % only account for past meals
                        time_from_meals = time_from_meals(time_from_meals>=0);
                        if isempty(time_from_meals)
                            time_from_meals = NaN;
                        end
                        % take closest meal
                        [dist,pos] = min(time_from_meals);
                        % meal related FP stats
                        FP_stats(p,1).time_of_last_meal(k) = meal_pos(pos);
                        FP_stats(p,1).meal_amount(k) = meals(meal_pos(pos));
                        FP_stats(p,1).since_last_meal(k) = dist*tu;
%                         % distance from fault
%                         fault_time = find(obj.label{p});
%                         FP_from_fault(k) = FP(k) - fault_time;
                    end
                end
                % plot FP stats histograms
                if plot_flag
                    plot_FP_stats(FP_stats)
                end
            end
            
        end
        
        
        function plot(obj, N, varargin)
            p = inputParser;
            addParameter(p,'snooze_alarms',0,@isnumeric)
            parse(p,varargin{:});
            snooze_alarms = p.Results.snooze_alarms;
            
            %%
            my_alarms = obj.get_alarms('snooze_alarms', snooze_alarms);
            my_alarm = logical(my_alarms{N})*1;
            
            my_label = obj.labels{N};
            
            my_t = obj.time;
            
            labels_pos = find(my_label);
            for k =1:length(labels_pos)
                x = labels_pos(k);
                TPtime = my_label(x);
                my_label(x:x+TPtime*60/5) = 1;
            end
            
            my_t = datetime(datestr(my_t));
            
            my_label = double(my_label);
            isTP = ((my_label .* my_alarm) > 0);
%             isFP = ((my_label - my_alarm) < 0);
            FP_x = obj.FP_pos{N};
            
            figure
            ax(1) = subplot(2,1,1);
            L = [];
            L(1) = plot(my_t, obj.t1d_data{N}.CGM, 'k', 'DisplayName', 'labels');
                        
            ax(2) = subplot(2,1,2);
            L = [];
            L(1) = plot(my_t, my_label, 'k', 'DisplayName', 'labels');
            hold on
            L(2) = stem(my_t(isTP), my_alarm(isTP), 'g', 'DisplayName', 'TP');
%             L(3) = stem(my_t(isFP), my_alarm(isFP), 'r', 'DisplayName', 'FP');
            L(3) = stem(my_t(FP_x), my_alarm(FP_x), 'r', 'DisplayName', 'FP');
            legend(L)
            xlim([my_t(1) my_t(end)])
            
            linkaxes(ax, 'x')
            
        end
        
        
        function fp_per_day(obj)
            % Plots distribution of FP/day among patients
            
            if isempty(fieldnames(obj.KPI))
                disp('did not perform evaluation of this experiment')
            else
                n_days = obj.time(end)-obj.time(1);
                figure('Color','w')
                boxplot(obj.single_cm.FP/n_days)
                title('FP/day')
            end
        end
        
        
        function my_alarms = get_alarms(obj,varargin)
            % Return alarms in double format
            % optional: silence during day or night
            % optional: set custom numeric value
            % optional: silence after calibrations
            % optional: snooze alarm for a duration
            
            p = inputParser;
            addParameter(p,'day_portion','all_day',@ischar)
            addParameter(p,'label',1,@isnumeric)
            addParameter(p,'removeCalibrations',0,@isnumeric)
            addParameter(p,'snooze_alarms',0,@isnumeric)
            parse(p,varargin{:});
            day_portion = p.Results.day_portion;
            label = p.Results.label;
            removeCalibrations = p.Results.removeCalibrations;
            snooze_alarms = p.Results.snooze_alarms;
            
            if snooze_alarms > 0
                obj =  obj.snooze_alarms(snooze_alarms);
            end
            
            % make a copy
            my_alarms = obj.alarms;
            
            % for each subject
            for N = 1:length(my_alarms)
                
                % put label
                A = my_alarms{N};
                A(A > 0) = A(A > 0) * label;
                my_alarms{N,1} = A;
                %                 my_alarms{N} = my_alarms{N}.*label;
                
                % get only night or day if requester
                switch day_portion
                    case 'all_day'
                        %
                    case 'night_only'
                        [isDay,~] = day_night_portions(obj.time);
                        %                         my_alarms{N}(isDay) = 0;
                        A(isDay) = 0;
                    case 'day_only'
                        [~,isNight] = day_night_portions(obj.time);
                        %                         my_alarms{N}(isNight) = 0;
                        A(isNight) = 0;
                end
                % remove calibrations portions if requested
                if removeCalibrations
                    calibration_h = [12 24]; % hours of the day
                    calibration_duration  = 30; % minutes
                    isCalibration = time_portions(obj.time, calibration_h, calibration_duration);
                    %                     my_alarms{N}(isCalibration) = 0;
                    A(isCalibration) = 0;
                end
                
                my_alarms{N,1} = A;
            end
        end
        
        
        function obj = snooze_alarms(obj, snooze_duration)
            % Apply snooze for specified duration after an alarm goes off
            %
            % snooze_duration: number of samples to snooze
            
            % for each experiment subject
            for N = 1:length(obj.alarms)
                alarm = obj.alarms{N};
                % snooze alarms
                y = rising_edge(alarm);
                pos = find(y);
                for k = pos(:)' % for every alarm
                    y(k+1:k+1+snooze_duration) = 0;
                end
                y = y(1:length(alarm)); % it can happen
                isON = boolean(y);
                % re-save in object
                obj.alarms{N} = isON.*alarm;
            end
        end
        
        
        function obj = snooze_meals(obj,snooze_duration)
            % Apply snooze after a meal is announced
            
            % for each experiment subject
            for N = 1:length(obj.alarms)
                alarm = obj.alarms{N};
                meals = obj.t1d_data{N}.meal;
                                
                % meal positions
                meals = rising_edge(meals);
                m_pos = find(meals);
                
                % snooze alarms after meals
                for k = 1:length(m_pos)
                    snooze_start = m_pos(k);
                    snooze_end = min(length(alarm), m_pos(k)+snooze_duration);
                    alarm(snooze_start:snooze_end) = 0;
                end
                
                % re-save in object
                obj.alarms{N} = alarm;
            end
        end
        
    end
end

%% utilities
function [isDay,isNight] = day_night_portions(time)
% detect portions of the day (day or night)
%
% time: time in datenum
% isDay: day indices
% isNight: night indices

day_start = 7; %h
day_end = 24; %h

isDay = ( mod(time,1) >= day_start/24 & mod(time,1) < day_end/24 );
isNight = ( mod(time,1) < day_start/24 | mod(time,1) >= day_end/24 );

end


function isPortion = time_portions(time, HH, MM)
% Return portions of time from HH:00 to HH:MM

time = datenum(time);

% set calibration interval (in datetime)
HH = mod(HH,24);
HH = HH/24;
MM = MM/60/24;
time_interval = [HH; HH + MM];

% remove calibration portions
for ind = 1:size(time_interval,2)
    
    cal_start = time_interval(1,ind);
    cal_end = time_interval(2,ind);
    
    isPortion = ( mod(time,1) >= cal_start & mod(time,1) <= cal_end );
end

end


function ynew = rising_edge(y)
% Apply rising edge trigger

ynew = zeros(length(y),1);

x = find(y);
if ~isempty(x)
    edges_= [x(1); x(find(diff(x) > 1) + 1)];
    ynew(edges_) = 1;
end

end

%%
function plot_FP_stats(FP_stats)

FP_per_patient = zeros(length(FP_stats),1);
for n = 1:length(FP_stats)
    FP_per_patient(n) = length(FP_stats(n).ToD);
end

ToD = cat(1,FP_stats.ToD);
time_since_meal = [FP_stats.since_last_meal]';

% FP_from_fault = [FP_stats.FP_from_fault];
% FP_from_fault(FP_from_fault<0)=inf;

% plot FP times
figure('Color','w')

% FP count
FP_hist_axes = subplot(3,1,1);
my_bins = 0.5:1:max(FP_per_patient) + 1.5;
histogram(FP_per_patient, 'BinEdges', my_bins);
xticks(1:max(FP_per_patient))
title('FP per patients');
xlim([0.5 max(FP_per_patient) + 0.5]);
ylabel('# patients');
xlabel('# FP')

% FP time of the day
subplot(3,1,2)
histogram(ToD,'BinEdges',0:0.5:24)
xticks(0:1:24)
title('FP times')
ylabel('# FP')
xlabel('ToD (h)')
xlim([0 24])


% FP time since meals
subplot(3,1,3)
my_bins = 0.5:1:max(FP_per_patient) + 1.5;
histogram(time_since_meal, 'BinEdges', my_bins)
xticks(0:1:max(time_since_meal))
title('FP times since last meal')
ylabel('# FP')
xlabel('time (h)')
xlim([0 12])

% % FP time from fault
% subplot(4,1,4)
% post_fault_duration=18; % [hours]
% post_fault_in_minutes=post_fault_duration*60/5;
% x=FP_from_fault<post_fault_in_minutes;
% histogram(x,'BinEdges',[-0.5:1:1.5])
% xticks([0,1])
% xticklabels({'away from fault',[num2str(post_fault_duration) ' hours after fault start']})
% ylabel('# FP')

% see below
setappdata(FP_hist_axes,'FP_count',FP_per_patient);
cursorMode = datacursormode(gcf);
set(cursorMode,'enable','on','UpdateFcn',@FP_callback);

end

function output_txt = FP_callback(~,event_obj)
% Callback function for "FP per patients" histogram

pos = get(event_obj,'Position');
bin_pos = pos(1);
bin_count = pos(2);

hAxes = get(get(event_obj,'Target'),'Parent');
FP_count = getappdata(hAxes,'FP_count');
elements_in_bin = find(FP_count == bin_pos);

el_str = strcat(num2str(elements_in_bin(:)'));
output_txt = {['bin: ',num2str(bin_pos)],['count: ',num2str(bin_count)],['subj: ',num2str(el_str)]};

end



