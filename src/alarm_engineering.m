function alarm_new = alarm_engineering(alarm,Ts)

% get only first alarm in consecutive alarms
y_first = rising_edge(alarm);

% snooze mechanism
minimum_snooze_time = 180/Ts;
% [isON,predicted_off_time] = alarm_snooze(alarm,minimum_snooze_time);
isON = alarm_snooze_simple(alarm,minimum_snooze_time);

% apply snooze to first alarms
% alarm_new = and(isON,y_first).*1;
alarm_new = and(isON,y_first) .* y_first;

end

function ynew = rising_edge(y)
% Apply rising edge trigger

ynew = zeros(length(y),1);
x = find(y);

if ~isempty(x)
    edges_ = [x(1); x(find(diff(x) > 1) + 1)];
%     ynew(edges_) = 1;
    ynew(edges_) = y(edges_);
end


end


function isON = alarm_snooze_simple(alarm_signal,snooze_time)

y = rising_edge(alarm_signal);

pos = find(y);
for k = pos(:)'
   y(k+1:k+1+snooze_time) = 0; 
end

y = y(1:length(alarm_signal));
isON = boolean(y);

% figure
% subplot(2,1,1)
% plot(alarm_signal)
% subplot(2,1,2)
% plot(y)

end

function [isON,predicted_off_time] = alarm_snooze(alarm_signal,minimum_snooze_time)
% Determine on/off status for system based on an alarm signal.
%
% alarm_signal: 1 or 0, depending whether the alarm is raised or not
% Ts: sampling time
%
% System is temporarily shut off for a default duration whenever an alarm is raised during on mode.
% For each alarm raised during off time, off time is prolonged by one sample.
%

if minimum_snooze_time>0
    %% ----------------- shutoff using default time -----------------
    
    % initialization
    predicted_off_time=zeros(size(alarm_signal));
    off_time=0;
    
    % for each sample of the alarm signal
    for j=1:length(alarm_signal)-1
        
        switch alarm_signal(j)
            
            case 1 % === alarm raised ===
                
                % if it's the first raised during on time jump to default off time, otherwise increase the predicted off time
                if off_time == 0
                    off_time=minimum_snooze_time;
                else
                    off_time=off_time+1;
                end
                
                
            case 0 % === alarm not raised ===
                
                % 1 sample on off time was spent, decrease but don't go below zero
                off_time=off_time-1;
                off_time=max(off_time,0);
                
        end
        
        predicted_off_time(j+1)=off_time;
        
    end
    
else
    
    %% ----------------- no shutoff at all -----------------
    
    predicted_off_time=zeros(size(alarm_signal));
    
end

%% on/off status in binary notation
isON = (predicted_off_time == 0);

end