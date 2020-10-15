function my_slope = moving_slope(y, TW, Ts)
% Get slope of y inside moving window via 1st grade polynomial fit
%
% y: input signal
% TW: moving window length (in hours)
% Ts: sampling time in minutes

my_slope = zeros(size(y));

% time window in number of samples
tw = TW*60/Ts;

for t = tw+1:length(y)
    
    a = y(t-tw:t);
    x = 0:Ts:length(a)*Ts-Ts;
    
    G = [x' ones(size(x'))];
    theta = G\a;
    my_slope(t) = theta(1);
    
% %     this is slower
%     [p] = polyfit(x',a,1);
%     my_slope(t) = p(1);
    
    % debug
    %     figure
    %     hold on
    %     plot(x,y,'o')
    %     plot(x,p(1)*x+p(2),'r')
    
end

end