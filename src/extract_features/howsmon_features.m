function [GFM, IFM, CGMslope] = howsmon_features(CGM, insulin, K, LW, SW, Ts)
%
% Extracts features used in
% Howsmon et al., "CGM enables detection of LISAs", Sensors 2017
%
% CGM: Continuous Glucose Measurement (mg/dL)
% insulin: injected insulin (U/min)
% K: decaying constant
% LW: long window (hours)
% SW: short window (hours)
% Ts: sampling time (minutes)

% Get GFM
[GFM, ~, ~] = gfm(CGM, LW, SW, Ts);
% correction
GFM = GFM/60;

% Get IFM
[PIE, ~] = twoCompSysD(insulin, K, Ts);
[IFM, ~, ~] = ifm(PIE, LW, SW, Ts);

% Get CGM slope in short window
CGMslope = moving_slope(CGM, SW, Ts);

end
