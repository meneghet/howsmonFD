function [IFM,PIE_LW,PIE_SW] = ifm(PIE, LW, SW, Ts)
%
% Calculates Insulin Fault Metric (IFM) as shown in
% Howsmon et al., "CGM enables detection of LISAs", Sensors 2017
%
% inputs
% - LW: length of long window for moving average [hours]
% - SW: length of short window for moving average [hours]
% - Ts: sampling time [minutes]
% outputs
% - IFM: insulin fault metric
% - PIE: plasma insulin estimation


% Moving average (short window)
SW_samples = SW*60/Ts;
PIE_SW = movmean(PIE,[SW_samples 0]);

% Moving average (long window)
LW_samples = LW*60/Ts;
PIE_LW = movmean(PIE,[LW_samples 0]);

% IFM
IFM = (PIE_SW./PIE_LW) - 1;

end
