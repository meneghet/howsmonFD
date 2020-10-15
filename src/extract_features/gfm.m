function [GFM, CGM_SW, CGM_LW] = gfm(CGM, LW, SW, Ts)
%
% Calculates Glucose Fault Metric (GFM) as shown in
% Howsmon et al., "CGM enables detection of LISAs", Sensors 2017
%
% inputs
% - CGM: continuous glucose monitor data [mg/dL]
% - LW: length of long window for moving average [hours]
% - SW: length of short window for moving average [hours]
% - Ts: sampling time [minutes]
% outputs
% - GFM: glucose fault metric
% - CGM_SW: short window moving average
% - CGM_LW: long window moving average


% Moving average
LW_samples = LW*60/Ts;
SW_samples = SW*60/Ts;
CGM_SW = movmean(CGM,[SW_samples 0]);
CGM_LW = movmean(CGM,[LW_samples 0]);

% AUC
AUC = (CGM_SW-CGM_LW) * Ts;

% GFM
GFM = zeros(length(CGM),1);
for k = 2:length(CGM)
    if (CGM_SW(k)/CGM_LW(k))>1
        GFM(k) = GFM(k-1)+AUC(k);
    else
        GFM(k) = 0;
    end
end

end
