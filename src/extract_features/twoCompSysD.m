function [X1, X2] = twoCompSysD(u, K, Ts, varargin)
%
% Discrete two-compartmental system
%
% my beautiful drawing below:
%
%        state1      state2
%  ----> ( X2 ) ---> ( X1 ) --->
%  u(t)          K           K
%
%

p = inputParser;
addOptional(p,'x0',[],@isnumeric);
parse(p,varargin{:});
x0=p.Results.x0;

%% define SS system
A = [1-K K; 0 1-K];
B = [0; 1];
C = [1 0];
D = 0;

% sys = ss(A,B,C,D,Ts,'TimeUnit','minutes');
sys = ss(A,B,C,D,1,'TimeUnit','minutes');
sys = d2d(sys,Ts);

%% filter

% time vector
t = 0:Ts:length(u)*Ts-Ts;

% filtering
if ~isempty(x0)
    [y, t, x] = lsim(sys, u, t, x0);
else
    [y, t, x] = lsim(sys, u, t);
end

% outputs
X1 = x(:,1);
X2 = x(:,2);

end

