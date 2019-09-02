function [trfnew] = shift_rft(trf, t_max)
%
% This version process time axis, do not change real data
% Using function [t_1stmax, wavef_1stmax] = shift_ampt(t, wavef, istest)
% to determine the time of free surface P t_1stmax of RZ receiver function
% Then use this value to shift RF time 
% ------ t_1stmax < 0, Peak of free surface P is at "-" of time axis
% ------ t_1stmax > 0, Peak of free surface P is at "+" of time axis
% All t of time axis change to t = t - t_max;
% Input:
% trf: time series
% t_max: t_1stmax
% Output:
% trfnew: shift time series
%

trfnew = trf - t_max;
