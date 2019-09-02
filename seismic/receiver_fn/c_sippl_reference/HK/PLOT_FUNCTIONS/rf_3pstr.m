function [t_0p1s, t_2p1s, t_1p2s] = rf_3pstr(rayp, vp, H, k)
%
% Calculate the travel time of phases 0p1s, 2p1s, and 1p2s from direct P
% according to a given Vp, H, and k for different ray parameters.
% Input:
% rayp: ray parameters (sec/deg). It could be an array.
% vp: if scalar, average crustal velocity
%     if vector, crustal velocity model of each layer
% H: if scalar, Moho depth
%    if vector, thickness of each layer 
% k: if scalar, Vp/Vs ratio
%    if vector, crustal Vp/Vs ratio of each layer
% If Vp, H, k are vectors, they must have the same length
% Output:
% t_0p1s: P-S travel time
% t_2p1s: 2p1s travel time
% t_1p2s: 1p2s travel time
% The three phases could be arrays if rayp is an array
% The time is calculated from direct P
%

r0 = 6371.227;

if isscalar(vp) && isscalar(H) && isscalar(k)
   p = rayp * 180.0 / (r0 - H) / pi;
   slow_p = sqrt(1 / vp / vp - p .* p);
   tp = H * slow_p;
   
   slow_s = sqrt(k * k / vp / vp - p .* p);
   ts = H * slow_s;   
   
   t_0p1s = ts - tp;
   t_2p1s = ts + tp;
   t_1p2s = 2 * ts;
   t_0p1s = tocol(t_0p1s);
   t_2p1s = tocol(t_2p1s);
   t_1p2s = tocol(t_1p2s);

else
   if length(vp) ~= length(H) || length(vp) ~= length(k) || length(H) ~= length(k)
      error('vp, H, k have no the same length');
   end
   
   vp = tocol(vp);
   H = tocol(H);
   k = tocol(k);

   t_0p1s = zeros(length(rayp), 1);
   t_2p1s = zeros(length(rayp), 1);
   t_1p2s = zeros(length(rayp), 1);

   for ii = 1:length(rayp)
      
      p = rayp(ii) .* 180.0 ./ (r0 - cumsum(H)) ./ pi;

      slow_p = sqrt(1 ./ vp ./ vp - p .* p);
      tp = sum(H .* slow_p);
      
      slow_s = sqrt(k .* k ./ vp ./ vp - p .* p);
      ts = sum(H .* slow_s);

      t_0p1s(ii) = ts - tp;
      t_2p1s(ii) = ts + tp;
      t_1p2s(ii) = 2 * ts;
   end
         
end   
