function [x, y] = ellipse(a, b, aphi, x0, y0, dummy)
%
% Calculate ellipse centered at [x0, y0] with angle of major axis aphi.
%

t = [0:360];

X = a * cosd(t);
Y = b * sind(t);

x = X * cosd(aphi) - Y * sind(aphi) + x0;
y = X * sind(aphi) + Y * cosd(aphi) + y0;

if nargin > 5
   semimajor_x1 = a * cos(aphi*pi/180);
   semimajor_y1 = a * sin(aphi*pi/180);
   semimajor_x2 = -a * cos(aphi*pi/180);
   semimajor_y2 = -a * sin(aphi*pi/180);

   semiminor_x1 = -b * sin(aphi*pi/180);
   semiminor_y1 = b * cos(aphi*pi/180);
   semiminor_x2 = b * sin(aphi*pi/180);
   semiminor_y2 = -b * cos(aphi*pi/180);

   [xltan, ii] = min(x);
   yltan = y(ii);
   [xrtan, ii] = max(x);
   yrtan = y(ii);

   [ybtan, ii] = min(y);
   xbtan = x(ii);
   [yutan, ii] = max(y);
   xutan = x(ii);

   figure(1);
   clf;
   plot(x, y);
   hold on;
   line([xltan, xltan], [ybtan, yutan]);
   line([xrtan, xrtan], [ybtan, yutan]);
   line([xltan, xrtan], [ybtan, ybtan]);
   line([xltan, xrtan], [yutan, yutan]);
   %line([semimajor_x1, semimajor_x2], [semimajor_y1, semimajor_y2]);
   %line([semiminor_x1, semiminor_x2], [semiminor_y1, semiminor_y2]);
   %axis equal;
   %set(gca, 'YDir','Reverse');
   %xlim([70 90]);
   %ylim([1.6 1.75]);
   grid on;
   box on;
end
