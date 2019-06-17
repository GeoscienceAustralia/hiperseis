function [time, data, headf, headi, headc] = rsac3(sacfilename,machineformat)

% Read all data of one single sac file and output the two column ascii data.
% Input:
% sacfilename: sac file name including dir
% machineformat: machine data format choice, consisten with fopen
% Output:
% time: time column
% data: data column
% headf, headi, headc: 3 types of sac header

% For Unix, open binary file as IEEE floating point with big-endian byte ordering,
% For Linux, open binary file as IEEE floating point with little-endian byte ordering
fid=fopen(sacfilename, 'r', machineformat);

if (fid==-1)
  disp('can not open input data file format, press CTRL-C to exit \n');
  pause;
end

headf = fread(fid, [5, 14], 'float32');
headi = fread(fid, [5, 8], 'int32');
headc = fread(fid, [24, 8], 'char');
headf = headf'; 
headi = headi'; 
headc = headc';

delta = headf(1, 1);
% --- Determine the accuracy of delta
ii = 1;
while(1)
   deltar = round(delta * 10^ii) / 10^ii;
   if abs(deltar - delta) < 1 / 10^(ii+1)
      ndigts = ii;
      break;
   else
      ii = ii + 1;
   end
end

b = headf(2, 1);
b = round(b * 10^ndigts) / 10^ndigts;
% e = headf(2, 2);
npts=headi(2, 5);

time = (b : delta : delta*(npts-1)+b)';

data = fread(fid, npts, 'float32');

fclose(fid);





