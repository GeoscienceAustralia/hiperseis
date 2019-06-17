function [machineformat] = sacmachine(sacfilename)

% Inspect machine format of sac file machine format: 'ieee-be' or 'ieee-le'
% Input:
% sacfilename: sac file name including dir
% output:
% machineformat: machine data format choice, consisten with fopen

% Open binary file as IEEE floating point with big-endian byte ordering,
% if SAC file is written in Unix.
% Open binary file as IEEE floating point with little-endian byte ordering,
% if SAC file is written in Linux.

% unused field in header:
% headi(6,3:5), headi(7,1:5), headi(8,5) = -12345
% headf(13,5), headf(14,1:5) = -12345.0

fid=fopen(sacfilename, 'r','ieee-be');
if (fid==-1)
  error(['can not open sacfile ', sacfilename]);
end

headf = fread(fid, [5, 14], 'float32');
headi = fread(fid, [5, 8], 'int32');
headc = fread(fid, [24, 8], 'char');
headf = headf'; headi = headi'; headc = headc';

if headi(7,1) == -12345 || headi(7,2) == -12345 || headi(7,3) == -12345 || headi(7,4) == -12345 || headi(7,5) == -12345 || ...
      headf(14,1) == -12345.0 || headf(14,2) == -12345.0 || headf(14,3) == -12345.0 || headf(14,4) == -12345.0 || headf(14,5) == -12345.0
   machineformat = 'ieee-be';
   fclose(fid);
   return;
else
   fclose(fid);

   fid=fopen(sacfilename, 'r','ieee-le');
   if (fid==-1)
      error(['can not open sacfile ', sacfilename]);
   end

   headf = fread(fid, [5, 14], 'float32');
   headi = fread(fid, [5, 8], 'int32');
   headc = fread(fid, [24, 8], 'char');
   headf = headf'; headi = headi'; headc = headc';

   if headi(7,1) == -12345 || headi(7,2) == -12345 || headi(7,3) == -12345 || headi(7,4) == -12345 || headi(7,5) == -12345 || ...
         headf(14,1) == -12345.0 || headf(14,2) == -12345.0 || headf(14,3) == -12345.0 || headf(14,4) == -12345.0 || headf(14,5) == -12345.0
      machineformat = 'ieee-le';
      fclose(fid);
      return;
   else
      machineformat = 'unknown';
      fclose(fid);
      return;
   end
   
end
   


