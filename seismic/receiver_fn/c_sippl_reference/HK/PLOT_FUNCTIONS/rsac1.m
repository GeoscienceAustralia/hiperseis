function [time, data, head] = rsac1(sacfilename,machineformat)
%
% Read all data of one single sac file and output the two column ascii data.
% Input:
% sacfilename: sac file name including dir
% machineformat: machine data format choice, consisten with fopen
% Output:
% time: time column
% data: data column
% head: sac header structure
% 

% For Unix, open binary file as IEEE floating point with big-endian byte ordering,
% For Linux, open binary file as IEEE floating point with little-endian byte ordering
%
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
npts=headi(2, 5);
e = delta * (npts - 1) + b;
e = round(e * 10^ndigts) / 10^ndigts;

%time = (b : delta : e)';
time = linspace(b, e, npts)';
data = fread(fid, npts, 'float32');

fclose(fid);

%----------- File Field --------------------
head.npts = headi(2,5);
head.nvhdr = headi(2,2);
head.b = b; %headf(2,1);
head.e = e; %headf(2,2);
head.delta = delta; %headf(1,1);
head.odelta = headf(1,5);
head.iftype = headi(4,1);
head.leven = headi(8,1);
head.idep = headi(4,2);
head.depmin = headf(1,2);
head.depmax = headf(1,3);
head.scale = headf(1,4);
head.depmen = headf(12,2);
head.nzyear = headi(1,1);
head.nzjday = headi(1,2);
head.nzhour = headi(1,3);
head.nzmin = headi(1,4);
head.nzsec = headi(1,5);
head.nzmsec = headi(2,1);
head.iztype = headi(4,3);
head.o = headf(2,3);
head.ko = strtrim(char(headc(2,9:16)));

%---------------Phase Picks Field --------------------
head.a = headf(2,4);
head.ka = strtrim(char(headc(2,17:24)));
head.f = headf(5,1);
head.kf = strtrim(char(headc(6,9:16)));
head.t0 = headf(3,1);
head.t1 = headf(3,2);
head.t2 = headf(3,3);
head.t3 = headf(3,4);
head.t4 = headf(3,5);
head.t5 = headf(4,1);
head.t6 = headf(4,2);
head.t7 = headf(4,3);
head.t8 = headf(4,4);
head.t9 = headf(4,5);
head.kt0 = strtrim(char(headc(3,1:8)));
head.kt1 = strtrim(char(headc(3,9:16)));
head.kt2 = strtrim(char(headc(3,17:24)));
head.kt3 = strtrim(char(headc(4,1:8)));
head.kt4 = strtrim(char(headc(4,9:16)));
head.kt5 = strtrim(char(headc(4,17:24)));
head.kt6 = strtrim(char(headc(5,1:8)));
head.kt7 = strtrim(char(headc(5,9:16)));
head.kt8 = strtrim(char(headc(5,17:24)));
head.kt9 = strtrim(char(headc(6,1:8)));

%------------------Instrument Field --------------------
head.kinst = strtrim(char(headc(8,17:24)));
head.iinst = headi(4,5);
head.resp0 = headf(5,2);
head.resp1 = headf(5,3);
head.resp2 = headf(5,4);
head.resp3 = headf(5,5);
head.resp4 = headf(6,1);
head.resp5 = headf(6,2);
head.resp6 = headf(6,3);
head.resp7 = headf(6,4);
head.resp8 = headf(6,5);
head.resp9 = headf(7,1);

%---------------Station Field --------------------
head.knetwk = strtrim(char(headc(8,1:8)));
head.kstnm = strtrim(char(headc(1,1:8)));
head.istreg = headi(5,1);
head.stla = headf(7,2);
head.stlo = headf(7,3);
head.stel = headf(7,4);
head.stdp = headf(7,5);
head.cmpaz = headf(12,3);
head.cmpinc = headf(12,4);
head.kcmpnm = strtrim(char(headc(7,17:24)));
head.lpspol = headi(8,2);

%--------------- Event Field --------------------
head.kevnm = strtrim(char(headc(1,9:24)));
head.ievreg = headi(5,2);
head.evla = headf(8,1);
head.evlo = headf(8,2);
head.evel = headf(8,3);
head.evdp = headf(8,4);
head.mag = headf(8,5);
head.imagtyp = headi(6,1);
head.imagsrc = headi(6,2);
head.ievtyp = headi(5,3);
head.norid = headi(2,3);
head.nevid = headi(2,4);
head.nwfid = headi(3,2);
head.khole = strtrim(char(headc(2,1:8)));
head.dist = headf(11,1);
head.az = headf(11,2);
head.baz = headf(11,3);
head.gcarc = headf(11,4);

%--------------------Miscellaneous Fields --------------------
head.internal_f2_5 = headf(2,5);
head.internal_f11_5 = headf(11,5);
head.internal_f12_1 = headf(12,1);
head.internal_i3_1 = headi(3,1);
head.user0 = headf(9,1);
head.user1 = headf(9,2);
head.user2 = headf(9,3);
head.user3 = headf(9,4);
head.user4 = headf(9,5);
head.user5 = headf(10,1);
head.user6 = headf(10,2);
head.user7 = headf(10,3);
head.user8 = headf(10,4);
head.user9 = headf(10,5);
head.lcalda = headi(8,4);
head.iqual = headi(5,4);
head.isynth = headi(5,5);
head.kdatrd = strtrim(char(headc(8,9:16))); 
head.lovrok = headi(8,3);
head.mxsize = headi(3,3);
head.nysize = headi(3,4);
head.xminimum = headf(12,5);
head.xmaximum = headf(13,1);
head.yminimum = headf(13,2);
head.ymaximum = headf(13,3);
head.unused_f13_4 = headf(13,4);
head.unused_f13_5 = headf(13,5);
head.unused_f14_1 = headf(14,1);
head.unused_f14_2 = headf(14,2);
head.unused_f14_3 = headf(14,3);
head.unused_f14_4 = headf(14,4);
head.unused_f14_5 = headf(14,5);
head.nxsize = headi(3,3);
head.nysize = headi(3,4);
head.unused_i3_5 = headi(3,5);
head.unused_i4_4 = headi(4,4);
head.unused_i6_3 = headi(6,3);
head.unused_i6_4 = headi(6,4);
head.unused_i6_5 = headi(6,5);
head.unused_i7_1 = headi(7,1);
head.unused_i7_2 = headi(7,2);
head.unused_i7_3 = headi(7,3);
head.unused_i7_4 = headi(7,4);
head.unused_i7_5 = headi(7,5);
head.unused_i8_5 = headi(8,5);
head.kuser0 = strtrim(char(headc(6,17:24)));
head.kuser1 = strtrim(char(headc(7,1:8)));
head.kuser2 = strtrim(char(headc(7,9:16)));
