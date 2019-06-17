function [nerr] = wsac1(sacfilename, head, data, machineformat)

% Write data into binary sac file.
% Input:
% sacfilename: sac file name including dir
% machineformat: machine data format choice, consisten with fopen
% head: structure for sac header
% data: waveform data
% if length(data) > npts (info in headi), only the first npts is written
% if length(data) < npts (info in headi), pad zeros upto npts into data and be written
% Output:
% nerr = 1 if no error

% For Unix, open binary file as IEEE floating point with big-endian byte ordering,
% For Linux, open binary file as IEEE floating point with little-endian byte ordering
% If machineformat is neglected, use native platform ordering

if nargin == 3
    fid=fopen(sacfilename, 'w');
elseif nargin == 4
    fid=fopen(sacfilename, 'w', machineformat);
end

if (fid==-1)
  error('can not open output data file format, press CTRL-C to exit \n');
end

%-------- Initialize headf, headi, headc
headf = ones(14, 5) .* -12345.0;
headi = ones(8, 5) .* -12345;
headc = ones(8, 24) .* abs(' ');
% headc(1, 1) = '-'; headc(1, 2) = '1'; headc(1, 3) = '2'; headc(1, 4) = '3'; headc(1, 5) = '4'; headc(1, 6) = '5';
% headc(1, 9) = '-'; headc(1, 10) = '1'; headc(1, 11) = '2'; headc(1, 12) = '3'; headc(1, 13) = '4'; headc(1, 14) = '5';
% for irow = 2:8
%    headc(irow, 1) = '-'; headc(irow, 2) = '1'; headc(irow, 3) = '2'; headc(irow, 4) = '3'; headc(irow, 5) = '4'; headc(irow, 6) = '5';
%    headc(irow, 9) = '-'; headc(irow, 10) = '1'; headc(irow, 11) = '2'; headc(irow, 12) = '3'; headc(irow, 13) = '4'; headc(irow, 14) = '5';
%    headc(irow, 17) = '-'; headc(irow, 18) = '1'; headc(irow, 19) = '2'; headc(irow, 20) = '3'; headc(irow, 21) = '4'; headc(irow, 22) = '5';
% end

%-------- Put head structure value to headf, headi, headc
headf(1,1) = head.delta;
headf(1,2) = head.depmin;
headf(1,3) = head.depmax;
headf(1,4) = head.scale;
headf(1,5) = head.odelta;
headf(2,1) = head.b;
headf(2,2) = head.e;
headf(2,3) = head.o;
headf(2,4) = head.a;
headf(2,5) = head.internal_f2_5;
headf(3,1) = head.t0;
headf(3,2) = head.t1;
headf(3,3) = head.t2;
headf(3,4) = head.t3;
headf(3,5) = head.t4;
headf(4,1) = head.t5;
headf(4,2) = head.t6;
headf(4,3) = head.t7;
headf(4,4) = head.t8;
headf(4,5) = head.t9;
headf(5,1) = head.f;
headf(5,2) = head.resp0;
headf(5,3) = head.resp1;
headf(5,4) = head.resp2;
headf(5,5) = head.resp3;
headf(6,1) = head.resp4;
headf(6,2) = head.resp5;
headf(6,3) = head.resp6;
headf(6,4) = head.resp7;
headf(6,5) = head.resp8;
headf(7,1) = head.resp9;
headf(7,2) = head.stla;
headf(7,3) = head.stlo;
headf(7,4) = head.stel;
headf(7,5) = head.stdp;
headf(8,1) = head.evla;
headf(8,2) = head.evlo;
headf(8,3) = head.evel;
headf(8,4) = head.evdp;
headf(8,5) = head.mag;
headf(9,1) = head.user0;
headf(9,2) = head.user1;
headf(9,3) = head.user2;
headf(9,4) = head.user3;
headf(9,5) = head.user4;
headf(10,1) = head.user5;
headf(10,2) = head.user6;
headf(10,3) = head.user7;
headf(10,4) = head.user8;
headf(10,5) = head.user9;
headf(11,1) = head.dist;
headf(11,2) = head.az;
headf(11,3) = head.baz;
headf(11,4) = head.gcarc;
headf(11,5) = head.internal_f11_5;
headf(12,1) = head.internal_f12_1;
headf(12,2) = head.depmen;
headf(12,3) = head.cmpaz;
headf(12,4) = head.cmpinc;
headf(12,5) = head.xminimum;
headf(13,1) = head.xmaximum;
headf(13,2) = head.yminimum;
headf(13,3) = head.ymaximum;
headf(13,4) = head.unused_f13_4;
headf(13,5) = head.unused_f13_5;
headf(14,1) = head.unused_f14_1;
headf(14,2) = head.unused_f14_2;
headf(14,3) = head.unused_f14_3;
headf(14,4) = head.unused_f14_4;
headf(14,5) = head.unused_f14_5;

headi(1,1) = head.nzyear;
headi(1,2) = head.nzjday;
headi(1,3) = head.nzhour;
headi(1,4) = head.nzmin;
headi(1,5) = head.nzsec;
headi(2,1) = head.nzmsec;
headi(2,2) = head.nvhdr;
headi(2,3) = head.norid;
headi(2,4) = head.nevid;
headi(2,5) = head.npts;
headi(3,1) = head.internal_i3_1;
headi(3,2) = head.nwfid;
headi(3,3) = head.nxsize;
headi(3,4) = head.nysize;
headi(3,5) = head.unused_i3_5;
headi(4,1) = head.iftype;
headi(4,2) = head.idep;
headi(4,3) = head.iztype;
headi(4,4) = head.unused_i4_4;
headi(4,5) = head.iinst;
headi(5,1) = head.istreg;
headi(5,2) = head.ievreg;
headi(5,3) = head.ievtyp;
headi(5,4) = head.iqual;
headi(5,5) = head.isynth;
headi(6,1) = head.imagtyp;
headi(6,2) = head.imagsrc;
headi(6,3) = head.unused_i6_3;
headi(6,4) = head.unused_i6_4;
headi(6,5) = head.unused_i6_5;
headi(7,1) = head.unused_i7_1;
headi(7,2) = head.unused_i7_2;
headi(7,3) = head.unused_i7_3;
headi(7,4) = head.unused_i7_4;
headi(7,5) = head.unused_i7_5;
headi(8,1) = head.leven;
headi(8,2) = head.lpspol;
headi(8,3) = head.lovrok;
headi(8,4) = head.lcalda;
headi(8,5) = head.unused_i8_5;

%---------------------------------------------------------------------------
nstr = length(head.kstnm);
if nstr <= 8
%    headc(1,1:8) = [32 32 32 32 32 32 32 32];
   headc(1,1:nstr) = head.kstnm;
else
   headc(1,1:8) = head.kstnm(1:8);
end
   
nstr = length(head.kevnm);
if nstr <= 16
   headc(1,9:9+nstr-1) = head.kevnm;
else
   headc(1,9:24) = head.kevnm(1:16);
end

%---------------------------------------------------------------------------
nstr = length(head.khole);
if nstr <= 8
   headc(2,1:nstr) = head.khole;
else
   headc(2,1:8) = head.khole(1:8);
end
   
nstr = length(head.ko);
if nstr <= 8
   headc(2,9:9+nstr-1) = head.ko;
else
   headc(2,9:16) = head.ko(1:8);
end

nstr = length(head.ka);
if nstr <= 8
   headc(2,17:17+nstr-1) = head.ka;
else
   headc(2,17:24) = head.ka(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.kt0);
if nstr <= 8
   headc(3,1:nstr) = head.kt0;
else
   headc(3,1:8) = head.kt0(1:8);
end

nstr = length(head.kt1);
if nstr <= 8
   headc(3,9:9+nstr-1) = head.kt1;
else
   headc(3,9:16) = head.kt1(1:8);
end

nstr = length(head.kt2);
if nstr <= 8
   headc(3,17:17+nstr-1) = head.kt2;
else
   headc(3,17:24) = head.kt2(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.kt3);
if nstr <= 8
   headc(4,1:nstr) = head.kt3;
else
   headc(4,1:8) = head.kt3(1:8);
end

nstr = length(head.kt4);
if nstr <= 8
   headc(4,9:9+nstr-1) = head.kt4;
else
   headc(4,9:16) = head.kt4(1:8);
end

nstr = length(head.kt5);
if nstr <= 8
   headc(4,17:17+nstr-1) = head.kt5;
else
   headc(4,17:24) = head.kt5(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.kt6);
if nstr <= 8
   headc(5,1:nstr) = head.kt6;
else
   headc(5,1:8) = head.kt6(1:8);
end

nstr = length(head.kt7);
if nstr <= 8
   headc(5,9:9+nstr-1) = head.kt7;
else
   headc(5,9:16) = head.kt7(1:8);
end

nstr = length(head.kt8);
if nstr <= 8
   headc(5,17:17+nstr-1) = head.kt8;
else
   headc(5,17:24) = head.kt8(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.kt9);
if nstr <= 8
   headc(6,1:nstr) = head.kt9;
else
   headc(6,1:8) = head.kt9(1:8);
end

nstr = length(head.kf);
if nstr <= 8
   headc(6,9:9+nstr-1) = head.kf;
else
   headc(6,9:16) = head.kf(1:8);
end

nstr = length(head.kuser0);
if nstr <= 8
   headc(6,17:17+nstr-1) = head.kuser0;
else
   headc(6,17:24) = head.kuser0(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.kuser1);
if nstr <= 8
   headc(7,1:nstr) = head.kuser1;
else
   headc(7,1:8) = head.kuser1(1:8);
end

nstr = length(head.kuser2);
if nstr <= 8
   headc(7,9:9+nstr-1) = head.kuser2;
else
   headc(7,9:16) = head.kuser2(1:8);
end

nstr = length(head.kcmpnm);
if nstr <= 8
   headc(7,17:17+nstr-1) = head.kcmpnm;
else
   headc(7,17:24) = head.kcmpnm(1:8);
end

%---------------------------------------------------------------------------
nstr = length(head.knetwk);
if nstr <= 8
   headc(8,1:nstr) = head.knetwk;
else
   headc(8,1:8) = head.knetwk(1:8);
end

nstr = length(head.kdatrd);
if nstr <= 8
   headc(8,9:9+nstr-1) = head.kdatrd;
else
   headc(8,9:16) = head.kdatrd(1:8);
end

nstr = length(head.kinst);
if nstr <= 8
   headc(8,17:17+nstr-1) = head.kinst;
else
   headc(8,17:24) = head.kinst(1:8);
end


%---------------------------------------------------------------------------
npts = headi(2,5);
if npts <= 0
   disp(['Error in npts to write in ',sacfilename]);
   fclose(fid);
   nerr = 5;
   return;
end

%------------- Write header -----------------------------------------------
headf = headf'; headi = headi'; headc = headc';

[chdf] = fwrite(fid, headf, 'float32');
if chdf ~= 14*5
	disp(['Error in writing float header in ',sacfilename]);
	fclose(fid);
	nerr = 1;
	return;
end
[chdi] = fwrite(fid, headi, 'int32');
if chdi ~= 8*5
	disp(['Error in writing integer header in ',sacfilename]);
	fclose(fid);
	nerr = 2;
	return;
end
[chdc] = fwrite(fid, headc, 'char');
if chdc ~= 8*24
	disp(['Error in writing char header in ',sacfilename]);
	fclose(fid);
	nerr = 3;
	return;
end

%------------- Write data -----------------------------------------------
ndata = length(data);
if ndata > npts   
   [cdata] = fwrite(fid, data(1:npts), 'float32');
   if cdata ~= npts
      fprintf('npts = %d, ndata = %d, but write %d data\n', npts, ndata, cdata);
      disp(['Error in writing data in ',sacfilename]);
      fclose(fid);
      nerr = 4;
      return;
   end
elseif ndata < npts
   data(ndata+1:npts) = 0.0;
   [cdata] = fwrite(fid, data, 'float32');
   if cdata ~= length(data)
      fprintf('npts = %d, ndata = %d, but write %d data\n', npts, ndata, cdata);
      disp(['Error in writing data in ',sacfilename]);
      fclose(fid);
      nerr = 4;
      return;
   end
else
   [cdata] = fwrite(fid, data, 'float32');
   if cdata ~= npts
      fprintf('npts = %d, ndata = %d, but write %d data\n', npts, ndata, cdata);
      disp(['Error in writing data in ',sacfilename]);
      fclose(fid);
      nerr = 4;
      return;
   end
end

fclose(fid);
nerr = 1;



