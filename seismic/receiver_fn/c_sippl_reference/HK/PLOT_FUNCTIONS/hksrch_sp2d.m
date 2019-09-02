%
% Plot the H-k search results
% 1. The contour map of stacked amplitudes in H-k domain
% 2. The Receiver Functions with 3 phases arrivals responding to optimal H and k
%

clear;

dirwrk = fullfile('..', 'Work', 'Sta_Files'); 
% reg = 'Reg_CenterChina';
% reg = 'Reg_EastChina';
% reg = 'Reg_NorthEast';
% reg = 'Reg_NorthChina';
% reg = 'Reg_NorthWest';
% reg = 'Reg_SouthChina';
reg = 'Reg_SouthWest';
% reg = 'Regional';
netcode = 'XZ';
stacode = 'LIZ';
vpdir = 'Vp_Tomo';
% vpdir = 'Vp_Crt20';
% vpdir = 'Vp_Crt51';
% vpdir = 'Vp_DSS';
bazdir = 'BAZ_a0.360';
gausdir = 'Gauss_2.5';
psdir = 'RF_PSH5';
dirsacs = fullfile('G:', reg, netcode, stacode, psdir, gausdir);
%dirsacs = fullfile('~', 'Seismograms', reg, netcode, stacode, psdir, gausdir);
ascf = ['sfs2d.', netcode, '.', stacode, '.x1'];
sumf = ['hks2d.', netcode, '.', stacode,'.x1.00'];
% ascf = ['sfr2d.b90x7.', stacode];
% sumf = ['hkr2d.b90x7.', stacode,'.00'];
nmfactor = 0.8;     % 0.5 or 0.8

%---- H-k summary file -----------------
file_sum = fullfile(dirwrk, netcode, stacode, vpdir, bazdir, sumf);

fid = fopen(file_sum, 'r');
while(1)
	tline = fgetl(fid);
	if ~ischar(tline), break, end;
   if ~isempty(strfind(tline, '%Selected average crust Vp'))
      tmpcell = textscan(tline, '%*[%]Selected average crust Vp = %f km/s', 1);
      vp_sel = tmpcell{1};
   end
   if ~isempty(strfind(tline, '%Vp/Vs ratio from '))
      tmpcell = textscan(tline, '%*[%]Vp/Vs ratio from %f to %f at step %f with %d samples', 1);
      rvpvslb = tmpcell{1};
      rvpvsub = tmpcell{2};
      rvpvsinc = tmpcell{3};
      nrvpvs = tmpcell{4};
   end
   if ~isempty(strfind(tline, '%Moho depth from '))
      tmpcell = textscan(tline, '%*[%]Moho depth from %f to %f at step %f with %d samples', 1);
      mhlb = tmpcell{1};
      mhub = tmpcell{2};
      mhinc = tmpcell{3};
      nmoho = tmpcell{4};
   end
   if ~isempty(strfind(tline, '%Vp/Vs ratio = '))
      tmpcell = textscan(tline, '%*[%]Vp/Vs ratio = %f +/- %f', 1);
      rvpvs1 = tmpcell{1};
      sigma_rv = tmpcell{2};
   end
   if ~isempty(strfind(tline, '%Moho depth = '))
      tmpcell = textscan(tline, '%*[%]Moho depth = %f +/- %f km', 1);
      moho1 = tmpcell{1};
      sigma_mh = tmpcell{2};
   end
   if ~isempty(strfind(tline, '%a =  '))
      tmpcell = textscan(tline, '%*[%]a =  %f, b =  %f, alpha = %f', 1);
      a = tmpcell{1};
      b = tmpcell{2};
      magld = tmpcell{3};
   end
end

fseek(fid, 0, -1);
tmpcell = textscan(fid,'%s %f %f %f %f %f %f %f%*[^\n]', 'commentstyle','%');
rfsac_hk = tmpcell{1};
rayp_hk = tmpcell{2};
t1 = tmpcell{3};
s1 = tmpcell{4};
t2 = tmpcell{5};
s2 = tmpcell{6};
t3 = tmpcell{7};
s3 = tmpcell{8};

clear tmpcell;
fclose(fid);

%---- H-k s-function --------------------
% [nrvpvs, rvpvs] = setnvec(rvpvslb, rvpvsub, rvpvsinc);
% [nmoho, moho] = setnvec(mhlb, mhub, mhinc);
rvpvs = linspace(rvpvslb, rvpvsub, nrvpvs);
moho = linspace(mhlb, mhub, nmoho);

[moho_mtr, rvpvs_mtr] = meshgrid(moho, rvpvs);

file_sfun = fullfile(dirwrk, netcode, stacode, vpdir, bazdir, ascf);
fid = fopen(file_sfun, 'r' ,'l');
sfunc_mtr = fread(fid, [nmoho, nrvpvs], 'float32');
sfunc_mtr = sfunc_mtr';
fclose(fid);

%---- 'AXX.rayp.list', need epc, baz and ray parameter etc ---------------
%---- 'AXX.list' is in fact not required anymore -----------
file_list = fullfile(dirwrk, netcode, stacode, [netcode, '.', stacode, '.rayps.list']);
fid = fopen(file_list, 'r');
tmp_cell = textscan(fid, '%s %f %f %f %s %f %f%*[^\n]', 'commentstyle','%');
fclose(fid);
fsacz_rp = tmp_cell{1};
epc_rp = tmp_cell{2};
az_rp = tmp_cell{3};
baz_rp = tmp_cell{4};
ph1st_rp = tmp_cell{5};
travt_rp = tmp_cell{6};
rayp_rp = tmp_cell{7};
clear tmp_cell;
bazr_rp = deg2rad(baz_rp);

fsacz_rpcom = cell(length(fsacz_rp), 1); 
for ii = 1:length(fsacz_rp)
   tmpstr = fliplr(fsacz_rp{ii});
   [sacext, tmpstr] = strtok(tmpstr, '.');
   [comp, tmpstr] = strtok(tmpstr, '.');
   comp = fliplr(comp);
   fsacz_rpcom{ii} = fliplr(tmpstr);
end

%---- 'AXX.sht', need time shift --------------------
file_sht = fullfile(dirwrk, netcode, stacode, [netcode, '.', stacode, '.sht']);
fid = fopen(file_sht, 'r');
tmp_cell = textscan(fid, '%s %f%*[^\n]', 'commentstyle','%');
rfsac_sh = tmp_cell{1};
sht_sh = tmp_cell{2};
clear tmp_cell;
fclose(fid);

rfsac_shcom = cell(length(rfsac_sh), 1);
for ii = 1:length(rfsac_sh)
   tmpstr = fliplr(rfsac_sh{ii});
   [sacext, tmpstr] = strtok(tmpstr, '.');
   [comp, tmpstr] = strtok(tmpstr, '.');
   comp = fliplr(comp);
   rfsac_shcom{ii} = fliplr(tmpstr);
end
rfsac_shcom = rfsac_shcom';

%---- Loop rfsac_hk
%---- Search fsacz_rp for ray epc and baz -----------------
%---- Search rfsac_sh for time shift -------
rfsac_hkcom = cell(length(rfsac_hk), 1);
epc = zeros(length(rfsac_hk), 1);
baz = zeros(length(rfsac_hk), 1); 
bazr = zeros(length(rfsac_hk), 1);
sht = zeros(length(rfsac_hk), 1);

for ii = 1:length(rfsac_hk)
	tmpstr = fliplr(rfsac_hk{ii});
	[sacext, tmpstr] = strtok(tmpstr, '.');
	[comp, tmpstr] = strtok(tmpstr, '.');
	comp = fliplr(comp);
	rfsac_hkcom{ii} = fliplr(tmpstr);

   loce = strcmp(rfsac_hkcom{ii}, fsacz_rpcom);
   if all(~loce)
      error(['No location info for ', rfsac_hkcom{ii}]);
   elseif sum(loce) > 1
      error(['Multiple location info for ', rfsac_hkcom{ii}]);
   end
   epc(ii) = epc_rp(loce);
   baz(ii) = baz_rp(loce);
   bazr(ii) = deg2rad(baz(ii));

   loch = strcmp(rfsac_hkcom{ii}, rfsac_shcom);
   if all(~loch)
      error(['No time shift info for ', rfsac_hk{jj}]);
   elseif sum(loch) > 1
      error(['Multiple time shift info for ', rfsac_hk{jj}]);
   end
   sht(ii) = sht_sh(loch);
end

%---- Locate position of optimal solution in H-k matrix -----------------
[rowii, coljj] = find(abs(rvpvs_mtr - rvpvs1) <= 0.0001 & abs(moho_mtr - moho1) <= 0.001);
if length(rowii) ~= 1 || length(coljj) ~= 1
   error('Error in determining rvpvs_max and moho_max');
end
rvpvs_max = rvpvs_mtr(rowii, coljj);
moho_max = moho_mtr(rowii, coljj);

rvpvs_lb = rvpvs_mtr(1, coljj);
rvpvs_ub = rvpvs_mtr(end, coljj);
moho_lb = moho_mtr(rowii, 1);
moho_ub = moho_mtr(rowii, end);

%---- Take Moho row at optimal rvpvs ---------------
sfun_in1 = sfunc_mtr(rowii, :) >= 0;
moho1 = moho_mtr(1, sfun_in1);
sfunc1 = sfunc_mtr(rowii, sfun_in1);
moho1 = [moho1(1), moho1, moho1(end)];
sfunc1 = [0, sfunc1, 0];

%---- Take rvpvs col at optimal moho ---------------
sfun_in2 = sfunc_mtr(:, coljj) >= 0;
rvpvs2 = rvpvs_mtr(sfun_in2, 1);
sfunc2 = sfunc_mtr(sfun_in2, coljj);
rvpvs2 = [rvpvs2(1); rvpvs2; rvpvs2(end)];
sfunc2 = [0; sfunc2; 0];

%
%dmh = moho - moho_max;
%drv = rvpvs - rvpvs_max;
%
%%theta = [reshape(moho_mtr, nrvpvs*nmoho, 1), reshape(rvpvs_mtr, nrvpvs*nmoho, 1)];
%%H_inv = cov(theta - ones(nrvpvs*nmoho,1)*[moho_max, rvpvs_max]);
%
%%dmh = moho - moho(jmhmax);
%rho_mh = mean(dmh.*dmh);
%%drv = rvpvs - rvpvs(irvmax);
%rho_rv = mean(drv.*drv);
%rho_rvmh = mean(drv.*dmh);
%
%magl = atan2(2*rho_rvmh*sigma_rv*sigma_mh, ...
%   sigma_mh*sigma_mh-sigma_rv*sigma_rv) / 2;
%
%a2 = sigma_mh*sigma_mh*sigma_rv*sigma_rv*(sin(magl)*sin(magl)-cos(magl)*cos(magl)) / ...
%   (sigma_mh*sigma_mh*sin(magl)*sin(magl) - sigma_rv*sigma_rv*cos(magl)*cos(magl));
%a = sqrt(abs(a2));
%b2 =  sigma_mh*sigma_mh*sigma_rv*sigma_rv*(cos(magl)*cos(magl)-sin(magl)*sin(magl)) / ...
%   (sigma_mh*sigma_mh*cos(magl)*cos(magl) - sigma_rv*sigma_rv*sin(magl)*sin(magl));
%b = sqrt(abs(b2));
%magld = magl * 180 / pi;
%
[x, y] = ellipse(a, b, magld, moho_max, rvpvs_max);

%---- Sort RFs in the order of ray parameter -----------------
[rayp_hk, sindex] = sort(rayp_hk);
rfsac_hk = rfsac_hk(sindex);
t1 = t1(sindex);
s1 = s1(sindex);
t2 = t2(sindex);
s2 = s2(sindex);
t3 = t3(sindex);
s3 = s3(sindex);
sht = sht(sindex);

tref = cell(1, length(rfsac_hk));
xref = cell(1, length(rfsac_hk));
pmax = zeros(length(rfsac_hk), 1);
t1loc = zeros(length(rfsac_hk), 1);
t2loc = zeros(length(rfsac_hk), 1);
t3loc = zeros(length(rfsac_hk), 1);
for ii = 1:length(rfsac_hk)
	fsep_loc = strfind(rfsac_hk{ii}, '/');
	if isempty(fsep_loc)
		filesac = fullfile(dirsacs,rfsac_hk{ii});
	else
		qname = rfsac_hk{ii}(1:fsep_loc-1);
		tmpname = rfsac_hk{ii}(fsep_loc+1:end);
		filesac = fullfile(dirsacs, qname, gausdir, tmpname);
	end
	machineformat = sacmachine(filesac);
	[tref{ii}, xref{ii}] = rsac3(filesac, machineformat);
   tref{ii} = shift_rft(tref{ii}, sht(ii));
	pmax(ii) = max(abs(xref{ii}));
% 	xref_x{ii} = (xref{ii}./ pmax(ii) + rayp(ii)*2)./2;
	xref{ii} = (xref{ii} + rayp_hk(ii)*nmfactor) ./ nmfactor;
   [~, t1loc(ii)] = min(abs(tref{ii} - t1(ii)));
   [~, t2loc(ii)] = min(abs(tref{ii} - t2(ii)));
   [~, t3loc(ii)] = min(abs(tref{ii} - t3(ii)));
end


%---- Plot ---------------------------------------------------
figure(13);
clf;
h111 = subplot('Position', [0.1, 0.1, 0.7, 0.7]);  %[left bottom width height]
[c1, ~] = contourf(moho_mtr, rvpvs_mtr, 10.^sfunc_mtr);
% h1 = pcolor(moho_mtr, rvpvs_mtr, 10.^sfunc_mtr);
% shading interp; %faceted; %interp;
%imagesc(moho_mtr(1,:), rvpvs_mtr(:,1), 10.^sfunc_mtr);
%colorbar('SouthOutside');
hold on;
%plot(moho_max, rvpvs_max, 'w.', 'MarkerSize', 8);
%plot(moho_max, rvpvs_max, 'wo', 'MarkerSize', 5, 'MarkerFaceColor','w');
plot(x, y, 'w', 'LineWidth',1);
hln1 = line([moho_lb, moho_ub], [rvpvs_max, rvpvs_max]);
set(hln1, 'Color', 'w');
hln2 = line([moho_max, moho_max], [rvpvs_lb, rvpvs_ub]);
set(hln2, 'Color', 'w');
% xlim([moho_lb, 50]);
xlabel('Moho (km)', 'FontWeight','Bold', 'FontSize', 12);
ylabel('Vp/Vs', 'FontWeight','Bold', 'FontSize', 12);
set(h111, 'FontWeight','Normal', 'FontSize', 12);
%title([stacode, ', Vp=', num2str(vp_sel),'km, Moho=', num2str(moho_max), 'km, Vp/Vs=', num2str(rvpvs_max)]); 

h112 = subplot('Position', [0.1, 0.87, 0.7, 0.1]);  %[left bottom width height]
fill(moho1, sfunc1, 'r', 'EdgeColor','none', 'LineStyle','none');
hold on;
plot(moho_mtr(1,:), sfunc_mtr(rowii, :));
set(h112, 'FontWeight','Normal', 'FontSize', 12);
xlim([moho_lb, moho_ub]);
% xlim([moho_lb, 50]);
ylim([-1.0, 1.0]);

h113 = subplot('Position', [0.87, 0.1, 0.1, 0.7]);  %[left bottom width height]
fill(sfunc2, rvpvs2, 'r',  'EdgeColor','none', 'LineStyle','none');
hold on;
plot(sfunc_mtr(:, coljj), rvpvs_mtr(:, 1));
set(h113, 'FontWeight','Normal', 'FontSize', 12);
xlim([-1.0, 1.0]);
ylim([rvpvs_lb, rvpvs_ub]);

ymin = min(rayp_hk)-0.3;
ymax = max(rayp_hk)+0.8;
figure(14);
clf;
ax = axes;
for ii = 1:length(rfsac_hk)
   %x_in = xref{ii} >= rayp(ii);
   xin1 = xref{ii}; %(x_in);
   tin1 = tref{ii}; %(x_in);
   tin1 = [tin1(1); tin1; tin1(end)];
   xin1 = [rayp_hk(ii); xin1; rayp_hk(ii)];
	fill(tin1, xin1, [0.8,0.8,0.8],  'EdgeColor','none', 'LineStyle','none');
	hold on;
end

for ii = 1:length(rfsac_hk)
   plot(tref{ii}, xref{ii});
   %text(tref{ii}(105), xref{ii}(1), [num2str(rayp(ii))], ...
	%	'HorizontalAlignment','left','VerticalAlignment','bottom',...
	%	'FontSize',7,'Color','k'); %'FontWeight','bold',

   plot(tref{ii}(t1loc(ii)), xref{ii}(t1loc(ii)), 'r.');
	plot(tref{ii}(t2loc(ii)), xref{ii}(t2loc(ii)), 'r.');
	plot(tref{ii}(t3loc(ii)), xref{ii}(t3loc(ii)), 'r.');
end
ylim([ymin ymax]);
xlim([-5 50]);
xlabel('Time (sec)', 'FontWeight','Bold', 'FontSize', 12); 
ylabel('Ray Parameter (sec/deg)', 'FontWeight','Bold', 'FontSize', 12);
%set(ax,'YTickLabel','');
set(ax, 'FontSize',12, 'FontWeight','Normal');
title([stacode, ', Vp=', num2str(vp_sel),'km, Moho=', num2str(moho_max), 'km, Vp/Vs=', num2str(rvpvs_max)]); 

% figure(10);
% clf;
% plot(rayp_hk, pmax', 'o');
% xlabel('Ray Parameter'); %, 'FontSize',8);
% ylabel('Max P amplitude'); %, 'FontSize',8);

figure(9);
clf;
h2 = mmpolar(bazr, epc, '', 'Style','compass', ...
   'RLimit',[0, 100], 'RTickAngle',90, 'RTickValue',[10:10:100], ...
   'RTickLabel',{'','','30^o','','50^o','','70^o','','90^o',''}, ...
   'RGridColor',[0.5 0.5 0.5], 'RTickLabelHalign','right', 'RTickLabelValign','bottom', ...
   'TGridColor',[0.5 0.5 0.5], 'TTickDelta',10, 'TTickSign','+');
set(h2, 'Color','b', 'LineStyle','none', 'Marker','.', 'MarkerSize',8);
hold on;
text(0, 0, stacode,...
   'HorizontalAlignment','left', 'VerticalAlignment','bottom',...
   'FontSize',8, 'Color','k', 'FontWeight','bold', 'BackgroundColor', 'w');
h1 = polar(0, 0);
set(h1, 'LineStyle','none', 'Marker','p', 'Color','r', 'MarkerFaceColor', 'y', 'MarkerSize', 12);

