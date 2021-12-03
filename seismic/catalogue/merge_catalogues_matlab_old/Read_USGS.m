%% Reads different catalogues and produces ascii files in a unified format
% =========================================================================
% Babak Hejrani, Geoscience Australia, January 2019
%
% This program reads GA's data and wties output ascii files
% with the following format (suggested by Alexei Gorbatov):
% #EV-ID Time Lon Lat Dep Mag Nph
% ST-ID Time Phase
% Update 28 Oct 2018: I prefer a new format, described below ...
% #EV ,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,ph ,ma
% ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
%
% =========================================================================

clear all; close all; clc; if ispc; sep='\'; else sep='/'; end % This code will work both in Linux and windows
comment = 0; s1 = datenum(2000,1,1,0,0,0)-datenum(2000,1,1,0,0,1); % One second in datenum format
% The top folder where all catalogues are located
tf = '/g/data1a/ha3/Passive/Events'; % This is top folder that all the data is located, This is my computer at ANU ...
% The location of different catalogues, folder names under the top folder
folds = {'USGS'};
% format of the files in each folder: ascii or xml
forms = {'/*.csv'};

% Take at on folder
folder = [tf,sep,folds{1},sep];
disp(['Working on folder: ',folder]);
disp([' ']);
d = dir([folder,forms{1}]); % Matlab 2016b and later
disp(['--> Number of files for: ',forms{1},' is: ',num2str(length(d))]);
if length(d) < 1
    disp(['--> No file with that format in this folder!'])
    keyboard
end
nuev=0; % number of events loaded so far ...
tic
for id =1:length(d) % Loop over files ...
    %file = [char(d(id).folder),sep,char(d(id).name)]; %NCI, MATLAB 2016b and later
    file = [folder,sep,char(d(id).name)]; % ANU Matlab 2014a
    
    %   1                         2        3       4    5     6     7   8    9   10   11  12         13                        14                             15          16              17         18        19    20    21              22
    %time                    ,latitude,longitude,depth,mag,magType,nst,gap, dmin,rms ,net,id        ,updated                 ,place                          ,type      ,horizontalError,depthError,magError,magNst,status,locationSource,magSource
    %2016-05-01T18:56:05.860Z,-20.7528,169.9387 ,62.56,5.1,mb     ,   ,111,1.923,1   ,us ,us10005cvg,2016-07-19T20:03:08.040Z,"150km SSE of Isangel, Vanuatu",earthquake,8.3            ,5.3       ,0.07    ,66,    reviewed,us,          us
    %2018-12-31T10:49:41.340Z,-31.7676,-69.3073 ,101.2,5.3,mww    ,   ,30 ,0.575,1.36,us ,us2000izex,2019-01-16T17:22:11.813Z,"49km S of Calingasta, Argentina",earthquake,6.4,4,0.055,32,reviewed,us,us
    %2018-12-31T14:10:25.200Z,37.4856 ,141.4451 ,43.26,5.1,mww    ,   ,106,0.645,0.67,us ,us2000izfu,2018-12-31T17:15:56.594Z,"39km E of Namie, Japan",earthquake,3.2,5.5,0.073,18,reviewed,us,us
    %2018-12-31T20:22:15.900Z,-17.5085,-174.67  ,170.91,5.1,mb    ,   ,51 ,4.534,1.38,us ,us2000izi1,2018-12-31T20:40:16.040Z,"145km NNW of Neiafu, Tonga",earthquake,7.8,6.9,0.04,200,reviewed,us,us
    %2018-12-31T20:45:31.020Z,15.2791 ,-46.3283 ,10   ,5  ,mb     ,   ,74 ,11.861,1.09,us,us2000izi5,2018-12-31T22:06:38.040Z,"Northern Mid-Atlantic Ridge",earthquake,10.4,1.9,0.039,209,reviewed,us,us
    
    fid = fopen(file,'r');
    header = fgetl(fid);
    while ~feof(fid)
        lin = fgetl(fid);
        nuev = nuev + 1;
        ind = strfind(lin,',');
        YYtest  = str2num(lin(1:4));
        if isempty(YYtest)
            YY(nuev,1) = -999;
        else
            YY(nuev,1) = YYtest;
        end
        MMtest  = str2num(lin(6:7));
        if isempty(MMtest)
            MM(nuev,1) = -999;
        else
            MM(nuev,1) = MMtest;
        end
        DDtest  = str2num(lin(9:10));
        if isempty(DDtest)
            DD(nuev,1) = -999;
        else
            DD(nuev,1) = DDtest;
        end
        hhtest  = str2num(lin(12:13));
        if isempty(hhtest)
            hh(nuev,1) = -999;
        else
            hh(nuev,1) = hhtest;
        end
        mmtest  = str2num(lin(15:16));
        if isempty(mmtest)
            mm(nuev,1) = -999;
        else
            mm(nuev,1) = mmtest;
        end
        sstest      = str2num(lin(18:23));
        if isempty(sstest)
            ss(nuev,1) = -999;
        else
            ss(nuev,1) = sstest;
        end
        la(nuev,1)  = str2num(lin(ind(1)+1:ind(2)-1));
        lo(nuev,1)  = str2num(lin(ind(2)+1:ind(3)-1)); 
        detest      = str2num(lin(ind(3)+1:ind(4)-1));
        if isempty(detest)
            de(nuev,1) = -999;
        else
            de(nuev,1) = detest;
        end
        matest      = str2num(lin(ind(4)+1:ind(5)-1));
        if isempty(matest)
            ma(nuev,1) = -999;
        else
            ma(nuev,1) = matest;
        end
        mtype(nuev,1)  = {lin(ind(5)+1:ind(6)-1)}; % ML, Ms, mb, Mw
        
        Nphtest     = str2num(lin(ind(6)+1:ind(7)-1));
        if isempty(Nphtest)
            Nph(nuev,1) = -999;
        else
            Nph(nuev,1) = Nphtest;
        end
        
        gaptest     = str2num(lin(ind(7)+1:ind(8)-1));
        if isempty(gaptest)
            gap(nuev,1) = -999;
        else
            gap(nuev,1) = gaptest;
        end
        
        ID(nuev,1)  = {lin(ind(11)+1:ind(12)-1)};
    end
    fclose(fid);
    
    % Indexes
    indml  = find(strcmpi(mtype,'ml'));
    indmb  = find(strcmpi(mtype,'mb'));
    indms  = find(strcmpi(mtype,'Ms'));
    indmw  = find(strcmpi(mtype,'mw'));
    indmww = find(strcmpi(mtype,'mww'));
    
    % Creates a variable of NaN or a cell of '', the later does not
    % work because we want to build the GG structure as a matrix
    % containing numbers ...
    %ml=NaN(size(mtype));
    %mb=NaN(size(mtype));
    %ms=NaN(size(mtype));
    %mw=NaN(size(mtype));
    % Creates a variables of -999
    ml = repmat(-999,size(mtype));
    mb = repmat(-999,size(mtype));
    ms = repmat(-999,size(mtype));
    mw = repmat(-999,size(mtype));
    
    % Re-fill the ml, mb, ms and mw using Indexes
    ml(indml)  = ma(indml);
    mb(indmb)  = ma(indmb);
    ms(indms)  = ma(indms);
    mw(indmw)  = ma(indmw);
    mw(indmww) = ma(indmww);
    
    % Gary Gibson catalogue is from 1857. We should cut it from
    % 20?? - ...
    % The analyst data start from 2012-06-21T05:18:40.822000Z
    %otdate = datenum([YY MM DD hh mm ss]);
    %A      = [otdate lo la de Nph ma];
    %[B,I]  = sortrows(A); % Gray's data set is sorted, This is not neccessary
    %CutFromHere = datenum([2000 01 01 00 00 00]); % datenum([2012 06 21 05 01 00.00]);
    %ind    = find(B(:,1)> CutFromHere);
    %GG.ev  = B(ind,:);
    %             1  2  3  4  5  6  7  8  9  10 11 12 13 14 15                  16
    USGS.ev   = [YY MM DD hh mm ss lo la de Nph mb ms ml mw -999*ones(size(YY)) gap];
    save USGS USGS;
    
    % Write the output with the following format:
    % #Reporter,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,Nph ,mb, ms, ml, ID
    % ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
    fid = fopen('USGS.csv','w+');
    for i = 1:length(YY)
        %            Rep     %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb     ms      ml      mw      ID     gap
        fprintf(fid,'#USGS , %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f,  %s,    %04.2f\n',...
                             YY(i),   MM(i),  DD(i),  hh(i),  mm(i),  ss(i),   lo(i),   la(i),   de(i),  Nph(i),  mb(i), ms(i),  ml(i),  mw(i),  char(ID(i)), gap(i));
    end
    fclose(fid);
    clear GG ev picks YY MM DD hh mm ss lo la de Nph ma ml mb ms mw gap
end
toc