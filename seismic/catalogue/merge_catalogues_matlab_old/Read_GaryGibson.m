%% Reads different catalogues and produces ascii files in a unified format
% =========================================================================
% Babak Hejrani, Geoscience Australia, 22 Oct 2018
%
% This program reads csv file containg the Gary Gibson catalogue and wties output ascii files
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
folds = {'GaryGibson'};
% format of the files in each folder: ascii or xml
forms = {'/*.csv'};

% For loop over events
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
    
    fid = fopen(file,'r');
    % s       s    f   s          f    f     f  f    f      f      s      f      f          f         f     f     s  f    f      f    f    f      f      f        f        f     f   f    f         f       s              s                s           s
    % 1       2    3   4          5    6     7  8    9      10     11     12     13         14        15    16    17 18   19     20   21   22     23     24       25       26    27  28   29        30      31             32               33          34
    %RunAuth,Type ,Ind,Ass      ,Year,Month,Day,Hour,Minute,Second,STCode,STCorr,LongitudeE,LatitudeS,Depth,ZCode,Mx,Mval,TimUnc,Smaj,Smin,SmajAz,DepUnc,Arrivals,SDResids,Sites,Gap,Gap2,NearestSG,Nearest,HDPlace       ,Origin (text)   ,YearF      ,Magnitude text
    %ADE    ,local,1  ,Preferred,1840,3    ,31 ,06  ,30    ,      ,UTC   ,0     ,138.6     ,34.9     ,7    ,     ,MP,2.7 ,      ,    ,    ,      ,      ,        ,        ,     ,   ,    ,         ,       ,"Adelaide, SA",1840-03-31  0630,1840.246642,MP 2.7 ADE
    %RunAuth,Type ,Ind,Ass      ,Year,Month,Day,Hour,Minute,Second,STCode,STCorr,LongitudeE,LatitudeS,Depth,ZCode,Mx,Mval,TimUnc,Smaj,Smin,SmajAz,DepUnc,Arrivals,SDResids,Sites,Gap,Gap2,NearestSG,Nearest,HDPlace                                    ,Origin (text)   ,YearF      ,Magnitude text
    %AUST   ,local,1  ,Preferred,2017,7    ,15 ,21  ,27    ,26    ,ACST  ,9.5   ,131.897   ,26.122   ,0    ,?    ,ML,3.2 ,2.7   ,1   ,1   ,2     ,22    ,13      ,3.87    ,13   ,79 ,    ,WRKA     ,3.429  ,"30 km NW of Ernabella, Musgrave Range, SA",2017-07-15  2127,2017.536696,ML 3.2 AUSTMs 3.5 AUST
    header = fgetl(fid);
    while ~feof(fid)
        lin = fgetl(fid);
        nuev = nuev + 1;
        ind = strfind(lin,',');
        YYtest  = str2num(lin(ind(4)+1:ind(5)-1));
        if isempty(YYtest)
            YY(nuev,1) = -999;
        else
            YY(nuev,1) = YYtest;
        end
        MMtest  = str2num(lin(ind(5)+1:ind(6)-1));
        if isempty(MMtest)
            MM(nuev,1) = -999;
        else
            MM(nuev,1) = MMtest;
        end
        DDtest  = str2num(lin(ind(6)+1:ind(7)-1));
        if isempty(DDtest)
            DD(nuev,1) = -999;
        else
            DD(nuev,1) = DDtest;
        end
        hhtest  = str2num(lin(ind(7)+1:ind(8)-1));
        if isempty(hhtest)
            hh(nuev,1) = -999;
        else
            hh(nuev,1) = hhtest;
        end
        mmtest  = str2num(lin(ind(8)+1:ind(9)-1));
        if isempty(mmtest)
            mm(nuev,1) = -999;
        else
            mm(nuev,1) = mmtest;
        end
        sstest      = str2num(lin(ind(9)+1:ind(10)-1));
        if isempty(sstest)
            ss(nuev,1) = -999;
        else
            ss(nuev,1) = sstest;
        end
        lo(nuev,1)  = str2num(lin(ind(12)+1:ind(13)-1));
        la(nuev,1)  = str2num(lin(ind(13)+1:ind(14)-1)); % This is a positive number! The southern hemisphere is usually negative latitude. 
        detest      = str2num(lin(ind(14)+1:ind(15)-1));
        if isempty(detest)
            de(nuev,1) = -999;
        else
            de(nuev,1) = detest;
        end
        matest      = str2num(lin(ind(17)+1:ind(18)-1));
        if isempty(matest)
            ma(nuev,1) = -999;
        else
            ma(nuev,1) = matest;
        end
        mtype(nuev,1)  = {lin(ind(16)+1:ind(17)-1)}; % ML, Ms, mb, Mw
        
        Nphtest     = str2num(lin(ind(23)+1:ind(24)-1));
        if isempty(Nphtest)
            Nph(nuev,1) = -999;
        else
            Nph(nuev,1) = Nphtest;
        end
        
        gaptest     = str2num(lin(ind(26)+1:ind(27)-1));
        if isempty(gaptest)
            gap(nuev,1) = -999;
        else
            gap(nuev,1) = gaptest;
        end
        % This is special csv: ","
        %s = textscan(lin,'%s %s %f %s %f %f %f %f %f %f %s %f %f %f %f','Delimiter',',');
    end
    fclose(fid);
    
    % Indexes
    indml = find(strcmpi(mtype,'ML'));
    indmb = find(strcmpi(mtype,'mb'));
    indms = find(strcmpi(mtype,'Ms'));
    indmw = find(strcmpi(mtype,'Mw'));
    
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
    ml(indml) = ma(indml);
    mb(indmb) = ma(indmb);
    ms(indms) = ma(indms);
    mw(indmw) = ma(indmw);
    
    % Gary Gibson catalogue is from 1857. We should cut it from
    % 20?? - ...
    % The analyst data start from 2012-06-21T05:18:40.822000Z
    %otdate = datenum([YY MM DD hh mm ss]);
    %A      = [otdate lo la de Nph ma];
    %[B,I]  = sortrows(A); % Gray's data set is sorted, This is not neccessary
    %CutFromHere = datenum([2000 01 01 00 00 00]); % datenum([2012 06 21 05 01 00.00]);
    %ind    = find(B(:,1)> CutFromHere);
    %GG.ev  = B(ind,:);
    %          1   2  3  4  5  6  7  8    9  10  11 12 13 14 15                  16
    GG.ev   = [YY MM DD hh mm ss lo -1*la de Nph mb ms ml mw -999*ones(size(YY)) gap];
    save GG GG;
    
    % Write the output with the following format:
    % #Reporter,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,Nph ,mb, ms, ml, ID
    % ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
    fid = fopen('GG.csv','w+');
    for i = 1:length(YY)
        %            Rep   %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID
        fprintf(fid,'#GG , %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f, -999, %04.2f\n',...
            YY(i),  MM(i),  DD(i),  hh(i),  mm(i),  ss(i),   lo(i),   la(i),   de(i),  Nph(i),  mb(i), ms(i),  ml(i),  mw(i),       gap(i));
    end
    fclose(fid);
    clear GG ev picks YY MM DD hh mm ss lo la de Nph ma ml mb ms mw gap
end
toc