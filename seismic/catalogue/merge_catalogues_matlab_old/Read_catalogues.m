%% Reads different catalogues and produces ascii files in a unified format
% =========================================================================
% Babak Hejrani, Geoscience Australia, 22 Oct 2018
%
% This program reads ascii and xml files and wties output ascii files 
% with the following format (suggested by Alexei Gorbatov):
% #EV-ID Time Lon Lat Dep Mag Nph
% ST-ID Time Phase
% Update 28 Oct 2018: I prefer a new format, described below ... 
% #EV ,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,ph ,ma
% ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
%
% Some weird cases in the GA data: '1088942.xml, Nph 56, but only 14
% readings! File: 1102369.xml seems to be in an inifinte loop!
% 
% Some weird cases in EHB: 202.xml has no picks! but he Nph is 607 :)
% 
% At the end all files will be merged together to form a new catalogue ... 
% 
% =========================================================================

clear all; close all; clc; if ispc; sep='\'; else sep='/'; end % This code will work both in Linux and windows
comment = 0; s1 = datenum(2000,1,1,0,0,0)-datenum(2000,1,1,0,0,1); % One second in datenum format
% The top folder where all catalogues are located
tf = 'C:\Users\u87316\OneDrive - Geoscience Australia\Babak\GA\BodyWaveTomo\Catalogue\'; % This is top folder that all the data is located, This is my computer at ANU ... 
% tf = '/home/babak/mega/GA/BodyWaveTomo/Catalogue/';
% tf = '/home/babak/GA/';
% The location of different catalogues, folder names under the top folder
folds = {'engdahl-events'...
    'analyst-reviewed-Alexei-csv'...
    'GaryGibson' ...
    'ISC-Arrivals'};
% format of the files in each folder: ascii or xml
forms = {'/*.xml' ...
    '/*.txt' ...
    '/*.csv' ...
    '/Arrivals_*'};

phnumisccsv = 0; % number of ISC events loaded so far ...

% For loop over events
for ifo =  2 % [2 1 4 3] % 1:length(folds)
    % Take at on folder
    folder = [tf,sep,folds{ifo},sep];
    disp(['Working on folder: ',folder]);
    disp([' ']);
    txt =  forms{ifo};
    d = dir([folder,txt]); % Matlab 2016b and later
    disp(['--> Number of files for: ',txt,' is: ',num2str(length(d))]);
    if length(d) < 1
        disp(['--> No file with that format in this folder!'])
        keyboard
    end
    nuev=0; % number of events loaded so far ...
    tic
    for id =1:length(d) % Loop over files ...
        %file = [char(d(id).folder),sep,char(d(id).name)]; %NCI, MATLAB 2016b and later
        file = [folder,sep,char(d(id).name)]; % ANU Matlab 2014a
        
        if strcmp(folds{ifo},'GaryGibson') % READ CSV
            
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
                la(nuev,1)  = str2num(lin(ind(13)+1:ind(14)-1));
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
            %          1   2  3  4  5  6  7  8  9  10 11 12 13 14 15                  16
            GG.ev   = [YY MM DD hh mm ss lo la de Nph mb ms ml mw -999*ones(size(YY)) gap]; 
            save -v7.3 GG73 GG; % GG has no event-id
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
                        
        elseif strcmp(folds{ifo},'analyst-reviewed-Alexei-csv')
            
            fid = fopen(file,'r');
            header = fgetl(fid);
            %{
            % Alexei's version one:
            %                 1  2  3  4  5  6  7   8   9   10  11  12  13 14 15 16 17 18 19 20
            ts=textscan(fid,'%f %f %s %s %s %f  %f  %f  %f  %f  %f  %f  %s %s %s %s %s %s %f %f');
            %   1         2           3      4        5         6          7        8       9      10      11     12   13     14    15          16      17      18           19        20
            %  evid     orid        Date(M/DD/YYYY) Time       lon         lat      depth   nass    mb      ms     ml   snet   sta  phase    arrival.time                   timeres   delta
            %  285245   400639    7/02/2010 (183)  4:59:51.297  139.3222   -4.4598   49.0825   13    5.34    3.87 -999.00 PS    JAY    P       7/02/2010 (183)  5:00:26.911   -1.023    2.384
            %}
            % Alexei's second version:
            %   1      2           3          4        5          6      7       8       9    10   11     12   13   14      15      16     17  18          19     20   21              22      23          24    25    26     27
            % evid   prefor  Date(M/DD/YYYY)  Day    time        lon     lat     depth   nass ndef review mb   ms   ml      timeres vmodel net arrival.sta iphase chan Date(M/DD/YYYY) Day    Time        snr   delta seaz    esaz
            % 285516 323023  7/18/2010       (199) 13:04:15.684 150.5645 -6.2684 74.4467 94   94   p      6.66 7.10 -999.00 -0.446  iasp91 IM  PMG         P      BHZ  7/18/2010       (199) 13:05:22.444 11221 4.613 47.33  226.87
            % evid   prefor     date          day    time        lon     lat     depth   nass ndef review mb   ms   ml      timeres vmodel net arrival.sta iphase chan arrival.time     daye  time        snr   delta  seaz   esaz
            %                 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
            ts=textscan(fid,'%f %f %s %s %s %f %f %f %f %f %s %f %f %f %f %s %s %s %s %s %s %s %s %f %f %f %f');
            fclose(fid);
            %{
            %Fast: Date Format: (M/DD/YYYY)
            % Strfind works, but then we need a for loop over the strfind outputs which is a structure
            %date1 = char(ts{3}); dind = strfind(date1,'/'); 
            %time1 = char(ts{5}); tind = strfind(time1,':');
            %YY = str2num(date1(:,end-3:end)); % str2num(date1(:,end-3:end));
            %MM = str2num(date1(:,end-6:end-5));
            %DD = str2num(date1(:,1:end-8));
            %hh = str2num(time1(:,end-5:end));
            %mm = str2num(time1(:,end-8:end-7));
            %ss = str2num(time1(:,1:end-10));
            %ot = datenum(YY,MM,DD,hh,mm,ss);
            %}
            % Slow
            %evymd = ts{3};
            %evtim = ts{5};
            space = repmat(' ',size(ts{3}));
            tic; ot = datenum([char(ts{3}),space,char(ts{5})]);
            [YY MM DD hh mm ss] = datevec(ot); toc;           
            num  = find(diff(ot));
            endp = [num;length(ot)]; % End of each event
            begp = [1;num+1]; % Begin of each event
            % Replace the magnitudes of -999.00 by NaN
            %mb   = ts{10}(begp); %indmb = find(mb == -999.00); mb(indmb)=NaN;
            %ms   = ts{11}(begp); %indms = find(ms == -999.00); ms(indms)=NaN;
            %ml   = ts{12}(begp); %indml = find(ml == -999.00); ml(indml)=NaN;
            mw = repmat(-999,size(ts{12}(begp)));
            %      %YY         MM      DD       hh        mm      ss       lo             la         de          ph       mb           ms           ml           mw ID
            %       1          2        3        4        5        6        7             8          9           10       11           12           13           14 15    
            GA.ev= [YY(begp) MM(begp) DD(begp) hh(begp) mm(begp) ss(begp) ts{6}(begp) ts{7}(begp) ts{8}(begp) ts{9}(begp) ts{12}(begp) ts{13}(begp) ts{14}(begp) mw ts{1}(begp)];
            %{
            %Fast
            %date1 = char(ts{16}); %dind = strfind(date1,'/'); % Strfind works, but then we need a for loop (slow)
            %time1 = char(ts{18}); %tind = strfind(time1,':');
            %tic;
            %YY = str2num(date1(:,end-3:end)); % str2num(date1(:,end-3:end));
            %MM = str2num(date1(:,end-6:end-5));
            %DD = str2num(date1(:,1:end-8));
            %hh = str2num(time1(:,end-5:end));
            %mm = str2num(time1(:,end-8:end-7));
            %ss = str2num(time1(:,1:end-10));
            %at = datenum(YY,MM,DD,hh,mm,ss);
            %toc;
            %}
            %Slow
            %stymd = ts{16};
            %sttim = ts{18};
            tic; at=datenum([char(ts{21}),space,char(ts{23})]);
            [YY MM DD hh mm ss] = datevec(at); toc;
            %st=ts{14};
            %ne=ts{13};
            %ph=ts{15};
            for iev = 1:length(begp)
                ind = begp(iev) : endp(iev);
                GA.picks(iev).at = [YY(ind) MM(ind) DD(ind) hh(ind) mm(ind) ss(ind)];
                GA.picks(iev).dist = ts{25}(ind);  
                GA.picks(iev).ch = char(ts{20}(ind));
                GA.picks(iev).ne = char(ts{17}(ind));
                GA.picks(iev).st = char(ts{18}(ind));
                GA.picks(iev).ph = char(ts{19}(ind));               
                g    = ts{26}(ind); % This is baz for this event
                g = sort(g); % sort them,
                g = [g; min(g)+360];% add 360 to the minimum
                diffg = max(diff(g)); % Take the maximum of differences, this is the maximum azimutham gap
                gap(iev,1)    = diffg;
            end
            GA.ev = [GA.ev, gap]; clear gap g diffg;
            % save -v7.3 GA73 GA; % Save the GA in GA.mat format ...
            save GA GA; clear ev picks;
            % #EV ,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,Nph ,ma
            % ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
            fid=fopen('GA.csv','w+');
            for i=1:length(GA.ev(:,1))
                % Write the event info
                %                  1       2       3        4      5       6        7       8        9        10       11    12      13      14       15     16
                %            Rep  %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID     Gap
                fprintf(fid,'#GA , %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f, %7.0f, %7.4f\n',...
                                  GA.ev(i,1:16)                                                                                            );
                for j=1:length(GA.picks(i).at(:,1))
                    fprintf(fid,'%s, ',GA.picks(i).st(j,:)); % stations
                    fprintf(fid,'%s, ',' '); % channel
                    fprintf(fid,'%s, ',' '); % location code 00, 10, ... 
                    fprintf(fid,'%s, ',GA.picks(i).ne(j,:)); % network
                    fprintf(fid,'%s, ',' '); % lon
                    fprintf(fid,'%s, ',' '); % lat
                    fprintf(fid,'%s, ',' '); % ele
                    fprintf(fid,'%s, ',GA.picks(i).ph(j,:)); % phase
                    fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',GA.picks(i).at(j,:)); % 
                    fprintf(fid,'%07.4f \n',GA.picks(i).dist(j));
                end
            end
            fclose(fid);
            clear GA YY MM DD hh mm ss ts mb ms ml at ot
            
        elseif strcmp(folds{ifo},'engdahl-events') % READ XMLs
            
            % Try 1:
            %{
            % try 1:
            fid = fopen(file,'r');
            n21 = 1;
            cl = 0;
            num= 0;
            while ~feof(fid)
                cl=cl+1;
                lin=fgetl(fid);
                if lin21 == 1
                    if length(lin) < 21
                        continue
                    end
                    test1 = strcmp(lin(1:21),'    <origin publicID=');
                    if test1
                        ev_info_collected = 0;
                        while ev_info_collected < 5
                            str = fgetl(fid); cl=cl+1;
                            otfound = strcmp(str,'      <time>'); % <time>
                            lafound = strcmp(str,'      <latitude>'); % <time>
                            lofound = strcmp(str,'      <longitude>'); % <time>
                            defound = strcmp(str,'      <depth>'); % <time>
                            if length(str) > 29
                                Nphfound = strcmp(str(1:30),'        <associatedPhaseCount>'); % <time>
                            else
                                Nphfound = 0;
                            end
                            if otfound
                                str = fgetl(fid); cl=cl+1; % <value>2005-09-16T09:00:51.000150Z</value>
                                larind = strfind(str,'>'); smaind = strfind(str,'<');
                                otstr  = str(larind(1)+1:smaind(2)-1);
                                otdate = datenum([otstr(1:10),' ',otstr(12:end-1)]);
                                ev_info_collected = ev_info_collected + 1;
                            end
                            if lafound
                                str = fgetl(fid); cl=cl+1; % <value>2005-09-16T09:00:51.000150Z</value>
                                larind = strfind(str,'>'); smaind = strfind(str,'<');
                                la = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                            end
                            if lofound
                                str = fgetl(fid); cl=cl+1; % <value>2005-09-16T09:00:51.000150Z</value>
                                larind = strfind(str,'>'); smaind = strfind(str,'<');
                                lo = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                            end
                            if defound
                                str = fgetl(fid); cl=cl+1; % <value>2005-09-16T09:00:51.000150Z</value>
                                larind = strfind(str,'>'); smaind = strfind(str,'<');
                                de = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                            end
                            if Nphfound
                                larind = strfind(str,'>'); smaind = strfind(str,'<');
                                Nph = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                            end
                        end
                        if ev_info_collected == 5
                            lin21 = 0;  % We have the lon lat dep ot ... Search for magnitude
                        end
                    else
                        continue
                    end
                end
                if length(lin) < 26
                    continue
                end
                test2 = strcmp(lin(1:26),'      <magnitude publicID=');
                if test2
                    num = num +1;
                    str = fgetl(fid); cl=cl+1; % <magnitude>
                    str = fgetl(fid); cl=cl+1; % <value
                    larind = strfind(str,'>'); smaind = strfind(str,'<');
                    ma(num)  = str2num(str(larind(1)+1:smaind(2)-1));
                else
                    continue
                end
            end
            fclose(fid);
            % If you reach here, you have found an event with magnitude
            if exist('otdate','var') == 1
                nuev = nuev + 1; % event counter
                ev(nuev,:) = [otdate lo la de Nph]; %
                if exist('ma','var') == 1
                    evm(nuev)  = {ma};
                else
                    evm(nuev)  = {NaN};
                end
            else
                disp(['Could not found origin time ... '])
                keyboard
            end
            clear otdate lo la de Nph ma
            %}
            
            % Try 2:
            if comment; disp(['-->      Working on file: ',char(d(id).name)]); end
            fid = fopen(file,'r'); % open the file
            ev_info_collected = 0;
            ev_distance_collected = 0;
            while ~feof(fid)
                str=fgetl(fid);
                Nphr = 0; % Number of phases read so far ...
                % Load the picks:
                while isempty(strfind(str,'<origin publicID=')) % We will not read all lines down to <event>, we will exit the loop if we reach <magnitude>
                    
                    while ~isempty(strfind(str,'<pick publicID'))
                        % This is a pick
                        while isempty(strfind(str,'</pick>'))
                            if ~isempty(strfind(str,'<time>'))
                                str=fgetl(fid); % This is <value>
                                if isempty(strfind(str,'<uncertainty>')) %
                                    % So this is value, go on ...
                                else
                                    % So this is <uncertainty>, read again
                                    str=fgetl(fid);
                                end
                                Nphr        = Nphr + 1; % One time is read
                                larind      = strfind(str,'>'); smaind = strfind(str,'<');
                                atstr       = str(larind(1)+1:smaind(2)-1);
                                at(Nphr,1)  = datenum([atstr(1:10),' ',atstr(12:end-1)]); % Phase arrival time
                            elseif ~isempty(strfind(str,'<waveformID'))
                                cind        = strfind(str,'"');
                                ne(Nphr,1)  = {str(cind(1)+1:cind(2)-1)};
                                st(Nphr,1)  = {str(cind(3)+1:cind(4)-1)};
                                ch(Nphr,1)  = {str(cind(5)+1:cind(6)-1)};
                            elseif ~isempty(strfind(str,'<phaseHint'))
                                phtest = strfind(str,'<phaseHint/>');
                                if isempty(phtest)
                                    larind      = strfind(str,'>'); smaind = strfind(str,'<');
                                    ph(Nphr,1)  = {str(larind(1)+1:smaind(2)-1)}; % Phase arrival time
                                else
                                    ph(Nphr,1)  = {'-'};
                                end
                                
                            end
                            str = fgetl(fid);
                        end
                        %{
                        %            <pick publicID="quakeml:ga.ga.gov.au/pick/4311408">
                        % 				<creationInfo>
                        % 					<creationTime>2018-05-22T21:33:56.186Z</creationTime>
                        % 					<author>dbp:lsaikal:181</author>
                        % 				</creationInfo>
                        % 				<backazimuth>
                        % 					<value>330.03</value>
                        % 				</backazimuth>
                        % 				<time>
                        % 					<value>2018-05-22T06:31:57.768Z</value>
                        % 				</time>
                        % 				<waveformID channelCode="BHE" locationCode="" networkCode="AU" stationCode="MULG">quakeml:ga.ga.gov.au/waveform/1270</waveformID>
                        % 				<evaluationStatus>preliminary</evaluationStatus>
                        % 				<phaseHint>S</phaseHint>
                        % 				<evaluationMode>automatic</evaluationMode>
                        % 			</pick>
                        %}
                    end
                    str=fgetl(fid);
                end
                
                % While1: Keep reading until you reach </origin>
                while isempty(strfind(str,'</origin>'))
                    
                    str=fgetl(fid);
                    if ~isempty(strfind(str,'<longitude>')) % length(lin) >= 27
                        sv = 1; % To start the saerch for the <value>
                        while sv == 1 % Enter the while loop
                            str = fgetl(fid); % read a line
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<'); % If yes, then extract the value
                                lo       = str2num(str(larind(1)+1:smaind(2)-1)); % store the value into longitude
                                ev_info_collected = ev_info_collected + 1; % one information is loaded
                                if comment; disp('---        Found lo'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<depth>')) % length(lin) >= 23
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                de  = str2num(str(larind(1)+1:smaind(2)-1)); % Depth in km. 
                                ev_info_collected = ev_info_collected + 1;
                                sv=0; % You found the value, get out
                                if comment; disp('---        Found dep'); end
                            end
                        end
                    elseif ~isempty(strfind(str,'<time>')) % length(lin) >= 22
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                otstr  = str(larind(1)+1:smaind(2)-1);
                                otdate = datenum([otstr(1:10),' ',otstr(12:end-1)]); % convert the origin time into datenum
                                ev_info_collected = ev_info_collected + 1;
                                if comment; disp('---        Found ot'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<latitude>')) % length(lin) >= 26
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                la  = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                                if comment; disp('---        Found la'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<associatedPhaseCount>')) 
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        Nph  = str2num(str(larind(1)+1:smaind(2)-1));
                        ev_info_collected = ev_info_collected + 1;
                        if comment; disp('---        Found Nph'); end
                        
                    elseif ~isempty(strfind(str,'<arrival>')) 
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<distance>')) % Is it <distance>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                ev_distance_collected = ev_distance_collected + 1;
                                evst_dist(ev_distance_collected,1) = str2num(str(larind(1)+1:smaind(2)-1)); % Depth in km
                                sv=0; % You found the <distance>, get out
                                if comment; disp('---        Found distance'); end
                            end
                        end
                    elseif ~isempty(strfind(str,'<magnitude publicID'))
                        while isempty(strfind(str,'<value>'))
                            str = fgetl(fid);
                        end
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        ma  = str2num(str(larind(1)+1:smaind(2)-1));
                        while isempty(strfind(str,'<type>'))
                            str = fgetl(fid);
                        end
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        mt  = str(larind(1)+1:smaind(2)-1);
                        ev_info_collected = ev_info_collected + 1; % This is 6
                        if comment; disp('---        Found Mag'); end
                        
                        while isempty(strfind(str,'<azimuthalGap>'))
                            str = fgetl(fid);
                        end
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        gap  = str2num(str(larind(1)+1:smaind(2)-1));
                        ev_info_collected = ev_info_collected + 1; % This is 7 
                        break
                    end
                end
                
                % Did we get the magnitude?0
                if ev_info_collected == 6 % Mag but no azimuthal coverage
                    % You found the magnitude, the str is '<value>1.2...'
                    nuev = nuev + 1; % event counter
                    ev(nuev,:) = [otdate lo la de Nph -999];
                    mvalu(nuev,1) = ma;
                    mtype(nuev,1) = {mt};
                    ID(nuev,1)    = str2num(char(d(id).name(1:end-4)));
                    %clear otdate lo la de Nph ma
                    ev_info_collected = 0;
                elseif ev_info_collected == 5 % no mag no az
                    keyboard % It seems that Endahl always has magnitude ...
                    % No Magnitude, the str is '<pick publicID'
                    nuev = nuev + 1; % event counter
                    ev(nuev,:) = [otdate lo la de Nph -999];
                    mvalu(nuev,1) = -999;
                    mtype(nuev,1) = {mt};
                    ID(nuev,1)    = str2num(char(d(id).name(1:end-4)));
                    %clear otdate lo la de Nph ma
                    ev_info_collected = 0;
                elseif ev_info_collected == 7 % both mag and az
                    % You found the magnitude, the str is '<value>1.2...'
                    nuev = nuev + 1; % event counter
                    ev(nuev,:) = [otdate lo la de Nph gap];
                    mvalu(nuev,1) = ma;
                    mtype(nuev,1) = {mt};
                    ID(nuev,1)    = str2num(char(d(id).name(1:end-4)));
                    %clear otdate lo la de Nph ma
                    ev_info_collected = 0;
                end
                
                % check if the same number of phases and distances were
                % read:
                if Nphr == ev_distance_collected
                    disp(['Number of phases read: ',num2str(Nphr),' '])
                    disp(['Number of distances read: ',num2str(ev_distance_collected),' '])
                else
                    keyboard
                end
                
                if Nphr == 0
                    disp(['Number of phases in file: ',char(d(id).name), ' is zero']);
                    disp(['Build empty variables to be stored in picks ... '])
                    at = [];
                    evst_dist = [];
                    ch='';
                    ne='';
                    st='';
                    ph='';
                else
                    if Nph == length(at)
                        % great, we expected to have Nph, and we read the same
                        % number of phases ...
                        if comment; disp(['Successfull, Number of phases read: ',num2str(Nph)]); end
                    else
                        if comment; disp(['This file must contain: ',num2str(Nph),' phase[s]']); end
                        if comment; disp(['But it contains: ',num2str(length(at))]); end
                        if Nph < length(at)
                            if comment; disp(['The program over-writes the number of phases']); end
                            Nph = length(at);
                        else
                            if comment; disp(['You have less  arrival information!']); end
                            Nph = length(at);
                        end
                    end
                end
                picks(nuev).at = datevec(at);
                % File: 10770.xml produces 3455 distance values ! 
                % I checked the 10770.xml, it seems that it actually has
                % more than 3000 phases and distances ...
                % maybe the problem starts when the next file is read.
                % 04 Dec 2018, It seems that I have solved the problem ...
                %  if length(dist) == 3455
                %              keyboard
                %  else
                picks(nuev).dist = evst_dist;
                %  end
                picks(nuev).ch = char(ch);
                picks(nuev).ne = char(ne);
                picks(nuev).st = char(st);
                picks(nuev).ph = char(ph);
                % If you reach here, you'r done, get out of the while ...
                break
            end
            fclose(fid);
            clear ot lo la de Nph ma at ch ne st ph gap evst_dist
            
            % Display ... 
            if mod(id,10) == 1
                disp([num2str(id),' files read: ']);
                disp(['Size of the ev: ',ByteSize(ev)]);
                disp(['Size of the picks: ',ByteSize(picks)]);
            end
            % Save to a file ... 
            if id == length(d)
                                
                % Figure out the magnitudes and their types:
                % Indexes
                indml = find(strcmpi(mtype,'ML'));
                indmb = find(strcmpi(mtype,'Mb'));
                indms = find(strcmpi(mtype,'Ms'));
                indmw = find(strcmpi(mtype,'Mw'));
                
                % Creates a variable of NaN or a cell of '', the later does not
                % work because we want to build the GG structure as a matrix
                % containing numbers ...
                % Creates a variables of -999
                ml = repmat(-999,size(mtype));
                mb = repmat(-999,size(mtype));
                ms = repmat(-999,size(mtype));
                mw = repmat(-999,size(mtype));
                
                % Re-fill the ml, mb, ms and mw using Indexes
                ml(indml) = mvalu(indml);
                mb(indmb) = mvalu(indmb);
                ms(indms) = mvalu(indms);
                mw(indmw) = mvalu(indmw);
                
                % [otdate lo la de Nph -999]
                A = [ev mb ms ml mw ID];
                [B,I] = sortrows(A,1);
                ev = [B(:,1:5), B(:,7:11), B(:,6)]; % Event information, sorted
                picks = picks(I); % pick information, pick information sorted
                clear A B I
                %                       10 11 12 13 14 15
                EHB.ev = [datevec(ev(:,1)) ev(:,2:end)];
                EHB.picks=picks; 
                clear ev picks
                % You have finished reading files ... Save them ... 
                %save -v7.3 EHB73 EHB;
                save EHB EHB;
                
                % Write the output with the following format:
                % #Reporter,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,Nph ,mb, ms, ml, ID
                % ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
%                 fid=fopen('EHB.csv','w+');
%                 for i=1:length(EHB.ev(:,1))
%                     % Write the event info
%                     %                  1       2       3        4      5       6        7       8        9        10       11    12      13      14       15     16
%                     %            Rep  %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID     Gap
%                     fprintf(fid,'#EHB, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f, %7.0f, %04.2f \n',...
%                         EHB.ev(i,1:16)                                                                                            );
%                     for j=1:length(EHB.picks(i).at(:,1))
%                         fprintf(fid,'%s, ',EHB.picks(i).st(j,:)); % stations
%                         fprintf(fid,'%s, ',EHB.picks(i).ch(j,:)); % channel
%                         fprintf(fid,'%s, ',' '); % location code 00, 10, ...
%                         fprintf(fid,'%s, ',' '); % network
%                         fprintf(fid,'%s, ',' '); % lon
%                         fprintf(fid,'%s, ',' '); % lat
%                         fprintf(fid,'%s, ',' '); % ele
%                         fprintf(fid,'%s, ',EHB.picks(i).ph(j,:)); % phase
%                         fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',EHB.picks(i).at(j,:)); % YY(j),MM(j),DD(j),hh(j),mm(j),ss(j));
%                         fprintf(fid,'%07.4f \n',EHB.picks(i).dist(j)); % YY(j),MM(j),DD(j),hh(j),mm(j),ss(j));
%                     end
%                 end
%                 fclose(fid);
                clear EHB YY MM DD hh mm ss ts mb ms ml at ot
            end
            
        elseif strcmp(folds{ifo},'ISC-Arrivals') % READ monthly csv files
            
            %DATA_TYPE ARRIVAL:ASSOCIATED CSV
            %--EVENT--|---------------------------------------------------------ARRIVAL DATA----------------------------------------------------------|---------------ORIGIN DATA (PRIME HYPOCENTRE)------------|---EVENT MAGNITUDE--
            %EVENTID  ,REPORTER ,STA  ,LAT     ,LON      ,ELEV   ,CHN,DIST  ,BAZ  ,ISCPHASE,REPPHASE,DATE      ,TIME       ,RES  ,TDEF,AMPLITUDE,PER  ,AUTHOR   ,DATE      ,TIME       ,LAT     ,LON      ,DEPTH,AUTHOR   ,TYPE  ,MAG
            %    1         2       3      4        5          6    7   8       9    10        11        12        13         14    15     16       17    18        19         20          21         22     23       24     25     26
            %608346729,         ,JNBK , 42.2800, 142.7527,  115.0,???,  0.17,187.1,Pn      ,P       ,2016-01-01,00:38:11.90,  0.5,True,         ,     ,ISC      ,2016-01-01,00:37:58.24, 42.4538, 142.7817, 92.0,      ISC,mb    , 3.9
            %608346729,         ,JSHD , 42.4098, 142.4667,   79.0,???,  0.24,259.4,Pn      ,P       ,2016-01-01,00:38:11.80,  0.2,True,         ,     ,ISC      ,2016-01-01,00:37:58.24, 42.4538, 142.7817, 92.0,      ISC,mb    , 3.9
            fid=fopen(file,'r'); %NCI
            str=fgetl(fid);
            while isempty(strfind(str,'EVENTID  ,'))
                str=fgetl(fid);
            end %            1   2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
            if comment; disp([str]); end
            ts=textscan(fid,'%f %s %s %f %f %f %s %f %f %s %s %s %s %f %s %f %f %s %s %s %f %f %f %s %s %f','Delimiter',',');
            fclose(fid);
            % Arrival times, Origing time, lon lat depth mag of the events
            space1 = repmat(' ',size(ts{12}));
            date1 = char(ts{12});
            time1 = char(ts{13});
            % This one takes 16 seconds for a file of size 55 MB
            tic;
            YYa = str2num(date1(:,1:4));
            MMa = str2num(date1(:,6:7));
            DDa = str2num(date1(:,9:10));
            hha = str2num(time1(:,1:2));
            mma = str2num(time1(:,4:5));
            ssa = str2num(time1(:,7:end));
            %at = datenum(YY,MM,DD,hh,mm,ss);
            toc;
            % This one takes 88 seconds for a file of size 55 MB
            %tic; at=datenum([char(ts{12}),space1,char(ts{13})]); toc;
            %ot = datenum([char(evymd),space,char(evtim)]); % takes a lot of time
            num = find(diff(ts{1}));
            endp = [num;length(ts{1})]; % End of each event
            begp = [1;num+1]; % Begin of each event
            
            % Figure out the magnitudes and their types:
            % Indexes
            m = cellstr(deblank(char(ts{25}(begp)))); % Convert to character, remove the space (deblank), convert back to cell
            indempty = find(strcmpi(m,'')); % Some of the ISC events don't have magnitude
            % Remove the events with no magnitude
            begp(indempty) = [];
            endp(indempty) = [];
            m   (indempty) = [];
            % Find the ones with Ml, mb ms and mw
            indmb = find(strcmpi(m,'mb'));
            indms = find(strcmpi(m,'MS'));
            indml = find(strcmpi(m,'ML'));
            indmw = find(strcmpi(m,'Mw'));
            
            % Creates a variable of NaN or a cell of '', the later does not
            % work because we want to build the GG structure as a matrix
            % containing numbers ...
            % Creates a variables of -999
            mb = repmat(-999,size(m));
            ms = repmat(-999,size(m));
            ml = repmat(-999,size(m));
            mw = repmat(-999,size(m));
            
            % Re-fill the ml, mb, ms and mw using Indexes
            mb(indmb) = ts{26}(begp(indmb));
            ms(indms) = ts{26}(begp(indms));
            ml(indml) = ts{26}(begp(indml));
            mw(indmw) = ts{26}(begp(indmw));
            
            %space2 = repmat(' ',size(begp));
            % This one is fast
            date1 = char(ts{19}(begp));
            time1 = char(ts{20}(begp));
            %tic; 
            YY = str2num(date1(:,1:4));
            MM = str2num(date1(:,6:7));
            DD = str2num(date1(:,9:10));
            hh = str2num(time1(:,1:2));
            mm = str2num(time1(:,4:5));
            ss = str2num(time1(:,7:end));
            ot = datenum(YY,MM,DD,hh,mm,ss);
            %toc;
            % This one is slow
            %tic; ot = datenum([char(ts{19}(begp)),space2,char(ts{20}(begp))]); toc;
            
            if id == 1 % The first file ...
                %    Origin   Lon        Lat          Dep        NumberOfPhases   Mag         EVENTID
                ev  = [ot ts{22}(begp) ts{21}(begp) ts{23}(begp) (endp-begp)+1 mb ms ml mw ts{1}(begp)];
                %evnumisccsv = evnumisccsv + length(begp); % The last row stored for the first file
            else
                ev1 = [ot ts{22}(begp) ts{21}(begp) ts{23}(begp) (endp-begp)+1 mb ms ml mw ts{1}(begp)];
                ev  = [ev ; ev1]; clear ev1;
            end
            
            ch=ts{7}; qind = find(strcmp(ch,'???')); ch(qind) = {'?'};
            for iev = 1:length(begp)
                ind = begp(iev) : endp(iev);
                phnumisccsv = phnumisccsv +1;
                picks(phnumisccsv).at = [YYa(ind) MMa(ind) DDa(ind) hha(ind) mma(ind) ssa(ind)];
                picks(phnumisccsv).dist = ts{8}(ind);
                picks(phnumisccsv).ch = deblank(char(ch(ind)));
                picks(phnumisccsv).st = deblank(char(ts{3}(ind)));
                picks(phnumisccsv).lo = ts{5}(ind);
                picks(phnumisccsv).la = ts{4}(ind);
                picks(phnumisccsv).el = ts{6}(ind);
                picks(phnumisccsv).ph = deblank(char(ts{10}(ind)));
                
                g    = ts{9}(ind); % This is baz for this event
                g = sort(g); % sort them, 
                g = [g; min(g)+360];% add 360 to the minimum
                diffg = max(diff(g)); % Take the maximum of differences, this is the maximum azimutham gap
                gap(phnumisccsv,1)    = diffg;
            end
            
            clear qind ch ml mb ms mw indml indmb indms indmw g diffg
            clear YY MM DD hh mm ss YYa MMa DDa hha mma ssa space1 time1 date1 ts
            
            if mod(id,10) == 1
                disp([num2str(id),' files read: ']);
                disp(['Size of the ev: ',ByteSize(ev)]);
                disp(['Size of the picks: ',ByteSize(picks)]);
            end
            
            if id == length(d)
                disp(['Reading ISC finished ... ']);
                if length(ev(:,1)) == length(picks) % This is a strange test ... 
                    disp(['  ---- > Number of events: ',num2str(length(picks))]);
                else
                    disp(['Problem in reading ISC ... ']);
                    keyboard
                end
                ev = [ev, gap];
                % Sort the ISC based on origin time
                [ev,I]    = sortrows(ev,1); % Sort the events by their origin time
                ISC.ev    = [datevec(ev(:,1)), ev(:,2:end)]; clear ev; % Store the event origin time in YYYY MM DD hh mm ss.ss, remove ev
                ISC.picks = picks(I); clear picks; % Sort the picks structure by the same order ... 
                
                save -v7.3 ISC73 ISC % save it as a MATLAB structure
                
                % Maybe store it in two 1GB files without the -v7.3 ... 
                
                fid=fopen('ISC.csv','w+'); % Write into an ascii file (csv file). 
                for i=1:length(ISC.ev(:,1))
                    % Write the event info
                    %                  1       2       3        4      5       6        7       8        9        10       11    12      13      14       15       16
                    %            Rep  %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID      Gap
                    fprintf(fid,'#ISC, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f, %7.0f, %04.2f\n',...
                                      ISC.ev(i,1:16));
                    for j=1:length(ISC.picks(i).at(:,1))
                        fprintf(fid,'%s, ',char(ISC.picks(i).st(j,:))); % stations
                        fprintf(fid,'%s, ',char(ISC.picks(i).ch(j,:))); % channel
                        fprintf(fid,'%s, ',' '); % location code 00, 10, ...
                        fprintf(fid,'%s, ',' '); % network
                        fprintf(fid,'%010.5f, ',ISC.picks(i).lo(j)); % lon
                        fprintf(fid,'%010.5f, ',ISC.picks(i).la(j)); % lat
                        fprintf(fid,'%010.5f, ',ISC.picks(i).el(j)); % ele
                        fprintf(fid,'%s, ',char(ISC.picks(i).ph(j,:))); % phase
                        fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',ISC.picks(i).at(j,:));                        
                        fprintf(fid,'%07.4f \n',ISC.picks(i).dist(j));
                    end
                end
                fclose(fid);
                clear ISC;
            end
            
        elseif strcmp(folds{ifo},'analyst-reviewed2018Aug') % READ XMLs Different format
            
            % =============================================================
            % FOR LATER 
            % =============================================================
            disp(['-->      Working on file: ',char(d(id).name)]);
            fid = fopen(file,'r'); % open the file
            %				<longitude> 27
            %				<depth> 23
            %				<time> 22
            %				<latitude> 26
            %                                        <usedPhaseCount> 56
            %					<associatedPhaseCount> 42
            %			<magnitude publicID= 32
            ev_info_collected = 0;
            num =0;
            cl =0 ;
            while ~feof(fid)
                str=fgetl(fid);
                cl=cl+1;
                % While1: Keep reading until you reach </origin>
                while isempty(strfind(str,'</origin>'))
                    str=fgetl(fid);
                    if ~isempty(strfind(str,'<longitude>')) % length(lin) >= 27
                        sv = 1; % To start the saerch for the <value>
                        while sv == 1 % Enter the while loop
                            str = fgetl(fid); % read a line
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<'); % If yes, then extract the value
                                lo       = str2num(str(larind(1)+1:smaind(2)-1)); % store the value into longitude
                                ev_info_collected = ev_info_collected + 1; % one information is loaded
                                if comment; disp('---        Found lo'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<depth>')) % length(lin) >= 23
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                de  = str2num(str(larind(1)+1:smaind(2)-1))/1000; % Depth in meters
                                ev_info_collected = ev_info_collected + 1;
                                sv=0; % You found the value, get out
                                if comment; disp('---        Found dep'); end
                            end
                        end
                    elseif ~isempty(strfind(str,'<time>')) % length(lin) >= 22
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                otstr  = str(larind(1)+1:smaind(2)-1);
                                otdate = datenum([otstr(1:10),' ',otstr(12:end-1)]); % convert the origin time into datenum
                                ev_info_collected = ev_info_collected + 1;
                                if comment; disp('---        Found ot'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<latitude>')) % length(lin) >= 26
                        sv = 1;
                        while sv == 1
                            str = fgetl(fid);
                            if ~isempty(strfind(str,'<value>')) % Is it <value>?
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                la  = str2num(str(larind(1)+1:smaind(2)-1));
                                ev_info_collected = ev_info_collected + 1;
                                if comment; disp('---        Found la'); end
                                sv=0; % You found the value, get out
                            end
                        end
                    elseif ~isempty(strfind(str,'<associatedPhaseCount>')) % length(lin) >= 42
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        Nph  = str2num(str(larind(1)+1:smaind(2)-1));
                        ev_info_collected = ev_info_collected + 1;
                        if comment; disp('---        Found Nph'); end
                    end
                    %                    if ev_info_collected == 5
                    %                         nuev = nuev + 1; % event counter
                    %                         GA.ev(nuev,:) = [otdate lo la de Nph];
                    %                         clear otdate lo la de Nph ma
                    %                         break % get out of the while loop
                    %                    end
                end
                if ev_info_collected == 5
                    % You found the ot lo ladep Nph, Search for the magnitude
                    
                else
                    disp('Reached </origin> ... ')
                    disp(['file ',char(d(id).name),' does not contain the 5 information: ot lo la de Nph'])
                    keyboard
                end
                
                % Before you reach '<pick publicID', magnitude might appear
                while isempty(strfind(str,'<pick publicID'))
                    % The 'pick publicID' hasn't appear yet...
                    str = fgetl(fid);
                    if ~isempty(strfind(str,'<magnitude publicID'))
                        while isempty(strfind(str,'<value>'))
                            str = fgetl(fid);
                        end
                        larind   = strfind(str,'>'); smaind = strfind(str,'<');
                        ma  = str2num(str(larind(1)+1:smaind(2)-1));
                        ev_info_collected = ev_info_collected + 1;
                        if comment; disp('---        Found Mag'); end
                        break
                    end
                end
                
                % Did we get the magnitude?
                if ev_info_collected == 6
                    % You found the magnitude, the str is '<value>1.2...'
                    nuev = nuev + 1; % event counter
                    GA.ev(nuev,:) = [otdate lo la de Nph ma];
                    %clear otdate lo la de Nph ma
                    ev_info_collected = 0;
                else % Or not?
                    % No Magnitude, the str is '<pick publicID'
                    ma = NaN;
                    nuev = nuev + 1; % event counter
                    GA.ev(nuev,:) = [otdate lo la de Nph ma];
                    %clear otdate lo la de Nph ma
                    ev_info_collected = 0;
                end
                
                Nphr = 0; % Number of phases read so far ...
                % Load the picks:
                while isempty(strfind(str,'</event>')) % We will not read all lines down to <event>, we will exit the loop if we reach <magnitude>
                    while ~isempty(strfind(str,'pick publicID'))
                        % This is a pick
                        while isempty(strfind(str,'</pick>'))
                            if ~isempty(strfind(str,'<time>'))
                                str=fgetl(fid); % This is <value>
                                if isempty(strfind(str,'<uncertainty>')) %
                                    % So this is value, go on ...
                                else
                                    % So this is <uncertainty>, read again
                                    str=fgetl(fid);
                                end
                                Nphr = Nphr + 1; % One time is read
                                larind    = strfind(str,'>'); smaind = strfind(str,'<');
                                atstr     = str(larind(1)+1:smaind(2)-1);
                                at(Nphr)  = datenum([atstr(1:10),' ',atstr(12:end-1)]); % Phase arrival time
                            elseif ~isempty(strfind(str,'<waveformID'))
                                cind = strfind(str,'"');
                                ch(Nphr)  = {str(cind(1)+1:cind(2)-1)};
                                lc(Nphr)  = {str(cind(3)+1:cind(4)-1)};
                                ne(Nphr)  = {str(cind(5)+1:cind(6)-1)};
                                st(Nphr)  = {str(cind(7)+1:cind(8)-1)};
                            elseif ~isempty(strfind(str,'<phaseHint>'))
                                larind   = strfind(str,'>'); smaind = strfind(str,'<');
                                ph(Nphr)  = {str(larind(1)+1:smaind(2)-1)}; % Phase arrival time
                            end
                            str = fgetl(fid);
                        end
                        %            <pick publicID="quakeml:ga.ga.gov.au/pick/4311408">
                        % 				<creationInfo>
                        % 					<creationTime>2018-05-22T21:33:56.186Z</creationTime>
                        % 					<author>dbp:lsaikal:181</author>
                        % 				</creationInfo>
                        % 				<backazimuth>
                        % 					<value>330.03</value>
                        % 				</backazimuth>
                        % 				<time>
                        % 					<value>2018-05-22T06:31:57.768Z</value>
                        % 				</time>
                        % 				<waveformID channelCode="BHE" locationCode="" networkCode="AU" stationCode="MULG">quakeml:ga.ga.gov.au/waveform/1270</waveformID>
                        % 				<evaluationStatus>preliminary</evaluationStatus>
                        % 				<phaseHint>S</phaseHint>
                        % 				<evaluationMode>automatic</evaluationMode>
                        % 			</pick>
                    end
                    str=fgetl(fid);
                end
                if Nph == length(at)
                    % great, we expected to have Nph, and we read the same
                    % number of phases ...
                    disp(['Successfull, Number of phases read: ',num2str(Nph)]);
                else
                    disp(['This file must contain: ',num2str(Nph),' phase[s]'])
                    disp(['But it contains: ',num2str(length(at))])
                    if Nph < length(at)
                        disp(['The program over-writes the number of phases']);
                        Nph = length(at);
                    else
                        disp(['You have less  arrival information!'])
                        Nph = length(at);
                    end
                end
                GA.picks(nuev).at = at;
                GA.picks(nuev).ch = ch;
                GA.picks(nuev).lc = lc;
                GA.picks(nuev).ne = ne;
                GA.picks(nuev).st = st;
                GA.picks(nuev).ph = ph;
                % If you reach here, you'r done, get out of the while ...
                break
            end
            fclose(fid);
            clear ot lo la de Nph ma at ch lc ne st ph
            
        elseif strcmp(folds{ifo},'isc_events') % READ monthly XMLs
            % =============================================================
            % FOR LATER 
            % =============================================================     
            
        end
    end
    toc
end