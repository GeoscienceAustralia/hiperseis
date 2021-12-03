%% Reads different catalogues and produces ascii files in a unified format
% =========================================================================
% Babak Hejrani, Geoscience Australia, 22 Oct 2018
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
folds = {'GeoscienceAustralia'};
% format of the files in each folder: ascii or xml
forms = {'/*.txt'};

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
    ot = datenum([char(ts{3}),space,char(ts{5})]);
    [YY MM DD hh mm ss] = datevec(ot);
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
    at=datenum([char(ts{21}),space,char(ts{23})]);
    [YY MM DD hh mm ss] = datevec(at); 
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
end
toc
