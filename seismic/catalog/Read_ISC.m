%% Reads ISC arrivals and produces ascii files in a new format format
% =========================================================================
% Babak Hejrani, Geoscience Australia, 14 January 2019
%
% This program reads ascii ISC-Arrivals and wties output ascii files
% with the following format (suggested by Alexei Gorbatov):
% #EV-ID Time Lon Lat Dep Mag Nph
% ST-ID Time Phase
% Update 28 Oct 2018: I prefer a new format, described below ...
% #EV ,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,ph ,ma
% ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
%
% =========================================================================

clear all; close all; clc; if ispc; sep='\'; else sep='/'; end % This code will work both in Linux and windows
comment = 1; s1 = datenum(2000,1,1,0,0,0)-datenum(2000,1,1,0,0,1); % One second in datenum format
% The top folder where all catalogues are located
tf = '/g/data1a/ha3/Passive/Events/';
% The location of different catalogues, folder names under the top folder
folds = {'ISC-Arrivals'};
% format of the files in each folder: ascii or xml
forms = {'/Arrivals_*'};


% The ISC data is from 1970 till today, Do you wanted all? Do you want a specifioc period of time? 
starting = 1993; % start
ending   = 2016; % end

% For loop over events
folder = [tf,sep,folds{1},sep]; % Take a list of files in the folder
disp(['Working on folder: ',folder]);
disp([' ']);
d = dir([folder,forms{1}]); % Matlab 2016b and later
disp(['--> Number of files for: ',forms{1},' is: ',num2str(length(d))]);
if length(d) < 1
    disp(['--> No file with that format in this folder!'])
    keyboard
end

% Found some files, Find the ones that are within [starting:ending].
for i = 1:length(d)
    ind = strfind(d(i).name,'_');
    yyy(i,1) = str2num(d(i).name(ind(2)+1:ind(2)+4));
end
yyy_ind = find(yyy >= starting & yyy <= ending);

nuev=0; % number of events loaded so far ...
phnumisccsv = 0; % number of ISC events loaded so far ...

for id =yyy_ind' % Loop over files ...
    %file = [char(d(id).folder),sep,char(d(id).name)]; %NCI, MATLAB 2016b and later
    file = [folder,sep,char(d(id).name)]; % ANU Matlab 2014a and later
    tic
    % READ monthly csv files
    %DATA_TYPE ARRIVAL:ASSOCIATED CSV
    %--EVENT--|---------------------------------------------------------ARRIVAL DATA----------------------------------------------------------|---------------ORIGIN DATA (PRIME HYPOCENTRE)------------|---EVENT MAGNITUDE--
    %EVENTID  ,REPORTER ,STA  ,LAT     ,LON      ,ELEV   ,CHN,DIST  ,BAZ  ,ISCPHASE,REPPHASE,DATE      ,TIME       ,RES  ,TDEF,AMPLITUDE,PER  ,AUTHOR   ,DATE      ,TIME       ,LAT     ,LON      ,DEPTH,AUTHOR   ,TYPE  ,MAG
    %    1         2       3      4        5          6    7   8       9    10        11        12        13         14    15     16       17    18        19         20          21         22     23       24     25     26
    %608346729,         ,JNBK , 42.2800, 142.7527,  115.0,???,  0.17,187.1,Pn      ,P       ,2016-01-01,00:38:11.90,  0.5,True,         ,     ,ISC      ,2016-01-01,00:37:58.24, 42.4538, 142.7817, 92.0,      ISC,mb    , 3.9
    %608346729,         ,JSHD , 42.4098, 142.4667,   79.0,???,  0.24,259.4,Pn      ,P       ,2016-01-01,00:38:11.80,  0.2,True,         ,     ,ISC      ,2016-01-01,00:37:58.24, 42.4538, 142.7817, 92.0,      ISC,mb    , 3.9
    fid=fopen(file,'r'); 
    str=fgetl(fid);
    while isempty(strfind(str,'EVENTID  ,'))
        str=fgetl(fid);
    end %            1   2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    if comment; disp(['File: ',char(d(id).name)]); end
    ts=textscan(fid,'%f %s %s %f %f %f %s %f %f %s %s %s %s %f %s %f %f %s %s %s %f %f %f %s %s %f','Delimiter',',');
    fclose(fid);
    % Arrival times, Origing time, lon lat depth mag of the events
    % space1 = repmat(' ',size(ts{12}));
    date1 = char(ts{12});
    time1 = char(ts{13});
    % This one takes 16 seconds for a file of size 55 MB
    %tic;
    YYa = str2num(date1(:,1:4));
    MMa = str2num(date1(:,6:7));
    DDa = str2num(date1(:,9:10));
    hha = str2num(time1(:,1:2));
    mma = str2num(time1(:,4:5));
    ssa = str2num(time1(:,7:end));
    
    %at = datenum(YY,MM,DD,hh,mm,ss);
    %toc;
    % This one takes 88 seconds for a file of size 55 MB
    %tic; at=datenum([char(ts{12}),space1,char(ts{13})]); toc;
    %ot = datenum([char(evymd),space,char(evtim)]); % takes a lot of time
    
    % Find the events: 
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
    
    if id == yyy_ind(1) % The first file ...
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
        
        % Check the size, I don't want to deal with the -v 7.3. The version
        % 7.3 can not be read in python. 
        bscat = ByteSize(ISC); % output is a string. 
        space_ind = strfind(bscat,' '); % Find the space. 
        bsv = str2num(bscat(1:space_ind-1)); % Size of the volume (a number)
        bsb = bscat(space_ind+1:end); % GB, Kb or TB?
        
        if strcmp(bsb,'Gb')
            if bsv > 1.75 % Is it larger than 1 Gb?
                % This is less than 2 Gb, divide it into 4. 
                lent  = length(ISC.ev(:,1));
                len4  = fix(lent/4); % divide the length by 4. 
                % save the first one                
                ISC1.ev    = ISC.ev(1:len4,:);
                ISC1.picks = ISC.picks(1:len4);
                save ISC1 ISC1; % save it as a MATLAB structure
                clear ISC1;
                % save the second one
                ISC2.ev    = ISC.ev(len4+1:2*len4,:);
                ISC2.picks = ISC.picks(len4+1:2*len4);
                save ISC2 ISC2; % save it as a MATLAB structure
                clear ISC2;
                % save the third one
                ISC3.ev    = ISC.ev(2*len4+1:3*len4,:);
                ISC3.picks = ISC.picks(2*len4+1:3*len4);
                save ISC3 ISC3; % save it as a MATLAB structure
                clear ISC3;                
                % save the forth one
                ISC4.ev    = ISC.ev(3*len4+1:end,:);
                ISC4.picks = ISC.picks(3*len4+1:end);
                save ISC4 ISC4; % save it as a MATLAB structure
                clear ISC4;
            end
        elseif strcmp(bsb,'Mb')
            save ISC ISC % save it as a MATLAB structure
        end
                
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
    toc
end