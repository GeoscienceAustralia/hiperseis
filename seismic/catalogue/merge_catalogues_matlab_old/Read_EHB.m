%% Reads different catalogues and produces ascii files in a unified format
% =========================================================================
% Babak Hejrani, Geoscience Australia, 14 January 2019
%
% This program reads ascii and xml files and wties output ascii files
% with the following format (suggested by Alexei Gorbatov):
% #EV-ID Time Lon Lat Dep Mag Nph
% ST-ID Time Phase
% Update 28 Oct 2018: I prefer a new format, described below ...
% #EV ,YY ,MM ,DD ,hh ,mm ,ss ,lo ,la ,de ,ph ,ma
% ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT
%
% Some weird cases in EHB: 202.xml has no picks! but he Nph is 607 :)
%
% =========================================================================

clear all; close all; clc; if ispc; sep='\'; else sep='/'; end % This code will work both in Linux and windows
comment = 0; s1 = datenum(2000,1,1,0,0,0)-datenum(2000,1,1,0,0,1); % One second in datenum format
% The top folder where all catalogues are located
tf = '/g/data1a/ha3/Passive/Events/'; % This is top folder that all the data is located,
% The location of different catalogues, folder names under the top folder
folds = {'engdahl-events'};
% format of the files in each folder: ascii or xml
forms = {'/*.xml'};

% Take a lif of xml files in the folder
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
    if comment; disp(['File: ',char(d(id).name)]); end
    fid = fopen(file,'r'); % open the file
    
    ev_info_collected = 0; % Info collected so far. 
    ev_distance_collected = 0; % Distance information collected so far. 
    
    
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
        fid=fopen('EHB.csv','w+');
        for i=1:length(EHB.ev(:,1))
            % Write the event info
            %                  1       2       3        4      5       6        7       8        9        10       11    12      13      14       15     16
            %            Rep  %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID     Gap
            fprintf(fid,'#EHB, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, %010.5f, %010.5f, %010.5f, %05.0f, %04.2f, %04.2f, %04.2f, %04.2f, %7.0f, %04.2f \n',...
                EHB.ev(i,1:16)                                                                                            );
            for j=1:length(EHB.picks(i).at(:,1))
                fprintf(fid,'%s, ',EHB.picks(i).st(j,:)); % stations
                fprintf(fid,'%s, ',EHB.picks(i).ch(j,:)); % channel
                fprintf(fid,'%s, ',' '); % location code 00, 10, ...
                fprintf(fid,'%s, ',' '); % network
                fprintf(fid,'%s, ',' '); % lon
                fprintf(fid,'%s, ',' '); % lat
                fprintf(fid,'%s, ',' '); % ele
                fprintf(fid,'%s, ',EHB.picks(i).ph(j,:)); % phase
                fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',EHB.picks(i).at(j,:)); % YY(j),MM(j),DD(j),hh(j),mm(j),ss(j));
                fprintf(fid,'%07.4f \n',EHB.picks(i).dist(j)); % YY(j),MM(j),DD(j),hh(j),mm(j),ss(j));
            end
        end
        fclose(fid);
        clear EHB YY MM DD hh mm ss ts mb ms ml at ot
    end
    
end
toc