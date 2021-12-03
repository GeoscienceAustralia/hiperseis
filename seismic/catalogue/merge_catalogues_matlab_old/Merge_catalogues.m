%% Merge 5 catalogues into one.
% =========================================================================
% Babak Hejrani, Geoscience Australia, 08 Nov 2018 - 5 Feb 2019
% 
% The priorities are:
% EHB : 2005/?? - 2009/01 => These are xml files, it is a pain to get the arrivals out
% ISC : 2009/01 - 2016/04 => I have downloaded the csv monthly arrival files
% USGS: 2016/04 - 2018/12 => Haven't looked at it yet ...
% GA  : 2010/07 - 2018/05 => These are xml files, it is a pain to get the arrivals out, Here even the Magnitude is complicated [See 1466768.xml and 1000544.xml]
% GG  : 1857    - 2017    => These are just hypocenters, It would be good to compare the two catalogues (GA VS GG)
% 
% =========================================================================
clear all; close all; clc; format compact;
one_sec=datenum(2000,1,1,0,0,1)-datenum(2000,1,1,0,0,0);
comment = 0;
distmin=50; % Minimum distance
GAhypoGG =0;
mmin = 5.4; % Minimum magnitude for telseismic events in ISC catalogue

azgmin = 180; % Inside Australian continent. 
lonAUS1 = 105; % Australian continent in ISC catalogue
lonAUS2 = 160;
latAUS1 = -40;
latAUS2 = -13;

m1min = 4.9;
lon1OCE1 = 015; % Oceanic areas in ISC catalogue
lon1OCE2 = 090;
lat1OCE1 = -60;
lat1OCE2 =  30;

m2min = 4.9;
lon2OCE1 = 090.001; % Oceanic areas in ISC catalogue
lon2OCE2 = 180;
lat2OCE1 = -60;
lat2OCE2 = -45;

%% LOAD DATA:
disp(['Loading the data: '])
tic

usgsd = dir('USGS.mat');
disp([' --> Looking for USGS: '])
if isempty(usgsd)
    disp(['        No USGS data!']);
else
    load USGS; bsUSGS = ByteSize(USGS)
end

gad = dir('GA.mat');
disp([' --> Looking for GA: '])
if isempty(gad)
    disp(['        No GA data!']);
else
    load GA; bsGA = ByteSize(GA)
end

ggd = dir('GG.mat');
disp([' --> Looking for Gary Gibson: '])
if isempty(ggd)
    disp(['        No Gary Gibson data!']);
else
    load GG; bsGG = ByteSize(GG)
end

iscd = dir('ISC*.mat');
disp([' --> Looking for ISC: '])
if isempty(iscd)
    disp(['        No ISC data!']);
else
    for i = 1:length(iscd)
        eval(['load ',char(iscd(i).name)]);
    end
    ISC.ev = [ISC1.ev; ISC2.ev; ISC3.ev; ISC4.ev];    
    table1 = struct2table(ISC1.picks); % Convert the ISC1 picks into table
    table2 = struct2table(ISC2.picks); % Convert the ISC2 picks into table
    table3 = struct2table(ISC3.picks); % Convert the ISC3 picks into table
    table4 = struct2table(ISC4.picks); % Convert the ISC4 picks into table
    table5 = [table1; table2; table3; table4]; clear table1 table2 table3 table4; % Merge them
    ISC.picks = table2struct(table5);
    clear ISC1 ISC2 ISC3 ISC4;
    bsISC = ByteSize(ISC)
end

ehbd = dir('EHB.mat');
disp([' --> Looking for EHB: ']);
if isempty(ehbd)
    disp(['        No EHB data!']);
else
    load EHB; bsEHB = ByteSize(EHB)
end

toc

disp(' =============================== ')
disp(['Total number of events loaded (a lot of overlap): ']);
length([ISC.ev(:,1);EHB.ev(:,1);GA.ev(:,1);GG.ev(:,1);USGS.ev(:,1)])
disp(' =============================== ')

%% Visualize
f1 = figure; f1.Color = 'w';
axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 110]);
% load coastlines;
% patchm(coastlat,coastlon,'w','FaceColor','none');

% Geoshow can map things
%hst = geoshow(stlat, stlon, 'DisplayType', 'Point', 'Marker', 'd', 'Color', 'red','MarkerFaceColor', 'red', 'MarkerSize', 4);
%hev = geoshow(evlat, evlon, 'DisplayType', 'Point', 'Marker', '.', 'Color', 'blue','MarkerFaceColor', 'blue', 'MarkerSize', 1);
% Scatterm is much better
a=5;
scatterm(ISC.ev(:,8),ISC.ev(:,7),a,'b','filled'); % cb=colorbar; cb.Label.String = 'Earthquake Depth [km]';
hold on;
scatterm(EHB.ev(:,8),EHB.ev(:,7),a,'r','filled');
scatterm(GA.ev(:,8),GA.ev(:,7),a,'g','filled');
scatterm(GG.ev(:,8),GG.ev(:,7),a,'k','filled');
scatterm(USGS.ev(:,8),USGS.ev(:,7),a,'c','filled');

% hst = scatterm(stlat,stlon,a,      'filled','MarkerFaceColor','red','Marker','d','SizeData',20);

%% Merge the catalogues
%       1       2       3        4      5       6        7       8        9        10       11    12      13      14       15
% Rep  %YY      MM      DD      hh      mm      ss      lo       la       de       ph       mb,   ms      ml,     mw,      ID
% 1  ,2  ,3  ,4  ,5  ,6  ,7  ,8  ,9 - 14
% ST ,CH ,LC ,NE ,LON,LAT,ELE,PH ,AT

% =====================  EHB -2009 % =====================================
% The first event in this catalogue is from 2005 9 16 07 28 39.00
% The last  event in this catalogue is from 2009 3 31 22 21 53.00
% 1x11633 struct array with fields:
%     at, 9 - 14
%     ch, 2
%     ne, 4
%     st, 1
%     ph, 8
% EHB.ev(:,9) = EHB.ev(:,9)*1000; % A correction to EHB. This was a mistake in an earlier version of the Read_catalogues.m
ev   = [EHB.ev, ones(length(EHB.ev(:,1)),1)];
picks= EHB.picks; clear EHB;
picks(1).lc=[]; % 3
picks(1).lo=[]; % 5
picks(1).la=[]; % 6
picks(1).el=[]; % 7
% visualize
f2 = figure; f2.Color = 'w';
a=5;
axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 110]);
scatterm(ev(:,7),ev(:,8),a,'b','filled');
title('EHB');

% ===================== ISC 2009 - 2016 ==================================
% Cut out events in ISC occurred before 2009 04 01
% 1x175057 struct array with fields:
%     at, 9 - 14
%     ch, 2
%     st, 1
%     lo, 5
%     la, 6
%     el, 7
%     ph, 8
ISC.picks(1).lc=[]; % 3
ISC.picks(1).ne=[]; % 4
ot    = datenum(ISC.ev(:,1:6));
mb    = ISC.ev(:,11);
ms    = ISC.ev(:,12);
ml    = ISC.ev(:,13);
mw    = ISC.ev(:,14);

lo    = ISC.ev(:,7);
la    = ISC.ev(:,8);
gap   = ISC.ev(:,16);
%el    = ISC.ev(:,9);

% Filter 1: Exclude the events at 2005 9 16 - 2009 04 01
indot = find(ot < datenum([2005 09 15 23 59 59]) | ot > datenum([2009 04 01 00 00 00])); 
% Filter 2: Inside Australia, take the events with azimuthal gap < 180
% Filter 3: In the Pacific and Indian ocean, take events with magnitude 5. Box 1: lon1OCE1 & lon1OCE2 & lat1OCE1 & lat1OCE2 
% Filter 4: In the Pacific and Indian ocean, take events with magnitude 5. Box 2: lon2OCE1 & lon2OCE2 & lat2OCE1 & lat2OCE2 
% Filter 5: For the rest of the catalogue take events with magnitude 5.5.
indlm = find( ( lo(indot)>lonAUS1 & lo(indot)<lonAUS2 & la(indot)>latAUS1 & la(indot)<latAUS2 & gap(indot)<azgmin ) | ...
              ( lo(indot)>lon1OCE1 & lo(indot)<lon1OCE2 & la(indot)>lat1OCE1 & la(indot)<lat1OCE2 & ( mb(indot)>m1min | ms(indot)>m1min | ml(indot)>m1min | mw(indot)>m1min ) ) | ...
              ( lo(indot)>lon2OCE1 & lo(indot)<lon2OCE2 & la(indot)>lat2OCE1 & la(indot)<lat2OCE2 & ( mb(indot)>m2min | ms(indot)>m2min | ml(indot)>m2min | mw(indot)>m2min ) ) | ...
              ( mb(indot)>mmin | ms(indot)>mmin | ml(indot)>mmin | mw(indot)>mmin ) );

% for comparison and plot: 
m1min_test = 5.5;
m2min_test = 5.5;
indlm_test = find( ( lo(indot)>lonAUS1 & lo(indot)<lonAUS2 & la(indot)>latAUS1 & la(indot)<latAUS2 & gap(indot)<azgmin ) | ...
              ( lo(indot)>lon1OCE1 & lo(indot)<lon1OCE2 & la(indot)>lat1OCE1 & la(indot)<lat1OCE2 & ( mb(indot)>m1min_test | ms(indot)>m1min_test | ml(indot)>m1min_test | mw(indot)>m1min_test ) ) | ...
              ( lo(indot)>lon2OCE1 & lo(indot)<lon2OCE2 & la(indot)>lat2OCE1 & la(indot)<lat2OCE2 & ( mb(indot)>m2min_test | ms(indot)>m2min_test | ml(indot)>m2min_test | mw(indot)>m2min_test ) ) | ...
              ( mb(indot)>mmin | ms(indot)>mmin | ml(indot)>mmin | mw(indot)>mmin ) );

% visualize
f3 = figure; f3.Color = 'w';
axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 110]);
scatterm(la(indot(indlm)),lo(indot(indlm)),a,'b','filled');
scatterm(la(indot(indlm_test)),lo(indot(indlm_test)),a,'r','filled');
title('ISC: Red M>5.4, Blue M is a relaxed for oceanic events')

clear lo la mb ms ml mw ot;

ev    = [ev; [ISC.ev(indot(indlm),:), 2*ones(length(indot(indlm)),1)] ]; % add the ISC's list of events to EHB
% Sort the events: 
otdum1 = datenum(ev(:,1:6));
[B,I] = sortrows(otdum1);

ev = ev(I,:);
clear otdum1 B;

table1= struct2table(picks); % Convert the EHB picks into table
table2= struct2table(ISC.picks(indot(indlm))); clear ISC; % Convert the picks of selected ISC events into table
table3= [table1; table2]; clear table1 table2; % Merge them
picks = table2struct(table3); clear table3; % Convert back to structure
picks = picks(I); % Sort the picks according to events
clear I; 

% ===================== USGS 2016 - 2018 =================================
otdum1 = datenum(USGS.ev(:,1:6));
[B,I] = sortrows(otdum1);
USGS.ev = USGS.ev(I,:); % The data is sorted. 

% Filter the data the way we filter the ISC data. 
mb    = USGS.ev(:,11);
ms    = USGS.ev(:,12);
ml    = USGS.ev(:,13);
mw    = USGS.ev(:,14);
lo    = USGS.ev(:,7);
la    = USGS.ev(:,8);
gap   = USGS.ev(:,16);
%el    = ISC.ev(:,9);

% Filter 2: Inside Australia, take the events with azimuthal gap < 180
% Filter 3: In the Pacific and Indian ocean, take events with magnitude 5. Box 1: lon1OCE1 & lon1OCE2 & lat1OCE1 & lat1OCE2 
% Filter 4: In the Pacific and Indian ocean, take events with magnitude 5. Box 2: lon2OCE1 & lon2OCE2 & lat2OCE1 & lat2OCE2 
% Filter 5: For the rest of the catalogue take events with magnitude 5.5.
indlm = find( ( lo>lonAUS1 & lo<lonAUS2 & la>latAUS1 & la<latAUS2 & gap<360 ) | ...
              ( lo>lon1OCE1 & lo<lon1OCE2 & la>lat1OCE1 & la<lat1OCE2 & ( mb>m1min | ms>m1min | ml>m1min | mw>m1min ) ) | ...
              ( lo>lon2OCE1 & lo<lon2OCE2 & la>lat2OCE1 & la<lat2OCE2 & ( mb>m2min | ms>m2min | ml>m2min | mw>m2min ) ) | ...
              ( mb>mmin | ms>mmin | ml>mmin | mw>mmin ) );

% for comparison and plot: 
m1min_test = 5.5;
m2min_test = 5.5;
indlm_test = find( ( lo>lonAUS1 & lo<lonAUS2 & la>latAUS1 & la<latAUS2 & gap<azgmin ) | ...
              ( lo>lon1OCE1 & lo<lon1OCE2 & la>lat1OCE1 & la<lat1OCE2 & ( mb>m1min_test | ms>m1min_test | ml>m1min_test | mw>m1min_test ) ) | ...
              ( lo>lon2OCE1 & lo<lon2OCE2 & la>lat2OCE1 & la<lat2OCE2 & ( mb>m2min_test | ms>m2min_test | ml>m2min_test | mw>m2min_test ) ) | ...
              ( mb>mmin | ms>mmin | ml>mmin | mw>mmin ) );

f4 = figure; f4.Color = 'w';
axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 110]);
scatterm(la(indlm)     ,lo(indlm)     ,a,'b','filled');
scatterm(la(indlm_test),lo(indlm_test),a,'r','filled');
title('USGS: Red M>5.4, Blue M is a relaxed for oceanic events')

lenev = length(ev);
ev    = [ev; [USGS.ev(indlm,:), 3*ones(length(indlm),1)] ]; % Add USGS to ISC - EHB: 
% preallocate the picks witht the same size of ev.
for i = lenev+1:length(ev)
    picks(i).at = [];
end

% ===================== GG ================================================
otGG = datenum(GG.ev(:,1:6));
GGind = find(otGG>datenum(1993,1,1,0,0,0)); % Exclude data before 1993
GG.ev = [GG.ev(GGind,:) 4*ones(length(GGind),1)];

% ===================== GA ================================================
% Note: Rakib reported an event with strange latitude ... See this file: GA-lat-lon.pdf
% figure; subplot(2,1,1); plot(GA.ev(:,8),'.'); ylabel('Latitude'); subplot(2,1,2); plot(GA.ev(:,7),'.'); ylabel('Longitude');
otGA = datenum(GA.ev(:,1:6));
GAind = find(otGA>datenum(1993,1,1,0,0,0)); % Exclude data before 1993
GA.ev = [GA.ev(GAind,:) 5*ones(length(GAind),1)];
GA.picks = GA.picks(GAind);

if GAhypoGG % Replace GA with GG.
    % Jan 2019: Replace GA with GG.
    % The Gary Gibson has priority.
    % So, I search for similarities within +/-5 s and then check the distance,
    % if the distance is less then 50 km, I consider them the same event and I
    % replace GA's hypocenter information with GG's hypocenter.
    % Comment on 11 Feb 2019: I tried to compare GG adn GA. The result is as
    % follow:
    % num
    % num =
    %         3787
    % num5
    % num5 =
    %    294
    % num5AUS
    % num5AUS =
    %         3489
    % num1LargeDist
    % num1LargeDist =
    %      0
    % num1LargeDistAUS
    % num1LargeDistAUS =
    %      4
    % num1LessDist
    % num1LessDist =
    %     23
    % num1LessDistAUS
    % num1LessDistAUS =
    %         1076
    % numMore1
    % numMore1 =
    %     27
    otdnGA = datenum(GA.ev(:,1:6));
    otdnGG = datenum(GG.ev(:,1:6));
    
    load_ind=[]; num=0;
    load_ind5=[]; num5=0;
    load_ind5AUS=[]; num5AUS=0;
    
    load_ind1LargeDist   = []; num1LargeDist   =0;
    load_ind1LargeDistAUS= []; num1LargeDistAUS=0;
    
    load_ind1LessDist    = []; num1LessDist    =0;
    load_ind1LessDistAUS = []; num1LessDistAUS =0;
    
    %load_indMore1        = [];
    numMore1        =0;
    
    for iev = 1:length(otdnGG)
        % Origin time and location of one event in the GA catalogue
        otdn1 = otdnGG(iev); % Origin time for one event in GA catalogue
        ll1   = GG.ev(iev,7:8); % Lon, lat for one event in GA catalogue
        % Determines whether the earthquake is inside Australia or not ...
        GGAUS = ( ll1(1)>lonAUS1 & ll1(1)<lonAUS2 & ll1(2)>latAUS1 & ll1(2)<latAUS2 );
        % Find the GA event in the compiled one: ev
        ind   = find(abs(otdnGA - otdn1)<5*one_sec); % Will result in 3: >1, =1, <1
        if length(ind) < 1 % you found a new event. Add this to the data ...
            if GGAUS
                if comment; disp(['IND<1; checked +/-05 s, This must be a new event and it is inside Australia ... ']); end
                num=num+1;
                evn(num,:) = GG.ev(iev,:);
                load_ind(num,:)           = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind; iev];
                num5AUS=num5AUS+1;
                load_ind5AUS(num5AUS,:)   = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind5AUS; iev];
            else
                if comment; disp(['IND<1; checked +/-05 s, This must be a new event and it is outside Australia ... ']); end
                num=num+1;
                evn(num,:) = GG.ev(iev,:);
                load_ind(num,:)           = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind; iev];
                num5=num5+1;
                load_ind5(num5,:)         = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind5; iev];
            end
        elseif length(ind) == 1 % There is one event with similar time on ISC catalogue
            % Check the location
            %               lon1       lon2              lat1       lat2
            %ll20 = (abs(ev(ind,7) - ll1(1))<0.2 & abs(ev(ind,8) - ll1(2))<0.2);
            [ar,az]=distance(GA.ev(ind,8),GA.ev(ind,7),ll1(2),ll1(1));
            km = deg2km(ar);
            if km > distmin
                if GGAUS
                    if comment; disp(['IND=1; But the distance is larger than: ',num2str(distmin),' This could be a new event inside Australia ... ']); end
                    num=num+1;
                    evn(num,:) = GG.ev(iev,:);
                    load_ind(num,:)                           = [iev otdn1 ll1 ind otdnGA(ind) km]; % [load_ind; iev];
                    num1LargeDistAUS=num1LargeDistAUS+1;
                    load_ind1LargeDistAUS(num1LargeDistAUS,:) = [iev otdn1 ll1 ind otdnGA(ind) km]; % [load_ind1LargeDistAUS; iev];
                else
                    if comment; disp(['IND=1; But the distance is larger than: ',num2str(distmin),' This could be a new event outside Australia ... ']); end
                    num1LargeDist=num1LargeDist+1;
                    load_ind1LargeDist(num1LargeDist,:)       = [iev otdn1 ll1 ind otdnGA(ind) km]; % [load_ind1LargeDist; iev];
                end
            else
                if GGAUS
                    if comment; disp(['IND=1; and the distance is less than: ',num2str(distmin),' This is not a new event inside Australia ... ']); end
                    num1LessDistAUS=num1LessDistAUS+1;
                    load_ind1LessDistAUS(num1LessDistAUS,:)   = [iev otdn1 ll1 ind otdnGA(ind) km]; % [load_ind1LessDistAUS; iev];
                else
                    if comment; disp(['IND=1; and the distance is less than: ',num2str(distmin),' This is not a new event outside Australia ... ']); end
                    num1LessDist=num1LessDist+1;
                    load_ind1LessDist(num1LessDist,:)         = [iev otdn1 ll1 ind otdnGA(ind) km]; % [load_ind1LessDist; iev];
                end
            end
        elseif length(ind) > 1
            % Several events were found ...
            %ll20 = find(abs(ev(ind,7) - ll1(1))<0.2 & abs(ev(ind,8) - ll1(2))<0.2);
            if comment; disp('More than one event matched ... '); end
            [ar,az]=distance(GA.ev(ind,8),GA.ev(ind,7),ll1(2),ll1(1));
            km = deg2km(ar);
            numMore1 = numMore1+1;
            load_indMore1(numMore1,:) = {{[iev otdn1 ll1]} {[ind otdnGA(ind) km]}}; % [load_indMore1; iev];
        end
    end
    disp([''])
    % The differences and similarities between GG and GA is described by the following variables:
    % load_ind=[]; num=0; % Accpeted events: 1 2 3 4 below
    disp(['num']); num
    % 1. load_ind5=[]; num5=0;                          % Event in GG with no similar ot (+/-5s) in GA, outside Australia
    % 2. load_ind5AUS=[]; num5AUS=0;                    % Event in GG with no similar ot (+/-5s) in GA, inside Australia
    disp(['num5']); num5
    disp(['num5AUS']); num5AUS
    % 3. load_ind1LargeDist   = []; num1LargeDist   =0; % Event in GG found 1 similar ot (+/-5s) in GA, Distance > mindist, outside Australia
    % 4. load_ind1LargeDistAUS= []; num1LargeDistAUS=0; % Event in GG found 1 similar ot (+/-5s) in GA, Distance > mindist, inside Australia
    disp(['num1LargeDist']); num1LargeDist
    disp(['num1LargeDistAUS']); num1LargeDistAUS
    % 5. load_ind1LessDist    = []; num1LessDist    =0; % Event in GG found 1 similar ot (+/-5s) in GA, Distance < mindist, outside Australia
    % 6. load_ind1LessDistAUS = []; num1LessDistAUS =0; % Event in GG found 1 similar ot (+/-5s) in GA, Distance < mindist, inside Australia
    disp(['num1LessDist']); num1LessDist
    disp(['num1LessDistAUS']); num1LessDistAUS
    % 7. load_indMore1        = []; numMore1        =0; % Event in GG found > 1 similar ot (+/-5s) in GA
    disp(['numMore1']); numMore1
        % Replace all GA's hypo with similar GG's hypo
    %                         1    2    34  5     6         7
    % load_ind1LessDistAUS = [iev otdn1 ll1 ind otdnGA(ind) km]
    % ind (5) is GA
    % iev (1) is GG    
    GA.ev(load_ind1LessDistAUS(:,5),[1:9,end]) = GG.ev(load_ind1LessDistAUS(:,1),[1:9,end]);
end

% Filter GG_GA for azimuth gap < 180 ...
disp(['Number of events in GA databse:']);
length(GA.ev(:,1))
indGAazgap = find(GA.ev(:,end-1) < azgmin);
disp(['Number of events with az gap < 180 in GA databse:']);
length(indGAazgap)
GA.ev    = [GA.ev(indGAazgap,:)];
GA.picks = GA.picks(indGAazgap);

%% Prepare for the large loop:
otdn  = datenum(ev(:,1:6));

load_ind=[]; num=0;
load_ind5=[]; num5=0;
load_ind5AUS=[]; num5AUS=0;

load_ind1LargeDist   = []; num1LargeDist   =0;
load_ind1LargeDistAUS= []; num1LargeDistAUS=0;

load_ind1LessDist    = []; num1LessDist    =0;
load_ind1LessDistAUS = []; num1LessDistAUS =0;

%load_indMore1        = [];
numMore1        =0;

otdnGA = datenum(GA.ev(:,1:6));

for iev = 1:length(otdnGA)
    % Origin time and location of one event in the GA catalogue
    otdn1 = otdnGA(iev); % Origin time for one event in GA catalogue
    ll1   = GA.ev(iev,7:8); % Lon, lat for one event in GA catalogue
    % Determines whether the earthquake is inside Australia or not ...
    GAAUS = ( ll1(1)>lonAUS1 & ll1(1)<lonAUS2 & ll1(2)>latAUS1 & ll1(2)<latAUS2 );
    %{
%         if GAAUS
%             % This event is in AUstralia, there is no need for further
%             % tests.. just save it ...
%             num=num+1;
%             evn(num,:) = GA.ev(iev,:);
%             load_ind = [load_ind; iev];
%         else
%             % This event is outside Australia ...
%             ind   = find(abs(otdn - otdn1)<10*one_sec);
%             if length(ind) < 1
%                 disp(['IND<1,AUS=0; +/-05 s => no event, +/-10 s => no event']);
%             else
%                 disp(['IND>1,AUS=0: +/-05 s => no event, +/-10 s => ',num2str(length(ind)),' event(s) ']);
%             end
%             keyboard
%         end
    %}
    % Find the GA event in the compiled one: ev
    ind   = find(abs(otdn - otdn1)<5*one_sec); % 3 options: >1, =1, <1
    if length(ind) < 1 % you found a new event. Add this to the data ...
        %{
% This would have been the case if we were excluding all ISC events from Australia, But we are not.
% But first check if it is inside Australia ...
        AUS = ( ll1(1)>lon1 & ll1(1)<lon2 & ll1(2)>lat1 & ll1(2)<lat2 );
        if AUS
            % This event is in AUstralia, there is no need for further
            % tests.. just save it ...
            num=num+1;
            evn(num,:) = GA.ev(iev,:);
            load_ind = [load_ind; iev];
        else
            % This event is outside Australia ...
            ind   = find(abs(otdn - otdn1)<10*one_sec);
            if length(ind) < 1
                disp(['IND<1,AUS=0; +/-05 s => no event, +/-10 s => no event']);
            else
                disp(['IND>1,AUS=0: +/-05 s => no event, +/-10 s => ',num2str(length(ind)),' event(s) ']);
            end
            keyboard
        end
        %}
        %{
        % No need for 10 second search ...
%         ind   = find(abs(otdn - otdn1)<10*one_sec);
%         if length(ind) < 1
%             if GAAUS
%                 disp(['IND<1; checked +/-05 s, checked +/-10 s, This must be a new event and it is inside Australia ... ']);
%                 num=num+1;
%                 evn(num,:) = GA.ev(iev,:);
%                 load_ind   = [load_ind; iev];
%                 load_ind5  = [load_ind5; iev];
%             else
%                 disp(['IND<1; checked +/-05 s, checked +/-10 s, This must be a new event and it is outside Australia ... ']);
%                 num=num+1;
%                 evn(num,:) = GA.ev(iev,:);
%                 load_ind   = [load_ind; iev];
%                 load_ind5AUS= [load_ind5AUS; iev];
%             end
%
%         else
%             disp(['IND>1; checked +/-05 s, For +/-10 s => ',num2str(length(ind)),' event(s) was found']);
%             if ind == 1
%                 indll20 = find(abs(ev(ind,7) - ll1(1))<0.2 & abs(ev(ind,8) - ll1(2))<0.2);
%                 if length(indll20) < 1
%                     disp(['indll20<1; The search within the 0.2 degree of the lon-lat of the event, finds nothing! ']);
%                     disp(['This is a new event ... ']);
%                     num=num+1;
%                     evn(num,:) = GA.ev(iev,:);
%                     load_ind   = [load_ind; iev];
%                 elseif length(indll20) == 1
%                     disp(['IND=1; The search for 0.2 degree around the lon-lat of the event, reveals Exactly one event ']);
%                 elseif length(indll20) > 1
%                     disp(['IND=1; The search for 0.2 degree around the lon-lat of the event, reveals more than one event! ']);
%                 end
%             elseif ind>1
%
%             end
%         end
        %}
        if GAAUS
            if comment; disp(['IND<1; checked +/-05 s, This must be a new event and it is inside Australia ... ']); end
            num=num+1;
            evn(num,:)     = GA.ev(iev,:);
            load_ind(num,:)           = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind; iev];
            num5AUS=num5AUS+1;
            load_ind5AUS(num5AUS,:)   = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind5AUS; iev];
        else
            if comment; disp(['IND<1; checked +/-05 s, This must be a new event and it is outside Australia ... ']); end
            num=num+1;
            evn(num,:) = GA.ev(iev,:);
            load_ind(num,:)           = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind; iev];
            num5=num5+1;
            load_ind5(num5,:)         = [iev otdn1 ll1 NaN NaN NaN]; % [load_ind5; iev];
        end
    elseif length(ind) == 1 % There is one event with similar time on ISC catalogue
        % Check the location
        %               lon1       lon2              lat1       lat2
        %ll20 = (abs(ev(ind,7) - ll1(1))<0.2 & abs(ev(ind,8) - ll1(2))<0.2);
        [ar,az]=distance(ev(ind,8),ev(ind,7),ll1(2),ll1(1));
        km = deg2km(ar);
        if km > distmin
            if GAAUS
                if comment; disp(['IND=1; But the distance is larger than: ',num2str(distmin),' This could be a new event inside Australia ... ']); end
                num=num+1;
                evn(num,:) = GA.ev(iev,:);
                load_ind(num,:)                           = [iev otdn1 ll1 ind otdn(ind) km]; % [load_ind; iev];
                num1LargeDistAUS=num1LargeDistAUS+1;
                load_ind1LargeDistAUS(num1LargeDistAUS,:) = [iev otdn1 ll1 ind otdn(ind) km]; % [load_ind1LargeDistAUS; iev];
            else
                if comment; disp(['IND=1; But the distance is larger than: ',num2str(distmin),' This could be a new event outside Australia ... ']); end
                num1LargeDist=num1LargeDist+1;
                load_ind1LargeDist(num1LargeDist,:)       = [iev otdn1 ll1 ind otdn(ind) km]; % [load_ind1LargeDist; iev];
            end
        else
            if GAAUS
                if comment; disp(['IND=1; and the distance is less than: ',num2str(distmin),' This is not a new event inside Australia ... ']); end
                num1LessDistAUS=num1LessDistAUS+1;
                load_ind1LessDistAUS(num1LessDistAUS,:)   = [iev otdn1 ll1 ind otdn(ind) km]; % [load_ind1LessDistAUS; iev];
            else
                if comment; disp(['IND=1; and the distance is less than: ',num2str(distmin),' This is not a new event outside Australia ... ']); end
                num1LessDist=num1LessDist+1;
                load_ind1LessDist(num1LessDist,:)         = [iev otdn1 ll1 ind otdn(ind) km]; % [load_ind1LessDist; iev];
            end
        end
    elseif length(ind) > 1
        % Several events were found ...
        %ll20 = find(abs(ev(ind,7) - ll1(1))<0.2 & abs(ev(ind,8) - ll1(2))<0.2);
        if comment; disp('More than one event matched ... '); end
        [ar,az]=distance(ev(ind,8),ev(ind,7),ll1(2),ll1(1));
        km = deg2km(ar);
        numMore1 = numMore1+1;
        load_indMore1(numMore1,:) = {{[iev otdn1 ll1]} {[ind otdn(ind) km]}}; % [load_indMore1; iev];
    end
end

%% The differences and similarities between GA and ISC+EHB is described by the following variables:
% load_ind=[]; num=0; % Accpeted events: 1 2 3 4 below
disp(['num']); num
% 1. load_ind5=[]; num5=0;                          % Event in GA with no similar ot (+/-5s) in EHB+ISC, outside Australia
% 2. load_ind5AUS=[]; num5AUS=0;                    % Event in GA with no similar ot (+/-5s) in EHB+ISC, inside Australia
disp(['num5']); num5
disp(['num5AUS']); num5AUS
% 3. load_ind1LargeDist   = []; num1LargeDist   =0; % Event in GA found 1 similar ot (+/-5s) in EHB+ISC, Distance > mindist, outside Australia
% 4. load_ind1LargeDistAUS= []; num1LargeDistAUS=0; % Event in GA found 1 similar ot (+/-5s) in EHB+ISC, Distance > mindist, inside Australia
disp(['num1LargeDist']); num1LargeDist
disp(['num1LargeDistAUS']); num1LargeDistAUS
% 5. load_ind1LessDist    = []; num1LessDist    =0; % Event in GA found 1 similar ot (+/-5s) in EHB+ISC, Distance < mindist, outside Australia
% 6. load_ind1LessDistAUS = []; num1LessDistAUS =0; % Event in GA found 1 similar ot (+/-5s) in EHB+ISC, Distance < mindist, inside Australia
disp(['num1LessDist']); num1LessDist
disp(['num1LessDistAUS']); num1LessDistAUS
% 7. load_indMore1        = []; numMore1        =0; % Event in GA found > 1 similar ot (+/-5s) in EHB+ISC
disp(['numMore1']); numMore1
% Manual analysis of the 3, 4, 5, 6 index:
% Number 3: % Event in GA found 1 similar ot (+/-5s) in EHB+ISC, Distance > mindist, outside Australia

%% ADD evn and new picks to ev and picks and filter the aftershocks:
ev = [ev; evn]; % The length must be replaced by azimuth later ...
clear evn;
table1= struct2table(picks); % Convert the EHB picks into table
picksGA = GA.picks(load_ind(:,1)); % Take the events you want from GA ...
% Add few fields
% picksGA(1).ch=[]; % 7
picksGA(1).lc=[]; % 3
picksGA(1).lo=[]; % 5
picksGA(1).la=[]; % 6
picksGA(1).el=[]; % 7
% Convert to table and add them together ...
table2= struct2table(picksGA); clear GA; % Convert the picks of selected ISC events into table
table3= [table1; table2]; clear table1 table2; % Merge them
picks = table2struct(table3); clear table3; % Convert back to structure

%% Visualize
%
f5 = figure; f5.Color = 'w';
% axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 180]);
% load coastlines;
% patchm(coastlat,coastlon,'w','FaceColor','none');
axesm('MapProjection','hammer','Grid','on','Frame','on','Origin',[0 110]);

% Geoshow can map things
%hst = geoshow(stlat, stlon, 'DisplayType', 'Point', 'Marker', 'd', 'Color', 'red','MarkerFaceColor', 'red', 'MarkerSize', 4);
%hev = geoshow(evlat, evlon, 'DisplayType', 'Point', 'Marker', '.', 'Color', 'blue','MarkerFaceColor', 'blue', 'MarkerSize', 1);
% Scatterm is much better
a=5;
scatterm(ev(:,8),ev(:,7),a,'r','filled'); %cb=colorbar; cb.Label.String = 'Earthquake Depth [km]';
load coast; % Matlab 2014b
patchm(lat,long,'k','FaceColor','none');

disp(['Number of events after magnitude and azimuthal gap filtering: ']);
length(ev(:,1))
% hst = scatterm(stlat,stlon,a,      'filled','MarkerFaceColor','red','Marker','d','SizeData',20);

% Sort the catalogue
otdum1 = datenum(ev(:,1:6));
[B,I] = sortrows(otdum1);

ev = ev(I,:);
picks = picks(I); % pick information sorted

% make sure the number of phases in 'ev' is the same as the number phases in 'picks'.  
for i=1:length(ev(:,1))
    ev(i,10) = length(picks(i).dist);
end

% Save depend on name
if GAhypoGG
    CATGG.ev    = ev;
    CATGG.picks = picks;
    save CATGG CATGG;
else
    CAT.ev      = ev;
    CAT.picks   = picks;
    save CAT CAT;
end

if GAhypoGG
    fid=fopen('CATGG.csv','w+');
else
    fid=fopen('CAT.csv','w+');
end
l = 0;
for i=1:length(ev(:,1))
    % Write the event info
    %                        1       2       3        4      5       6        7     8       9       10     11     12     13     14     15     16
    %                       %YY      MM      DD      hh      mm      ss      lo     la      de      ph     mb,    ms     ml,    mw,    ID     GAP    EVnumber
    if ev(i,17) == 1
        fprintf(fid,'#EHB, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 2
        fprintf(fid,'#ISC, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 3
        fprintf(fid,'#USGS, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 4
        fprintf(fid,'#GG, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 5
        fprintf(fid,'#GA, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    end
    l = l+1; % line number
    mod1000 = mod(l,1000);
    if mod1000 == 1
        disp(['Event number: ',num2str(i)])
    end
    if isempty(picks(i).at)
        
    else
        for j=1:length(picks(i).at(:,1))
            fprintf(fid,'%s, ',picks(i).st(j,:)); % 1. stations
            if isempty(picks(i).ch)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).ch(j,:)); % 2. channel
            end
            if isempty(picks(i).lc)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).lc(j,:)); % 3. location code 00, 10, ...
            end
            if isempty(picks(i).ne)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).ne(j,:)); % 4. network
            end
            if isempty(picks(i).lo)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).lo(j)); % 5. lon
            end
            if isempty(picks(i).la)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).la(j)); % 6. lat
            end
            if isempty(picks(i).el)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).el(j)); % 7. ele
            end
            fprintf(fid,'%s, ',picks(i).ph(j,:)); % 8. phase
            fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',picks(i).at(j,:)); % 9. 10. 11. 12. 13. 14. Arrival Time
            fprintf(fid,'%7.3f \n',picks(i).dist(j)); % 15. distance
            l=l+1;
        end
    end
end
fclose(fid);