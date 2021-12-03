%% Download ISC data in monthly file for all EQs with Magntiude > 3
clear all; close all; clc;
% first text
t1='http://www.isc.ac.uk/cgi-bin/web-db-v4?iscreview=&out_format=CSV&ttime=on&ttres=on&tdef=on&phaselist=&sta_list=&stnsearch=GLOBAL&stn_ctr_lat=&stn_ctr_lon=&stn_radius=&max_stn_dist_units=deg&stn_top_lat=&stn_bot_lat=&stn_left_lon=&stn_right_lon=&stn_srn=&stn_grn=&bot_lat=&top_lat=&left_lon=&right_lon=&ctr_lat=&ctr_lon=&radius=&max_dist_units=deg&searchshape=GLOBAL&srn=&grn=&start_year=';
% start year here ... 
% Then second text
t2='&start_month=';
% start month here ...
t3='&start_day=01&start_time=00%3A00%3A00&end_year=';
% end year here ...
t4='&end_month=';
% end moths here ...
t5='&end_day=01&end_time=00%3A00%3A00&min_dep=&max_dep=&min_mag=3&max_mag=&req_mag_type=Any&req_mag_agcy=ISC&include_links=on&request=STNARRIVALS';
% http://www.isc.ac.uk/cgi-bin/web-db-v4?iscreview=&out_format=CSV&ttime=on&ttres=on&tdef=on&phaselist=&sta_list=&stnsearch=GLOBAL&stn_ctr_lat=&stn_ctr_lon=&stn_radius=&max_stn_dist_units=deg&stn_top_lat=&stn_bot_lat=&stn_left_lon=&stn_right_lon=&stn_srn=&stn_grn=&bot_lat=&top_lat=&left_lon=&right_lon=&ctr_lat=&ctr_lon=&radius=&max_dist_units=deg&searchshape=GLOBAL&srn=&grn=&start_year=2015&start_month=10&start_day=01&start_time=00%3A00%3A00&end_year=2016&end_month=11&end_day=01&end_time=00%3A00%3A00&min_dep=&max_dep=&min_mag=4&max_mag=&req_mag_type=Any&req_mag_agcy=ISC&include_links=on&request=STNARRIVALS
sy = 'start_year=';
sm = 'start_month=';
ey = 'end_year=';
em = 'end_month=';
syind = strfind(text,sy);
smind = strfind(text,sm);
eyind = strfind(text,ey);
emind = strfind(text,em);

%% Output 
file = 'Download_ISC_monthly_mag3';
fid=fopen(file,'w+');
years = 2016:-1:1970;
for iy = years
    systr=num2str(iy);
    if iy == 2016
        months = 1:3;
    else
        months = 1:11; 
    end
    for im = months
        smstr1=num2str(im);
        smstr2=num2str(im,'%02.0f');
        if im == 11
            eystr=num2str(iy+1);
            emstr1=num2str(1);
            emstr2=num2str(1,'%02.0f');
        else
            eystr=num2str(iy);
            emstr1=num2str(im+1);
            emstr2=num2str(im+1,'%02.0f');
        end
        disp(['Building a text for: ',systr,'/',smstr2,' ',eystr,'/',emstr2])
        text = [t1,systr,t2,smstr1,t3,eystr,t4,emstr1,t5];
        name = ['Arricals_ISC_',systr,'-',smstr2,'_',eystr,'-',emstr2];
        fprintf(fid,'%s\n',['curl -L -o ',name,' "',text,'"']);
    end
end
fclose(fid);