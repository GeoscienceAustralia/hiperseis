mkdir TEMP

# cd into the TEMP dir and do the following
cd TEMP
wget ftp://hazards.cr.usgs.gov/CWBQuery/EdgeCWBRelease.tar.gz

tar -zxvf EdgeCWBRelease.tar.gz
tar -xvf bin_release.tar
tar -xvf scripts_release.tar

mv bin ~vdl
chmod g+rw Jars/*.jar
cp -p Jars/*.jar ~vdl/bin/
chmod -R 755 ~vdl/bin
chmod g+rw ~vdl/bin/*.jar

cp -r scripts ~vdl/

cp .bash* ~vdl/

cd ~vdl

# make sure java class lib has jars
ln -s  /home/vdl/bin bin/lib

# the scripts/installCWBRelease.bash has issues. as a reference only
dirs="bin log LOG EDGEMOM config SEEDLINK EW log/MSLOG log/CWBLOG log/q330";
# create any missing directories
for dir in $dirs; do
if [ ! -e ${dir} ]; then
        echo "Create ${dir} directory"
        mkdir ${dir}
fi
done

tar -xvf TEMP/dbconn.tar
#.dbconn/dbconn.conf
#.dbconn/dbconn_mysql_localhost_3306.conf
#.dbconn/dbconn_mysql_localhost.conf

# creat databases - require type in mysql user info: root, passwd, y
cd scripts/INSTALL/
expect -c "spawn bash ./makedb.bash; expect -re \"+ read acct\"; send \"\n\"; expect -re \"+ read -s rpw\"; send \"CWBr00tpwd\"; expect -re \"config(y/n)?\"; send \"y\";"
#bash ./makedb.bash   # if run success, this will create a few db

expect -c "spawn bash createketchum.bash; expect -re \"Enter password:\"; send \"CWBr00tpwd\";"
# type the mysqlroot passpwd

expect =c "spawn bash createvdlro.bash; expect -re \"Enter password:\"; send \"CWBr00tpwd\";"
# type the mysqlroot passpwd

cd ~vdl/

echo "
#/home/ketchum/anss.prop by Anss via Util.saveProperties() cp=/home/vdl/bin/Anss.jar
#Fri Apr 02 21:43:54 GMT+00:00 2010
PrinterCommand=lpr -P shaky2
SMTPServer=136.177.7.24
PrinterFile=anss.lpt
StatusServer=localhost
DBServer=localhost/3306:\anss\:mysql\:anss
StartY=112.0
StartX=0.0
RTSServer=localhost
HostIP=136.177.36.20
SessionFile=SESSION.OUT" >> anss.prop

echo "
#chandisp.prop /home/vdl
#Mon Jul 28 21:04:00 GMT+00:00 2008
PrinterFile=anss.lpt
SessionFile=SESSION.OUT
Comp=[BEHS]HZ Only
Sort=0
RTSServer=localhost
PrinterCommand=lpr -P shaky2
StartECY=12.0
StartECX=10.0
HostIP=localhost
SNWGroupID=0
StatusServer=localhost
MySQLServer=localhost
DBServer=localhost/3306\:anss\:mysql\:anss
SoundOff=false
StartY=25.0
StartX=91.0
Refresh=1
popupY=24
popupX=47
SMTPServer=136.177.7.24
AlertMode=false" >> chandisp.prop

echo "
#/home/ketchum/edgecon.prop by EdgeConfig via Util.saveProperties() cp=/home/vdl/bin/EdgeConfig.jar
#Wed Dec 26 21:45:00 GMT+00:00 2012
AllowInstance=y
PrinterCommand=print
DBServer=localhost/3306\:edge\:mysql\:edge
SMTPServer=136.177.7.24
PrinterFile=edgecon.lpt
MetaDBServer=localhost/3306\:metadata\:mysql\:metadata
MySQLServer=localhost
StatusDBServer=localhost/3306\:status\:mysql\:status
StartY=0.0
StartX=383.0
SSLOff=true
SessionFile=SESSION.OUT" >> edgecon.prop

echo "
#Tue Aug 08 03:59:03 UTC 2017
#
#/home/vdl/edge_${HOSTNAME}.prop
Node=2
ndatapath=1
nday=100
datapath=/datas/2015/
DBServer=localhost/3306:edge:mysql:edge
MetaDBServer=localhost/3306:metadata:mysql:metadata
StatusDBServer=localhost/3306:status:mysql:status
StatusServer=localhost
logfilepath=log/
SMTPServer=mailx
emailTO=anss.alarm@usgs.gov
SMTPFrom=${HOSTNAME}-vdl@usgs.gov
daysize=10000
extendsize=2000
ebqsize=2000
AlarmIP=localhost
instanceconfig=true" >> edge_$HOSTNAME.prop

echo "
#Tue Aug 08 03:59:03 UTC 2017
#
#/home/vdl/edge_${HOSTNAME}.prop
Node=2
ndatapath=1
nday=100
datapath=/datas/2015/
DBServer=localhost/3306:edge:mysql:edge
MetaDBServer=localhost/3306:metadata:mysql:metadata
StatusDBServer=localhost/3306:status:mysql:status
StatusServer=localhost
logfilepath=log/
SMTPServer=mailx
emailTO=anss.alarm@usgs.gov
SMTPFrom=${HOSTNAME}-vdl@usgs.gov
daysize=10000
extendsize=2000
ebqsize=2000
AlarmIP=localhost
instanceconfig=true" >> edge.prop

echo "
StartX=501.0
PrinterFile=edgecon.lpt
MySQLServer=localhost
DBServer=localhost/3306\:edge\:mysql\:edge
MetaDBServer=localhost/3306\:metadata\:mysql\:metadata
StatusDBServer=localhost/3306\:status\:mysql\:status
SMTPServer=136.177.7.24
HostIP=136.177.30.35
SessionFile=SESSION.OUT
PrinterCommand=print
StartY=233.0" >> metagui.prop

touch msread.prop

echo "
cwbip=localhost
metadataserver=137.227.230.1
MySQLServer=localhost
DBServer=localhost/3306\:edge\:mysql\:edge
queryport=2061" >> query.prop

echo "
ebqsize=20000
KetchumOverride=
StatusServer=localhost
MySQLServer=localhost
SMTPFrom=ip-172-31-21-121-vdl<ketchum@usgs.gov>
emailTo=niket.chhajed@ga.gov.au
extendsize=4000
logfilepath=log/
AlarmIP=localhost
ndatapath=1
nday=100
datapath=/datas/2015/" >> queryserver.prop

echo "
#
# Edgemomsetup file for - 2#0
#
#
# Base edgemom thread, must be in every config
 Mom:EdgeMom:-empty >>edgemom
# For loading historical archived data (usually from DMC) only from this host
DMCLoad:MiniSeedServer:-p 7900-noudpchan -nohydra >>dmcload
# Maintain channel information for access by threads
echn:EdgeChannelServer:-empty >>echn
#Start Output infrastructure to make rings for RRP, LISS, EWExport
#OInf:OutputInfrastructure:-dbg >>oi" >> "EDGEMOM/edgemom_2#0.setup"

echo "
#Written for querymominstance 2Q0
 Mom:EdgeMom:-empty
# The EdgeChannelServer is used for configuration parameters from channel table (FilterPicker)
 Echn:EdgeChannelServer:-empty >>echnqm
QS:gov.usgs.cwbquery.EdgeQueryServer:-p 2061 -mdsip 137.227.224.97 -mdsport 2052 -allowrestricted >>queryserver" >> "EDGEMOM/querymom_2Q0.setup"

echo "* * * * * bash scripts/chkCWB alarm ^ 350 -alarm -nocmd -snw -noudpchan >>LOG/chkcwb.log 2>&1" > crn.tmp
echo "* * * * * bash scripts/chkCWB edgemom 2#0 350 -f EDGEMOM/edgemom_2#0.setup >>LOG/chkcwb.log 2>&1" >> crn.tmp
echo "* * * * * bash scripts/chkCWB querymom 2Q0 450 -max -f EDGEMOM/querymom_2Q0.setup >>LOG/chkcwb.log 2>&1" >> crn.tmp

crontab crn.tmp


