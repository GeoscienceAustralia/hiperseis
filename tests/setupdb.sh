#!/usr/bin/env bash
DBFLAG="mysql://sysop:sysop@localhost/seiscomp3"
SQLITEINIT=${SEISCOMP_ROOT}/share/db/sqlite3.sql
PBDB=$1
INVENTORY=$2
CONFIG=$3

if [[ $# -lt 3 ]]; then
    echo $0: usage: ./setupdb.sh pbdb inventory config
    exit 1
fi

echo "Retrieving inventory: $INVENTORY"
scxmldump -f -I -d ${DBFLAG}  > ${INVENTORY}
echo "Retrieving configuration: $CONFIG"
scxmldump -f -C -d ${DBFLAG}  > ${CONFIG}
echo "Initializing sqlite database: $PBDB"
if [ -f "${PBDB}" ]; then
    rm "${PBDB}"
fi
sqlite3 -batch -init ${SQLITEINIT} ${PBDB} .exit
echo "Populating sqlite database: $PBDB"
scdb --plugins dbsqlite3 -d sqlite3://${PBDB} -i ${INVENTORY}
scdb --plugins dbsqlite3 -d sqlite3://${PBDB} -i ${CONFIG}
