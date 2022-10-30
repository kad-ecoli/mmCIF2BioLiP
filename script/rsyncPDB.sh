#!/bin/sh
FILE=`readlink -e $0`
BINDIR=`dirname $FILE`
ROOTDIR=`dirname $BINDIR`

MIRRORDIR=$ROOTDIR/pdb                 # your top level rsync directory
LOGFILE=$BINDIR/logs               # file for storing logs
RSYNC=rsync                             # location of local rsync

##########################################################################################
#
#        YOU MUST UNCOMMENT YOUR CHOICE OF SERVER AND CORRESPONDING PORT BELOW
#
SERVER=rsync.wwpdb.org::ftp                                   # RCSB PDB server name
PORT=33444                                                    # port RCSB PDB server is using
#
#SERVER=rsync.ebi.ac.uk::pub/databases/rcsb/pdb-remediated     # PDBe server name
#PORT=873                                                      # port PDBe server is using
#
#SERVER=pdb.protein.osaka-u.ac.jp::ftp                         # PDBj server name
#PORT=873                                                      # port PDBj server is using
#
##########################################################################################


mkdir -p $MIRRORDIR/data/structures/divided/mmCIF
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/data/structures/divided/mmCIF/ $MIRRORDIR/data/structures/divided/mmCIF > $LOGFILE 2>/dev/null

mkdir -p $MIRRORDIR/derived_data/index
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/derived_data/index/resolu.idx $MIRRORDIR/derived_data/index/ >> $LOGFILE 2>/dev/null
${RSYNC} -rlpt -v -z --delete --port=$PORT ${SERVER}/derived_data/index/compound.idx $MIRRORDIR/derived_data/index/ >> $LOGFILE 2>/dev/null
