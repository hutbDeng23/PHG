#!/bin/bash
#
# $Id: opendx.sh,v 1.3 2006/02/17 07:14:26 zlb Exp $

if [ $# -ne 1 ]; then
    echo 1>&2 "Usage: $0 DX_filename"
    exit 1
fi

if ! test -r "$1"; then
    echo 1>&2 "Error: cannot open DX file \"$1\"!"
    exit 1
fi

. `dirname $0`/functions

dxfile=`canonicalize "$1"`
netfile=`canonicalize "${0%\.sh}.net"`
	echo script=$0 netfile=$netfile
tmpfile=/tmp/opendx-$$.net
export DXDATA=`dirname $dxfile`
trap "/bin/rm -f $tmpfile" 0 1 2 3 15

(
    oldname=`grep '^main_FileSelector_.*' "$netfile" \
		| sed -e s'/.*\"\(.*\)\";/\1/g'`
    sed -e 's!\"'$oldname'\"!\"'$dxfile'\"!' \
	-e 's!\"'`basename $oldname`'\"!\"'`basename $dxfile`'\"!' "$netfile"
#-------------------------
###    cat
###) | dx -script
#-------------------------
) > $tmpfile
dx -image $tmpfile
#-------------------------
