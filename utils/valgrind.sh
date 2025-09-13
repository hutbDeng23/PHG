#!/bin/bash
#
# This script executes a prog under valgrind
#
# $Id: valgrind.sh,v 1.19 2021/07/13 02:15:50 zlb Exp $

export PHG_OPTIONS="$PHG_OPTIONS -mem_sampling_rate=0 +use_perfctr"

if test $# -lt 2; then
    echo "Usage: $0 nprocs prog arguments"
    exit 1
fi

. `dirname $0`/functions

prog=`canonicalize $2`
nprocs=$1
shift 2

supp=`dirname $0`/valgrind.supp
suppressions="--gen-suppressions=no --suppressions=`canonicalize $supp`"
curdir=`pwd`
out="valgrind-output"

cat <<END >wrapper_valgrind.sh
#!/bin/bash
#
export GLIBCPP_FORCE_NEW=1	# gcc >= 3.2.2
export GLIBCXX_FORCE_NEW=1	# gcc >= 3.4
if ! type valgrind >/dev/null 2>&1; then
    echo 1>&2 "Cannot execute valgrind, abort.!"
    exit 1
fi
if valgrind --tool=memcheck true >/dev/null 2>&1; then
    # new version
    opts="--tool=memcheck"
else
    # old version
    opts=""
fi
logfile="${curdir}/${out}.\`hostname -s\`.\$\$"
if valgrind \$opts --log-file="\$logfile" true >/dev/null 2>&1; then
    # new version
    opts="\$opts --log-file=\$logfile"
else
    # old version
    opts="\$opts --logfile=\$logfile"
fi
/bin/rm -f "\$logfile".* "\$logfile"
yes 2>/dev/null | valgrind \$opts --leak-check=yes --show-reachable=yes \
	--run-libc-freeres=yes --error-limit=no $suppressions $prog "\$@"
exit \$?
END
chmod a+x wrapper_valgrind.sh

if false; then
    # for MPICH-1 using P4 group file
    echo "localhost 0 `pwd`/wrapper_valgrind.sh" >p4file
    i=1
    while [ $i -lt $nprocs ]; do
	echo "localhost 1 `pwd`/wrapper_valgrind.sh" >>p4file
	i=$((i+1))
    done
    trap "/bin/rm -f p4file wrapper_valgrind.sh" 0 1 2 3 15
    /bin/rm -f valgrind[0-9]*.*[0-9]*
    if test $nprocs -eq 1; then
	bash wrapper_valgrind.sh "$@"
    else
	bash wrapper_valgrind.sh -p4pg p4file "$@"
    fi
else
    # use mpirun
    trap "/bin/rm -f wrapper_valgrind.sh" 0 1 2 3 15
    /bin/rm -f ${out}.*[0-9]*
    if type mpijob >/dev/null 2>&1; then
	mpijob -np $nprocs ./wrapper_valgrind.sh "$@"
    elif type mpirun >/dev/null 2>&1; then
	mpirun -np $nprocs ./wrapper_valgrind.sh "$@"
    elif type mpiexec >/dev/null 2>&1; then
	mpiexec -n $nprocs ./wrapper_valgrind.sh "$@"
    else
	bash wrapper_valgrind.sh "$@"
    fi
fi
