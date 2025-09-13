#!/bin/bash
#
# This script runs all test cases in quad_test3-triangle.c
#
# $Id: quad_test3-triangle-testall.sh,v 1.14 2021/02/08 08:14:51 zlb Exp $

test -r ~/.profile && . ~/.profile

test -z "$MAX_LEVEL" && MAX_LEVEL=5	# max level
test -z "$TOL" && TOL=1e-10		# tolerance

set -e

(cd ../src; make -s lib)

for ((case=1; case<100; case++)); do
    rm -f quad_test3-triangle quad_test3-triangle.o
    make -s USER_CFLAGS="-DTEST_CASE=$case" quad_test3-triangle >/dev/null 2>&1
    for ((level=0; level<=MAX_LEVEL; level++)); do
	for ((type=-1; type<=1; type++)); do
	    tmp="none dot cross"
	    test "$type" -ne "0" && tmp="none"
	    for proj in $tmp; do
		###echo 1>&2 "Running \"./quad_test3-triangle -type $type -proj $proj -qi_subdiv_level0 $level $@\""
		res=`./quad_test3-triangle \
			-type $type -proj $proj \
			-qi_subdiv_level0 $level "$@" 2>/dev/null \
			| grep error || exit 1`
		error="`echo $res | awk '{print \$7}'`"
		check="`echo "$error < $TOL" \
			| sed -e 's/e[+-]0* / /g' -e 's/e/*10^/g' | bc -l`"
		# Term colors: red=31, green=32, yellow=33, blue=34, restore=0
		test "$check" -eq 0 && /bin/echo -n "[01;31m"
		/bin/echo -n "Case $case, subdiv $level, "
		echo $type $proj $error \
			| awk '{printf "type %2d, proj %5s, error %0.3e", \
				 $1, $2, $3}'
		if test "$check" -eq 0; then
		    echo $TOL "[01;0m" \
			| awk '{printf " [failed (tol=%0.3g)%s]", $1, $2}'
		fi
		echo ""
	    done
	done
    done
done
