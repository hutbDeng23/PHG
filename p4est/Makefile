# $Id: Makefile,v 1.18 2022/07/30 11:01:29 zlb Exp $
#
# For P4EST make variables like P4EST_CPPFLAGS, P4EST_LDFLAGS and P4EST_LIBS
#	include $(strip ${P4EST_DIR})/etc/Makefile.p4est.mk
# with P4EST_DIR denoting the p4est installation directory.

include ../Makefile.inc
include Makefile.dep

CPPFLAGS	:= ${CPPFLAGS} ${BUILD_CPPFLAGS} -DPHG_TO_P4EST
CFLAGS		:= ${CFLAGS} ${BUILD_CFLAGS} ${USER_CFLAGS}
CXXFLAGS	:= ${CXXFLAGS} ${BUILD_CXXFLAGS} ${USER_CXXFLAGS}

default:
all:

lib:
	(cd ../src; $(MAKE) lib)

clean:
	-/bin/rm -f *.o libphg2p4est.*
distclean: clean

.PHONY: default all lib clean distclean
