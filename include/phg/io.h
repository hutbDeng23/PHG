/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: io.h,v 1.4 2021/12/08 02:48:52 zlb Exp $ */

#ifndef PHG_IO_H

#ifdef __cplusplus
extern "C" {
#endif

BOOLEAN phgIOOpen(MPI_Comm comm, const char *fn);
int phgIOPrintf(MPI_Comm comm, char *fmt, ...);
void phgIOWriteRoot(MPI_Comm comm, void *buf, size_t size, int count);
void phgIOWrite(MPI_Comm comm, void *buffer, size_t size, size_t count,
		size_t nglobal, BTYPE flags[], INT map[], BOOLEAN map_flag);
void phgIOClose(void);
const char *phgIOAddExtension(const char *fn, const char *ext);

void *phgIOAllocBuffer(size_t size);
void phgIOFreeBuffer(void);
void phgIOResetBuffer(void);
void phgIOBigEndianAppend(void *buf, size_t size, int count);

void phgCloseInputFile_(FILE *fp);
FILE *phgOpenInputFile_(const char *filename);

#ifdef __cplusplus
}
#endif

#define PHG_IO_H
#endif
