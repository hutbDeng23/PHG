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

/* Parallel I/O functions
 *
 * $Id: io.c,v 1.22 2021/12/08 02:48:52 zlb Exp $
 */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>	/* PATH_MAX */

static struct {
    MPI_Comm	comm;
    int		rank, nprocs;
} dummy_g, *g = &dummy_g;

#if USE_MPI
# define DUMMY_g	g->comm = comm; MPI_Comm_rank(comm, &g->rank); \
			MPI_Comm_size(comm, &g->nprocs);
#else	/* USE_MPI */
# define DUMMY_g	g->comm = comm; g->rank = 0; g->nprocs = 1;
#endif

static BOOLEAN use_mpiio = 0;
static BOOLEAN mpiio_atomicity = 0;

#if USE_MPI
static struct {		/* struct for sorting maps (L2Gmap_vert) */
    INT gindex, lindex;
} *map_sort, *ms0, *ms1;

static int
map_sort_compare(const void *m0, const void *m1)
{
    INT i;
    ms0 = (void *)m0; ms1 = (void *)m1;
    return (i = ms0->gindex - ms1->gindex) > 0 ? 1 : (i < 0 ? -1 : 0);
}

#if USE_MPIIO
static int
merge_blocks(int n, int *blocklens, int *indices)
{
    int *blocklens0 = blocklens, *p = indices;

    if (n <= 0)
	return n;

    while (--n > 0) {
	if (p[1] - p[0] == 1) {
	    ++(*blocklens);
	    ++p;
	}
	else {
	    ++blocklens;
	    *(++indices) = *(++p);
	}
    }

    return blocklens - blocklens0 + 1;
}
#endif	/* USE_MPIIO */
#endif	/* USE_MPI */

#if USE_MPIIO

/* Use MPI-2 parallel I/O functions */

static MPI_File fh = MPI_FILE_NULL;

static BOOLEAN
phgIOOpenMPI(MPI_Comm comm, const char *fn)
{
    DUMMY_g

    MPI_File_open(g->comm, (char *)fn,
			MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (fh == MPI_FILE_NULL) {
	return FALSE;
    }
    MPI_File_set_atomicity(fh, (int)mpiio_atomicity);
    MPI_File_set_size(fh, 0);	/* truncate initial file size to 0 */
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    return TRUE;
}

static int
phgIOPrintfMPI(MPI_Comm comm, char *fmt, va_list ap)
/* serial fprintf analog to an MPI file: at first rank 0 writes the data,
   then new File View set to end of file */
{
    char s[1024];
    int ret = 0;
    MPI_Offset offset;
    MPI_Status status;

    DUMMY_g

    if (g->rank == 0) {
	ret = vsprintf(s, fmt, ap);
	MPI_File_write(fh, s, strlen(s), MPI_BYTE, &status);
	MPI_File_get_position(fh, &offset);
    }
    MPI_Bcast(&offset, sizeof(offset), MPI_BYTE, 0,
				g == NULL ? phgComm : g->comm);
#if 0
    MPI_File_sync(fh);
    MPI_Barrier(g->comm);
    MPI_File_sync(fh);
#endif
    MPI_File_seek(fh, offset, MPI_SEEK_SET);
    MPI_File_seek_shared(fh, offset, MPI_SEEK_SET);

    return ret;
}

static void
phgIOWriteRootMPI(MPI_Comm comm, void *buf, size_t size, int count)
{
    MPI_Offset offset;
    MPI_Status status;

    DUMMY_g

    if (g->rank == 0) {
        MPI_File_write(fh, buf, size * count, MPI_BYTE, &status);
        MPI_File_get_position(fh, &offset);
    }
    MPI_Bcast(&offset, sizeof(offset), MPI_BYTE, 0,
				g == NULL ? phgComm : g->comm);
#if 0
    MPI_File_sync(fh);
    MPI_Barrier(g->comm);
    MPI_File_sync(fh);
#endif
    MPI_File_seek(fh, offset, MPI_SEEK_SET);
    MPI_File_seek_shared(fh, offset, MPI_SEEK_SET);

    return;
}

static void
phgIOWriteMPI(MPI_Comm comm, void *buffer, size_t size, size_t count,
	 	size_t nglobal, BTYPE flags[], INT map[], BOOLEAN map_flag)
/* Parallel data write.
 *
 * Each process has 'count' entries of size 'size' stored in 'buffer[]',
 * 'flags[i]' indicates whether entry i is owned by current process,
 * 'map[i]' is the global id of entry i 
 *
 * If flags == NULL then all entries are owned by current process
 * If map == NULL then map[i] == i + (sum of count of prev procs) is assumed.
 *
 * If map_flag == FALSE the flags and map are NULL on all processes,
 * otherwise they are only NULL on processes with count==0.
 */
{
    MPI_Datatype type, type0;
    MPI_Status status;
    MPI_Offset offset, disp;
    INT i, j, k;
    int *blocklens, *indices, jj, kk;

    DUMMY_g

    MPI_File_get_position(fh, &offset);
    MPI_File_get_byte_offset(fh, offset, &disp);
    if (g->nprocs > 1 && map_flag) {
	map_sort = phgAlloc(count * sizeof(*map_sort));
	if (count > 0) {
	    /* sort map[] (blocks must have increasing offsets in filetype) */
	    for (i = 0; i < count; i++) {
		map_sort[i].gindex = map[i];
		map_sort[i].lindex = i;
	    }
	    qsort(map_sort, count, sizeof(*map_sort), map_sort_compare);
	}
	blocklens = phgAlloc(count * sizeof(*blocklens));
	indices = phgAlloc(count * sizeof(*indices));
	j = 0;
	for (i = 0; i < count; i++) {
	    k = map_sort[i].lindex;
	    if (flags != NULL && !(flags[k] & OWNER)) continue;
	    blocklens[j] = 1;
	    indices[j++] = map_sort[i].gindex;
	}
	phgDebug((3, "Before merging contiguous blocks: %d blocks\n", j));
	j = merge_blocks(j, blocklens, indices); /* merge contiguous blocks */
	phgDebug((3, "After merging contiguous blocks:  %d blocks\n", j));
	MPI_Type_contiguous(size, MPI_BYTE, &type0);
	MPI_Type_commit(&type0);
	MPI_Type_indexed(j, blocklens, indices, type0, &type);
	MPI_Type_commit(&type);
	MPI_File_set_view(fh, disp, type0, type, "native", MPI_INFO_NULL);
	MPI_Type_free(&type);
	j = 0;
	for (i = 0; i < count; i++) {
	    k = map_sort[i].lindex;
	    if (flags != NULL && !(flags[k] & OWNER)) continue;
	    blocklens[j] = 1;
	    indices[j++] = k;
	}
	j = merge_blocks(j, blocklens, indices); /* merge contiguous blocks */
	MPI_Type_indexed(j, blocklens, indices, type0, &type);
	MPI_Type_commit(&type);
	MPI_File_write_all(fh, buffer, 1, type, &status);
	MPI_Type_free(&type);
	MPI_Type_free(&type0);
	MPI_File_set_view(fh, disp + nglobal * size, MPI_BYTE, MPI_BYTE,
		"native", MPI_INFO_NULL);
	phgFree(blocklens);
	phgFree(indices);
	phgFree(map_sort);
    }
    else {
#if 0	/* FIXME: using write_ordered gets an "io No such file or directory"
	   warning when closing the file on non root process with MPICH 1.2.7
	   ch_p4 */
	MPI_Type_contiguous(count * size, MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_File_write_ordered(fh, buffer, 1, type, &status);
	MPI_File_seek(fh, offset + nglobal * size, MPI_SEEK_SET);
	MPI_File_seek_shared(fh, offset + nglobal * size, MPI_SEEK_SET);
	MPI_Type_free(&type);
#else
	j = count;
	MPI_Scan(&j, &k, 1, MPI_INT, MPI_SUM, g == NULL ? phgComm : g->comm);
	k -= j;		/* k = location */
	j = 1;		/* j = blocksize */
	MPI_Type_contiguous(size, MPI_BYTE, &type0);
	MPI_Type_commit(&type0);
	jj = j;
	kk = k;
	MPI_Type_indexed(1, &jj, &kk, type0, &type);
	MPI_Type_commit(&type);
	MPI_File_set_view(fh, disp, type0, type, "native", MPI_INFO_NULL);
	MPI_Type_free(&type);
	MPI_File_write_all(fh, buffer, count, type0, &status);
	MPI_Type_free(&type0);
	MPI_File_set_view(fh, disp + nglobal * size, MPI_BYTE, MPI_BYTE,
		"native", MPI_INFO_NULL);
#endif
    }
}

static void
phgIOCloseMPI()
{
    MPI_File_close(&fh);
    fh = MPI_FILE_NULL;
}

#endif	/* USE_MPIIO */

/*------------------------- Sequential I/O functions ------------------------*/

static FILE *fp = NULL;

static BOOLEAN
phgIOOpen_(MPI_Comm comm, const char *fn)
{
    DUMMY_g

    if (g->rank > 0)
	return TRUE;
    return (fp = fopen(fn, "wb")) != NULL;
}

static int
phgIOPrintf_(MPI_Comm comm, char *fmt, va_list ap)
/* serial fprintf analog to an MPI file: at first rank 0 writes the data,
   then new File View set to end of file */
{
    DUMMY_g

    if (g->rank > 0)
	return 0; 

    return vfprintf(fp, fmt, ap);
}

static void
phgIOWriteRoot_(MPI_Comm comm, void *buf, size_t size, int count)
{
    DUMMY_g

    if (g->rank > 0)
        return; 

    if (fwrite(buf, size, count, fp) != count)
	phgError(1, "%s: fwrite error, abort.\n", __func__);
    return;
}

static void
phgIOWrite_(MPI_Comm comm, void *buffer, size_t size, size_t count,
		size_t nglobal, BTYPE flags[], INT map[], BOOLEAN map_flag)
/* Parallel data write.
 *
 * Each process has 'count' entries of size 'size' stored in 'buffer[]',
 * 'flags[i]' indicates whether entry i is owned by current process,
 * 'map[i]' is the global id of entry i 
 *
 * If flags == NULL then all entries are owned by current process
 * If map == NULL then map[i] == i + (sum of count of prev procs) is assumed.
 */
{
    DUMMY_g

#if USE_MPI
    if (g->nprocs > 1) {
#if 0
	/* old algorithm: each process sends its buffer to proc 0, and proc 0
	 * writes received data segments to file.
	 *
	 * FIXME: this algorithm is extremely slow on LSSC-3 (because of non
	 * contiguous writes) and may produce wrong output if nprocs is large
	 * (>= 32) */
	MPI_Datatype type, type0;
	MPI_Status status;
	MPI_Request req;
	INT i, j, k;
	int *blocklens, *indices, count_max, pass;
	long pos = 0;
	unsigned char *buf, *ptr;

	if (g->rank == 0)
	    pos = ftell(fp);

	pass = count;
	MPI_Allreduce(&pass, &count_max, 1, MPI_INT, MPI_MAX, g->comm);
	buf = phgAlloc(size * count_max);
	blocklens = phgAlloc(2 * (count_max + 1) * sizeof(*blocklens));
	indices = blocklens + count_max + 1;
	map_sort = phgAlloc(count * sizeof(*map_sort));

	j = 0;
	if (!map_flag) {
	    int tmp = count, map0;
	    MPI_Scan(&tmp, &map0, 1, MPI_INT, MPI_SUM, g->comm);
	    map0 -= tmp;
	    for (i = 0; i < count; i++) {
		map_sort[i].gindex = i + map0;
		map_sort[i].lindex = i;
	    }
	}
	else if (count > 0) {
	    /* sort map[] (blocks must have increasing offsets in filetype) */
	    for (i = 0; i < count; i++) {
		map_sort[i].gindex = map[i];
		map_sort[i].lindex = i;
	    }
	    qsort(map_sort, count, sizeof(*map_sort), map_sort_compare);
	}

	j = 0;
	for (i = 0; i < count; i++) {
	    k = map_sort[i].lindex;
	    if (flags != NULL && !(flags[k] & OWNER)) continue;
	    blocklens[j] = 1;
	    indices[j++] = k;
	}
	j = merge_blocks(j, blocklens, indices);
	MPI_Type_contiguous(size, MPI_BYTE, &type0);
	MPI_Type_commit(&type0);
	MPI_Type_indexed(j, blocklens, indices, type0, &type);
	MPI_Type_commit(&type);

	j = 0;
	for (i = 0; i < count; i++) {
	    k = map_sort[i].lindex;
	    if (flags != NULL && !(flags[k] & OWNER)) continue;
	    blocklens[j] = 1;
	    indices[j++] = map_sort[i].gindex;
	}
	count = merge_blocks(j, blocklens, indices);
	blocklens[count] = -1;
	indices[count] = j;

	for (pass = 0; pass < g->nprocs; pass++) {
	    if (g->rank != 0 && g->rank != pass)
		continue;

	    if (pass != 0) {
		/* send blocklens/indices to proc 0 */
		if (g->rank == pass)
		    MPI_Send(blocklens, 2*(count_max+1) * sizeof(*blocklens),
			 MPI_BYTE, 0, 11, g->comm);
		else
		    MPI_Recv(blocklens, 2*(count_max+1) * sizeof(*blocklens),
			 MPI_BYTE, pass, 11, g->comm, &status);
	    }

	    /* send data to process 0 */
	    if (g->rank == pass) {
		MPI_Isend(buffer, 1, type, 0, 22, g->comm, &req);
	    }
	    if (g->rank == 0) {
		/* get count */
		j = 0;
		while (blocklens[j] >= 0)
		    j++;
		MPI_Recv(buf, indices[j], type0, pass, 22, g->comm, &status);
	    }
	    if (g->rank == pass)
		MPI_Wait(&req, &status);

	    if (g->rank != 0)
		continue;

	    /* write data */
	    ptr = buf;
	    for (j = 0; blocklens[j] >= 0; j++) {
		fseek(fp, pos + indices[j] * size, SEEK_SET);
		if (fwrite(ptr, size, blocklens[j], fp) != blocklens[j])
		    phgError(1, "%s: fwrite error, abort.\n", __func__);
		ptr += blocklens[j] * size;
	    }
	}	/* for (pass ...) */

	MPI_Type_free(&type);
	MPI_Type_free(&type0);

	phgFree(blocklens);
	phgFree(map_sort);
	phgFree(buf);

	if (g->rank == 0)
	    fseek(fp, pos + nglobal * size, SEEK_SET);
#else	/* 0 */
	/* new algorithm: collect and write data by truncs */
	static size_t TRUNK_SIZE = (16 * 1024 * 1024);	/* 16MB */
	size_t size0 = sizeof(map_sort->gindex);
	INT trunc_count = ((size_t)TRUNK_SIZE / size);
	unsigned char *sbuffer, *rbuffer, *p;	/* the file buffer */
	INT i, j, m, n, n0, ralloc;
	int rank, np, tag;
	MPI_Datatype type;
	MPI_Status status;
	MPI_Comm comm;

	rbuffer = NULL;
	ralloc = 0;
	if (g->rank == 0)
	    sbuffer = phgAlloc(trunc_count * size);
	else
	    sbuffer = phgAlloc(trunc_count * (size + size0));

	/* sort map[], which avoids repeated scanning of the arrays */
	map_sort = phgAlloc(count * sizeof(*map_sort));

	if (!map_flag) {
	    int tmp = count, map0;
	    MPI_Scan(&tmp, &map0, 1, MPI_INT, MPI_SUM, g->comm);
	    map0 -= tmp;
	    for (i = 0; i < count; i++) {
		map_sort[i].gindex = i + map0;
		map_sort[i].lindex = i;
	    }
	}
	else if (count > 0) {
	    /* sort map[] (blocks must have increasing offsets in filetype) */
	    for (i = 0; i < count; i++) {
		map_sort[i].gindex = map[i];
		map_sort[i].lindex = i;
	    }
	    qsort(map_sort, count, sizeof(*map_sort), map_sort_compare);
	}

	MPI_Comm_dup(g->comm, &comm);
	MPI_Type_contiguous(size + size0, MPI_BYTE, &type);
	MPI_Type_commit(&type);

	i = 0;
	tag = 0;
	for (n0 = 0; n0 < nglobal; n0 += trunc_count, tag = 1 - tag) {
	    while (i < count && map_sort[i].gindex < n0)
		i++;
	    p = sbuffer;
	    m = 0;
	    for (; i < count && map_sort[i].gindex < n0 + trunc_count; i++) {
		if (flags != NULL && !(flags[map_sort[i].lindex] & OWNER)) {
		    /* not owner of the entry, skip */
		    continue;
		}
		if (g->rank == 0) {
		    /* root process: copy data to sbuffer[] */
		    memcpy(sbuffer + size * (map_sort[i].gindex - n0),
			   buffer + size * map_sort[i].lindex,
			   size);
		    m++;
		    continue;
		}
		memcpy(p, buffer + size * map_sort[i].lindex, size);
		p += size;
		/* append gindex to buffer */
		memcpy(p, &map_sort[i].gindex, size0);
		p += size0;
	    }
	    if (g->rank > 0) {
		/* send data to root process */
		n = (INT)((p - sbuffer) / (size + size0));
		MPI_Send(sbuffer, n, type, 0, tag, comm);
		continue;
	    }
	    /* root process: receive data from all other processes */
	    for (np = 1; np < g->nprocs; np++) {
		int tmp;
#if 0
		/* More efficient, but less reliable */
		MPI_Probe(MPI_ANY_SOURCE, tag, comm, &status);
#else
		MPI_Probe(np, tag, comm, &status);
#endif
		rank = status.MPI_SOURCE;
		MPI_Get_count(&status, type, &tmp);
		n = tmp;
		assert(n <= trunc_count);
		if (ralloc < n) {
		    phgFree(rbuffer);
		    rbuffer = phgAlloc((ralloc = n) * (size + size0));
		}
		MPI_Recv(rbuffer, n, type, rank, tag, comm, &status);
		/* process received buffer */
		p = rbuffer;
		for (j = 0; j < n; j++) {
		    INT k;
		    memcpy(&k, p + size, size0);
		    k -= n0;
		    assert(k >= 0 && k < nglobal - n0 && k < trunc_count);
		    memcpy(sbuffer + k * size, p, size);
		    p += size + size0;
		}
		m += n;
	    }
	    if (!(m == trunc_count ||
		   (n0 + trunc_count > nglobal && m == nglobal - n0))) {
		phgError(1, "%s:%d: unexpected.\n", __FILE__, __LINE__);
	    }
	    /* write sbuffer */
	    if (fwrite(sbuffer, size, m, fp) != m)
		phgError(1, "%s: fwrite error, abort.\n", __func__);
	}

	MPI_Type_free(&type);
	MPI_Comm_free(&comm);

	phgFree(map_sort);
	phgFree(sbuffer);
	phgFree(rbuffer);
#endif	/* 0 */
	return;
    }
#endif	/* USE_MPI */

    assert(count == nglobal);
    if (fwrite(buffer, size, count, fp) != count)
	phgError(1, "%s: fwrite error, abort.\n", __func__);;
}

static void
phgIOClose_(void)
{
    if (fp != NULL) {
	fclose(fp);
	fp = NULL;
    }
}

/*---------------------- End of sequential I/O functions --------------------*/

const char *
phgIOAddExtension(const char *fn, const char *ext)
/* adds '.ext' to the given filename if the latter does not contain
 * an extension, returns the address of the new filename which points either
 * to 'fn' or a static buffer */
{
    static char new_name[PATH_MAX];
    const char *p;

    if (ext == NULL)
	return fn;
    if ((p = strrchr(fn, '/')) == NULL)
	p = fn;
    if ((p = strchr(p, '.')) != NULL)
	return fn;
    sprintf(new_name, "%s.%s", fn, ext);
    return new_name;
}

/*---------------------------- The wrapper functions ------------------------*/

BOOLEAN
phgIOOpen(MPI_Comm comm, const char *fn)
{
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	initialized = TRUE;
	phgOptionsRegisterTitle("\n[MPI I/O options]", "\n", "mpi");
	phgOptionsRegisterNoArg("-mpi_io",
				"Use MPI parallel I/O functions when available",
				&use_mpiio);
	phgOptionsRegisterNoArg("-mpi_io_atomicity",
				"Set MPI parallel I/O atomicity",
				&mpiio_atomicity);
	return TRUE;
    }

#if USE_MPIIO
    if (use_mpiio)
	return phgIOOpenMPI(comm, fn);
    else
#endif	/* USE_MPIIO */
    return phgIOOpen_(comm, fn);
}

int
phgIOPrintf(MPI_Comm comm, char *fmt, ...)
{
    va_list ap;
    int ret;
    int (*print)(MPI_Comm comm, char *fmt, va_list ap);

#if USE_MPIIO
    if (use_mpiio)
	print = phgIOPrintfMPI;
    else
#endif	/* USE_MPIIO */
    print = phgIOPrintf_;

    va_start(ap, fmt);
    ret = print(comm, fmt, ap);
    va_end(ap);

    return ret;
}

void
phgIOWriteRoot(MPI_Comm comm, void *buf, size_t size, int count)
{
#if USE_MPIIO
    if (use_mpiio)
	phgIOWriteRootMPI(comm, buf, size, count);
    else
#endif	/* USE_MPIIO */
    phgIOWriteRoot_(comm, buf, size, count);

    return;
}

void
phgIOWrite(MPI_Comm comm, void *buffer, size_t size, size_t count,
		size_t nglobal, BTYPE flags[], INT map[], BOOLEAN map_flag)
{
#if USE_MPIIO
    if (use_mpiio)
	phgIOWriteMPI(comm, buffer, size, count, nglobal, flags, map, map_flag);
    else
#endif	/* USE_MPIIO */
    phgIOWrite_(comm, buffer, size, count, nglobal, flags, map, map_flag);
}

void
phgIOClose(void)
{
#if USE_MPIIO
    if (use_mpiio)
	phgIOCloseMPI();
    else
#endif	/* USE_MPIIO */
    phgIOClose_();
}

/*-------------------- output buffer and managing functions -----------------*/

static unsigned char *buffer = NULL;
static size_t buffer_pos, buffer_size;

static void *
swap_endian(void *buff, size_t n)
{
#if PHG_BIGENDIAN
    return buff;
#else	/* PHG_BIGENDIAN */
    static unsigned char buff0[8];
    unsigned char *p0, *p1;

    if (n <= 1)
	return buff;

    assert(sizeof(buff0) >= n);

    p0 = buff0 + n;
    p1 = buff;
    while (p0 > buff0) *(--p0) = *(p1++);
    return buff0;
#endif	/* PHG_BIGENDIAN */
}

static void
bigendian_copy(void *buf, size_t size, size_t count, size_t loc)
{
    size_t i;

    assert(buffer != NULL && loc + size * count <= buffer_size);
    for (i = 0; i < count; i++) {
	memcpy(buffer + loc, swap_endian(buf, size), size);
	buf += size;
	loc += size;
    }
}

void *
phgIOAllocBuffer(size_t size)
{
    buffer_pos = 0;
    buffer_size = size;
    buffer = phgAlloc(size);
    return buffer;
}

void
phgIOFreeBuffer(void)
{
    phgFree(buffer);
    buffer = NULL;
    buffer_pos = buffer_size = 0;
}

void
phgIOResetBuffer(void)
{
    buffer_pos = 0;
}

void
phgIOBigEndianAppend(void *buf, size_t size, int count)
{
    size_t size0 = size, count0 = count;
    bigendian_copy(buf, size0, count0, buffer_pos);
    buffer_pos += size0 * count0;
}

/*------------------ END output buffer and managing functions ---------------*/

#if defined(GZIP_PROG) || defined(BZIP2_PROG)
#define MAX_FILES 16
static struct {
    FILE *fp;
    BOOLEAN is_pipe;
} open_files[MAX_FILES];
#endif

void
phgCloseInputFile_(FILE *fp)
{
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    int i;
    BOOLEAN is_pipe = FALSE;
    for (i = 0; i < MAX_FILES; i++) {
	if (open_files[i].fp == fp) {
	    is_pipe = open_files[i].is_pipe;
	    open_files[i].fp = NULL;
	    break;
	}
    }
    (is_pipe ? pclose : fclose)(fp);
#else
    if (fp != NULL) fclose(fp);
#endif
}

FILE *
phgOpenInputFile_(const char *filename)
{
    FILE *fp;
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    char cmd[PATH_MAX + 16];
    int len = strlen(filename);
    BOOLEAN is_pipe = FALSE;
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	int i;
	initialized = TRUE;
	for (i = 0; i < MAX_FILES; i++)
	    open_files[i].fp = NULL;
    }
#endif

#ifdef GZIP_PROG
    if (len > 3 && !phgMemcmp(filename + len - 3, ".gz", 3)) {
	if ((fp = fopen(filename, "r")) != NULL) {
	    fclose(fp);
	    phgInfo(1, "gzipped file.\n");
	    sprintf(cmd, GZIP_PROG " -dc %s 2>/dev/null", filename);
	    fp = popen(cmd, "r");
	    is_pipe = TRUE;
	}
    } else
#endif
#ifdef BZIP2_PROG
    if (len > 4 && !phgMemcmp(filename + len - 4, ".bz2", 4)) {
	if ((fp = fopen(filename, "r")) != NULL) {
	    fclose(fp);
	    phgInfo(1, "bzip2 file.\n");
	    sprintf(cmd, BZIP2_PROG " -dc %s 2>/dev/null", filename);
	    fp = popen(cmd, "r");
	    is_pipe = TRUE;
	}
    } else
#endif
    {
	fp = fopen(filename, "r");
    }
    if (fp == NULL)
	phgWarning("can't open file \"%s\"!\n", filename);
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    else {
	int i;
	for (i = 0; i < MAX_FILES; i++)
	    if (open_files[i].fp == NULL) break;
	if (i >= MAX_FILES) {
	    phgWarning("too many open files.\n");
	    (is_pipe ? pclose : fclose)(fp);
	    fp = NULL;
	}
	else {
	    open_files[i].fp = fp;
	    open_files[i].is_pipe = is_pipe;
	}
    }
#endif

    return fp;
}
