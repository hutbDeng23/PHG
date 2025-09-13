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

/* $Id: checkpoint.c,v 1.10 2022/06/24 02:42:26 zlb Exp $ */

#include "phg.h"
#include "phg/checkpoint.h"
#include <stdarg.h>
#include <sys/stat.h>
#include <unistd.h>	/* access() */
#include <string.h>	/* strcmp() */	

static void 
save_profile(GRID *g, const char *dirname, int ndof, DOF **dofs, INT datasize)
{
    char  fn[100];
    FILE *fp;
   
    sprintf(fn,"%s/profile.p%04d", dirname, g->rank);

    phgInfo(1, "* save data to %s\n", fn);
    if ((fp = fopen(fn, "w")) == NULL) {
        phgWarning("write profile data %s failed!\n", fn);
    } else {
        int i;
        fprintf(fp,"nprocs %d\nnelem %"dFMT"\nserial_no %d\nndofs %d\n%"dFMT"\n",
             g->nprocs, g->nelem_global, g->serial_no, ndof, datasize);

        for(i=0;i<ndof;i++){
           fprintf(fp,"%s\n%"dFMT" %"dFMT"\n", dofs[i]->name,
		(INT)DofDataCount(dofs[i]), (INT)DofDataCountGlobal(dofs[i]));
        }
        fclose(fp);
    }

    return;
}

static void 
save_grid_dofs(GRID *g, const char *dirname, int ndof, DOF **dofs)
{
    int i;
    char  fn[100];
    FILE *fp;
    INT n;

    /* TODO: save Hierarchical grid data */
    sprintf(fn,"%s/dofdata.p%04d", dirname,g->rank);
    if ((fp = fopen(fn, "w")) == NULL) {
        phgWarning("write dof data %s failed!\n", fn);
    } else {
        for(i=0; i<ndof; i++)
        {
            n = DofDataCount(dofs[i]);
            fwrite(&n, sizeof(INT), 1, fp);
            fwrite(dofs[i]->data, sizeof(FLOAT), n, fp);
        }
        fclose(fp);
    }

    return;
}

static void 
save_userdata(GRID *g, const char *dirname, BYTE *userdata, INT datasize)
{
    char  fn[100];
    FILE *fp;
   
    if(userdata==NULL) 
         return;

    sprintf(fn,"%s/userdata.p%04d", dirname,g->rank);

    phgInfo(1, "* save user data to %s\n", fn);
    if ((fp = fopen(fn, "w")) == NULL) {
        phgWarning("write user data %s failed!\n", fn);
    } else {
        fwrite(&datasize, sizeof(INT), 1, fp);
        fwrite(userdata, sizeof(BYTE), datasize, fp);
        fclose(fp);
    }
    return;
}

void
phgCheckpointSave(GRID *g, const char *dirname, BOOLEAN distribute,
         INT datasize, BYTE *userdata, DOF *dof, ...)
/* Save the grid info., DOFs' data and user defined date into files for program resume
 *
 * when distribute == TRUE, the data is stored in seperated files with suffix ".p$rank"
 *      distribute == FALSE, all the process store the data in a single file. (TODO)
 *
 *      */
{
    int ndof;
    DOF **dofs;
    va_list ap;
   
    if(g==NULL){
       return;
    }

    dofs = phgAlloc(256 * sizeof(*dofs));

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    if(access(dirname, 0)==-1){
#ifndef __WIN__
	mkdir(dirname,S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH);
#else	/* !defined(__WIN__) */
	mkdir(dirname);
#endif	/* !defined(__WIN__) */
    }

    save_profile(g,dirname, ndof, dofs, datasize);
    save_grid_dofs(g,dirname,ndof,dofs);
    save_userdata(g, dirname, userdata, datasize);

    phgFree(dofs);

    return ;
}

static void 
check_profile(GRID *g, const char *dirname, int ndof, DOF **dofs, INT datasize)
{
    char  fn[100];
    FILE *fp;
   
    sprintf(fn,"%s/profile.p%04d", dirname, g->rank);

    if ((fp = fopen(fn, "r")) == NULL) {
        phgWarning("write profile data %s failed!\n", fn);
    } else {
        int i;
        int nprocs, serial_no , no_dof;
        char name[128];
        INT nelem_global, datacount, datacountglobal;

        if (fscanf(fp,"nprocs %d\nnelem %"dFMT"\nserial_no %d\nndofs %d\n%"dFMT"\n",
              &nprocs, &nelem_global, &serial_no, &no_dof, &datacount) != 5)
	    phgError(1,"%s:%d: error reading file \"%s\".\n",
						__FILE__, __LINE__, fn);
        if(nprocs != g->nprocs||
           nelem_global != g->nelem_global||
           serial_no != g->serial_no ||
           no_dof != ndof ||
           datacount != datasize)
               phgError(1,"*Mismatch: nprocs %d/%d nelem %"dFMT"/%"dFMT" serial_no %d/%d ndofs %d/%d userdatasize %"dFMT"/%"dFMT"\n", g->nprocs, nprocs, g->nelem_global,nelem_global, g->serial_no, serial_no, ndof, no_dof, datasize, datacount);

        for(i=0;i<ndof;i++){
           if (fscanf(fp, "%s\n%"dFMT" %"dFMT"\n",
			name, &datacount, &datacountglobal) != 3)
		phgError(1,"%s:%d: error reading file \"%s\".\n",
			__FILE__, __LINE__, fn);
           if(datacount != DofDataCount(dofs[i])||
              datacountglobal != DofDataCountGlobal(dofs[i])||
              strcmp(name, dofs[i]->name)!=0)
                phgError(1,"*Mismatch: dofdatacount %"dFMT"/%"dFMT" %"dFMT"/%"dFMT" dofname %s <=> %s\n", datacount, DofDataCount(dofs[i]), datacountglobal, DofDataCountGlobal(dofs[i]), name, dofs[i]->name);
        }
        fclose(fp);
    }

    return;
}


static void 
load_grid_dofs(GRID *g, const char *dirname, int ndof, DOF **dofs)
{
    
    int i;
    char  fn[100];
    FILE *fp;
    INT n;
 
    sprintf(fn,"%s/dofdata.p%04d", dirname,g->rank);
    if ((fp = fopen(fn, "r")) == NULL) {
        phgWarning("read dof data %s failed!\n", fn);
    } else {
        INT size;
        for(i=0; i<ndof; i++)
        {
            n = DofDataCount(dofs[i]);
            if (fread(&size, sizeof(INT), 1, fp) != 1)
		phgError(1,"%s:%d error reading dof data from \"%s\".\n",
				__FILE__, __LINE__, fn);
            if(size!=n || (size=fread(dofs[i]->data,sizeof(FLOAT),size,fp)) != n)
                phgError(1,"*[rank %d] Mismatch: dof[%s]'s data size: %d, load size: %d", g->rank, n, size);
        }
        fclose(fp);
    }
    return;
}

static int 
load_userdata(GRID *g, const char *dirname, BYTE *userdata, INT userdatasize)
{
    char  fn[100];
    FILE *fp;
   
    if(userdata==NULL) 
       return -1;

    sprintf(fn,"%s/userdata.p%04d", dirname,g->rank);

    if ((fp = fopen(fn, "r")) == NULL) {
        phgWarning("%s:%d: reading userdata file \"%s\" failed!\n",
						__FILE__, __LINE__, fn);
	return -1;
    } else {
        INT dd;
        if (fread(&dd, sizeof(INT), 1, fp) != 1) {
	    phgWarning("%s:%d: error reading file \"$s\".\n",
					__FILE__, __LINE__, fn);
	    return -1;
	}
        if(dd != userdatasize ||
		(dd = fread(userdata, sizeof(BYTE), dd, fp)) != userdatasize) {
             phgWarning("%s: %"dFMT" byte userdata read, which should be %"dFMT"\n",
			__FILE__, __LINE__, dd, userdatasize);
	    return -1;
	}
        fclose(fp);
    }
    return 0;
}

void
phgCheckpointRestore(GRID *g, const char *dirname, BOOLEAN distribute,
           INT userdatasize, BYTE *userdataout, DOF *dof, ...)
/* 
 * Load the grid info., DOFs' data and user defined date from files 
 *
 */
{
    int ndof;
    DOF **dofs;
    va_list ap;
    /*BYTE *userdatein;*/

    if(g==NULL){
       return;
    }

    dofs = phgAlloc(256 * sizeof(*dofs));

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    if(access(dirname, 0)==-1){
       phgError(1,"There is no restore diectory %s\n", dirname);
    }

    check_profile(g,dirname, ndof, dofs, userdatasize);
    load_grid_dofs(g,dirname,ndof,dofs);
    load_userdata(g,dirname,userdataout, userdatasize);

    phgFree(dofs);

    return;
}

