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

#include "vtkPhgData.h"

#include <vtkByteSwap.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCellArray.h>
#include <vtkFieldData.h>
#include <vtkDataSetAttributes.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#ifndef __cplusplus
#define __cplusplus
#endif

#include "phg.h"

/*vtkCxxRevisionMacro(vtkPhgData, "$Revision: 1.4 $");*/
vtkStandardNewMacro(vtkPhgData);

vtkPhgData::vtkPhgData()
{
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::New();
    this->SetOutput(output);
    output->ReleaseData();
    output->Delete();

    this->ScalarsName = NULL;
    this->VectorsName = NULL;
}

vtkPhgData::~vtkPhgData()
{
}
	
vtkUnstructuredGrid *vtkPhgData::GetOutput()
{
    return this->GetOutput(0);
}

vtkUnstructuredGrid* vtkPhgData::GetOutput(int idx)
{
    return vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(idx));
}

void vtkPhgData::SetOutput(vtkUnstructuredGrid *output)
{
    this->GetExecutive()->SetOutputData(0, output);
}

int vtkPhgData::RequestUpdateExtent(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  int piece, numPieces;

  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  if (piece < 0 || piece >= numPieces) return 0;

  return 1;
}

static vtkIdType *cells_ptr;
static unsigned int *data_ptr;
static int *types;
static int global_index, skip1, read2, skip3;
static int callback_type;

static BOOLEAN
phg_callback CB_ARGS(s)
{
    switch (callback_type) {
	case 0:		/* vertices and type of the cell */
	    if (skip1 > 0) {skip1--; return TRUE;}
	    if (read2 <= 0) return FALSE;
	    // cells->InsertNextCell(4, s->verts);
	    *(types++) = 10;
	    *(cells_ptr++) = 4;
	    *(cells_ptr++) = s->verts[0];
	    *(cells_ptr++) = s->verts[1];
	    *(cells_ptr++) = s->verts[2];
	    *(cells_ptr++) = s->verts[3];
	    read2--;
	    return read2 > 0 ? TRUE : FALSE;
	case 1:		/* element index */
	    *(data_ptr++) = global_index++;
	    break;
	case 2:		/* refinement type */
	    *(data_ptr++) = s->type;
	    break;
	case 3:		/* element mark */
	    *(data_ptr++) = s->mark;
	    break;
    }

    return TRUE;
}

#if USE_MPI
extern "C" GRID *get_submesh(GRID *g);
#endif
extern GRID *grid_;

int vtkPhgData::RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
	outInfo->Get(vtkDataObject::DATA_OBJECT()));
    if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER())
	> 0) return 1;

#if USE_MPI
    GRID *g = get_submesh(grid_);
#else
    GRID *g = grid_;
#endif
    int i, *types_save;
    int ncells, size, piece, numPieces/*, ghostLevel*/;
    vtkCellArray *cells;
    vtkDataSetAttributes *attr = (vtkDataSetAttributes*)(output->GetCellData());
    vtkUnsignedIntArray *data;

    if (g == NULL) {
	cerr << "vtkPhgData: grid undefined.\n";
	return 1;
    }

    /* vertices */
    vtkPoints *points=vtkPoints::New();
    points->SetDataTypeToFloat();
    points->SetNumberOfPoints(g->nvert);
    for (i = 0; i < g->nvert; i++) points->SetPoint(i, g->verts[i]);
    output->SetPoints(points);
    points->Delete();

    /* cells and cell_types */
    ncells = g->nleaf;
    size = 5 * ncells;
    piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    //output->GetUpdateExtent(piece, numPieces, ghostLevel);
    skip1 = piece * ncells / numPieces;
    read2 = ((piece+1) * ncells / numPieces) - skip1;
    skip3 = ncells - skip1 - read2;
    //phgInfo(0, "skip1=%d, read2=%d, skip3=%d\n", skip1, read2, skip3);

    cells = vtkCellArray::New();
    types = new int[read2];
    cells_ptr = cells->WritePointer(ncells, size);
    types_save = types;
    callback_type = 0;
    phgTraverseElements(g, phg_callback);
    //phgInfo(0, "Number of cells: %d\n", cells->GetNumberOfCells());
    output->SetCells(types_save, cells);
    cells->Delete();
    delete [] types_save;
    char *sname = this->ScalarsName;

    if (sname == NULL || !strcmp(sname, "element_index")) {
	// global index
	data = vtkUnsignedIntArray::New();
	data->SetName("element_index");
	data->SetNumberOfComponents(1);
	data->SetNumberOfTuples(ncells);
	data_ptr=((vtkUnsignedIntArray *)data)->WritePointer(0, ncells);
	global_index = 0;
	callback_type = 1;
	phgTraverseElements(g, phg_callback);
	attr->SetScalars(data);
	data->Delete();
    } else if (!strcmp(sname, "refinement_type")) {
	// refinement type
	// FIXME: seems VTK doesn't like scalars of type char?
	data = vtkUnsignedIntArray::New();
	data->SetName("refinement_type");
	data->SetNumberOfComponents(1);
	data->SetNumberOfTuples(ncells);
	data_ptr=((vtkUnsignedIntArray *)data)->WritePointer(0, ncells);
	callback_type = 2;
	phgTraverseElements(g, phg_callback);
	attr->SetScalars(data);
	data->Delete();
    } else if (!strcmp(sname, "element_mark")) {
	// element mark
	data = vtkUnsignedIntArray::New();
	data->SetName("element_mark");
	data->SetNumberOfComponents(1);
	data->SetNumberOfTuples(ncells);
	data_ptr=((vtkUnsignedIntArray *)data)->WritePointer(0, ncells);
	callback_type = 3;
	phgTraverseElements(g, phg_callback);
	attr->SetScalars(data);
	data->Delete();
    }

    if (g != grid_) phgFreeGrid(&g);

    return 1;
}

int vtkPhgData::FillOutputPortInformation(int, vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}
