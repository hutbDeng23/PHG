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

/* $Id: vtkPhgData.h,v 1.3 2016/06/24 08:53:33 zlb Exp $ */

/* Note: vtkDataReader is used in place of vtkSource, because we
 * need the Get/SetScalars/Vectors/etc functions. */

#include <vtkDataReader.h>

class vtkUnstructuredGrid;

class vtkPhgData : public vtkDataReader
{
public:
  static vtkPhgData *New();
  /*vtkTypeRevisionMacro*/vtkTypeMacro(vtkPhgData, vtkDataReader);
  vtkUnstructuredGrid *GetOutput();
  vtkUnstructuredGrid *GetOutput(int idx);
  void SetOutput(vtkUnstructuredGrid *output);

protected:
  vtkPhgData();
  ~vtkPhgData();
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                                  vtkInformationVector *);
  int RequestData(vtkInformation *, vtkInformationVector **,
                                  vtkInformationVector *);
  int FillOutputPortInformation(int, vtkInformation*);
};
