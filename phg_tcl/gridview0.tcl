#!/usr/bin/vtk
#
# vtk script for displaying an unstructured grid
#
# Sample file: /usr/share/doc/vtk-data-4.2/Data/tetraMesh.vtk (vtk-data)
#
# $Id: gridview0.tcl,v 1.2 2008/01/06 14:56:26 zlb Exp $

package require vtk
package require vtkinteraction

vtkRenderWindow renWin
renWin SetSize 800 392
vtkRenderer ren
ren SetViewport 0.01 0.01 0.495 0.99
ren SetBackground 0 0 0
renWin AddRenderer ren
vtkRenderer ren1
ren1 SetViewport 0.505 0.01 0.99 0.99
ren1 SetBackground 0 0 0
renWin AddRenderer ren1

vtkUnstructuredGridReader reader
  reader SetFileName [lindex $argv 0]

vtkDataSetMapper mapper
  mapper SetInput [reader GetOutput]

vtkActor actor
  actor SetMapper mapper
ren AddViewProp actor

vtkActor actor1
  actor1 SetMapper mapper
  [actor1 GetProperty] BackfaceCullingOn
  [actor1 GetProperty] SetDiffuseColor .2 1 1
  [actor1 GetProperty] SetSpecular .4
ren1 AddViewProp actor1

# Interaction
vtkInteractorStyleTrackballCamera style
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin
  iren SetInteractorStyle style
  iren Initialize

# FIXME: (only on pc111) VTK does not respond to window events if
#	 "renWin Render" is called before "iren Initialize" (vtk 4.2 and 4.5.0)
ren ResetCamera
ren1 ResetCamera
renWin Render
wm withdraw .
