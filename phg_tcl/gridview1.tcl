#!/usr/bin/vtk
#
# This script splits an unstructured grid into tetrahedra, and
# defines an actor for each tetrahedron.
#
# $Id: gridview1.tcl,v 1.1 2006/02/11 04:12:37 zlb Exp $

package require vtk
package require vtkinteraction

set width 500
set height 500

# set the name of data file
if {$argc == 0} {
    set file "output.vtk"
} else {
    set file [lindex $argv 0]
}

vtkUnstructuredGridReader reader
  reader SetFileName "$file"
  reader SetScalarsName "element_mark"
  puts "Number of scalars: [reader GetNumberOfScalarsInFile]"
  
vtkTetra tetra
  [tetra GetPointIds] SetNumberOfIds 4
  [tetra GetPointIds] SetId 0 0
  [tetra GetPointIds] SetId 1 1
  [tetra GetPointIds] SetId 2 2
  [tetra GetPointIds] SetId 3 3

for {set i 0} {$i < 8} {incr i} {
  case $i in {
    {0} {set r 1; set g 0; set b 0}
    {1} {set r 0; set g 1; set b 0}
    {2} {set r 0; set g 0; set b 1}
    {3} {set r 0; set g 0; set b 1}
    {4} {set r 0; set g 1; set b 1}
    {5} {set r 1; set g 0; set b 1}
    {6} {set r 1; set g 1; set b 0}
    {7} {set r 0.5; set g 0.5; set b 0.5}
  }
  vtkProperty prop($i)
    prop($i) SetRepresentationToSurface
    prop($i) BackfaceCullingOn
    prop($i) SetAmbient 0.5
    prop($i) SetAmbientColor $r $g $b
    prop($i) SetDiffuse 0.5
    prop($i) SetDiffuseColor $r $g $b
    prop($i) SetSpecular 0.5
    prop($i) SetSpecularColor $r $g $b
    prop($i) SetOpacity 0.5
}

vtkRenderer ren
  ren SetBackground 1 1 1
  
    reader Update
    set nCells [[reader GetOutput] GetNumberOfCells]
    for {set i 0} {$i < $nCells} {incr i} {
	# Note: vtkPoints::GetPoints returns the proc 'vtkTemp3'
	# we need to explicitly copy the points
	set pts [[[reader GetOutput] GetCell $i] GetPoints]
	vtkPoints points($i)
	  points($i) SetNumberOfPoints 4
	  for {set j  0} {$j < 4} {incr j} {
	    set p [$pts GetPoint $j]
	    points($i) InsertPoint $j [lindex $p 0] [lindex $p 1] [lindex $p 2]
	  }
	vtkUnstructuredGrid elem($i)
	  elem($i) Allocate 1 1
	  elem($i) InsertNextCell [tetra GetCellType] [tetra GetPointIds]
	  elem($i) SetPoints points($i)
	vtkDataSetMapper mapper($i)
	  mapper($i) SetInput elem($i)
	  mapper($i) ScalarVisibilityOn
	vtkActor actor($i)
	  actor($i) SetMapper mapper($i)
	  set da [[[reader GetOutput] GetCellData] GetScalars]
	  set j [$da GetTuple1 $i]
	  actor($i) SetProperty prop([expr $j%8])
	ren AddActor actor($i)
    }

    #puts "Number of actors = [[ren GetActors] GetNumberOfItems]."
    #for {set i 0} {$i < $nCells} {incr i} {
    #	 puts "  Element ($i) --- vertices:"
    #	 puts "    [[[elem($i) GetCell 0] GetPoints] GetPoint 0]"
    #	 puts "    [[[elem($i) GetCell 0] GetPoints] GetPoint 1]"
    #	 puts "    [[[elem($i) GetCell 0] GetPoints] GetPoint 2]"
    #	 puts "    [[[elem($i) GetCell 0] GetPoints] GetPoint 3]"
    #}
  
vtkRenderWindow renWin
  renWin AddRenderer ren
  renWin SetSize $width $height

#vtkInteractorStyleTrackballCamera style
vtkInteractorStyleTrackballActor style
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin
  iren SetInteractorStyle style
  iren Initialize
 
# prevent the tk window from showing up then start the event loop
wm withdraw .
