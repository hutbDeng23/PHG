#!/usr/bin/vtk
#
# vtk script for displaying an unstructured grid
#
# Sample file: /usr/share/doc/vtk-data-4.2/Data/tetraMesh.vtk (vtk-data)
#
# $Id: gridview.tcl,v 1.1 2006/02/11 04:12:38 zlb Exp $

package require vtk
package require vtkinteraction

set width 500
set height 500

# set the name of data file
if {$argc == 0} {
    puts "\nUsage: gridview.tcl vtk_file\n"
    exit
}
set file [lindex $argv 0]

vtkUnstructuredGridReader reader
  reader SetFileName "$file"
  #reader SetScalarsName "element_info"
  set n [reader GetNumberOfScalarsInFile]
  puts "Number of scalars: $n"
  for {set i 0} {$i < $n} {incr i} {
    puts "    [reader GetScalarsNameInFile $i]"
  }
  set n [reader GetNumberOfVectorsInFile]
  puts "Number of vectors: $n"
  for {set i 0} {$i < $n} {incr i} {
    puts "    [reader GetVectorsNameInFile $i]"
  }
  
  	  # retrieve scalars data array
  	  reader Update
          set da "[[[reader GetOutput] GetCellData] GetScalars]"
	  puts "Number of tuples     = [$da GetNumberOfTuples]"
	  puts "Number of components = [$da GetNumberOfComponents]"
	  puts "Data type size       = [$da GetDataTypeSize]"

# as in the example showMesh.tcl we will show only the elements
# or elements faces on the surface of the objects
# see the comments in showMesh.tcl

vtkGeometryFilter geomFilter
  geomFilter SetInput [reader GetOutput]

vtkPolyDataMapper elementsMapper
  elementsMapper SetInput [geomFilter GetOutput]
  elementsMapper ScalarVisibilityOff

vtkActor elementsActor
  elementsActor SetMapper elementsMapper
  [elementsActor GetProperty] SetColor 0 0 1
  #[elementsActor GetProperty] SetDiffuse 1
  #[elementsActor GetProperty] SetAmbient 1
  ###[elementsActor GetProperty] SetInterpolationToFlat
  #[elementsActor GetProperty] SetOpacity 0.8

vtkExtractEdges edgesFilter
  edgesFilter SetInput [geomFilter GetOutput]

vtkPolyDataMapper edgesMapper
  edgesMapper SetInput [edgesFilter GetOutput]
  edgesMapper ScalarVisibilityOff
  edgesMapper SetResolveCoincidentTopologyToPolygonOffset

vtkActor edgesActor
  edgesActor SetMapper edgesMapper

set edgesProp [edgesActor GetProperty]
  $edgesProp SetColor 1 0 0
  $edgesProp SetDiffuse 1
  $edgesProp SetAmbient 1
  $edgesProp SetLineWidth 1
  
vtkRenderer ren
  ren SetBackground 1 1 1
  
# Generate data arrays containing point ids
vtkIdFilter ids
  ids SetInput [reader GetOutput]
  ids PointIdsOn
  ids CellIdsOn

#--- Point labels
if {0} {
  vtkSelectVisiblePoints visPoints
    visPoints SetInput [ids GetOutput]
    visPoints SetRenderer ren
    visPoints SetSelection 1 $width 1 $height
    visPoints SelectionWindowOn
  vtkLabeledDataMapper pno
    pno SetInput [ids GetOutput]
    pno SetLabelFormat "%g"
    pno SetLabelModeToLabelFieldData
    [pno GetLabelTextProperty] ShadowOff
    [pno GetLabelTextProperty] SetFontSize 20
    [pno GetLabelTextProperty] BoldOff
  vtkActor2D pointLabels
    pointLabels SetMapper pno
    [pointLabels GetProperty] SetColor 1 0 0
  ren AddActor pointLabels
}

#--- Cell labels
if {0} {
  # Find the centers of the cells specified by the array ids
  vtkCellCenters cc
    cc SetInput [ids GetOutput]
  # Create labels for cells
  vtkSelectVisiblePoints visCells
    visCells SetInput [cc GetOutput]
    visCells SetRenderer ren
    ###visCells SetSelection $xmin $xmax $ymin $ymax $zmin $zmax
    # You can excersise difference of using switching 
    # selection on and off. 
    # Switching selection off makes all cell numbers visible
    # regardless of the fact that if they into selection
    # frame or are obscured by other elements or are inside volumes
    visCells SelectionWindowOn
  vtkLabeledDataMapper cno
    cno SetInput [visCells GetOutput]
    cno SetLabelFormat "%g"
    cno SetLabelModeToLabelFieldData
    [cno GetLabelTextProperty] ShadowOff
    [cno GetLabelTextProperty] SetFontSize 20
    [cno GetLabelTextProperty] BoldOff
  vtkActor2D cellLabels
    cellLabels SetMapper cno
    eval [cellLabels GetProperty] SetColor 0 1 0
  ren AddActor cellLabels
}

ren AddActor elementsActor
ren AddActor edgesActor

vtkRenderWindow renWin
  renWin AddRenderer ren
  renWin SetSize $width $height

vtkInteractorStyleTrackballCamera style
vtkRenderWindowInteractor iren
  iren SetInteractorStyle style
  iren SetRenderWindow renWin
  iren Initialize
  
# prevent the tk window from showing up then start the event loop
wm withdraw .
