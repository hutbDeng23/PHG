#!/usr/bin/vtk
#
# Volume rendering
#
# $Id: gridview2.tcl,v 1.2 2010/09/05 14:39:30 zlb Exp $

catch {package require vtk}
catch {package require vtkinteraction}

set width 700
set height 700

vtkUnstructuredGridReader reader
  reader SetFileName "output.vtk"
  reader SetScalarsName "element_index"
vtkRenderer ren
  ren SetBackground 1 1 1
vtkRenderWindow renWin
  renWin AddRenderer ren
  renWin SetSize $width $height

vtkInteractorStyleTrackballCamera style
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin
  iren SetInteractorStyle style

    vtkVersion v
    if {[v GetVTKMajorVersion] < 4 ||
	([v GetVTKMajorVersion] == 4 && [v GetVTKMinorVersion] < 5)} {
	puts "\nError: \"update_grid_v\" requires VTK version >= 4.5.0."
	puts "Please use \"update_grid_1\" or \"update_grid_n\" instead."
	puts "Abort."
	exit
    }

    vtkPiecewiseFunction opacity
	opacity AddPoint 0.0     0.85
	opacity AddPoint 65536.0 0.85

    vtkColorTransferFunction color
	color AddRGBPoint 0.0 0.0 0.0 0.0
	color AddRGBPoint 1.0 0.0 0.0 1.0
	color AddRGBPoint 2.0 0.0 1.0 0.0
	color AddRGBPoint 3.0 0.0 1.0 1.0
	color AddRGBPoint 4.0 1.0 0.0 0.0
	color AddRGBPoint 5.0 1.0 0.0 1.0
	color AddRGBPoint 6.0 1.0 1.0 0.0
	color AddRGBPoint 7.0 1.0 1.0 1.0
    vtkVolumeProperty vproperty
	vproperty SetColor color
	vproperty SetScalarOpacity opacity
	vproperty ShadeOn
	vproperty SetInterpolationTypeToLinear
    vtkUnstructuredGridVolumeRayCastMapper mapper
	mapper SetInput [reader GetOutput]
    vtkVolume volume
	volume SetMapper mapper
	volume SetProperty vproperty
    ren AddVolume volume

wm withdraw .
renWin Render
