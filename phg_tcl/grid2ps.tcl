#!./phg_tcl
#
# $Id: grid2ps.tcl,v 1.1 2006/02/11 04:12:37 zlb Exp $

#encoding system euc-cn

catch {package require vtk}
catch {package require vtkinteraction}

option clear
#option add "*Font" {systemfont 12 bold}
option add "*Font" {systemfont 16 bold}

if {$argc != 1} {
    puts "\nUsage: grid2ps.tcl mesh_file\n"
    exit
}
set filename [lindex $argv 0]

phgInit
if {$filename != ""} {phgImport $filename}
vtkPhgData phgGrid
  phgGrid SetScalarsName "element_index"

vtkProperty property
  property SetRepresentationToSurface
  #property SetRepresentationToWireframe

vtkRenderer ren
  ren SetBackground 1 1 1
  ren SetLightFollowCamera 1

vtkRenderWindow renWin
  renWin AddRenderer ren
  #renWin SetSize 1000 1000

vtkInteractorStyleTrackballCamera style
vtkRenderWindowInteractor iren
  iren SetRenderWindow renWin
  iren SetInteractorStyle style

vtkLookupTable lut
set c 1
if {$c == 1} {
    # Normal colors
    lut SetHueRange		0.2 0.8
    #lut SetSaturationRange	0.0 1.0
    #lut SetValueRange		0.5 1.0
    #lut SetAlphaRange		0.5 0.5 
} elseif {$c == 2} {
    # No color
    lut SetHueRange		0 0
    lut SetSaturationRange	0 0
    lut SetValueRange		1 1
    lut SetAlphaRange		0 0
} elseif {$c == 3} {
    # Grey scale
    lut SetHueRange		0.2 0.3
    lut SetSaturationRange	0 0.5
    lut SetValueRange		0.75 1.0
    #lut SetAlphaRange		0 0
}

vtkDataSetMapper mapper
  mapper SetInput [phgGrid GetOutput]
  mapper ScalarVisibilityOn
  mapper SetScalarModeToUseCellData
  mapper SetLookupTable lut
vtkActor actor
  actor SetMapper mapper
  actor SetProperty property
 
vtkGeometryFilter geom
  geom SetInput [phgGrid GetOutput]
vtkExtractEdges edges
  edges SetInput [geom GetOutput]
vtkPolyDataMapper mapper1
  mapper1 SetInput [edges GetOutput]
  mapper1 ScalarVisibilityOff
  mapper1 SetResolveCoincidentTopologyToPolygonOffset
vtkActor actor1
  actor1 SetMapper mapper1
  [actor1 GetProperty] SetColor 1 0 0
  #vtkProperty prop1
  #  prop1 SetOpacity 1
  #actor1 SetBackfaceProperty prop1
  #[actor1 GetProperty] BackfaceCullingOn

ren AddActor actor
ren AddActor actor1

# Draw axes
vtkTextProperty tprop
    tprop SetColor 0 0 1
    tprop ShadowOff
vtkCubeAxesActor2D axes
if {1} {
    axes SetInput [phgGrid GetOutput]
    axes SetFlyModeToOuterEdges
} else {
    axes SetProp actor
    axes SetFlyModeToClosestTriad
    axes ScalingOff
}
axes SetCamera [ren GetActiveCamera]
axes SetLabelFormat "%0.2g"
axes SetFontFactor 1.5
axes SetAxisTitleTextProperty tprop
axes SetAxisLabelTextProperty tprop
[axes GetProperty] SetColor 0 0 1
[axes GetProperty] SetLineWidth 2
#ren AddViewProp axes

proc ps_export {} {
	#{{Text Files}       {.txt}        }
	#{{TCL Scripts}      {.tcl}        }
	#{{C Source Files}   {.c}      TEXT}
	#{{GIF Files}        {.gif}        }
	#{{GIF Files}        {}        GIFF}
    set types {
	{{PostScript files}	{.eps *.ps}	}
	{{All Files}		{*}		}
    }

    regsub -all "\\..*\$" [phgGetFileName] "" prefix
    regsub -all "/\[^/\]*\$" $prefix "" dir
    if {$dir == $prefix} {
	set dir "."
    } else {
	regsub -all "^$dir/" $prefix "" prefix
    } 
    #set filename [tk_getSaveFile -filetypes $types -initialfile ${prefix}.eps \
    #	-initialdir $dir -title "PostScript export" -defaultextension ""]
    set filename ${prefix}.eps

    if {$filename == ""} {return}

    catch {ps Delete}
    mapper1 SetResolveCoincidentTopology 0
    vtkGL2PSExporter ps
	ps SetInput renWin
	if {[regexp "\\.ps\$" $filename]} {
	    ps SetFileFormatToPS
	    regsub -all "\\.ps\$" $filename "" prefix
	} else {
	    ps SetFileFormatToEPS
	    regsub -all "\\.eps\$" $filename "" prefix
	}
	ps SetFilePrefix [phgConvertFileName $prefix]
	ps SetCompress 0
	ps SetDrawBackground 0
	ps SetSimpleLineOffset 0
	ps SetPS3Shading 0
	ps SetGlobalLineWidthFactor 1.0
	ps SetSortToBSP
	ps Write
	puts "File \"$filename\" created."
    mapper1 SetResolveCoincidentTopologyToPolygonOffset
    return
}

proc set_view {sx sy sz vx vy vz} {
    set camera [ren GetActiveCamera]
    scan [$camera GetPosition] "%f %f %f" x y z
    set x [expr $sx * abs($x)]
    set y [expr $sy * abs($y)]
    set z [expr $sz * abs($z)]
    $camera SetFocalPoint 0 0 0
    #$camera SetPosition $x $y $z
    $camera SetPosition $sx $sy $sz
    $camera SetViewUp $vx $vy $vz
    ren ResetCamera
    renWin Render
}

phgGrid Update
set range [[[[phgGrid GetOutput] GetCellData] GetScalars] GetRange]
mapper SetScalarRange [lindex $range 0] [lindex $range 1]
ren AddActor actor
renWin Render
set_view  -1  -0.5 -0.5   -0.5 -0.5 -1
ps_export

exit
