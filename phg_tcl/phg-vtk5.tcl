#!./phg_tcl
#
# $Id: phg-vtk5.tcl,v 1.1 2014/12/27 07:03:34 zlb Exp $

#encoding system euc-cn

catch {package require vtk}
catch {package require vtkinteraction}

option clear
###option add "*Font" {systemfont 10 bold}
option add "*Font" {systemfont 10}
option add "*BorderWidth" 1
option add "*activeBorderWidth" 1
option add "*cursor" top_left_arrow
option add "*activeForeground" purple

if {$argc != 1} {
    puts "Usage: phg.tcl filename"
    set filename ""
} else {
    set filename [lindex $argv 0]
}

phgInit
if {$filename != ""} {phgImport $filename}
vtkPhgData phgGrid
  phgGrid SetScalarsName "element_index"

set nprocs [phgGetNProcs]

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

set statusline ""

proc update_status {} {
    global statusline
    set statusline [phgShowStatus]
}

set current_cursor 0

proc set_cursor {type} {
    global current_cursor
    set current_cursor $type
    if {$type == 0} {set cursor top_left_arrow} else {set cursor watch}
    .vtkw			configure -cursor $cursor
    .ctrl			configure -cursor $cursor
    .ctrl.logo			configure -cursor $cursor
    .ctrl.display_menu		configure -cursor $cursor
    .ctrl.refine_menu		configure -cursor $cursor
    .ctrl.color_menu		configure -cursor $cursor
    .ctrl.faces			configure -cursor $cursor
    .ctrl.refine.all		configure -cursor $cursor
    .ctrl.refine.random		configure -cursor $cursor
    .ctrl.refine.percent	configure -cursor $cursor
    .ctrl.refine.surface	configure -cursor $cursor
    .ctrl.refine.refine		configure -cursor $cursor
    .ctrl.refine		configure -cursor $cursor
    .ctrl.partition		configure -cursor $cursor
    .status			configure -cursor $cursor
    .menubar			configure -cursor $cursor
    update
}

set cursor_stack(0)	0
set cursor_sp		0

proc push_cursor {} {
    global cursor_stack cursor_sp current_cursor
    incr cursor_sp
    set cursor_stack($cursor_sp) $current_cursor
}

proc pop_cursor {} {
    global cursor_stack cursor_sp
    set_cursor $cursor_stack($cursor_sp)
    if {$cursor_sp > 0} {set cursor_sp [expr $cursor_sp - 1]}
}

set draw_faces 1
proc redraw {} {
  global draw_faces filename
  if {$filename == ""} {return}
  push_cursor
  set_cursor 1
  phgGrid Update
  if {[phgGrid GetScalarsName] == "refinement_type"} {
    mapper SetScalarRange 1.0 3.0
  } else {
    set range [[[[phgGrid GetOutput] GetCellData] GetScalars] GetRange]
    mapper SetScalarRange [lindex $range 0] [lindex $range 1]
  }
  if {$draw_faces} {ren AddActor actor} else {catch {ren RemoveActor actor}}
  renWin Render
  pop_cursor
}

proc load_grid {} {
    global display_menu refine_menu nprocs filename
    set types {
	{{Any mesh files}	{*.dat *.dat.gz *.dat.bz2 \
				 *.mesh *.mesh.gz *.mesh.bz2 \
				}}
	{{ALBERT files}		{*.dat *.dat.gz *.dat.bz2}}
	{{Medit mesh files}	{*.mesh *.mesh.gz *.mesh.bz2}}
	{{All Files}		{*}}
    }
    set filename [phgGetFileName]
    regsub -all "/\[^/\]*\$" $filename "" dir
    if {$dir == $filename} {
	set dir "."
    } else {
	regsub -all "^$dir/" $filename "" filename
    } 
    set filename [tk_getOpenFile -filetypes $types -initialfile $filename \
	-initialdir $dir -title "Load grid"]
    if {$filename == ""} {return}
    push_cursor
    set_cursor 1
    phgFreeGrid
#    for {set i 0} {$i < $nprocs} {incr i} {
#	$display_menu entryconf $i -state disabled
#	$refine_menu entryconf [expr $i + 2] -state disabled
#    }
    .ctrl.display_menu configure -state disabled
    .ctrl.refine_menu configure -state disabled
    # FIXME: how to set current menu item to 0 below?
    $display_menu activate 0
    $refine_menu activate 0
    phgImport $filename
    set filename [phgGetFileName]
    ren ResetCamera
    redraw
    update_status
    .ctrl.partition configure -text "Partition"
    pop_cursor
}

proc albert_export {} {
    set types {
	{{ALBERT files}		{*.dat}		}
	{{All Files}		{*}		}
    }
    regsub -all "\\..*\$" [phgGetFileName] "" prefix
    regsub -all "/\[^/\]*\$" $prefix "" dir
    if {$dir == $prefix} {
	set dir "."
    } else {
	regsub -all "^$dir/" $prefix "" prefix
    } 
    set filename [tk_getSaveFile -filetypes $types -initialfile "$prefix.dat" \
	-initialdir $dir -title "ALBERT Export"]
    if {$filename == ""} {return}
    push_cursor
    set_cursor 1
    phgExportALBERT $filename
    pop_cursor
}

proc medit_export {} {
    set types {
	{{Medit mesh files}	{*.mesh}	}
	{{All Files}		{*}		}
    }
    regsub -all "\\..*\$" [phgGetFileName] "" prefix
    regsub -all "/\[^/\]*\$" $prefix "" dir
    if {$dir == $prefix} {
	set dir "."
    } else {
	regsub -all "^$dir/" $prefix "" prefix
    } 
    set filename [tk_getSaveFile -filetypes $types -initialfile "$prefix.mesh" \
	-initialdir $dir -title "Medit Export"]
    if {$filename == ""} {return}
    push_cursor
    set_cursor 1
    phgExportMedit $filename
    pop_cursor
}

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
    set filename [tk_getSaveFile -filetypes $types -initialfile ${prefix}.eps \
    	-initialdir $dir -title "PostScript export" -defaultextension ""]

    if {$filename == ""} {return}

    push_cursor
    set_cursor 1
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
	ps SetGlobalLineWidthFactor 0.5
	ps SetSortToBSP
	ps Write
	puts "File \"$filename\" created."
    mapper1 SetResolveCoincidentTopologyToPolygonOffset
    pop_cursor
    return
}

proc vtk_export {} {
    set types {
	{{VTK files}		{*.vtk}		}
	{{All Files}		{*}		}
    }
    regsub -all "\\..*\$" [phgGetFileName] "" filename
    regsub -all "/\[^/\]*\$" $filename "" dir
    if {$dir == $filename} {
	set dir "."
    } else {
	regsub -all "^$dir/" $filename "" filename
    } 
    set filename "${filename}.vtk"
    set filename [tk_getSaveFile -filetypes $types -initialfile $filename \
	-initialdir $dir -title "VTK export"]
    if {$filename == ""} {return}
    push_cursor
    set_cursor 1
    catch {writer Delete}
    vtkUnstructuredGridWriter writer
      writer SetInput [phgGrid GetOutput]
      writer SetFileName [phgConvertFileName $filename]
      writer Write
    puts "File \"$filename\" created."
    pop_cursor
}

proc partition {} {
    global color_menu display_menu refine_menu nprocs
    push_cursor
    set_cursor 1
    if {[phgPartitionGrid] == 0} {pop_cursor; return}
    phgCheckConformity
    redraw
    if {[phgGetNSubmeshes] > 1} {
	.ctrl.display_menu configure -state normal
	.ctrl.refine_menu configure -state normal
	.ctrl.partition configure -text "Repartition"
    }
    update_status
    pop_cursor

## Dialog for getting number of subdomains
#    toplevel .nparts
#    frame .nparts.input -relief flat -width 100
#    label .nparts.input.prompt -text "Number of subgrids" -relief flat
#    entry .nparts.input.value -width 16 -takefocus 1 -relief sunken -width 4 \
#    	-validate key -validatecommand {
#    		set a -1
#		catch {set a [expr %S + 0]}
#		return [expr $a >= 0 && $a <= 9]
#	}
#    .nparts.input.value insert end "2"
#    #.nparts.input.value sel range 0 [string length [.nparts.input.value get]]
#
#    frame .nparts.buttons -relief flat
#    button .nparts.buttons.ok -text "OK" -underline 0 -default active \
#        -takefocus 0 -command {
#	    set nparts 0
#	    catch {set nparts [expr [.nparts.input.value get] + 0]}
#	    if {$nparts > 0} {
#		phgPartitionGrid $nparts
#		$color_menu invoke 3
#		destroy .nparts
#	    } else {
#		tk_messageBox -message "Invalid value!" -type ok
#	    }
#	}
#    button .nparts.buttons.cancel -text "Cancel" -underline 0 -takefocus 0 \
#	-command {
#	    destroy .nparts
#	}
#    pack .nparts.input.prompt -side left
#    pack .nparts.input.value -side left
#    pack .nparts.input -side top -pady 10
#
#    pack .nparts.buttons.ok -side left -padx 20
#    pack .nparts.buttons.cancel -side left -padx 20
#    pack .nparts.buttons -side top -pady 10
#
#    set w .nparts
#    bind .nparts <Return> {
#	.nparts.buttons.ok flash; .nparts.buttons.ok invoke}
#    bind .nparts o {
#	.nparts.buttons.ok flash; .nparts.buttons.ok invoke}
#    bind .nparts <Escape> {
#	.nparts.buttons.cancel flash; .nparts.buttons.cancel invoke}
#    bind .nparts c {
#	.nparts.buttons.cancel flash; .nparts.buttons.cancel invoke}
#    focus .nparts.input.value
}

proc confirm_quit {} {
    #set answer [tk_messageBox -type okcancel -message "Really quit?"]
	set answer "ok"
    if {$answer == "ok"} {
	phgFinalize
	exit
	#::vtk::cb_exit
    }
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

#------------------ menus, toolbars, etc. ----------------------

# The menubar
menu .menubar -type menubar -relief groove
  .menubar add cascade -menu .menubar.file  -label "File" -underline 0
  .menubar add cascade -menu .menubar.view  -label "View" -underline 0
  .menubar add cascade -menu .menubar.help  -label "Help" -underline 0
  menu .menubar.file
    .menubar.file add comm -lab "Load ..." -underline 0 -comm {load_grid}
    .menubar.file add separator
    .menubar.file add comm -lab "Export ALBERT ..." -underline 0 \
	-comm {albert_export}
    .menubar.file add comm -lab "Export Medit ..." -underline 0 \
	-comm {medit_export}
    .menubar.file add comm -lab "Export VTK ..." -underline 7 -comm {vtk_export}
    .menubar.file add comm -lab "Export EPS ..." -underline 7 -comm {ps_export}
    .menubar.file add separator
    .menubar.file add comm -lab "Quit" -underline 0 -comm {confirm_quit}
  menu .menubar.view
    .menubar.view add comm -lab "Show axes" -comm {ren AddViewProp axes}
    .menubar.view add comm -lab "Hide axes" -comm {ren RemoveViewProp axes}
    .menubar.view add separator
    .menubar.view add comm -lab "+X" -comm {set_view  1  0  0   0 0 1}
    .menubar.view add comm -lab "-X" -comm {set_view -1  0  0   0 0 1}
    .menubar.view add comm -lab "+Y" -comm {set_view  0  1  0   1 0 0}
    .menubar.view add comm -lab "-Y" -comm {set_view  0 -1  0   1 0 0}
    .menubar.view add comm -lab "+Z" -comm {set_view  0  0  1   0 1 0}
    .menubar.view add comm -lab "-Z" -comm {set_view  0  0 -1   0 1 0}
  menu .menubar.help
. configure -menu menubar
pack .menubar -fill x -padx 0 -pady 0

#
# put the VTK drawing widget in a frame in order to add a border for it
frame .vtkw -relief groove -borderwidth 2
###wm protocol . WM_DELETE_WINDOW ::vtk::cb_exit
wm geometry . 900x600
set vtkw [vtkTkRenderWidget .vtkw.ren -rw renWin]
###::vtk::bind_tk_render_widget $vtkw
pack $vtkw -side top -fill both -expand yes
pack .vtkw -side left -fill both -expand yes

frame .ctrl -relief groove -borderwidth 0

set color_menu [tk_optionMenu .ctrl.color_menu color_dummy \
    "Color by index" "Color by type" "Color by mark" "Color off"]
    .ctrl.color_menu configure -width 16 -anchor c -direction below
    $color_menu entryconf 0 -comm {
	mapper ScalarVisibilityOn
	phgGrid SetScalarsName "element_index"; redraw}
    $color_menu entryconf 1 -comm {
	mapper ScalarVisibilityOn
	phgGrid SetScalarsName "refinement_type"; redraw}
    $color_menu entryconf 2 -comm {
	mapper ScalarVisibilityOn
	phgGrid SetScalarsName "element_mark"; redraw}
    $color_menu entryconf 3 -comm {
	mapper ScalarVisibilityOff; redraw}

set tmp "set display_menu \[tk_optionMenu .ctrl.display_menu display_dummy"
for {set i 0} {$i < $nprocs} {incr i} {set tmp "$tmp \"Display submesh $i\""}
set tmp "$tmp\]"
eval $tmp
.ctrl.display_menu configure -width 20 -anchor w -state disabled
for {set i 0} {$i < $nprocs} {incr i} {
    set tmp "$display_menu entryconf $i -comm {phgSelectDisplay $i; redraw}"
    eval $tmp
#    $display_menu entryconf $i -state disabled
}

set tmp "set refine_menu \[tk_optionMenu .ctrl.refine_menu refine_dummy"
set tmp "$tmp \"Refine all submeshes\""
for {set i 0} {$i < $nprocs} {incr i} {set tmp "$tmp \"Refine submesh $i\""}
set tmp "$tmp\]"
eval $tmp
.ctrl.refine_menu configure -width 20 -anchor w -state disabled
$refine_menu insert 1 separator
$refine_menu entryconf 0 -comm {phgSelectRefine -1}
for {set i 0} {$i < $nprocs} {incr i} {
    set tmp "$refine_menu entryconf [expr $i + 2] -comm {phgSelectRefine $i}"
    eval $tmp
#    $refine_menu entryconf [expr $i + 2] -state disabled
}

checkbutton .ctrl.faces -text "Draw element faces" -relief flat -width 23 \
	-anchor w -variable draw_faces -command {redraw}

set levels 1
frame .ctrl.levels -relief groove
label .ctrl.levels.text -text "\"Refine all\"/\"Coarsen\" level" -padx 5 -pady 5
frame .ctrl.levels.group -relief flat
label .ctrl.levels.group.number -textvariable levels -padx 5 \
	 -background white -relief groove
button .ctrl.levels.group.plus -text "+" -padx 5 -pady 0 -command {
    if {$levels < 6} {
	set levels [expr $levels + 1]
    }
}
button .ctrl.levels.group.minus  -text "-" -padx 5 -pady 0 -command {
    if {$levels > 1} {
	set levels [expr $levels - 1]
    }
}

frame .ctrl.refine -relief groove
checkbutton .ctrl.refine.all    -text "Refine all" \
     -relief flat -width 23 -anchor w \
    -command {
	set refine "phgRefineAllElements"
	pack forget .ctrl.refine.percent
	.ctrl.refine.all select
	.ctrl.refine.surface deselect
	.ctrl.refine.random deselect}
checkbutton .ctrl.refine.random -text "Refine random" \
     -relief flat -width 23 -anchor w \
    -command {
	set refine "phgRefineRandomElements"
	pack .ctrl.refine.percent -after .ctrl.refine.random \
		-side top -padx 10 -pady 0 -fill none
	.ctrl.refine.all  deselect
	.ctrl.refine.surface deselect
	.ctrl.refine.random select}
checkbutton .ctrl.refine.surface  -text "Refine surface" \
     -relief flat -width 23 -anchor w \
    -command {
	set refine "phgRefineSurface"
	pack forget .ctrl.refine.percent
	.ctrl.refine.all deselect
	.ctrl.refine.surface select
	.ctrl.refine.random deselect}
set percent 25
proc percent_update {percent} {
    .ctrl.refine.percent configure -label "Percentage: %$percent"
}
scale .ctrl.refine.percent -length 180 -orient h -variable percent \
    -digits 2 -from 0 -to 100 -relief flat -showvalue 0 \
    -command {percent_update}

button .ctrl.refine.refine -text "Refine"  \
    -command {
	push_cursor
	set_cursor 1
	case $refine in {
	    {phgRefineAllElements}	{$refine "$levels"}
	    {phgRefineRandomElements}	{$refine "$percent%"}
	    {phgRefineSurface}		{$refine}
	    {*}				{puts "unexpected, abort."; exit 1}
	}
	phgCheckConformity
	phgGrid Update
	mapper SetScalarRange 0 [expr [[mapper GetInput] GetNumberOfCells]-1]
	redraw
	update_status
	pop_cursor
    }

button .ctrl.coarsen -text "Coarsen"  \
    -command {
	push_cursor
	set_cursor 1
	phgCoarsenGrid "-$levels"
	phgCheckConformity
	phgGrid Update
	mapper SetScalarRange 0 [expr [[mapper GetInput] GetNumberOfCells] - 1]
	redraw
	update_status
	pop_cursor
    }

button .ctrl.partition -text "Partition" -command {partition}

if {[encoding system] == "euc-cn"} {
    button .ctrl.quit -text "ÍË³ö" -cursor pirate -command {confirm_quit}
    #button .ctrl.quit -text [encoding convertfrom euc-cn "ÍË³ö"] \
    #  -cursor pirate -command {confirm_quit}
} else {
    button .ctrl.quit -text "Quit" -cursor pirate -command {confirm_quit}
}

image create photo phg_logo -file "[file dirname $argv0]/phg-logo.gif"
label .ctrl.logo -image phg_logo -padx 0 -relief flat
pack .ctrl.logo -side top -pady 5

pack .ctrl.display_menu
pack .ctrl.refine_menu

pack .ctrl.color_menu
pack .ctrl.faces		-side top

pack .ctrl.levels.text	-side top
pack .ctrl.levels.group.minus	-side left -padx 5 -pady 0 -fill none
pack .ctrl.levels.group.number	-side left -padx 10 -pady 0 -fill none
pack .ctrl.levels.group.plus	-side left -padx 5 -pady 0 -fill none
pack .ctrl.levels.group	-side top
pack .ctrl.levels	-side top -padx 5 -pady 5 -fill none

pack .ctrl.refine.all		-side top -padx 10 -pady 0 -fill none
pack .ctrl.refine.surface	-side top -padx 10 -pady 5 -fill none
pack .ctrl.refine.random	-side top -padx 10 -pady 0 -fill none
#pack .ctrl.refine.percent	-side top -padx 10 -pady 0 -fill none
pack .ctrl.refine.refine	-side top -padx 45 -pady 5 -fill none
pack .ctrl.refine		-side top -pady 5

pack .ctrl.coarsen		-side top -padx 90 -pady 5 -fill none

pack .ctrl.partition		-side top -padx 10 -pady 5 -fill none

pack .ctrl.quit			-side bottom -padx 10 -pady 5 -fill none

pack .ctrl			-side left -fill y

# status line
label .status -textvariable statusline -anchor w -takefocus 0 -relief sunken \
	-foreground #600060
pack .status -after .menubar -side bottom -fill x -expand no
update_status

ren ResetCamera
redraw
.ctrl.refine.all invoke
.ctrl.faces select

tkwait window .
