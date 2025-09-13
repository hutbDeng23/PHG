#!/usr/bin/vtk
#
# vtk script for writing a VTK data file in ASCII format
#
# $Id: ASCII.tcl,v 1.2 2006/02/22 14:12:39 zlb Exp $

package require vtk
package require vtkinteraction

# set the name of data file
if {$argc != 2} {
    puts "\nUsage: ASCII.tcl input.vtk output.vtk\n"
    exit
}
set infile [lindex $argv 0]
set outfile [lindex $argv 1]

vtkDataSetReader reader
  reader SetFileName "$infile"
  #reader SetDataByteOrderToBigEndian
  

vtkDataSetWriter writer
  writer SetInputConnection [reader GetOutputPort]
  writer SetFileName "$outfile"
  writer SetFileTypeToASCII
  writer Write

exit 0
