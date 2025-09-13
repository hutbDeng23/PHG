#!./phg_tcl
#
# $Id: phg0.tcl,v 1.1 2006/02/11 04:12:37 zlb Exp $

if {$argc != 1} {
    puts "Usage: phg.tcl filename"
    exit 1
} else {
    set file [lindex $argv 0]
}

phgInit
phgImport $file
for {set i 0} {$i < 5} {incr i} {phgRefineAllElements}
phgPartitionGrid
for {set i 0} {$i < 10} {incr i} {phgRefineRandomElements "10%"}

#phgExportALBERT "output.dat"

phgFinalize

exit
