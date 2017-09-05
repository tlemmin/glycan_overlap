
proc TMAlign {ref mobile} {
    $ref writepdb ref.pdb
    $mobile writepdb mobile.pdb
    exec        /usr/local/bin/TMalign_f mobile.pdb ref.pdb -m matrix.dat
    set M [read_matrix "matrix.dat"]
    file delete matrix.dat
    file delete ref.pdb
    file delete mobile.pdb
    return $M
}


proc read_matrix {matrix} {
 set f [open $matrix r]
 set data [split [read $f] "\n"]
 close $f
 set M [list [list [lindex $data {2 2}] [lindex $data {2 3}] [lindex $data {2 4}] [lindex $data {2 1}]]\
 [list [lindex $data {3 2}] [lindex $data {3 3}] [lindex $data {3 4}] [lindex $data {3 1}]] \
  [list [lindex $data {4 2}] [lindex $data {4 3}] [lindex $data {4 4}] [lindex $data {4 1}]] \
    [list 0.0 0.0 0.0 1.0]]
 return $M
}

