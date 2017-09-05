variable OLPATH [file dirname [file normalize [info script]]]

source ${OLPATH}/TMAlign.tcl

proc excluded_volume {id sel_glycan} {
  set pdbs [glob ${OLPATH}/*.pdb]
  set reference_idx none
  set idx []
  foreach pdb $pdbs {
    lappend idx [mol new $pdb]
    if {[string compare $pdb ${folder_ab}/reference_HIV.pdb]==0} {
      set reference_idx  [lindex $idx end]
    }
  }
  set glycan_seg [lsort -unique [[atomselect $id $sel_glycan] get segname]]
  array set idx_glycans [build_enum $glycan_seg]
  set contacts []
  foreach k $glycan_seg {
    lappend contacts 0
  }
  set glycans [atomselect $id $sel_glycan]
  set radius 3
  foreach p {1 2 3} {
    set M [TMAlign [atomselect $id "segname \"${p}P.*\""] [atomselect $reference_idx "protein"]]
    foreach i $idx {
      set Ab [atomselect $i all] 
      $Ab move $M
      if {$i !=$reference_idx } {
        set name [lindex [split [lindex [split [molinfo $i get name] "/"] end] "."] 0]
        set f [open "glycans_${name}_${p}.dat" w]
        puts $f $glycan_seg
        set num_steps [molinfo $id get numframes]
        for {set frame 0} {$frame < $num_steps} {incr frame} {
          $glycans frame $frame
          set glycans_idx [lindex [measure contacts $radius $glycans $Ab] 0]
          set Abcontacts $contacts
          if {[llength $glycans_idx] > 0} {
            set ss [atomselect $id "index $glycans_idx"]
            set contact_glycans [lsort -unique [$ss get {segname resid name}]]
            foreach cg $contact_glycans {
              lassign $cg segn resid n
              set v_init [lindex $Abcontacts $idx_glycans($segn)]
              set nbrcontact [expr $v_init + 1]
              lset Abcontacts $idx_glycans($segn) $nbrcontact
            }
            $ss delete
          }
          puts $f $Abcontacts
        }
        $Ab delete
        close $f
      }
    }
  }
}


proc build_enum {keys} {
  set i 0
  foreach k $keys {
    set enum($k) $i
    incr i
  }
  return [array get enum]
}

