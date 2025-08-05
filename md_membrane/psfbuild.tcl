package require psfgen
psfgen_logfile "cnt_mem_preparation.log"

topology /data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_5_45_none.rtf
topology /data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_10_45_none.rtf
topology /data01/genbiolab/mdanh/data/CHARMMFF/cnt_ff/final/multi_wall/cnt_15_45_oh.rtf
topology /data01/genbiolab/mdanh/data/CHARMMFF/toppar/top_all36_lipid.rtf

segment CNT1 {
    pdb cnt_5_45_none_centered.pdb
}

coordpdb cnt_5_45_none_centered.pdb CNT1

segment CNT2 {
    pdb cnt_10_45_none_centered_chainB.pdb
}

coordpdb cnt_10_45_none_centered_chainB.pdb CNT2

segment CNT3 {
    pdb cnt_15_45_oh_centered_chainC.pdb
}

coordpdb cnt_15_45_oh_centered_chainC.pdb CNT3


segment MEMB {
    pdb cleaned_mem.pdb    
}

coordpdb cleaned_mem.pdb MEMB

guesscoord
writepsf cnt_mem.psf
writepdb cnt_mem.pdb

psfgen_logfile close
quit