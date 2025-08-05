#load MDanalysis and VMD before running this script

python /data01/genbiolab/mdanh/data/MD_scripts/md_membrane/cnt_membrane_preparation.py \
    cnt_6_45_oh.pdb \
    DPPC_140.pdb \
    -o cnt_6_45_oh_DPPC_140 \
    -t cnt_6_45_oh.rtf /data01/genbiolab/mdanh/data/CHARMMFF/toppar/top_all36_lipid.rtf \
    --xsc_file /data01/genbiolab/mdanh/data/simulation/membrane/MEM_DPPC/new/DPPC_140/namd/step6.6_equilibration.xsc \
    --rotation_x 90 --rotation_y 0 --rotation_z 0 \
    --radius_factor 1.0 \
    --buffer_distance 0.8 \
    --symmetric_removal \
    --salt_concentration 0.15 \
    --constraint_force 2.5