mol new {/home/dx/Desktop/MD_general/Dock/GAPDH/1_APO/1_3PG/dock_prep/GAPDH_WT_APO-001_receptor_3PG_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
animate style Loop
set lst111 {dock_prep/GAPDH_WT_APO-002_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-004_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-005_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-007_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-012_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-013_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-014_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-015_receptor_3PG_docked.pdb dock_prep/GAPDH_WT_APO-016_receptor_3PG_docked.pdb}

foreach i $lst111 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 0 }
