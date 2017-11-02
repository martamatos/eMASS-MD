mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-002_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst102 {dock_prep/GAPDH_WT_APO-004_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-016_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-024_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-026_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-027_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-032_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-036_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-038_receptor_G3P_docked.pdb }

foreach i $lst102 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 0 }
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 0
mol modselect 2 0 resname G3H
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 0 ColorID 0
mol drawframes 0 2 {0:1000}
mol showrep 0 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-005_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst112 {dock_prep/GAPDH_WT_APO-012_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-033_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-034_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-048_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-049_receptor_G3P_docked.pdb }

foreach i $lst112 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 1 }
mol addrep 1
mol modselect 1 1 protein
mol modstyle 1 1 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 1
mol modselect 2 1 resname G3H
mol modstyle 2 1 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 1 ColorID 1
mol drawframes 1 2 {0:1000}
mol showrep 1 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-013_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst111 {dock_prep/GAPDH_WT_APO-017_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-042_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-045_receptor_G3P_docked.pdb }

foreach i $lst111 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 2 }
mol addrep 2
mol modselect 1 2 protein
mol modstyle 1 2 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 2
mol modselect 2 2 resname G3H
mol modstyle 2 2 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 2 ColorID 2
mol drawframes 2 2 {0:1000}
mol showrep 2 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-006_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst311 {dock_prep/GAPDH_WT_APO-008_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-010_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-035_receptor_G3P_docked.pdb }

foreach i $lst311 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 3 }
mol addrep 3
mol modselect 1 3 protein
mol modstyle 1 3 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 3
mol modselect 2 3 resname G3H
mol modstyle 2 3 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 3 ColorID 3
mol drawframes 3 2 {0:1000}
mol showrep 3 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-001_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst101 {dock_prep/GAPDH_WT_APO-030_receptor_G3P_docked.pdb dock_prep/GAPDH_WT_APO-046_receptor_G3P_docked.pdb }

foreach i $lst101 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 4 }
mol addrep 4
mol modselect 1 4 protein
mol modstyle 1 4 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 4
mol modselect 2 4 resname G3H
mol modstyle 2 4 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 4 ColorID 4
mol drawframes 4 2 {0:1000}
mol showrep 4 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-007_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst211 {dock_prep/GAPDH_WT_APO-015_receptor_G3P_docked.pdb }

foreach i $lst211 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 5 }
mol addrep 5
mol modselect 1 5 protein
mol modstyle 1 5 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 5
mol modselect 2 5 resname G3H
mol modstyle 2 5 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 5 ColorID 5
mol drawframes 5 2 {0:1000}
mol showrep 5 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-023_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst201 {dock_prep/GAPDH_WT_APO-025_receptor_G3P_docked.pdb }

foreach i $lst201 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 6 }
mol addrep 6
mol modselect 1 6 protein
mol modstyle 1 6 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 6
mol modselect 2 6 resname G3H
mol modstyle 2 6 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 6 ColorID 6
mol drawframes 6 2 {0:1000}
mol showrep 6 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-050_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst034 {}

foreach i $lst034 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 7 }
mol addrep 7
mol modselect 1 7 protein
mol modstyle 1 7 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 7
mol modselect 2 7 resname G3H
mol modstyle 2 7 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 7 ColorID 7
mol drawframes 7 2 {0:1000}
mol showrep 7 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-001_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst100 {}

foreach i $lst100 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 8 }
mol addrep 8
mol modselect 1 8 protein
mol modstyle 1 8 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 8
mol modselect 2 8 resname G3H
mol modstyle 2 8 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 8 ColorID 8
mol drawframes 8 2 {0:1000}
mol showrep 8 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-041_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst103 {}

foreach i $lst103 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 9 }
mol addrep 9
mol modselect 1 9 protein
mol modstyle 1 9 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 9
mol modselect 2 9 resname G3H
mol modstyle 2 9 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 9 ColorID 9
mol drawframes 9 2 {0:1000}
mol showrep 9 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-037_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst113 {}

foreach i $lst113 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 10 }
mol addrep 10
mol modselect 1 10 protein
mol modstyle 1 10 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 10
mol modselect 2 10 resname G3H
mol modstyle 2 10 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 10 ColorID 10
mol drawframes 10 2 {0:1000}
mol showrep 10 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-009_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst013 {}

foreach i $lst013 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 11 }
mol addrep 11
mol modselect 1 11 protein
mol modstyle 1 11 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 11
mol modselect 2 11 resname G3H
mol modstyle 2 11 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 11 ColorID 11
mol drawframes 11 2 {0:1000}
mol showrep 11 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-031_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst012 {}

foreach i $lst012 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 12 }
mol addrep 12
mol modselect 1 12 protein
mol modstyle 1 12 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 12
mol modselect 2 12 resname G3H
mol modstyle 2 12 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 12 ColorID 12
mol drawframes 12 2 {0:1000}
mol showrep 12 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-040_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst002 {}

foreach i $lst002 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 13 }
mol addrep 13
mol modselect 1 13 protein
mol modstyle 1 13 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 13
mol modselect 2 13 resname G3H
mol modstyle 2 13 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 13 ColorID 13
mol drawframes 13 2 {0:1000}
mol showrep 13 0 0


mol new {test_files/test_cluster_dock/dock_prep/GAPDH_WT_APO-044_receptor_G3P_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
set lst325 {}

foreach i $lst325 { mol addfile "$i" type {pdb} first 0 last -1 step 10 waitfor 1 14 }
mol addrep 14
mol modselect 1 14 protein
mol modstyle 1 14 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 14
mol modselect 2 14 resname G3H
mol modstyle 2 14 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 14 ColorID 14
mol drawframes 14 2 {0:1000}
mol showrep 14 0 0


mol new {test_files/test_cluster_dock/1DC4_chainA.pdb} type {pdb} first 0 last -1 step 1 waitfor 1
mol addrep 15
mol modselect 1 15 protein
mol modstyle 1 15 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 15
mol modselect 2 15 resname G3H
mol modstyle 2 15 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 15 ColorID 15
mol drawframes 15 2 {0:1000}
mol showrep 15 0 0


