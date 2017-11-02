mol new {test_files/G6PD_WT_APO_G6P_docked.prmtop} type {parm7} first 0 last -1 step 1 waitfor 1
mol addfile {test_files/G6PD_WT_APO_G6P_docked_cl111_10frames.dcd} type {dcd} first 0 last -1 step 1 waitfor 1 0
animate style Once
mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name
mol addrep 0
mol modselect 2 0 resname G6P
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor Name
mol showrep 0 0 0
