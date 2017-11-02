mol new {/home/mrama/Desktop/MD/TALB/0_Crystalline_structures/1UCW_wt_S7P_watbox.pdb} type {pdb} first 0 last -1 step 1 waitfor 1

animate style Once

mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 0
mol modselect 2 0 residue 130 32 241
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 0 ResID

mol addrep 0
mol modselect 3 0 resname S7P
mol modstyle 3 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 3 0 Name


mol showrep 0 0 0
