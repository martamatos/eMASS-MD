mol new {/home/mrama/Desktop/MD/ENO_1E9I/0_Crystal_structures/1e9i_subD_noSO4-HOH.pdb} type {pdb} first 0 last -1 step 1 waitfor 1

mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 0
mol modselect 2 0 resname MG
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor Name 


mol showrep 0 0 0
