mol new {/home/mrama/Desktop/MD/GAPDH/0_Crystal_structures/1DC4_chainA.pdb} type {pdb} first 0 last -1 step 1 waitfor 1

animate style Once

mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 0
mol modselect 2 0 residue 148 147 149 207 208 230 175
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 0 ColorID 7

mol addrep 0
mol modselect 3 0 resname G3H
mol modstyle 3 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 3 0 ColorID 8


mol showrep 0 0 0
