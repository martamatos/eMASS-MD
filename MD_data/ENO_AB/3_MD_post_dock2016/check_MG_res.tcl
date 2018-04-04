mol new {/home/mrama/Desktop/MD/eMASS-MD_complete_data/MD_data/ENO_AB/2_Dock/1_MG/2_PEP/dock_prep/ENO_AB_WT_MG-001_receptor_PEP_docked.pdb} type {pdb} first 0 last -1 step 1 waitfor 1


mol addrep 0
mol modselect 1 0 protein
mol modstyle 1 0 NewCartoon 0.300000 10.000000 4.100000 0
mol modcolor Name

mol addrep 0
mol modselect 2 0 resname MG
mol modstyle 2 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 2 0 ColorID 0
mol drawframes 0 2 {0:1000}
mol showrep 0 0 0

mol addrep 0
mol modselect 3 0 within 5 of resname MG
mol modstyle 3 0 CPK 1.000000 0.300000 12.000000 12.000000
mol modcolor 3 0 Name
mol drawframes 0 3 {0:1000}
mol showrep 0 0 0


label add Bonds 0/6407 0/4254
label add Bonds 0/6407 0/4655
label add Bonds 0/6407 0/3570

label add Bonds 0/6408 0/3570
label add Bonds 0/6408 0/4266
label add Bonds 0/6408 0/2409
