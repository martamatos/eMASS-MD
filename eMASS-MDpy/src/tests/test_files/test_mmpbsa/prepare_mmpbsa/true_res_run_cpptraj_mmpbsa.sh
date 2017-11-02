cd /home/marta/G6PD/3_MD_post_dock/2_NADP/2_6PGL/cl111
$AMBERHOME/bin/cpptraj -i test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_complex_solv.in -p G6PD_WT_NADP_6PGL_docked.prmtop
$AMBERHOME/bin/cpptraj -i test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_complex_vac.in -p G6PD_WT_NADP_6PGL_docked.prmtop
$AMBERHOME/bin/cpptraj -i test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_ligand_vac.in -p G6PD_WT_NADP_6PGL_docked.prmtop
$AMBERHOME/bin/cpptraj -i test_files/test_mmpbsa/prepare_mmpbsa/get_traj_frame_receptor_vac.in -p G6PD_WT_NADP_6PGL_docked.prmtop
