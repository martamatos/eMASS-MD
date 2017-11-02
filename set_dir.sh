enzyme=ENO

main_folder_local_escaped="\/home\/mrama\/Dropbox\/PhD_stuff\/Projects\/MD\/eMASS-MD\/enzyme_models"
main_folder_cluster_escaped="\/zhome\/89\/0\/90554\/kinetics_fit"

main_folder_local=/home/mrama/Dropbox/PhD_stuff/Projects/MD/eMASS-MD/enzyme_models
main_folder_cluster=/zhome/89/0/90554/kinetics_fit

enzyme_folder=${enzyme}
enzyme_subfolder=${enzyme}_param_scan
enzyme_folder_cluster=${enzyme}_param_scan

cd $main_folder_cluster/$enzyme_folder/$enzyme_subfolder/data

rename 's/psoParameters.*.txt/psoParameters.txt/' *.txt
rename 's/lmaParameters.*.txt/lmaParameters.txt/' *.txt


sed -i -- "s/$main_folder_local\/$enzyme_folder\/$enzyme_subfolder\/input/$main_folder_cluster\/$enzyme_folder_cluster\/data/g" *.dat psoParameters.txt lmaParameters.txt
