# eMASS-MDpy
All Python 2.7 code used for estimating binding energies differences from MD to integrate with enzyme-level kinetic models.
Together with the data in eMASS-MD, this code can be used to reproduce all results plots in chapter 4 of my PhD thesis. 

### Overview

#### cluster_dock
Includes scripts to cluster the docking results and get a representative frame for each cluster.

#### kinetics_integration
Includes a script to prepare and run enzyme-level kinetic model fittings on an HPC cluster, and several scripts to analyse the results of integrating the MD results with enzyme-level kinetic models. 

#### mmpbsa
Includes 1) scripts to prepare an mmpbsa run, from generating the necessary pdb files to generating the scripts that will be executed to start mmpbsa; and 2) scripts to harvest the results and process them.

#### parsing
Includes random scripts to fix issues i've been stumbling on along the way, e.g. adding the ligand name in a .mol2 file.

#### process_traj
Includes scripts to analyze a given MD trajectory, it needs both MDAnalysis and mdtraj.

#### tests
Includes unit tests. To run these you need to change `tests_dir` in the respective `__init__.py`file.


#### run_MD_analyses
Contains custom scripts to analyze MD and docking results - needs to be adapted to run in different systems.

#### run_kinetic_analyses
Contains custom scripts to analyze the results of integrating MD with enzyme-level kinetic models - can be run in any system by changing "base_folder".


### Requirements:
 - Jinja2==2.8
 - matplotlib==1.5.1
 - [MDAnalysis==0.14.0](http://www.mdanalysis.org/)
 - numpy==1.10.4
 - pandas==0.18.0
 - PyYAML==3.11
 - scipy==0.17.0
 - seaborn==0.8.1
 - [mdtraj==1.6.1](http://mdtraj.org/1.6.2/)


### Reproducing results:

To reproduce the graphs in my PhD thesis, you only need to install the requirements in requirements_plots.txt with pip.
For more instructions see the general eMASS-MD README.

To use all scripts, you need to install the requirements in requirements_all.txt with pip.

If you are using a virtualenv, you might need to update pip before installing the requirements (otherwise scipy will complain it can't find lapack and the likes). To do so run the following command inside the virtualenv:

`python -m pip install --upgrade pip`

