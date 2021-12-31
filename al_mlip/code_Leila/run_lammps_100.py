#!/bin/bash

###########################################################################################
## Import modules
import numpy as np; import os; import inspect; import sys;
import getpass
def lammps_run(main_folder, file_name, lammps_input):
    cmnd1 = 'cd ' + main_folder + '/' + file_name + '/' + ';/home/leila/Downloads/lammps-stable/lammps-7Aug19/src/lmp_mpi<' + lammps_input
    print(cmnd1)
    os.system(cmnd1)



stilt = '100';

############# Provide directories for simulations
### Provide the name of the directory where the input file will be placed.
if (stilt == '100'):
	stilt_dir = '/symm_tilts_100';
if (stilt == '110'):
	stilt_dir = '/symm_tilts_110';
if (stilt == '111'):
	stilt_dir = '/symm_tilts_111';
if ((getpass.getuser() == 'patala') or (getpass.getuser() == 'srikanthpatala')):
    dir_path = ('/Users/' + getpass.getuser() + '/Dropbox/NCSU_Research/Repos/GBPy_new/');
    lmp_inp_dir = '/Users/'+getpass.getuser()+'/Desktop/symm_tilts'+stilt_dir;
	### Provide the name of the directory where the potential file can be found.
    pot_file_dir = '/Users/'+getpass.getuser()+'/Dropbox/NCSU_Research/Repos/GBpy_new/GBpy/bicrystal_gen/potentials'
    pot_file_name = 'Al99.eam.alloy';
if (getpass.getuser() == 'leila'):
    dir_path = '/home/leila/Leila_sndhard/codes/';
    lmp_inp_dir = '/home/leila/Desktop/symm_tilts'+stilt_dir;
	### Provide the name of the directory where the potential file can be found.
    pot_file_dir = '/home/leila/Leila_sndhard/codes/GBpy/bicrystal_gen/potentials'
    pot_file_name = 'in.EAM.dfs' 

sys.path.append(dir_path)

###########################################################################################


filepath = 'symm_tilts_'+stilt+'.txt';

nsz1 = 100000;
sig_num = np.zeros((nsz1,1));
sig_id = np.zeros((nsz1,1));
bpn_po1 = np.zeros((nsz1,3));

st100_list = [];

with open(filepath) as fp:
	line = fp.readline()
	cnt = 0;
	while line:
		s1 = line.split();
		fname = s1[1];
		st100_list.append(fname);

		s2 = fname.split('_');
		sig_num[cnt] = int(s2[1][1:]);
		sig_id[cnt] = int(s2[2]);
		bpn_po1[cnt,:] = np.array([int(s2[4]), int(s2[5]), int(s2[6])]);
		line = fp.readline();
		cnt += 1;


ngb = np.shape(sig_num)[0];

# ct1 = 0;
for ct1 in range(ngb):
	file_name = st100_list[ct1]
	lammps_input= 'in.'+ file_name + '_1'
	lammps_run(lmp_inp_dir,file_name, lammps_input )


sys.path.remove(dir_path)