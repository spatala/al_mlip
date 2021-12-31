import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpinterfaces.utils import print_exception
from ase import Atoms
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.lammps.data import LammpsData
import numpy as np
from glob import glob
from ase.eos import EquationOfState
import pandas as pd
import subprocess
import json 
from ase.units import kJ

def rmse(e1,e2):
    return np.sqrt(np.sum([(e1[n]-e)**2 for n,e in enumerate(e2)])/len(e1))

def generate_xyz_training(na,e,fo=None,t='cfg',dest=None,name=None):
    if t=='cfg':
       k = {26.981539:'Al',26.981:'Al'}
       p = read(na,':')[0]
       species = [k[i] for i in p.get_masses()]
       lenA = len([a for a in species if a=='Al'])
       X = Atoms('Al{0}'.format(lenA), pbc=True, cell=p.cell,positions=p.positions)
       fx = p.get_array('fx')
       fy = p.get_array('fy')
       fz = p.get_array('fz')
       force = np.array([[fx[n],fy[n],fz[n]] for n,f in enumerate(fx)])
       sp = SinglePointCalculator(X,energy=e,forces=force)
       X.set_calculator(sp)
       write(na.replace('.cfg','with_energy_force.xyz'),X)
    elif t=='struct':
       p = AseAtomsAdaptor.get_atoms(na)
       sp = SinglePointCalculator(p,energy=e,forces=fo)
       p.set_calculator(sp)
       write(dest+os.sep+name,p)

def call_ext(cmd_string,osc=False):
    print("calling {}".format(cmd_string.split(' ')[0]))
    if osc:
        returncode=os.system(cmd_string)
    else:
        returncode = subprocess.call(cmd_string,shell=True)
    if returncode:
          print("ERROR in {0}".format(cmd_string))
    else:
          print("{} command exited successfully".format(cmd_string.split(' ')[0]))


def run_lammps(pattern,infile):
  for d in glob(pattern):
    os.chdir(d)
    try:
       call_ext('{0} < {1} > nohup.out'.format(lmp_path,infile))
       with open('log.lammps') as rf:
         lins = rf.readlines()
         index_top = min([n for n,l in enumerate(lins) if 'Step' in l])
         index_end = max([n for n,l in enumerate(lins) if 'Loop' in l])
         energies = [float(p.split(' ')[-2]) for p in lins[index_top+1:index_end]]
       rf.close()
       l_a = []
       l_b = []
       l_c = []
       for ef,f in enumerate(['dump.'+str(i) +'.cfg' for i in np.sort([int(f.split('.')[-2]) for f in glob('*.cfg')])]):
         generate_xyz_training(f,energies[ef])
    except:
       print(print_exception())
    os.chdir('../')

def sp_fmt(num):
    return str(np.format_float_scientific(num,precision=4))

def convert_vol_to_latt(v,t='conv'):
    return (4*v)**0.333333

def rmse(e1,e2):
    return np.sqrt(np.sum([(e1[n]-e)**2 for n,e in enumerate(e2)])/len(e1)) 
    
def calculate_eos_plot_from_lammps(jrecs,pfile):
    volumes = []
    energies = []
    dft_energies = []
    pot_forces = []
    dft_forces = []
    for j in jrecs:
        config_nam = j['dft_input_set']['system_name']
        config_latt = j['efs_data']['config_lattice']
        config_specs = j['efs_data']['config_species']
        config_pos = j['efs_data']['config_positions']
        config_forces = j['efs_data']['forces']
        dft_forces += [np.linalg.norm(cf) for cf in config_forces]
        config_energy = j['efs_data']['config_energy']
        lammps_dir = 'LammpsTest_{}'.format(config_nam)
        os.mkdir(lammps_dir)
        s = Structure(lattice=config_latt,coords=config_pos,species=config_specs,coords_are_cartesian=True)
        LammpsData.from_structure(s,atom_style='atomic').write_file(lammps_dir + os.sep + 'run.data')
        call_ext('cp {0} {1}'.format(pfile,lammps_dir))
        if 'GAP' in pfile:
            pot = 'GAP_2b_SOAP.xml'
            name_tag = glob("{}.*".format(pot))[0].split('.')[-1][:-1]
            call_ext("sed -i 's#placetag#{0}#' {1}".format(name_tag,lammps_dir+os.sep+'in.GAP.test'))
        run_lammps(lammps_dir,pfile)
        try:
           con = read(lammps_dir+os.sep+'dump.0with_energy_force.xyz',':')[0]
           energy = con.get_potential_energy()
           force = con.get_array('forces')
           pot_forces += [np.linalg.norm(pf) for pf in force]
           volumes.append(s.volume/len(s))
           energies.append(energy/len(s))
           dft_energies.append(config_energy/len(s))
        except:
           print (print_exception())
    Eos_choices = ['murnaghan','birch','birchmurnaghan']
    data_ev = pd.DataFrame({'Energies':energies,'DFT_Energies':dft_energies,'Volumes':volumes}).sort_values(by=['Volumes'])
    fit_data_ev = data_ev
    for ec in Eos_choices:
       EosC = EquationOfState(volumes=fit_data_ev['Volumes'],energies=fit_data_ev['Energies'],eos=ec)
       EosD = EquationOfState(volumes=fit_data_ev['Volumes'],energies=fit_data_ev['DFT_Energies'],eos=ec)
       e0_pot, e0_err_pot, v0_pot, v0_err_pot, B_pot, B_err_pot, BP_pot, BP_err_pot,E_pred_pot,E_rms_pot = EosC.fit()
       print ("{0} EOS {1}\n".format(pfile,ec), "E0 = {0}+-{1}\n".format(sp_fmt(e0_pot), sp_fmt(e0_err_pot)),\
             "L0 = {0}+-{1}\n".format(sp_fmt(convert_vol_to_latt(v0_pot)), sp_fmt(convert_vol_to_latt(v0_err_pot))),\
             "B = {0}+-{1}\n".format(sp_fmt(B_pot/kJ * 1.0e24), sp_fmt(B_err_pot/kJ * 1.0e24)),\
             "BP = {0}+-{1}\n".format(sp_fmt(BP_pot), sp_fmt(BP_err_pot)),\
             "Fit RMSE {}\n".format(sp_fmt(E_rms_pot)),\
             "Config Energy RMSE = {} \n".format(sp_fmt(rmse(energies,dft_energies))),\
             "Forces RMSE = {} \n".format(sp_fmt(rmse(dft_forces,pot_forces))))
       e0_dft, e0_err_dft, v0_dft, v0_err_dft, B_dft, B_err_dft, BP_dft, BP_err_dft,E_pred_dft,E_rms_dft = EosD.fit()
       print ("DFT EOS {1}\n".format(pfile,ec), "E0 = {0}+-{1}\n".format(sp_fmt(e0_dft), sp_fmt(e0_err_dft)),\
             "L0 = {0}+-{1}\n".format(sp_fmt(convert_vol_to_latt(v0_dft)), sp_fmt(convert_vol_to_latt(v0_err_dft))),\
             "B = {0}+-{1}\n".format(sp_fmt(B_dft/kJ * 1.0e24), sp_fmt(B_err_dft/kJ * 1.0e24)),\
             "BP = {0}+-{1}\n".format(sp_fmt(BP_dft), sp_fmt(BP_err_dft)),\
             "Fit RMSE {}\n".format(sp_fmt(E_rms_dft)))
       #print (str(np.format_float_scientific(E_rms,precision=4)))
       #print (str(np.format_float_scientific((v0*4)**0.333333,precision=4)))
       plt.scatter(fit_data_ev['Volumes'],fit_data_ev['Energies'],label="{}".format(pfile.split('.')[1]),color='blue')
       plt.plot(fit_data_ev['Volumes'],E_pred_pot,label="{0} EOS, RMSE = {1}".format(ec, sp_fmt(E_rms_pot)),color='blue')

       plt.scatter(fit_data_ev['Volumes'],fit_data_ev['DFT_Energies'],label="{}".format('DFT'),color='orange')
       plt.plot(fit_data_ev['Volumes'],E_pred_dft,label="{0} EOS, RMSE = {1}".format(ec, sp_fmt(E_rms_dft)),color='orange')
       
       plt.axvline(v0_pot,label='{0} Equi. Lattice= {1}+-{2}'.format(pfile.split('.')[1],sp_fmt(convert_vol_to_latt(v0_pot)),sp_fmt(convert_vol_to_latt(v0_err_pot))),linestyle=':',color='blue')
       plt.axvline(v0_dft,label='{0} Equi. Lattice= {1}+-{2}'.format('DFT',sp_fmt(convert_vol_to_latt(v0_dft)),sp_fmt(convert_vol_to_latt(v0_err_dft))),linestyle=':',color='orange')
       plt.xlabel('Volume  A^3 /atom')
       plt.ylabel('Energy eV/atom')
       plt.legend()
       plt.savefig('Al_EOS_{0}_{1}.png'.format(ec,pfile))
       plt.close()

def convert_to_atomicrex(struct,forces,cnum,dest='Atomicrex'):
  try:
    pos_file = struct.to(fmt='poscar').split('\n')
    pos_string = [l+'\n' for l in pos_file]
    P = pos_string[0:5] + [str(len(struct))+'\n','Direct\n']
    print (len(pos_string[8:]),len(forces))
    for n,ps in enumerate(pos_string[8:]):
        if ps!='\n':
             P += ' '.join([ps.replace(' Al\n','')]+[str(fi) for fi in forces[n]]) + '\n'
    nam = '{0}/POSCAR_Atomicrex_{1}'.format(dest,cnum.replace('.','_'))
    print ('Writing {}'.format(nam))
    with open(nam,'w') as wf:
         for li in P:
             wf.write(li)
    wf.close()
    with open('{}/structures.xml'.format(dest),'a') as wa:
         lf = ["  <user-structure id='structure_{}'>\n".format(cnum),"    <poscar-file>{}</poscar-file>\n".format(nam.replace(dest+'/','')),\
               "    <properties>\n","      <atomic-forces fit='true' output-all='true'/>\n",\
               "    </properties>\n","  </user-structure>\n\n"]
         for l in lf:
               wa.write(l)
    wa.close()
  except:
    print (print_exception())
    

def convert_to_mtp_cfg(struct,energy,forces,stresses):
    config_step = ['BEGIN_CFG\n',' Size\n','{}\n'.format(len(struct))]
    latt = ['Supercell\n'] + ['\t'.join([str(l) for l in stl])+'\n' for stl in struct.lattice.matrix]
    for la in latt:
        config_step += la
    config_step.append('\t'.join([' AtomData:  id type','cartes_x','cartes_y','cartes_z','fx','fy','fz\n']))
    for n,c in enumerate(struct.cart_coords):
        config_step.append('\t'.join([str(n+1),'0',str(c[0]),str(c[1]),str(c[2]),\
                          str(forces[n][0]),str(forces[n][1]),str(forces[n][2])+'\n']))
    config_step.append(' Energy\n')
    config_step.append('\t{}\n'.format(str(energy)))
    config_step.append('\t'.join([' Stress: xx','yy','zz','xy','yz','xz\n']))
    print (stresses)
    if len(stresses)==3:
       config_step.append('\t'.join([str(sts) for sts in [stresses[0][0],stresses[1][1],stresses[2][2],stresses[0][1],stresses[1][2],stresses[0][2]]])+'\n')
    else:
       config_step.append('\t'.join([str(sts) for sts in stresses])+'\n')
    config_step += [' Feature   EFS_by       VASP\n','END_CFG\n\n']
    with open('TrainMTP.cfg','a') as fw:
         for cs in config_step:
                fw.write(cs)
    fw.close()

def convert_to_xyz(struct,energy,forces,stresses):
    ats = AseAtomsAdaptor.get_atoms(struct)
    sp = SinglePointCalculator(ats,energy=energy,forces=forces)
    ats.set_calculator(sp)
    return ats

def fit_eam(jrecs,dest):
    if not os.path.exists(dest):
       os.mkdir(dest)
    call_ext('cp Atomicrex/*.xml {}'.format(dest))
    call_ext('cp {0}/structures_orig.xml {0}/structures.xml'.format(dest))
    for j in jrecs:
        config_nam = '_'.join([j['dft_input_set']['system_name'],j['unique_id_timestamp']])
        config_latt = j['efs_data']['config_lattice']
        config_specs = j['efs_data']['config_species']
        config_pos = j['efs_data']['config_positions']
        config_forces = j['efs_data']['forces']
        s = Structure(lattice=config_latt,coords=config_pos,species=config_specs,coords_are_cartesian=True)
        convert_to_atomicrex(s,config_forces,config_nam.replace('.','_'),dest)
    
    os.chdir('AtomicrexNew')
    with open('structures.xml','a') as fw:
         fw.write('</group>')
    fw.close()
    call_ext(atomicrex_fit)
    os.chdir('../')

def fit_mtp_gap(jrecs):
    call_ext('rm TrainMTP.cfg TrainGAP.xyz')
    gap_xyz = []
    for j in jrecs:
        config_nam = '_'.join([j['dft_input_set']['system_name'],j['unique_id_timestamp']])
        config_latt = j['efs_data']['config_lattice']
        config_specs = j['efs_data']['config_species']
        config_pos = j['efs_data']['config_positions']
        config_forces = j['efs_data']['forces']
        config_energy = j['efs_data']['config_energy']
        config_stress = j['efs_data']['stresses']
        s = Structure(lattice=config_latt,coords=config_pos,species=config_specs,coords_are_cartesian=True)
        convert_to_mtp_cfg(s,config_energy,config_forces,config_stress)
        gap_xyz.append(convert_to_xyz(s,config_energy,config_forces,config_stress))
    call_ext('{0}'.format(mtp_fit))
    write('TrainGAP.xyz',gap_xyz)
    call_ext('{0}'.format(gap_fit))


if __name__ == '__main__':
    # read json db
    # filter out jrecs containing FCC in system name
    lmp_path = "mpirun -n 4 /home/jgabriel/Software/lammps/src/lmp_mpi"
    mtp_fit = "mpirun -n 1 /home/jgabriel/Software/MTP_Lammps/mlip_joshua/make/mlp run mlip_fit.ini"
    atomicrex_fit = "/home/jgabriel/Software/atomicrex/build/atomicrex main.xml"
    gap_path = "/home/jgabriel/Software/QUIP_GAP/Latest/QUIP/gfortran_bin/gap_fit "
    i1 = 0.5
    i2 = 0.05
    pot_name = "GAP_2b_SOAP.xml"
    gap_fit = gap_path + ' e0={Al:-0.11083469} sparsify_only_no_fit=F sparse_jitter=1e-08 energy_parameter_name=energy force_parameter_name=forces do_copy_at_file=F sparse_separate_file=T gp_file='+'{}'.format(pot_name)+' at_file=TrainGAP.xyz default_sigma={0.008 0.4 0 0} gap={distance_2b cutoff=6 theta_fac=0.5 covariance_type=ard_se n_sparse=30 sparse_method=cur_points delta='+'{}'.format(i1)+' : soap cutoff=6 atom_sigma=0.5 l_max=6 n_max=6 covariance_type=dot_product n_sparse=30 sparse_method=cur_points zeta=4.0 delta='+'{}'.format(i2)+'}'
    call_ext('rm -rf LammpsTest* AtomicrexNew* *.xml*')
    js = json.load(open('TrainingTest_17_March_2020_14_58_23_429598.json'))
    # json.load(open('TrainingTest_17_March_2020_12_24_31_930863.json')) 
    #json.load(open('TrainingTest_17_March_2020_12_00_03_846344.json'))
    #json.load(open('TrainingTest_17_March_2020_11_44_21_552591.json'))#json.load(open('TrainingTest_17_March_2020_11_34_24_488034.json'))#json.load(open('TrainingTest_16_March_2020_14_57_07_854652.json'))
    #jsel =js
    jsel_fit = #[j for j in js if '2x1x1' not in j['dft_input_set']['system_name']]
    jsel_plot = [j for j in js if 'Volume' in j['dft_input_set']['system_name']]   

    fit_eam(jsel_fit,'AtomicrexNew')
    calculate_eos_plot_from_lammps(jsel_plot,'in.EAM.test')
    fit_mtp_gap(jsel_fit)
    call_ext('rm -rf LammpsTest*')
    calculate_eos_plot_from_lammps(jsel_fit,'in.MTP.test')
    call_ext('rm -rf LammpsTest*')
    calculate_eos_plot_from_lammps(jsel_plot,'in.GAP.test')
    
