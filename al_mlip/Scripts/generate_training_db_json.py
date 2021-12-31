"""
the goal is to create a training database
that contains keys: dft_input_set: {Pymatgen VIS + vasp binary long term, now {INCAR Cutoff, XC, kpra, poscar.comment},
efs_data: {energy, forces and stresses}, soap_vector:[]
unique_id_timestamp
"""
from glob import glob
import sys
import os
import subprocess
from mpinterfaces.utils import print_exception, jobs_from_file
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.ase  import AseAtomsAdaptor
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
import numpy as np
import time
from datetime import datetime
from mpinterfaces.default_logger import get_default_logger
from ase.io import read,write
import json
import re

def convert_to_atomicrex(struct,forces,cnum,dest='Atomicrex'):
  try:
    pos_file = struct.to(fmt='poscar').split('\n')
    pos_string = [l+'\n' for l in pos_file]
    #pprint (pos_string[8:])
    P = pos_string[0:5] + [str(len(struct))+'\n','Direct\n']
    #print (pos_string[8:])
    print (len(pos_string[8:]),len(forces))
    #sys.exit()
    for n,ps in enumerate(pos_string[8:]):
        if ps!='\n':
             P += ' '.join([ps.replace(' Al\n','')]+[str(fi) for fi in forces[n]]) + '\n'
    nam = '{0}/POSCAR_Atomicrex_{1}'.format(dest,cnum)
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

def convert_to_mtp_cfg(outcar_path):
    """
    run mtp converter
    """
    cfg_name = outcar_path.split('/')[-2]+'.cfg'
    os.system("{0} {1} {2}".format(mtp_binary_convert,outcar_path,cfg_name))
    

def extract_vis_efs(mpint_job=None,job_dir=None,cfg_file=None,cfg_set=None):
    """
    extract INCAR Cutoff, XC, kpra, poscar.comment, OUTCAR
    """
    stat = 1
    vis_dict = {key:None for key in ['plane_wave_cutoff','kppra','system_name','ExchangeCorrelation']}
    efs_dict = {key:None for key in ['config_energy','forces','stresses','config_lattice','config_positions','config_species']}
    vis_dict['ExchangeCorrelation'] = 'PBE'
    if mpint_job:
       vis_dict['plane_wave_cutoff'] = mpint_job.vis.as_dict()['incar']['ENCUT']
       jdict = mpint_job.vis.as_dict()
       kp_mesh = jdict['kpoints']['kpoints'][0]
       #print ("Here", jdict)
       pos_length  = len(jdict['poscar']['structure'])
       vis_dict['kppra'] = kp_mesh[0] * kp_mesh[1] * kp_mesh[2] * pos_length
       pos_comment = jdict['poscar']['comment']
       vis_dict['system_name'] = pos_comment

       job_path = mpint_job.parent_job_dir + os.sep + mpint_job.job_dir

       try:
           config_energy, forces, stresses, config_struct = extract_efs_from_vasprun(job_path)
       except:
           stat = 0 
           print (print_exception())
    elif job_dir:
       inc = Incar.from_file(job_dir + os.sep + 'INCAR')
       vis_dict["plane_wave_cutoff"] = inc["ENCUT"]
       kpts = Kpoints.from_file(job_dir + os.sep + 'KPOINTS').as_dict()
       kp_mesh = kpts['kpoints'][0]
       pos  = Poscar.from_file(job_dir + os.sep + 'POSCAR')
       pos_length = len(pos.structure)
       vis_dict['kppra'] = kp_mesh[0] * kp_mesh[1] * kp_mesh[2] * pos_length
       vis_dict['system_name'] = pos.comment
       try:
           config_energy, forces, stresses, config_struct = extract_efs_from_vasprun(job_dir)
       except:
           stat = 0
           print (print_exception())

    elif cfg_file:
       vis_dict['plane_wave_cutoff'] = 550
       vis_dict['kppra'] = 15*15*15*4
       la_set, sp_set, po_set, ns, config_energy, forces, stresses = cfg_set
       config_struct = Structure(la_set, sp_set, po_set, coords_are_cartesian=True)
       vis_dict['system_name'] = ns

    efs_dict['config_energy'] = config_energy
    efs_dict['forces'] = forces
    efs_dict['stresses'] = stresses
    efs_dict['config_lattice'] = [list(m) for m in config_struct.lattice.matrix]
    efs_dict['config_species'] = [str(cs) for cs in config_struct.species]
    efs_dict['config_positions'] = [list(c) for c in config_struct.cart_coords]

    #if not cfg_file:
    #   outs = job_path + os.sep + 'OUTCAR'
    #   convert_to_mtp_cfg(outs)
    soap_vectors = run_SOAP_calc(cutoff,config_struct)
    unique_id = datetime.now().strftime("%d_%B_%Y_%H_%M_%S_%f")
    #convert_to_atomicrex(config_struct,forces,unique_id)

    if stat:
       return vis_dict, efs_dict, unique_id, soap_vectors
    

def extract_efs_from_vasprun(jpath):
    try:
       vasp = Vasprun(jpath + os.sep + 'vasprun.xml').ionic_steps[0]
       config_energy = vasp['e_wo_entrp']
       forces = vasp['forces']
       stresses = vasp['stress']
       config_struct = vasp['structure']
       return config_energy, forces, stresses, config_struct
    except:
       print (print_exception())

def run_soap_unique(struct,cutoff,all_svecs):
    run_SOAP_calc(cutoff,struct=struct)
    counts = []
    with open('soap_calc.out','r') as d:
        lins = d.readlines()
        for nums,l in enumerate(lins):
            add = True
            soap_vec = [float(f) for f in re.findall("-?[\d.]+(?:E-?\d+)?",l) if f!='000']
            if all_svecs:
                for asv in all_svecs:
                    if Dissimilarity2(soap_vec,asv) < eps:
                         add = False
                         break
                if add:
                   all_svecs.append(soap_vec)
                   counts.append(nums)
            else:
                all_svecs.append(soap_vec)
    d.close()
    return counts

def run_SOAP_calc(cutoff,struct=None,atm=None):
    soap_vecs = []
    call_ext('rm SOAP_structure.xyz* soap_calc.out')
    if struct:
        write("SOAP_structure.xyz",AseAtomsAdaptor.get_atoms(struct))
    elif atm:
        write("SOAP_structure.xyz",atm)
    at = "SOAP_structure.xyz"
    call_ext(quip_path+' descriptor_str="soap cutoff={}'.format(cutoff)+' l_max=12 n_max=12 atom_sigma=0.5 n_species=1 cutoff_transition_width=.5 species_Z={13} Z=13 normalize=False" atoms_filename='+\
              '{}'.format(at)+' | grep "DESC" >> soap_calc.out')
    with open('soap_calc.out','r') as d:
        lins = d.readlines()
        for nums,l in enumerate(lins):
            add = True
            soap_vec = [float(f) for f in re.findall("-?[\d.]+(?:E-?\d+)?",l) if f!='000']
            soap_vecs.append(soap_vec)
    return soap_vecs

def call_ext(cmd_string,osc=False):
    log.info("calling {}".format(cmd_string.split(' ')[0]))
    if osc:
        returncode=os.system(cmd_string)
    else:
        returncode = subprocess.call(cmd_string,shell=True)
    if returncode:
          log.info("ERROR in {0}".format(cmd_string))
    else:
          log.info("{} command exited successfully".format(cmd_string.split(' ')[0]))


def convert_mtp_to_poscar(fset=None):
    names = []
    lattice_set = []
    species_set = []
    position_set = []
    stress_set = []
    force_set = []
    energy_set = []
    with open(fset,'r') as rf:
             lins = rf.readlines()
    rf.close()
    config_inits = [nl for nl,l in enumerate(lins) if 'BEGIN_CFG' in l]
    size_inits = [nl+1 for nl,l in enumerate(lins) if 'Size' in l] # line numbers of lines in lins immediately after word Size is found
    cell_inits = [nl+1 for nl,l in enumerate(lins) if 'SuperCell' in l or 'Supercell' in l] # line numbers of lins after word SuperCell is found
    atom_inits = [nl+1 for nl,l in enumerate(lins) if 'AtomData' in l] # after AtomData is found
    energy_inits = [nl+1 for nl,l in enumerate(lins) if 'Energy' in l]
    stress_inits = [nl+1 for nl,l in enumerate(lins) if 'Stress' in l]
    data_individuals = {'Name':[],'ForceMTP':[],'ForceDFT':[]}
    for nci,ci in enumerate(config_inits):
             n_atoms = int(re.findall("-?[\d.]+(?:e-?\d+)?",lins[size_inits[nci]])[0])
             print (nci,n_atoms)
             latt = [[float(p) for p in re.findall("-?[\d.]+(?:e-?\d+)?",ls)] for ls in  lins[cell_inits[nci]:cell_inits[nci]+3]]
             latt_A = np.linalg.norm(latt[0])
             latt_B = np.linalg.norm(latt[1])
             latt_C = np.linalg.norm(latt[2])
             pos_forces = [[float(p) for p in re.findall("-?[\d.]+(?:e-?\d+)?",ls)] for ls in  lins[atom_inits[nci]:atom_inits[nci]+n_atoms]]
             posits = [ps[2:5] for ps in pos_forces] # positions
             force_its = [np.linalg.norm(np.array(ps[5:])) for ps in pos_forces] # magnitude of the forces
             force_arr = [ps[5:] for ps in pos_forces]
             energy_its = float(re.findall("-?[\d.]+(?:e-?\d+)?",lins[energy_inits[nci]])[0])
             stress_its = [float(s) for s in re.findall("-?[\d.]+(?:e-?\d+)?",lins[stress_inits[nci]])]
             # check wrapping correction for Lammps input data file creation
             wrapped_corr_positions = []
             latt_dims = {0:latt_A,1:latt_B,2:latt_C}
             for pc in posits:
                 new_coord = []
                 for d in [0,1,2]:
                    # deal with each coordinate dimension
                    if pc[d] > latt_dims[d]:
                       new_coord.append(pc[d]-latt_dims[d])
                    elif 0 <= pc[d] <= latt_dims[d]:
                       new_coord.append(pc[d])
                    elif pc[d] < 0.0:
                       new_coord.append(pc[d]+latt_dims[d])
                 wrapped_corr_positions.append(np.array(new_coord))
             lattice_set.append(latt)
             species_set.append(['Al' for c in wrapped_corr_positions])
             position_set.append(wrapped_corr_positions)
             energy_set.append(energy_its)
             force_set.append(force_arr)
             stress_set.append(stress_its)
             names.append("FCC_{}".format(str(latt_A)))
    return lattice_set, species_set, position_set, names, energy_set, force_set, stress_set


if __name__ == '__main__':
    """
    1. generate DB json from the vasp directories from glob or mpinterfaces json
       -  init db as a list (or) read from existing json
       -  init template dictionary with the keys
       -  glob or mpinterfaces extract, add information to dict
       -  generate SOAP vector
       -  append to db list
       -  after all directories, write out as json
    2. convert DB record into atomicrex format in the atomicrex folder (extract pymatgen structure, forces from DB record)
    3. convert DB record to mtp cfg from path to OUTCAR
    # add unique dateimt stamp to json saving
    # add support for job dir and mtp cfg to gather training data
    """
    quip_path = '/home/jgabriel/Software/QUIP_GAP/Latest/QUIP/gfortran_bin/quip '
    mtp_binary_convert = 'mpirun -n 1 /home/jgabriel/Software/MTP_Lammps/mlip_joshua/make/mlp convert-cfg --input-format=vasp-outcar '
    js = jobs_from_file('vasp_TrainFCC_Longer.json') + jobs_from_file('vasp_TrainSlabsAgain.json')# for jf in ['vasp_TrainFCC.json','vasp_TrainSlabs.json']] #jobs_from_file('vasp_TrainRelax.json') + jobs_from_file('vasp_index3_SlabsTrainStatic.json')
    jdirlist = []
    jmtp = 'TrainBulk.cfg'
    log = get_default_logger("{}".format("Test_DB"))
    cutoff = 6
    db_records = []
    call_ext('rm Atomicrex/POSCAR*')
    call_ext('cp Atomicrex/structures_orig.xml Atomicrex/structures.xml')
    # first consider TrainBulk.cfg
    #l_set, s_set, p_set, nams, e_set, f_set, st_set = convert_mtp_to_poscar(jmtp)
    #for nl, ms in enumerate(l_set):
    #    db_rec = {key:None for key in ['dft_input_set','efs_data','unique_id_timestamp','soap_vector','soap_vector_input']}
    #    try: 
    #          vis, efs, uid, soap_vec = extract_vis_efs(cfg_file=nl, cfg_set=[N[nl] for N in [l_set, s_set, p_set, nams, e_set, f_set, st_set] ])
    #          stat_code = 1
    #    except:
    #          stat_code = 0
    #          print (print_exception())
    #    if stat_code:
    #        db_rec['dft_input_set'] = vis
    #        db_rec['efs_data'] = efs
    #        db_rec['unique_id_timestamp'] = uid
    #        db_rec['soap_vector_input'] = 'descriptor_str="soap cutoff={}'.format(cutoff)+' l_max=12 n_max=12 atom_sigma=0.5 n_species=1 cutoff_transition_width=.5 species_Z={13} Z=13 normalize=False"'
    #        db_rec['soap_vector'] = soap_vec
    #        print ("DB REC being added")
    #        db_records.append(db_rec)

    # second  consider mpint jobs
    for n,j in enumerate(js):
        db_rec = {key:None for key in ['dft_input_set','efs_data','unique_id_timestamp','soap_vector','soap_vector_input']}
        try: 
            vis, efs, uid, soap_vec= extract_vis_efs(mpint_job=j)
            stat_code = 1
        except:
            stat_code = 0
            print (print_exception())
        if stat_code:
            db_rec['dft_input_set'] = vis
            db_rec['efs_data'] = efs
            db_rec['unique_id_timestamp'] = uid
            db_rec['soap_vector_input'] = 'descriptor_str="soap cutoff={}'.format(cutoff)+' l_max=12 n_max=12 atom_sigma=0.5 n_species=1 cutoff_transition_width=.5 species_Z={13} Z=13 normalize=False"'
            db_rec['soap_vector'] = soap_vec
            print ("DB REC being added")
            db_records.append(db_rec)


    # third consider individual directories
    #for n,j in enumerate(jdirlist):
    #    db_rec = {key:None for key in ['dft_input_set','efs_data','unique_id_timestamp','soap_vector','soap_vector_input']}
    #    try: 
    #        vis, efs, uid, soap_vec= extract_vasp_vis_efs(job_dir=j)
    #        stat_code = 1
    #    except:
    #        stat_code = 0
    #        print (print_exception())
    #    if stat_code:
    #        db_rec['dft_input_set'] = vis
    #        db_rec['efs_data'] = efs
    #        db_rec['unique_id_timestamp'] = uid
    #        db_rec['soap_vector_input'] = 'descriptor_str="soap cutoff={}'.format(cutoff)+' l_max=12 n_max=12 atom_sigma=0.5 n_species=1 cutoff_transition_width=.5 species_Z={13} Z=13 normalize=False"'
    #        db_rec['soap_vector'] = soap_vec
    #        print ("DB REC being added")
    #        db_records.append(db_rec)


    # db communication to pymongo mlab, dump to json document for now
    with open('TrainingTest_{}.json'.format(datetime.now().strftime("%d_%B_%Y_%H_%M_%S_%f")),'w') as jw:
             json.dump(db_records, jw)
    #call_ext('cd Atomicrex/')
    #call_ext('/home/jgabriel/Software/atomicrex/build/atomicrex main.xml')    
