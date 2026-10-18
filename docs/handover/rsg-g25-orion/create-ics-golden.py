import matplotlib
matplotlib.use('Agg')  # GOLDEN: headless backend (no X on compute nodes)
from stellar_ics.single_star_set import SetInitialConditions
from stellar_ics.single_star import GiantStarMesaController
from stellar_ics.eoswrapper import IdealEOS
from stellar_ics.stellar_model import StellarModel
import stellar_ics.grid as grid

from loadmodules import gadget_readsnapname, gadget_add_grids, gadget_write_ics, gadget_add_grid, gadget_add_box

import shutil
import argparse
import sys
import module
import os
import numpy as np

from const import msol,rsol
import matplotlib.pyplot as plt
import matplotlib as mpl
import mesa
import pickle
from stellar_ics.tools import get_background_grid_from_mesa_surface


# ---- GOLDEN: command-line parameters (defaults = original create-ics.py values) ----
parser = argparse.ArgumentParser(description='Create AREPO single-star IC from a MESA profile')
parser.add_argument('--name', required=True, help='model name, e.g. RSGg_golden4_p407')
parser.add_argument('--profile', required=True, help='MESA profile (text, or pickle with --pickle); relative to ic/ or absolute')
parser.add_argument('--masslimitfactor', type=float, default=1e-4, help='HEALPix grid lower_mass_limit (orig 1e-4)')
parser.add_argument('--grid-ndir', type=int, default=128, help='grid_ndir (orig 128)')
parser.add_argument('--npart', type=float, default=3e5, help='number of particles (orig 3e5)')
parser.add_argument('--fraction-cut', type=float, default=0.03, help='core cut radius / R (orig 0.03)')
parser.add_argument('--pickle', action='store_true', help='profile is a pickled mesa.py object')
cli = parser.parse_args()

# GOLDEN: all relative paths below (EOS5_data, Config.sh, submit_*.sh, <name>/, ../<name>/, ../<name>_relax/)
# are relative to the ic/ directory -> always run from the directory containing this script.
ic_dir = os.path.dirname(os.path.abspath(__file__))
mesa_profile_file = os.path.abspath(cli.profile)
os.chdir(ic_dir)
print('working dir:', os.getcwd())
print('profile    :', mesa_profile_file)

name = cli.name
mesa_in_pickle_format = cli.pickle
arepo_dir = '/orion/u/jinma/AREPO/Arepo'  # GOLDEN: only substituted into ../<name>/<name>-run.cmd (+json); not used to build/run anything
cluster = 'haswell'
email = 'name.surname@mpa-garching.mpg.org' # useless at the moment
number_processes = 240 

npart = cli.npart  # GOLDEN (orig 3e5)
fraction_cut = cli.fraction_cut # GOLDEN (orig 0.03) # 0.1 # 0.05 #
grid_rho = 1e-17 # 2e-14 #
grid_t = 2e2 # 2000.0 #
masslimitfactor = cli.masslimitfactor # GOLDEN (orig 1e-4) # 1e-6 # 3e-4 #1e-3
isotopes = ['h1','he4','c12','n14','o16']#['h1','he4','b11','c12','c13','n13','n14','n15','o15','o16','o17','f18','ne19','ne20','ne21','ne22','na22','na23','mg23','mg24','mg25','mg26','al25','al26','al27','si28','si29','si30','p29','p30','p31','s31','s32','s33','cl33','cl34','cl35','ar36','ar37','ar38','ar39','k39','ca40','sc43','ti44','v47','cr48','mn51','fe52','fe56','co55','ni56','ni58','ni59']
grid_ndir = cli.grid_ndir # GOLDEN (orig 128) # 32 # 

grid_u = grid_xnuc = None
renorm_isotopes = True

if mesa_in_pickle_format:
	mesa_pickle = open(mesa_profile_file, 'rb')
	mesa_profile = pickle.load(mesa_pickle)
else:
	mesa_profile = mesa.mesa_profile(sldir=os.path.dirname(mesa_profile_file), num_type='explicit', give_filename=os.path.basename(mesa_profile_file))  # GOLDEN: works with absolute paths

grid_rho, grid_u, grid_xnuc = get_background_grid_from_mesa_surface(mesa_profile, isotopes, renorm_isotopes, grid_rho=grid_rho, grid_t=grid_t)

print('grid_xnuc')
print(grid_xnuc)

models = {
    name: dict(
                        ic_controller=GiantStarMesaController,
                        profile=mesa_profile,
                        npart=npart,
                        grid_ndir=grid_ndir,
                        arepo_dir=arepo_dir,
                        fraction_cut=fraction_cut,
                        n=3,
                        use_isotopes=True,
                        renorm_isotopes=renorm_isotopes,
                        isotopes=isotopes,
                        factor_gridtoradius=2, # does not do anything at the moment
                        grid=grid.HEALPixGridFactory(grid_rho=grid_rho, grid_u=grid_u, grid_xnuc=grid_xnuc, mode='combined', lower_mass_limit=masslimitfactor, factor_gridtoradius=6.),
                        verbose=1,
                        )
}



model_set_opal = SetInitialConditions.for_opal_run(models)
model_set_opal.create_ics(cluster=cluster, email=email,
                          number_processes=number_processes)
                          
s = gadget_readsnapname( name+"/rhoin-"+name+".ic", hdf5=True)
data = s.data
print(np.min(data['mass']))
print(np.min(data['rho']))
param_file = "../"+name+"/params"
with open(param_file,"r") as f:
    lines = f.readlines()
    for line in lines:
        if 'BoxSize' in line:
            linel = line.strip().split()
            print(linel[1])
            boxsize = float(linel[1])
data['count'] = len(data['mass'])
data['boxsize'] = boxsize
r = s.r()
xnuc = data['xnuc'][r[:-1] > 0.5*boxsize][-1]
#print(xnuc)
gadget_add_box( name+"/rhoin-"+name+".ic", name+"/rhoin-"+name+".ic", boxsize*2, 64, boxdensity=grid_rho, boxenergy=grid_u, xnuc=xnuc, offx = 0.5 * boxsize*2 - 0.5 * data['boxsize'], offy = 0.5 * boxsize*2 - 0.5 * data['boxsize'], offz = 0.5 * boxsize*2 - 0.5 * data['boxsize'], skipBfld=True, writeDensity=True, double=True, format='hdf5', input_hdf5=True )
gadget_add_box( name+"/rhoin-"+name+".ic", name+"/rhoin-"+name+".ic", boxsize*10, 32, boxdensity=grid_rho, boxenergy=grid_u, xnuc=xnuc, offx = 0.5 * boxsize*10 - 0.5 * boxsize*2, offy = 0.5 * boxsize*10 - 0.5 * boxsize*2, offz = 0.5 * boxsize*10 - 0.5 * boxsize*2, skipBfld=True, writeDensity=True, double=True, format='hdf5', input_hdf5=True )
gadget_add_box( name+"/rhoin-"+name+".ic", name+"/rhoin-"+name+".ic", boxsize*50, 32, boxdensity=grid_rho, boxenergy=grid_u, xnuc=xnuc, offx = 0.5 * boxsize*50 - 0.5 * boxsize*10, offy = 0.5 * boxsize*50 - 0.5 * boxsize*10, offz = 0.5 * boxsize*50 - 0.5 * boxsize*10, skipBfld=True, writeDensity=True, double=True, format='hdf5', input_hdf5=True )
#gadget_add_grids( data, [boxsize*2,boxsize*10,boxsize*100], 32*2, xnuc=xnuc )
#gadget_add_grid( data, boxsize*2, 32, xnuc=xnuc )
#gadget_write_ics( name+"/rhoin-"+name+".ic", data, double=True, format="hdf5" )
#s = gadget_readsnapname( name+"/rhoin-"+name+".ic", hdf5=True)
#data = s.data
#data['rho'][data['rho'] == 0.0] = grid_rho
#data['mass'][data['mass'] == 0.0] = grid_rho
#data['u'][data['u'] == 0.0] = grid_u
#print('check mass')
#print(data['mass'])
#print('check type')
#print(data['type'])
#print(data['type'].shape)
#gadget_write_ics( name+"/rhoin-"+name+".ic", data, double=True, format="hdf5" )

#make a stage directory
modeldirectory_path = "../"+name+"_relax"
if not os.path.exists(modeldirectory_path):
    command = 'mkdir -p ' + modeldirectory_path
    os.system(command)

#update parameter file
#param_file = "../"+name+"/params"
value_array = ['rhoin-'+name+'.ic', 'EOS5_data','a09_z0.02_x0.7_ma.data']
name_array = ["InitCondFile", "EosOpalTable","OpacityOpalTable"]
shutil.copy(param_file, modeldirectory_path+"/param.txt")
module.update_param_file(modeldirectory_path+"/param.txt",name_array,value_array)

#copy all files to the newly made directory
shutil.copy(mesa_profile_file, modeldirectory_path)
shutil.copy(name+"/rhoin-"+name+".ic.hdf5", modeldirectory_path)
shutil.copy("Config.sh", modeldirectory_path)
shutil.copy("EOS5_data", modeldirectory_path)
shutil.copy("a09_z0.02_x0.7_ma.data", modeldirectory_path)
shutil.copy("submit_mpa.sh", modeldirectory_path)
shutil.copy("submit_raven.sh", modeldirectory_path)
shutil.copy("submit_nhr.sh", modeldirectory_path)
shutil.copy("submit_orion.sh", modeldirectory_path)
print ("move files to ", modeldirectory_path)

## plot comparisons of pressure and internal energy
### WARNING ##
## this part may give latex errors and give no output. I have found, it works with tap_latex + tap_python. Just load those, install arepo-snap-utils via --user if you really want the plots as Manuel Kramer envisioned them.
#
#ic_mesa = StellarModel.from_mesa_profile(mesa_profile)
#ic_data = model_set_opal.ics[name]
#
##plt.style.use('mesa2.mplstyle')
#
#width = 10
##f, axs = plt.subplots(1, 2, sharex=True, figsize=(width, width*0.6/2))
#fig = plt.figure(figsize=(12,7))
#ax = plt.subplot()
#
#plt.loglog(ic_data.stellar_model.r/rsol, ic_data.stellar_model.rho,
#              label='opal model')
#plt.loglog(ic_data.grid.shells['rad']/rsol, ic_data.grid.shells['rho'], 'x',
#              label='opal shells')
#plt.loglog(ic_mesa.r/rsol, ic_mesa.rho, label='mesa')
#ax.set_xlabel('r [$\\rm Rj_\\odot$]')
#ax.set_ylabel('$\\rm \\rho [g\/cm^{-3}]$')
#plt.legend()
#
##axs[1].loglog(ic_data.stellar_model.r/rsol, ic_data.stellar_model.u,
##              label='ideal model')
##axs[1].loglog(ic_data.grid.shells['rad']/rsol, ic_data.grid.shells['u'], 'x',
##              label='ideal shells')
##axs[1].loglog(ic_data.stellar_model.r/rsol, ic_data.stellar_model.u,
##              label='opal model')
##axs[1].loglog(ic_data.grid.shells['rad']/rsol, ic_data.grid.shells['u'], 'x',
##              label='opal shells')
##axs[1].loglog(ic_mesa.r/rsol, ic_mesa.u, label='mesa')
###axs[1].set_xlabel('r [$\\rm R_\\odot$]')
##axs[1].set_ylabel('u [erg]')
##axs[1].legend()
#
##plt.set_ylim(1e1, 1e23)
##plt.set_ylim(1e11, 1e17)
#
##f.tight_layout()
#fig.savefig('comparison.png')
#
########
##############3
###########################################################3

