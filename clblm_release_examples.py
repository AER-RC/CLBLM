#clblm_release_examples.py
# copy this script to your clblm directory

import os, sys, glob
import numpy as np
import subprocess as sub
import shutil
import glob

# Set the CLBLM paths
clblm_path=os.getcwd()
print('cwd', clblm_path)
scratch_path=clblm_path+'/scratch'
example_path=clblm_path+'/release_examples/'
data_path=clblm_path+'/clblm_data/'

# array of examples to run
test_Cases = ['AERI','IASI_US_Standard_PWV', 'Solar', 'Jacobian', 'LW_Flux']
# Aeri convolved (with XS and ILS applied), Jacobian H2O, LW_Flux_10-2000 Lamberian emissivity=1, Solar NRL 3-component
out_string = ['AERI_example', 'IASI_example', 'Solar_example', 'jacobian_example','LW_flux_example']




for Case, out in zip(test_Cases, out_string):
   print ('case', Case)

# first delete the scratch directory to avoid any problems
   if os.path.exists(scratch_path):
      shutil.rmtree(scratch_path)
# then remake it
   os.mkdir(scratch_path)

# link to the correct line file
   if os.path.islink(data_path+'spectroscopy/TAPE3'):
      os.remove(data_path+'spectroscopy/TAPE3')
   if Case == "Solar":
      src= clblm_path+'/release_examples/TAPE3_files/TAPE3_aer_v_3.8.1_solar'
   else:
      src= clblm_path+'/release_examples/TAPE3_files/TAPE3_aer_v_3.8.1_ex_ir'
   os.symlink(src,data_path+'spectroscopy/TAPE3')


#for each case, setup the input files
#Cp example scene
   shutil.copyfile(example_path+Case+'/scenes_'+Case+'.nc',clblm_path+'/user_archive/scene_files/scenes.nc')

#Cp example json
   shutil.copyfile(example_path+Case+'/clblm_config.'+Case+'.json',clblm_path+'/clblm_config.json')

#Cp emissivity/reflectivity files, if needed
   if Case == "IASI_US_Standard_PWV":
      shutil.copyfile(example_path+Case+'/EMISSIVITY',clblm_path+'/EMISSIVITY')
      shutil.copyfile(example_path+Case+'/REFLECTIVITY',clblm_path+'/REFLECTIVITY')
      print('copy reflectivity and emissivity files')

#Cp solar irradiance files, if needed
   if Case == "Solar":
      print('copy solar irradiance file')
      shutil.copyfile(example_path+Case+'/SOLAR.RAD.nc',data_path+'solar_irradiance/SOLAR.RAD.nc')

#for each case, Run clblm
   print('Running CLBLM')
   sub.call(["clblm"])

#Cp output to Case directory
   out_path=clblm_path+'/clblm_out/'
   search_string = '%s/'+out+'*.nc'
   nc_outs = sorted(glob.glob(search_string % out_path))
   if len(nc_outs) > 0:
      for file in nc_outs:
         out_2 = file.replace(out,out+'_run')
         shutil.move(file,out_2)
         shutil.move(out_2,example_path+Case)
   if len(nc_outs) == 0:
      print()
      print('CLBLM run failed. See error message(s) above.')
      exit()

print()
print('All example case runs are complete!')
print()