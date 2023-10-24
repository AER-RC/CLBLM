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
Case = ['AERI','IASI_US_Standard_PWV', 'Solar', 'Jacobian', 'LW_Flux']
# Aeri convolved (with XS and ILS applied), Jacobian H2O, LW_Flux_10-2000 Lamberian emissivity=1, Solar NRL 3-component
out = ['AERI_example', 'IASI_example', 'Solar_example', 'Jacobian_example','LW_Flux_example']


index=0

for y in Case:
   print ('case', Case[index])

# first delete the scratch directory to avoid any problems
   shutil.rmtree(scratch_path)
# then remake it
   os.mkdir(scratch_path)

# link to the correct line file
   os.remove(data_path+'spectroscopy/TAPE3')
   if Case[index] == "Solar":
      src= clblm_path+'/release_examples/TAPE3_files/TAPE3_aer_v_3.8.1_solar'
   else:
      src= clblm_path+'/release_xamples/TAPE3_files//TAPE3_LNFL_v3.2_AER_v3.8.1'
   os.symlink(src,data_path+'spectroscopy/TAPE3')

#for each case, setup the input files
#Cp example scene
   shutil.copyfile(example_path+Case[index]+'/scenes_'+Case[index]+'.nc',clblm_path+'/user_archive/scene_files/scenes.nc')

#Cp example json
   shutil.copyfile(example_path+Case[index]+'/clblm_config.'+Case[index]+'.json',clblm_path+'/clblm_config.json')

#Cp emissivity/reflectivity files, if needed
   if Case[index] == "IASI_US_Standard_PWV":
      shutil.copyfile(example_path+Case[index]+'/EMISSIVITY',clblm_path+'/EMISSIVITY')
      shutil.copyfile(example_path+Case[index]+'/REFLECTIVITY',clblm_path+'/REFLECTIVITY')
      print('copy reflectivity and emissivity files')

#Cp solar irradiance files, if needed
   if Case[index] == "Solar":
      print('copy solar irradiance file')
      shutil.copyfile(example_path+Case[index]+'/SOLAR.RAD.nc',data_path+'solar_irradiance/SOLAR.RAD.nc')

#for each case, Run clblm
   print('Running CLBLM')
   sub.call(["clblm"])

#Cp output to Case directory
   out_path=clblm_path+'/clblm_out/'+out[index]+'*.nc'
   for file in glob.glob(out_path):
      shutil.move(file,example_path+Case[index])


   index=index+1




