#! /bin/bash

set -x 

echo $MODULE_VERSION
source /opt/Modules/3.2.6/init/csh
module purge

module load dot
set intel_ver = 16.0.3.210
module load intel-fc/$intel_ver
module load intel-cc/$intel_ver
module load openmpi/1.10.2

module load netcdf/4.3.3.1
module load hdf5/1.8.14

make -f makefile
