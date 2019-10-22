#! /bin/csh

set echo

echo $MODULE_VERSION
source /opt/cray/pe/modules/default/init/csh

module swap PrgEnv-cray PrgEnv-intel
module load cray-netcdf
module load cray-hdf5

make -f makefile.magnus
