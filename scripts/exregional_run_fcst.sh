#!/bin/bash

#
#-----------------------------------------------------------------------
#
# Source the variable definitions file and the bash utility functions.
#
#-----------------------------------------------------------------------
#
. ${GLOBAL_VAR_DEFNS_FP}
. $USHDIR/source_util_funcs.sh
#
#-----------------------------------------------------------------------
#
# Source other necessary files.
#
#-----------------------------------------------------------------------
#
. $USHDIR/create_model_configure_file.sh
. $USHDIR/create_diag_table_file.sh
#
#-----------------------------------------------------------------------
#
# Save current shell options (in a global array).  Then set new options
# for this script/function.
#
#-----------------------------------------------------------------------
#
{ save_shell_opts; set -u +x; } > /dev/null 2>&1
#
#-----------------------------------------------------------------------
#
# Get the full path to the file in which this script/function is located 
# (scrfunc_fp), the name of that file (scrfunc_fn), and the directory in
# which the file is located (scrfunc_dir).
#
#-----------------------------------------------------------------------
#
scrfunc_fp=$( readlink -f "${BASH_SOURCE[0]}" )
scrfunc_fn=$( basename "${scrfunc_fp}" )
scrfunc_dir=$( dirname "${scrfunc_fp}" )
#
#-----------------------------------------------------------------------
#
# Print message indicating entry into script.
#
#-----------------------------------------------------------------------
#
print_info_msg "
========================================================================
Entering script:  \"${scrfunc_fn}\"
In directory:     \"${scrfunc_dir}\"

This is the ex-script for the task that runs a forecast with FV3 for the
specified cycle.
========================================================================"
#
#-----------------------------------------------------------------------
#
# Specify the set of valid argument names for this script/function.  
# Then process the arguments provided to this script/function (which 
# should consist of a set of name-value pairs of the form arg1="value1",
# etc).
#
#-----------------------------------------------------------------------
#
valid_args=( \
"cdate" \
"cycle_type" \
"cycle_dir" \
"ensmem_indx" \
"slash_ensmem_subdir" \
)
process_args valid_args "$@"
#
#-----------------------------------------------------------------------
#
# For debugging purposes, print out values of arguments passed to this
# script.  Note that these will be printed out only if VERBOSE is set to
# TRUE.
#
#-----------------------------------------------------------------------
#
print_input_args valid_args
#
#-----------------------------------------------------------------------
#
# Load modules.
#
#-----------------------------------------------------------------------
#
case $MACHINE in

  "WCOSS_CRAY")
    ulimit -s unlimited
    ulimit -a

    if [ ${PE_MEMBER01} -gt 24 ];then
      APRUN="aprun -b -j1 -n${PE_MEMBER01} -N24 -d1 -cc depth"
    else
      APRUN="aprun -b -j1 -n24 -N24 -d1 -cc depth"
    fi
    ;;

  "WCOSS_DELL_P3")
    ulimit -s unlimited
    ulimit -a
    APRUN="mpirun -l -np ${PE_MEMBER01}"
    HOMEfv3=/gpfs/dell6/emc/modeling/noscrub/emc.campara/fv3lamda/regional_workflow
    module use ${HOMEfv3}/modulefiles/wcoss_dell_p3.5
    module load regional
    module load fv3
    APRUN="mpirun"
    ;;

  "HERA")
    ulimit -s unlimited
    ulimit -a
    APRUN="srun"
    OMP_NUM_THREADS=4
    ;;

  "ORION")
    ulimit -s unlimited
    ulimit -a
    APRUN="srun"
    ;;

  "JET")
    ulimit -s unlimited
    ulimit -a
    APRUN="srun"
    OMP_NUM_THREADS=2
    ;;

  "ODIN")
    module list
    ulimit -s unlimited
    ulimit -a
    APRUN="srun -n ${PE_MEMBER01}"
    ;;

  "CHEYENNE")
    module list
    nprocs=$(( NNODES_RUN_FCST*PPN_RUN_FCST ))
    APRUN="mpirun -np $nprocs"
    ;;

  "STAMPEDE")
    module list
    APRUN="ibrun -np ${PE_MEMBER01}"
    ;;

  *)
    print_err_msg_exit "\
Run command has not been specified for this machine:
  MACHINE = \"$MACHINE\"
  APRUN = \"$APRUN\""
    ;;

esac
#
#-----------------------------------------------------------------------
#
# Set the forecast run directory.
#
#-----------------------------------------------------------------------
#
run_dir="${cycle_dir}${slash_ensmem_subdir}"
#
#-----------------------------------------------------------------------
#
# Create links in the INPUT subdirectory of the current run directory to
# the grid and (filtered) orography files.
#
#-----------------------------------------------------------------------
#
print_info_msg "$VERBOSE" "
Creating links in the INPUT subdirectory of the current run directory to
the grid and (filtered) orography files ..."


set -x
# Create links to fix files in the FIXLAM directory.


cd_vrfy ${run_dir}/INPUT

relative_or_null=""
if [ "${RUN_TASK_MAKE_GRID}" = "TRUE" ] && [ "${MACHINE}" != "WCOSS_CRAY" ]; then
  relative_or_null="--relative"
fi

# Symlink to mosaic file with a completely different name.
#target="${FIXLAM}/${CRES}${DOT_OR_USCORE}mosaic.halo${NH4}.nc"   # Should this point to this halo4 file or a halo3 file???
#target="${FIXLAM}/${CRES}${DOT_OR_USCORE}mosaic.halo${NH3}.nc"   # Should this point to this halo4 file or a halo3 file???
#for EMC version
target="${FIXLAM}/${CRES}${DOT_OR_USCORE}mosaic.nc"   # Should this point to this halo4 file or a halo3 file???
symlink="grid_spec.nc"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"$target\""
fi

## Symlink to halo-3 grid file with "halo3" stripped from name.
#target="${FIXLAM}/${CRES}${DOT_OR_USCORE}grid.tile${TILE_RGNL}.halo${NH3}.nc"
#if [ "${RUN_TASK_MAKE_SFC_CLIMO}" = "TRUE" ] && \
#   [ "${GRID_GEN_METHOD}" = "GFDLgrid" ] && \
#   [ "${GFDLgrid_USE_GFDLgrid_RES_IN_FILENAMES}" = "FALSE" ]; then
#  symlink="C${GFDLgrid_RES}${DOT_OR_USCORE}grid.tile${TILE_RGNL}.nc"
#else
#  symlink="${CRES}${DOT_OR_USCORE}grid.tile${TILE_RGNL}.nc"
#fi

# Symlink to halo-3 grid file with "halo3" stripped from name.
mosaic_fn="grid_spec.nc"
grid_fn=$( get_charvar_from_netcdf "${mosaic_fn}" "gridfiles" )

target="${FIXLAM}/${grid_fn}"
symlink="${grid_fn}"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"$target\""
fi

# Symlink to halo-4 grid file with "${CRES}_" stripped from name.
#
# If this link is not created, then the code hangs with an error message
# like this:
#
#   check netcdf status=           2
#  NetCDF error No such file or directory
# Stopped
#
# Note that even though the message says "Stopped", the task still con-
# sumes core-hours.
#
target="${FIXLAM}/${CRES}${DOT_OR_USCORE}grid.tile${TILE_RGNL}.halo${NH4}.nc"
symlink="grid.tile${TILE_RGNL}.halo${NH4}.nc"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"$target\""
fi



relative_or_null=""
if [ "${RUN_TASK_MAKE_OROG}" = "TRUE" ] && [ "${MACHINE}" != "WCOSS_CRAY" ] ; then
  relative_or_null="--relative"
fi

# Symlink to halo-0 orography file with "${CRES}_" and "halo0" stripped from name.
target="${FIXLAM}/${CRES}${DOT_OR_USCORE}oro_data.tile${TILE_RGNL}.halo${NH0}.nc"
symlink="oro_data.nc"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"$target\""
fi


#
# Symlink to halo-4 orography file with "${CRES}_" stripped from name.
#
# If this link is not created, then the code hangs with an error message
# like this:
#
#   check netcdf status=           2
#  NetCDF error No such file or directory
# Stopped
#
# Note that even though the message says "Stopped", the task still con-
# sumes core-hours.
#
target="${FIXLAM}/${CRES}${DOT_OR_USCORE}oro_data.tile${TILE_RGNL}.halo${NH4}.nc"
symlink="oro_data.tile${TILE_RGNL}.halo${NH4}.nc"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"$target\""
fi

#
# If using the FV3_HRRR or FV3_RAP physics suites, there are two files 
# (that contain statistics of the orography) that are needed by the gravity 
# wave drag parameterization in that suite.  Below, create symlinks to these 
# files in the run directory.  Note that the symlinks must have specific names 
# that the FV3 model is hardcoded to recognize, and those are the names 
# we use below.
#
if [ "${CCPP_PHYS_SUITE}" = "FV3_HRRR" ] || \
   [ "${CCPP_PHYS_SUITE}" = "FV3_RAP" ]; then


  fileids=( "ss" "ls" )
  for fileid in "${fileids[@]}"; do
    target="${FIXLAM}/${CRES}${DOT_OR_USCORE}oro_data_${fileid}.tile${TILE_RGNL}.halo${NH0}.nc"
    symlink="oro_data_${fileid}.nc"
    if [ -f "${target}" ]; then
      ln_vrfy -sf ${relative_or_null} $target $symlink
    else
      print_err_msg_exit "\
Cannot create symlink because target does not exist:
  target = \"${target}\"
  symlink = \"${symlink}\""
    fi
  done

fi


#
#-----------------------------------------------------------------------
#
# The FV3 model looks for the following files in the INPUT subdirectory
# of the run directory:
#
#   gfs_data.nc
#   sfc_data.nc
#   gfs_bndy*.nc
#   gfs_ctrl.nc
#
# Some of these files (gfs_ctrl.nc, gfs_bndy*.nc) already exist, but
# others do not.  Thus, create links with these names to the appropriate
# files (in this case the initial condition and surface files only).
#
#-----------------------------------------------------------------------
#
print_info_msg "$VERBOSE" "
Creating links with names that FV3 looks for in the INPUT subdirectory
of the current run directory (run_dir), where
  run_dir = \"${run_dir}\"
..."

BKTYPE=1    # cold start using INPUT
if [ -r ${run_dir}/INPUT/fv_tracer.res.tile1.nc ]; then
  BKTYPE=0  # cycling using RESTART
fi
print_info_msg "$VERBOSE" "
The forecast has BKTYPE $BKTYPE (1:cold start ; 0 cycling)"

cd_vrfy ${run_dir}/INPUT
#ln_vrfy -sf gfs_data.tile${TILE_RGNL}.halo${NH0}.nc gfs_data.nc
#ln_vrfy -sf sfc_data.tile${TILE_RGNL}.halo${NH0}.nc sfc_data.nc

relative_or_null=""

if [ ${BKTYPE} -eq 1 ]; then
  target="gfs_data.tile${TILE_RGNL}.halo${NH0}.nc"
else
  target="fv_core.res.tile1.nc"
fi
symlink="gfs_data.nc"
if [ -f "${target}" ]; then
  ln_vrfy -sf ${relative_or_null} $target $symlink
else
  print_err_msg_exit "\
  Cannot create symlink because target does not exist:
  target = \"$target\""
fi

if [ ${BKTYPE} -eq 1 ]; then
  target="sfc_data.tile${TILE_RGNL}.halo${NH0}.nc"
  symlink="sfc_data.nc"
  if [ -f "${target}" ]; then
    ln_vrfy -sf ${relative_or_null} $target $symlink
  else
    print_err_msg_exit "\
    Cannot create symlink because target does not exist:
    target = \"$target\""
  fi
else
  if [ -f "sfc_data.nc" ]; then
    print_info_msg "$VERBOSE" "
    sfc_data.nc is available at INPUT directory"
  else
    print_err_msg_exit "\
    sfc_data.nc is not available for cycling"
  fi
fi

#
#-----------------------------------------------------------------------
#
# Create links in the current run directory to fixed (i.e. static) files
# in the FIXam directory.  These links have names that are set to the
# names of files that the forecast model expects to exist in the current
# working directory when the forecast model executable is called (and
# that is just the run directory).
#
#-----------------------------------------------------------------------
#
set -x
cd_vrfy ${run_dir}


print_info_msg "$VERBOSE" "
Creating links in the current run directory (run_dir) to fixed (i.e.
static) files in the FIXam directory:
  FIXam = \"${FIXam}\"
  run_dir = \"${run_dir}\""

FV3_VER=EMC

#-----------------------------------------------------------------------
#  EMC version fix files
#-----------------------------------------------------------------------

if [ "${FV3_VER}" == "EMC" ]; then
  FIXam=/gpfs/dell6/emc/modeling/noscrub/Shun.Liu/rrfs/fix/fix_am_emc
  #---------------------------------------------- 
  # Copy all the necessary fix files
  #---------------------------------------------- 
  cp $FIXam/global_solarconstant_noaa_an.txt            solarconstant_noaa_an.txt
  cp $FIXam/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77  global_o3prdlos.f77
  cp $FIXam/global_h2o_pltc.f77                         global_h2oprdlos.f77
  cp $FIXam/global_sfc_emissivity_idx.txt               sfc_emissivity_idx.txt
  cp $FIXam/global_co2historicaldata_glob.txt           co2historicaldata_glob.txt
  cp $FIXam/co2monthlycyc.txt                           co2monthlycyc.txt
  cp $FIXam/global_climaeropac_global.txt               aerosol.dat
  
  cp $FIXam/global_glacier.2x2.grb .
  cp $FIXam/global_maxice.2x2.grb .
  cp $FIXam/RTGSST.1982.2012.monthly.clim.grb .
  cp $FIXam/global_snoclim.1.875.grb .
  cp $FIXam/CFSR.SEAICE.1982.2012.monthly.clim.grb .
  cp $FIXam/global_soilmgldas.t1534.3072.1536.grb .
  cp $FIXam/seaice_newland.grb .
  cp $FIXam/global_shdmin.0.144x0.144.grb .
  cp $FIXam/global_shdmax.0.144x0.144.grb .
  
  cp $FIXam/merra2_fix/aeroclim* .
  cp $FIXam/merra2_fix/optics* .
  
  ln -sf $FIXLAM/C3359.maximum_snow_albedo.tile7.halo0.nc C3359.maximum_snow_albedo.tile1.nc
  ln -sf $FIXLAM/C3359.snowfree_albedo.tile7.halo0.nc C3359.snowfree_albedo.tile1.nc
  ln -sf $FIXLAM/C3359.slope_type.tile7.halo0.nc C3359.slope_type.tile1.nc
  ln -sf $FIXLAM/C3359.soil_type.tile7.halo0.nc C3359.soil_type.tile1.nc
  ln -sf $FIXLAM/C3359.vegetation_type.tile7.halo0.nc C3359.vegetation_type.tile1.nc
  ln -sf $FIXLAM/C3359.vegetation_greenness.tile7.halo0.nc C3359.vegetation_greenness.tile1.nc
  ln -sf $FIXLAM/C3359.substrate_temperature.tile7.halo0.nc C3359.substrate_temperature.tile1.nc
  ln -sf $FIXLAM/C3359.facsf.tile7.halo0.nc C3359.facsf.tile1.nc
  
  FIXco2=$FIXam/fix_co2_proj
  for file in `ls $FIXco2/global_co2historicaldata* ` ; do
    cp $file $(echo $(basename $file) |sed -e "s/global_//g")
  done

  #---------------------------------------------- 
  # Copy tile data and orography for regional
  #---------------------------------------------- 
  ntiles=1
  tile=7
  CASE=${CRES}
  cp $FIXLAM/${CASE}_grid.tile${tile}.halo3.nc INPUT/.
  cp $FIXLAM/${CASE}_grid.tile${tile}.halo4.nc INPUT/.
  cp $FIXLAM/${CASE}_oro_data.tile${tile}.halo0.nc INPUT/.
  cp $FIXLAM/${CASE}_oro_data_ls.tile${tile}.halo0.nc INPUT/.
  cp $FIXLAM/${CASE}_oro_data_ss.tile${tile}.halo0.nc INPUT/.
  cp $FIXLAM/${CASE}_oro_data.tile${tile}.halo4.nc INPUT/.
  cp $FIXLAM/${CASE}_mosaic.nc INPUT/.
  
  cd INPUT
  ln -sf ${CASE}_mosaic.nc grid_spec.nc
  ln -sf ${CASE}_grid.tile7.halo3.nc ${CASE}_grid.tile7.nc
  ln -sf ${CASE}_grid.tile7.halo4.nc grid.tile7.halo4.nc
  ln -sf ${CASE}_oro_data.tile7.halo0.nc oro_data.nc
  ln -sf ${CASE}_oro_data.tile7.halo4.nc oro_data.tile7.halo4.nc
  ln -sf ${CASE}_oro_data_ls.tile7.halo0.nc oro_data_ls.nc
  ln -sf ${CASE}_oro_data_ss.tile7.halo0.nc oro_data_ss.nc
  if [ $BKTYPE = 1 ]; then
   # cold start
    ln -sf sfc_data.tile7.halo0.nc sfc_data.nc
    ln -sf gfs_data.tile7.halo0.nc gfs_data.nc
  fi
  cd ..
  
  #-------------------------------------------------------------------
  # Copy or set up files data_table, diag_table, field_table,
  #   input.nml, input_nest02.nml, model_configure, and nems.configure
  #-------------------------------------------------------------------
  
  PARMfv3=/gpfs/dell6/emc/modeling/noscrub/Shun.Liu/rrfs/fix/parm
  MPSUITE=thompson

  CCPP=${CCPP:-"true"}
  
  if [ $MPSUITE = thompson ] ; then
  CCPP_SUITE=${CCPP_SUITE:-"FV3_GFS_v15_thompson_mynn_lam3km"}
   if [ $BKTYPE = 1 ]; then
   # cold start
   cp ${PARMfv3}/thompson/input_sar_firstguess.nml input.nml 
   cp ${PARMfv3}/model_configure_sar_firstguess.tmp model_configure.tmp
   elif [ $BKTYPE = 0 -a ${cycle_type} = spinup ] ; then
   #cp ${PARMfv3}/thompson/input_sar_da_hourly.nml input.nml
   cp ${PARMfv3}/thompson/input_sar_da_hourly.nml_no_gsi_bdy input.nml
   cp ${PARMfv3}/model_configure_sar_da_hourly.tmp model_configure.tmp
   else
   cp ${PARMfv3}/thompson/input_sar_da.nml_no_gsi_bdy input.nml
   cp ${PARMfv3}/model_configure_sar.tmp_writecomp model_configure.tmp
   fi
  cp ${PARMfv3}/thompson/diag_table.tmp .
  cp ${PARMfv3}/thompson/field_table .
  cp $PARMfv3/thompson/CCN_ACTIVATE.BIN                          CCN_ACTIVATE.BIN
  cp $PARMfv3/thompson/freezeH2O.dat                             freezeH2O.dat
  cp $PARMfv3/thompson/qr_acr_qgV2.dat                           qr_acr_qgV2.dat
  cp $PARMfv3/thompson/qr_acr_qsV2.dat                           qr_acr_qsV2.dat
  cp $PARMfv3/thompson/thompson_tables_precomp.sl                thompson_tables_precomp.sl
  fi
  
  cp ${PARMfv3}/fd_nems.yaml .
  cp ${PARMfv3}/data_table .
  cp ${PARMfv3}/nems.configure .

  # processing for inserting GSI into bndy files
  mkdir create_expanded_restart_files_for_DA
  cd create_expanded_restart_files_for_DA
  cp ../field_table .
  cp ../input.nml .
  HOMEfv3=/gpfs/dell6/emc/modeling/noscrub/Shun.Liu/fv3lamda/regional_workflow
  cp $HOMEfv3/regional_da_imbalance/create_expanded_restart_files_for_DA.x .
  ./create_expanded_restart_files_for_DA.x
  mv fv_core.res.tile1_new.nc ../RESTART/.
  mv fv_tracer.res.tile1_new.nc ../RESTART/.
  cd ..
  
  module use /gpfs/dell2/emc/modeling/noscrub/emc.nemspara/soft/modulefiles
  module load esmf/8.1.0bs21
  
  export OMP_NUM_THREADS=2
  nodes=45
  ncnode=40
  let nctsk=ncnode/OMP_NUM_THREADS
  let ntasks=896
  echo nctsk = $nctsk and ntasks = $ntasks
  
  CYCLEtm06=$cdate
  yr=`echo $CYCLEtm06 | cut -c1-4`
  mn=`echo $CYCLEtm06 | cut -c5-6`
  dy=`echo $CYCLEtm06 | cut -c7-8`
  hr=`echo $CYCLEtm06 | cut -c9-10`
  
  if [ $BKTYPE = 1 ]; then
   # cold start
cat > temp << !
${yr}${mn}${dy}.${hr}Z
$yr $mn $dy $hr 0 0
!
  else
cat > temp << !
${yr}${mn}${dy}.${hr}Z.${CASE}.32bit.non-hydro
$yr $mn $dy $hr 0 0
!
  fi
  
  cat temp diag_table.tmp > diag_table


  FCST_LEN_HRS_thiscycle=${FCST_LEN_HRS_CYCLES[${hr}]}

  if [ $BKTYPE = 1 ]; then
  NHRSguess=01	#-- Forecast length for 1st guess generation (now 1-h with GDAS coldstart)
  NFCSTHRS=$NHRSguess	#-- Forecast length for 1st guess generation (now 1-h with GDAS coldstart)
  cat model_configure.tmp | sed s/NTASKS/$ntasks/ | sed s/YR/$yr/ | \
      sed s/MN/$mn/ | sed s/DY/$dy/ | sed s/H_R/$hr/ | \
      sed s/NHRS/$NFCSTHRS/ | sed s/NTHRD/$OMP_NUM_THREADS/ | \
      sed s/NCNODE/$ncnode/  >  model_configure
  elif [ $BKTYPE = 0 -a ${cycle_type} = spinup ] ; then
  NFCSTHRS=01	#-- Forecast length for spinup cycle (now 1-h with GDAS coldstart)
  NRST=01
  cat model_configure.tmp | sed s/NTASKS/$ntasks/ | sed s/YR/$yr/ | \
      sed s/MN/$mn/ | sed s/DY/$dy/ | sed s/H_R/$hr/ | \
      sed s/NHRS/$NFCSTHRS/ | sed s/NTHRD/$OMP_NUM_THREADS/ | \
      sed s/NCNODE/$ncnode/ | sed s/NRESTART/$NRST/ >  model_configure
  else
  #NFCSTHRS=12	#-- Forecast length for spinup cycle (now 1-h with GDAS coldstart)
  NFCSTHRS=$FCST_LEN_HRS_thiscycle  #-- Forecast length for spinup cycle (now 1-h with GDAS coldstart)
  NRST=1
  cat model_configure.tmp | sed s/NTASKS/$ntasks/ | sed s/YR/$yr/ | \
      sed s/MN/$mn/ | sed s/DY/$dy/ | sed s/H_R/$hr/ | \
      sed s/NHRS/$NFCSTHRS/ | sed s/NTHRD/$OMP_NUM_THREADS/ | \
      sed s/NCNODE/$ncnode/ | sed s/NRESTART/$NRST/ >  model_configure
  fi

  #For hyperthreading
  export OMP_NUM_THREADS=4
  
  #Copy UPP parm files for inline post
  cp ${PARMfv3}/post_itag                 ./itag
  cp ${PARMfv3}/nam_micro_lookup.dat      ./eta_micro_lookup.dat
  cp ${PARMfv3}/postxconfig-NT-fv3lam.txt ./postxconfig-NT.txt
  cp ${PARMfv3}/postxconfig-NT-fv3lam.txt ./postxconfig-NT_FH00.txt
  cp ${PARMfv3}/params_grib2_tbl_new      ./params_grib2_tbl_new
  
# export pgm=regional_forecast.x
# . prep_step
  
  #----------------------------------------- 
  # Run the forecast
  #-----------------------------------------
# export pgm=regional_forecast.x
# . prep_step
  
# startmsg
# ${APRUNC} $EXECfv3/regional_forecast.x >$pgmout 2>err
export KMP_STACKSIZE=1024m
export KMP_AFFINITY=disabled
#export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=4
#export OMP_STACKSIZE=2048m
export SENDECF=NO
  EXECfv3=/gpfs/dell6/emc/modeling/noscrub/Shun.Liu/fv3lamda/regional_workflow/exec
  ${APRUN} $EXECfv3/regional_forecast.x >pgmout 2>err
  export err=$?
# export err=$?;err_chk
  
  # copy GRIB2 files
  domain=conus

####################
# calculate tmmark
####################
if [ ${cycle_type} = spinup ]; then
tindx="0 1"
for i in $tindx
do
  cyc_spinstart=${CYCL_HRS_SPINSTART[${i}]}
  cyc_prodstart=${CYCL_HRS_PRODSTART[${i}]}

  if [ ${cyc_prodstart} = 00 ]; then
    cyc_prodstart=24
  fi

  if [ $hr -ge $cyc_spinstart -a $hr -le $cyc_prodstart ]; then
    let tm=cyc_prodstart-hr
    tmmark=tm$(printf "%02d" $tm)
    cyc=${CYCL_HRS_PRODSTART[${i}]}
    echo this_cyc: $hr spinup_start:$cyc_spinstart prod_start:$cyc_prodstart tmmark:$tmmark
    break
  else
    echo this is not spinup cycle, no TM value
  fi
done
else
tmmark=tm00
cyc=$hr
fi
  
  fhour="00 01"
  
 COMOUT=$COMOUT_BASEDIR
 mkdir -p $COMOUT/${cyc}
 domain="COUNS"
 #cyc=run_fcst_spinup_2021081811.log
 #for fhr in $fhour
 
 for ((i=0;i<${NFCSTHRS};i++))
 do
 fhr=$(printf "%02d" $i)
 mv PRSLEV.GrbF${fhr} ${COMOUT}/${cyc}/fv3lam.t${cyc}z.${domain}.f${fhr}.${tmmark}.grib2
 mv NATLEV.GrbF${fhr} ${COMOUT}/${cyc}/fv3lam.t${cyc}z.${domain}.natlev.f${fhr}.${tmmark}.grib2
 done

fi
  exit




#-----------------------------------------------------------------------
# end EMC version fix files
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#  GSL version fix files
#-----------------------------------------------------------------------
if [ "${FV3_VER}" == "GSL"]; then

relative_or_null=""
if [ "${RUN_ENVIR}" != "nco" ] && [ "${MACHINE}" != "WCOSS_CRAY" ] ; then
  relative_or_null="--relative"
fi

regex_search="^[ ]*([^| ]+)[ ]*[|][ ]*([^| ]+)[ ]*$"
num_symlinks=${#CYCLEDIR_LINKS_TO_FIXam_FILES_MAPPING[@]}
for (( i=0; i<${num_symlinks}; i++ )); do

  mapping="${CYCLEDIR_LINKS_TO_FIXam_FILES_MAPPING[$i]}"
  symlink=$( printf "%s\n" "$mapping" | \
             sed -n -r -e "s/${regex_search}/\1/p" )
  target=$( printf "%s\n" "$mapping" | \
            sed -n -r -e "s/${regex_search}/\2/p" )

  symlink="${run_dir}/$symlink"
  target="$FIXam/$target"
  if [ -f "${target}" ]; then
    ln_vrfy -sf ${relative_or_null} $target $symlink
  else
    print_err_msg_exit "\
  Cannot create symlink because target does not exist:
    target = \"$target\""
  fi

done

fi
#-----------------------------------------------------------------------
# end GSL version fix files
#-----------------------------------------------------------------------

if [ "${FV3_VER}" == "GSL"]; then
#
#-----------------------------------------------------------------------
#
# If running this cycle/ensemble member combination more than once (e.g.
# using rocotoboot), remove any time stamp file that may exist from the
# previous attempt.
#
#-----------------------------------------------------------------------
#
cd_vrfy ${run_dir}
rm_vrfy -f time_stamp.out
#
#-----------------------------------------------------------------------
#
# Create links in the current run directory to cycle-independent (and
# ensemble-member-independent) model input files in the main experiment
# directory.
#
#-----------------------------------------------------------------------
#
print_info_msg "$VERBOSE" "
Creating links in the current run directory to cycle-independent model
input files in the main experiment directory..."

relative_or_null=""
if [ "${RUN_ENVIR}" != "nco" ] && [ "${MACHINE}" != "WCOSS_CRAY" ] ; then
  relative_or_null="--relative"
fi

ln_vrfy -sf ${relative_or_null} ${DATA_TABLE_FP} ${run_dir}
ln_vrfy -sf ${relative_or_null} ${FIELD_TABLE_FP} ${run_dir}
ln_vrfy -sf ${relative_or_null} ${NEMS_CONFIG_FP} ${run_dir}

if [ "${DO_ENSEMBLE}" = TRUE ]; then
  ln_vrfy -sf ${relative_or_null} "${FV3_NML_ENSMEM_FPS[$(( 10#${ensmem_indx}-1 ))]}" ${run_dir}/${FV3_NML_FN}
else
   if [ ${BKTYPE} -eq 0 ]; then
    # cycling, using namelist for cycling forecast
    ln_vrfy -sf ${relative_or_null} ${FV3_NML_RESTART_FP} ${run_dir}/input.nml
  else
    # cold start, using namelist for cold start
    ln_vrfy -sf ${relative_or_null} ${FV3_NML_FP} ${run_dir}
  fi
fi
#
#-----------------------------------------------------------------------
#
# Call the function that creates the model configuration file within each
# cycle directory.
#
#-----------------------------------------------------------------------
#
create_model_configure_file \
  cdate="$cdate" \
  cycle_type="$cycle_type" \
  nthreads=${OMP_NUM_THREADS:-1} \
  run_dir="${run_dir}" || print_err_msg_exit "\
Call to function to create a model configuration file for the current
cycle's (cdate) run directory (run_dir) failed:
  cdate = \"${cdate}\"
  run_dir = \"${run_dir}\""

#
#-----------------------------------------------------------------------
#
# Call the function that creates the model configuration file within each
# cycle directory.
#
#-----------------------------------------------------------------------
#
create_diag_table_file \
  run_dir="${run_dir}" || print_err_msg_exit "\
  Call to function to create a diag table file for the current.
cycle's (cdate) run directory (run_dir) failed:
  run_dir = \"${run_dir}\""

#
#-----------------------------------------------------------------------
#
# If running ensemble forecasts, create a link to the cycle-specific
# diagnostic tables file in the cycle directory.  Note that this link
# should not be made if not running ensemble forecasts because in that
# case, the cycle directory is the run directory (and we would be creating
# a symlink with the name of a file that already exists).
#
#-----------------------------------------------------------------------
#
if [ "${DO_ENSEMBLE}" = "TRUE" ]; then
  if [ "${MACHINE}" = "WCOSS_CRAY" ]; then
    relative_or_null=""
  else
    relative_or_null="--relative"
  fi
  diag_table_fp="${cycle_dir}/${DIAG_TABLE_FN}"
  ln_vrfy -sf ${relative_or_null} ${diag_table_fp} ${run_dir}
fi
#
#-----------------------------------------------------------------------
#
# Set and export variables.
#
#-----------------------------------------------------------------------
#
export KMP_AFFINITY=scatter
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} #Needs to be 1 for dynamic build of CCPP with GFDL fast physics, was 2 before.
export OMP_STACKSIZE=1024m

#
#-----------------------------------------------------------------------
#
# If INPUT/phy_data.nc exists, convert it from NetCDF4 to NetCDF3
# (happens for cycled runs, not cold-started)
#
#-----------------------------------------------------------------------
#
cd INPUT
if [[ -f phy_data.nc ]] ; then
  rm -f phy_data.nc3 phy_data.nc4
  cp -fp phy_data.nc phy_data.nc4
  if ( ! time ( module purge ; module load intel szip hdf5 netcdf nco ; module list ; set -x ; ncks -3 --64 phy_data.nc4 phy_data.nc3) ) ; then
    mv -f phy_data.nc4 phy_data.nc
    rm -f phy_data.nc3
    echo "NetCDF 3=>4 conversion failed. :-( Continuing with NetCDF 4 data."
  else
    mv -f phy_data.nc3 phy_data.nc
  fi
fi
cd ..
#
#-----------------------------------------------------------------------
#
# Run the FV3-LAM model.  Note that we have to launch the forecast from
# the current cycle's directory because the FV3 executable will look for
# input files in the current directory.  Since those files have been
# staged in the cycle directory, the current directory must be the cycle
# directory (which it already is).
#
#-----------------------------------------------------------------------
#
$APRUN ${FV3_EXEC_FP} || print_err_msg_exit "\
Call to executable to run FV3-LAM forecast returned with nonzero exit
code."
#
#-----------------------------------------------------------------------
#
# Print message indicating successful completion of script.
#
#-----------------------------------------------------------------------
#
print_info_msg "
========================================================================
FV3 forecast completed successfully!!!

Exiting script:  \"${scrfunc_fn}\"
In directory:    \"${scrfunc_dir}\"
========================================================================"
#===================
# end of GSL version
#===================
fi
#
#-----------------------------------------------------------------------
#
# Restore the shell options saved at the beginning of this script/func-
# tion.
#
#-----------------------------------------------------------------------
#
{ restore_shell_opts; } > /dev/null 2>&1

