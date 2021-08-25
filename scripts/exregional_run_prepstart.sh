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
# Save current shell options (in a global array).  Then set new options
# for this script/function.
#
#-----------------------------------------------------------------------
#
{ save_shell_opts; set -u -x; } > /dev/null 2>&1
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

This is the ex-script for the task that runs a analysis with FV3 for the
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
valid_args=( "cycle_dir" "cycle_type" "modelinputdir" "lbcs_root" "fg_root")
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
#-----------------------------------------------------------------------
#
# Extract from CDATE the starting year, month, day, and hour of the
# forecast.  These are needed below for various operations.
#
#-----------------------------------------------------------------------
#
START_DATE=$(echo "${CDATE}" | sed 's/\([[:digit:]]\{2\}\)$/ \1/')

YYYYMMDDHH=$(date +%Y%m%d%H -d "${START_DATE}")
JJJ=$(date +%j -d "${START_DATE}")

YYYY=${YYYYMMDDHH:0:4}
MM=${YYYYMMDDHH:4:2}
DD=${YYYYMMDDHH:6:2}
HH=${YYYYMMDDHH:8:2}
YYYYMMDD=${YYYYMMDDHH:0:8}
#
#-----------------------------------------------------------------------
#
# go to INPUT directory.
# prepare initial conditions for 
#     cold start if BKTYPE=1 
#     warm start if BKTYPE=0
#     spinupcyc + warm start if BKTYPE=2
#       the previous 6 cycles are searched to find the restart files
#       valid at this time from the closet previous cycle.
#
#-----------------------------------------------------------------------

BKTYPE=0
if [ ${cycle_type} == "spinup" ]; then
   echo "spin up cycle"
  for cyc_start in ${CYCL_HRS_SPINSTART[@]}; do
    if [ ${HH} -eq ${cyc_start} ]; then
      BKTYPE=1
    fi
  done
else
  echo " product cycle"
  for cyc_start in ${CYCL_HRS_PRODSTART[@]}; do
    if [ ${HH} -eq ${cyc_start} ]; then
      if [ ${DO_SPINUP} == "true" ]; then
        BKTYPE=2   # using 1-h forecast from spinup cycle
      else
        BKTYPE=1
      fi
    fi
  done
fi

cd_vrfy ${modelinputdir}

if [ ${BKTYPE} -eq 1 ] ; then  # cold start, use prepare cold strat initial files from ics
    bkpath=${cycle_dir}/ics
    if [ -r "${bkpath}/gfs_data.tile7.halo0.nc" ]; then
      cp_vrfy ${bkpath}/gfs_bndy.tile7.000.nc gfs_bndy.tile7.000.nc        
      cp_vrfy ${bkpath}/gfs_ctrl.nc gfs_ctrl.nc        
      cp_vrfy ${bkpath}/gfs_data.tile7.halo0.nc gfs_data.tile7.halo0.nc        
      cp_vrfy ${bkpath}/sfc_data.tile7.halo0.nc sfc_data.tile7.halo0.nc        
      print_info_msg "$VERBOSE" "cold start from $bkpath"
    else
      print_err_msg_exit "Error: cannot find cold start initial condition from : ${bkpath}"
    fi
else

# background directory: 
#    1. warm start from spinup cycle in product cycle. 
#    2. warm start in spinup cycle. 
  if [ ${cycle_type} == "spinup" ] || [ ${BKTYPE} -eq 2 ]; then
     fg_restart_dirname=fcst_fv3lam_spinup
  else
     fg_restart_dirname=fcst_fv3lam
  fi

  YYYYMMDDHHmInterv=$( date +%Y%m%d%H -d "${START_DATE} ${DA_CYCLE_INTERV} hours ago" )
  bkpath=${fg_root}/${YYYYMMDDHHmInterv}/${fg_restart_dirname}/RESTART  # cycling, use background from RESTART

#   let us figure out which backgound is available
#
#   the restart file from FV3 has a name like: ${YYYYMMDD}.${HH}0000.fv_core.res.tile1.nc
#   But the restart files for the forecast length has a name like: fv_core.res.tile1.nc
#   So the defination of restart_prefix needs a "." at the end.
#
  restart_prefix="${YYYYMMDD}.${HH}0000."

# Shun add for EMC version
# if [ ${cycle_type} == "spinup" ]; then
  if [ ${cycle_type} == "spinup" ] || [ ${BKTYPE} -eq 2 ]; then
  restart_prefix=""
  fi

  n=${DA_CYCLE_INTERV}
  while [[ $n -le 6 ]] ; do
    checkfile=${bkpath}/${restart_prefix}fv_core.res.tile1.nc
    if [ -r "${checkfile}" ]; then
      print_info_msg "$VERBOSE" "Found ${checkfile}; Use it as background for analysis "
      break
    else
      n=$((n + ${DA_CYCLE_INTERV}))
      YYYYMMDDHHmInterv=$( date +%Y%m%d%H -d "${START_DATE} ${n} hours ago" )
      bkpath=${fg_root}/${YYYYMMDDHHmInterv}/${fg_restart_dirname}/RESTART  # cycling, use background from RESTART
      if [ ${n} -eq ${FCST_LEN_HRS_SPINUP} ] && [ ${cycle_type} == "spinup" ]; then
        restart_prefix=""
      fi
      print_info_msg "$VERBOSE" "Trying this path: ${bkpath}"
    fi
  done
#
  checkfile=${bkpath}/${restart_prefix}fv_core.res.tile1.nc
  if [ -r "${checkfile}" ]; then
    cp_vrfy ${bkpath}/${restart_prefix}fv_core.res.tile1.nc       fv_core.res.tile1.nc 
    cp_vrfy ${bkpath}/${restart_prefix}fv_tracer.res.tile1.nc     fv_tracer.res.tile1.nc
    cp_vrfy ${bkpath}/${restart_prefix}sfc_data.nc                sfc_data.nc 
    cp_vrfy ${bkpath}/${restart_prefix}coupler.res                coupler.res
    cp_vrfy ${bkpath}/${restart_prefix}fv_core.res.nc             fv_core.res.nc
    cp_vrfy ${bkpath}/${restart_prefix}fv_srf_wnd.res.tile1.nc    fv_srf_wnd.res.tile1.nc
    cp_vrfy ${bkpath}/${restart_prefix}phy_data.nc                phy_data.nc
    cp_vrfy ${fg_root}/${YYYYMMDDHHmInterv}/${fg_restart_dirname}/INPUT/gfs_ctrl.nc  gfs_ctrl.nc
  else
    print_err_msg_exit "Error: cannot find background: ${checkfile}"
  fi
fi

#-----------------------------------------------------------------------
#
# go to INPUT directory.
# prepare boundary conditions:
#       the previous 12 cycles are searched to find the boundary files
#       that can cover the forecast length.
#       The 0-h boundary is copied and others are linked.
#
#-----------------------------------------------------------------------

if [ "${NET}" = "RTMA" ]; then
    #find a bdry file last modified before current cycle time and size > 100M 
    #to make sure it exists and was written out completely. 
    TIME1HAGO=$(date -d "${START_DATE}" +"%Y-%m-%d %H:%M:%S")
    bdryfile0=${lbcs_root}/$(cd $lbcs_root;find . -name "gfs_bndy.tile7.000.nc" ! -newermt "$TIME1HAGO" -size +100M | xargs ls -1rt |tail -n 1)
    bdryfile1=$(echo $bdryfile0 | sed -e "s/gfs_bndy.tile7.000.nc/gfs_bndy.tile7.001.nc/")
    ln_vrfy -snf ${bdryfile0} .
    ln_vrfy -snf ${bdryfile1} .

else
  num_fhrs=( "${#FCST_LEN_HRS_CYCLES[@]}" )
  ihh=$( expr ${HH} + 0 )
  if [ ${num_fhrs} -gt ${ihh} ]; then
     FCST_LEN_HRS_thiscycle=${FCST_LEN_HRS_CYCLES[${ihh}]}
  else
     FCST_LEN_HRS_thiscycle=${FCST_LEN_HRS}
  fi
  if [ ${cycle_type} == "spinup" ]; then
     FCST_LEN_HRS_thiscycle=${FCST_LEN_HRS_SPINUP}
  fi 
  
  print_info_msg "$VERBOSE" " The forecast length for cycle (\"${HH}\") is
                 ( \"${FCST_LEN_HRS_thiscycle}\") "

#   let us figure out which boundary file is available
  bndy_prefix=gfs_bndy.tile7
  n=${EXTRN_MDL_LBCS_SEARCH_OFFSET_HRS}
  end_search_hr=$(( 12 + ${EXTRN_MDL_LBCS_SEARCH_OFFSET_HRS} ))
  YYYYMMDDHHmInterv=$(date +%Y%m%d%H -d "${START_DATE} ${n} hours ago")
  lbcs_path=${lbcs_root}/${YYYYMMDDHHmInterv}/lbcs
  while [[ $n -le ${end_search_hr} ]] ; do
    last_bdy_time=$(( n + ${FCST_LEN_HRS_thiscycle} ))
    last_bdy=$(printf %3.3i $last_bdy_time)
    checkfile=${lbcs_path}/${bndy_prefix}.${last_bdy}.nc
    if [ -r "${checkfile}" ]; then
      print_info_msg "$VERBOSE" "Found ${checkfile}; Use it as boundary for forecast "
      break
    else
      n=$((n + 1))
      YYYYMMDDHHmInterv=$(date +%Y%m%d%H -d "${START_DATE} ${n} hours ago")
      lbcs_path=${lbcs_root}/${YYYYMMDDHHmInterv}/lbcs
    fi
  done
#
  relative_or_null="--relative"
  nb=1
  if [ -r "${checkfile}" ]; then
    while [ $nb -le ${FCST_LEN_HRS_thiscycle} ]
    do
      bdy_time=$(( ${n} + ${nb} ))
      this_bdy=$(printf %3.3i $bdy_time)
      local_bdy=$(printf %3.3i $nb)

      if [ -f "${lbcs_path}/${bndy_prefix}.${this_bdy}.nc" ]; then
        ln_vrfy -sf ${relative_or_null} ${lbcs_path}/${bndy_prefix}.${this_bdy}.nc ${bndy_prefix}.${local_bdy}.nc
      fi

      nb=$((nb + 1))
    done
# check 0-h boundary condition
    if [ ! -f "${bndy_prefix}.000.nc" ]; then
      this_bdy=$(printf %3.3i ${n})
      cp_vrfy ${lbcs_path}/${bndy_prefix}.${this_bdy}.nc ${bndy_prefix}.000.nc 
    fi
  else
    print_err_msg_exit "Error: cannot find boundary file: ${checkfile}"
  fi

fi 
#
#-----------------------------------------------------------------------
#
# Print message indicating successful completion of script.
#
#-----------------------------------------------------------------------
#
print_info_msg "
========================================================================
Prepare start completed successfully!!!

Exiting script:  \"${scrfunc_fn}\"
In directory:    \"${scrfunc_dir}\"
========================================================================"
#
#-----------------------------------------------------------------------
#
# Restore the shell options saved at the beginning of this script/func-
# tion.
#
#-----------------------------------------------------------------------
#
{ restore_shell_opts; } > /dev/null 2>&1

