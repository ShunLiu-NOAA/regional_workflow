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
valid_args=( "cycle_dir" "cycle_type" "analworkdir" )
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
#
"WCOSS_C" | "WCOSS")
#
  module load NCO/4.7.0
  module list
  ulimit -s unlimited
  ulimit -a
  APRUN="mpirun -l -np ${PE_MEMBER01}"
  ;;
#
"WCOSS_DELL_P3")
#
  module load NCO/4.7.0
  module list
  ulimit -s unlimited
  ulimit -a
  #APRUN="mpirun -l -np ${PE_MEMBER01}"
  #APRUN="mpirun -l -np 384"
  HOMEfv3=/gpfs/dell2/emc/modeling/noscrub/emc.campara/fv3lamdax/regional_workflow
  module use ${HOMEfv3}/modulefiles/wcoss_dell_p3.5 
  module load regional
  APRUN="mpirun"
  ;;
#
"THEIA")
#
  ulimit -s unlimited
  ulimit -a
  np=${SLURM_NTASKS}
  APRUN="mpirun -np ${np}"
  ;;
#
"HERA")
  ulimit -s unlimited
  ulimit -a
  export OMP_NUM_THREADS=1
  export OMP_STACKSIZE=300M
  APRUN="srun"
  ;;
#
"ORION")
  ulimit -s unlimited
  ulimit -a
  export OMP_NUM_THREADS=1
  export OMP_STACKSIZE=1024M
  APRUN="srun"
  ;;
#
"JET")
  ulimit -s unlimited
  ulimit -a
  APRUN="srun"
  ;;
#
"ODIN")
#
  module list

  ulimit -s unlimited
  ulimit -a
  APRUN="srun -n ${PE_MEMBER01}"
  ;;
#
esac
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
# go to working directory.
# define fix and background path
#
#-----------------------------------------------------------------------

cd_vrfy ${analworkdir}

fixgriddir=$FIX_GSI/${PREDEF_GRID_NAME}
if [ ${cycle_type} == "spinup" ]; then
  bkpath=${cycle_dir}/fcst_fv3lam_spinup/INPUT
else
  bkpath=${cycle_dir}/fcst_fv3lam/INPUT
fi
# decide background type
if [ -r "${bkpath}/phy_data.nc" ]; then
  BKTYPE=0              # warm start
else
  BKTYPE=1              # cold start
fi

print_info_msg "$VERBOSE" "FIX_GSI is $FIX_GSI"
print_info_msg "$VERBOSE" "fixgriddir is $fixgriddir"
print_info_msg "$VERBOSE" "default bkpath is $bkpath"
print_info_msg "$VERBOSE" "background type is is $BKTYPE"


#-----------------------------------------------------------------------
#
# Make a list of the latest GFS EnKF ensemble
#
#-----------------------------------------------------------------------

stampcycle=$(date -d "${START_DATE}" +%s)
minHourDiff=100
loops="006"    # or 009s for GFSv15
ens_type="nc"  # or nemsio for GFSv15
foundens="false"
cat "no ens found" >> filelist03

case $MACHINE in

"WCOSS_C" | "WCOSS" | "WCOSS_DELL_P3")

  for loop in $loops; do
    for timelist in $(ls ${ENKF_FCST}/enkfgdas.*/*/atmos/mem080/gdas*.atmf${loop}.${ens_type}); do
      availtimeyyyymmdd=$(echo ${timelist} | cut -d'/' -f9 | cut -c 10-17)
      availtimehh=$(echo ${timelist} | cut -d'/' -f10)
      availtime=${availtimeyyyymmdd}${availtimehh}
      avail_time=$(echo "${availtime}" | sed 's/\([[:digit:]]\{2\}\)$/ \1/')
      avail_time=$(date -d "${avail_time}")

      stamp_avail=$(date -d "${avail_time} ${loop} hours" +%s)

      hourDiff=$(echo "($stampcycle - $stamp_avail) / (60 * 60 )" | bc);
      if [[ ${stampcycle} -lt ${stamp_avail} ]]; then
         hourDiff=$(echo "($stamp_avail - $stampcycle) / (60 * 60 )" | bc);
      fi

      if [[ ${hourDiff} -lt ${minHourDiff} ]]; then
         minHourDiff=${hourDiff}
         enkfcstname=gdas.t${availtimehh}z.atmf${loop}
         eyyyymmdd=$(echo ${availtime} | cut -c1-8)
         ehh=$(echo ${availtime} | cut -c9-10)
         foundens="true"
      fi
    done
  done

  if [ ${foundens} = "true" ]
  then
    ls ${ENKF_FCST}/enkfgdas.${eyyyymmdd}/${ehh}/atmos/mem???/${enkfcstname}.nc > filelist03
  fi

  ;;
"JET")

  for loop in $loops; do
    for timelist in $(ls ${ENKF_FCST}/*.gdas.t*z.atmf${loop}.mem080.${ens_type}); do
      availtimeyy=$(basename ${timelist} | cut -c 1-2)
      availtimeyyyy=20${availtimeyy}
      availtimejjj=$(basename ${timelist} | cut -c 3-5)
      availtimemm=$(date -d "${availtimeyyyy}0101 +$(( 10#${availtimejjj} - 1 )) days" +%m)
      availtimedd=$(date -d "${availtimeyyyy}0101 +$(( 10#${availtimejjj} - 1 )) days" +%d)
      availtimehh=$(basename ${timelist} | cut -c 6-7)
      availtime=${availtimeyyyy}${availtimemm}${availtimedd}${availtimehh}
      avail_time=$(echo "${availtime}" | sed 's/\([[:digit:]]\{2\}\)$/ \1/')
      avail_time=$(date -d "${avail_time}")

      stamp_avail=$(date -d "${avail_time} ${loop} hours" +%s)

      hourDiff=$(echo "($stampcycle - $stamp_avail) / (60 * 60 )" | bc);
      if [[ ${stampcycle} -lt ${stamp_avail} ]]; then
         hourDiff=$(echo "($stamp_avail - $stampcycle) / (60 * 60 )" | bc);
      fi

      if [[ ${hourDiff} -lt ${minHourDiff} ]]; then
         minHourDiff=${hourDiff}
         enkfcstname=${availtimeyy}${availtimejjj}${availtimehh}00.gdas.t${availtimehh}z.atmf${loop}
         foundens="true"
      fi
    done
  done

  if [ $foundens = "true" ]; then
    ls ${ENKF_FCST}/${enkfcstname}.mem0??.${ens_type} >> filelist03
  fi

esac

#
#-----------------------------------------------------------------------
#
# set default values for namelist
#
#-----------------------------------------------------------------------

cloudanalysistype=0
ifsatbufr=.false.
ifsoilnudge=.false.
ifhyb=.false.

# Determine if hybrid option is available
memname='atmf009'
nummem=$(more filelist03 | wc -l)
nummem=$((nummem - 3 ))
if [[ ${nummem} -eq 80 ]]; then
  print_info_msg "$VERBOSE" "Do hybrid with ${memname}"
  ifhyb=.true.
  print_info_msg "$VERBOSE" " Cycle ${YYYYMMDDHH}: GSI hybrid uses ${memname} with n_ens=${nummem}" 
fi

#
#-----------------------------------------------------------------------
#
# link or copy background and grib configuration files
#
#  Using ncks to add phis (terrain) into cold start input background. 
#           it is better to change GSI to use the terrain from fix file.
#  Adding radar_tten array to fv3_tracer. Should remove this after add this array in
#           radar_tten converting code.
#-----------------------------------------------------------------------
dothis="true"

if [ $dothis = "true" ]; then
#cp_vrfy ${fixgriddir}/fv3_akbk                     fv3_akbk
#cp_vrfy ${fixgriddir}/fv3_grid_spec                fv3_grid_spec

if [ ${BKTYPE} -eq 1 ]; then  # cold start uses background from INPUT
  cp_vrfy ${fixgriddir}/fix_temp_halo3/fv_core.res.nc     fv3_akbk
  cp_vrfy ${fixgriddir}/fix_temp_halo3/grid_spec.nc       fv3_grid_spec
  cp_vrfy ${bkpath}/gfs_data.tile7.halo0.nc        gfs_data.tile7.halo0.nc_b
  #ncks -A -v  phis ${fixgriddir}/phis.nc           gfs_data.tile7.halo0.nc_b
  ncks -A -v  phis ${fixgriddir}/fv3_dynvars        gfs_data.tile7.halo0.nc_b

  cp_vrfy ${bkpath}/sfc_data.tile7.halo0.nc        fv3_sfcdata
  cp_vrfy gfs_data.tile7.halo0.nc_b                fv3_dynvars
# cp /gpfs/dell2/ptmp/emc.campara/fv3lamdax/fv3lamdax.20211019/12/guess.tm06/sfc_data.tile7.nc fv3_sfcdata
# cp /gpfs/dell2/ptmp/emc.campara/fv3lamdax/fv3lamdax.20211019/12/guess.tm06/gfs_data.tile7.nc fv3_dynvars
  ln_vrfy -s fv3_dynvars                           fv3_tracer

  fv3lam_bg_opt=1
else                          # cycle uses background from restart
  cp_vrfy ${fixgriddir}/fix_temp_halo3/fv_core.res.nc     fv3_akbk
  cp_vrfy ${fixgriddir}/fix_temp_halo3/grid_spec.nc       fv3_grid_spec

  cp_vrfy  ${bkpath}/fv_core.res.tile1.nc             fv3_dynvars
  cp_vrfy  ${bkpath}/fv_tracer.res.tile1.nc           fv3_tracer
  cp_vrfy  ${bkpath}/sfc_data.nc                      fv3_sfcdata
  fv3lam_bg_opt=0
fi

# update times in coupler.res to current cycle time
cp_vrfy ${fixgriddir}/fv3_coupler.res          coupler.res
sed -i "s/yyyy/${YYYY}/" coupler.res
sed -i "s/mm/${MM}/"     coupler.res
sed -i "s/dd/${DD}/"     coupler.res
sed -i "s/hh/${HH}/"     coupler.res

fi

#
#-----------------------------------------------------------------------
#
# link observation files
# copy observation files to working directory 
#
#-----------------------------------------------------------------------
obs_source=rap
if [[ ${HH} -eq '00' || ${HH} -eq '12' ]]; then
  obs_source=rap_e
fi

case $MACHINE in

"WCOSS_C" | "WCOSS" | "WCOSS_DELL_P3")
   obsfileprefix=${obs_source}
   obspath_tmp=${OBSPATH}/${obs_source}.${YYYYMMDD}
  ;;
"JET" | "HERA")
   obsfileprefix=${YYYYMMDDHH}.${obs_source}
   obspath_tmp=${OBSPATH}
  ;;
"ORION" )
   obs_source=rap
   #obsfileprefix=${YYYYMMDDHH}.${obs_source}               # rap observation from JET.
   obsfileprefix=${obs_source}.${YYYYMMDD}/${obs_source}    # observation from operation.
   obspath_tmp=${OBSPATH}
  ;;
*)
   obsfileprefix=${obs_source}
   obspath_tmp=${OBSPATH}
esac


obs_files_source[0]=${obspath_tmp}/${obsfileprefix}.t${HH}z.prepbufr.tm00
obs_files_target[0]=prepbufr

obs_files_source[1]=${obspath_tmp}/${obsfileprefix}.t${HH}z.satwnd.tm00.bufr_d
obs_files_target[1]=satwndbufr

obs_files_source[2]=${obspath_tmp}/${obsfileprefix}.t${HH}z.nexrad.tm00.bufr_d
obs_files_target[2]=l2rwbufr

obs_number=${#obs_files_source[@]}
for (( i=0; i<${obs_number}; i++ ));
do
  obs_file=${obs_files_source[$i]}
  obs_file_t=${obs_files_target[$i]}
  if [ -r "${obs_file}" ]; then
    ln -s "${obs_file}" "${obs_file_t}"
  else
    print_info_msg "$VERBOSE" "Warning: ${obs_file} does not exist!"
  fi
done

####################
# calculate tmmark
####################
hr=$HH
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

 tm=`echo $tmmark | cut -c3-4`
 thistime=`$NDATE +${tm} ${CDATE}`

 CYCrun=`echo $thistime | cut -c9-10`

 COMOUT=$COMOUT_BASEDIR/${thistime}
 mkdir -p $COMOUT

 COMINgfs=/gpfs/dell1/nco/ops/com/gfs/prod
 COMINgdas=/gpfs/dell1/nco/ops/com/gfs/prod
 COMINnam=/gpfs/dell1/nco/ops/com/nam/prod
 COMINrap=/gpfs/hps/nco/ops/com/rap/prod
 COMINpararap=/gpfs/hps/nco/ops/com/rap/para
 COMINbias=$COMOUT
 COMINrtma=/gpfs/dell2/nco/ops/com/rtma/prod

#-----------------------------------------------------------------------
#
# Create links to fix files in the FIXgsi directory.
# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (regional only)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)
#
#-----------------------------------------------------------------------


FV3_VER=EMC
#if [ "${FV3_VER}" == "EMC" ]; then
PARMfv3=/gpfs/dell6/emc/modeling/noscrub/emc.campara/Shun.Liu/rrfs/fix/parm_RRFS_NA_3km
EXECfv3=/gpfs/dell6/emc/modeling/noscrub/emc.campara/Shun.Liu/rrfs/ufs-srweather-app/bin
export fixgsi=/gpfs/dell6/emc/modeling/noscrub/Eric.Rogers/fv3lam_for_dellp3.5/sorc/regional_gsi.fd/fix
export fixcrtm=/usrx/local/nceplibs/dev/NCEPLIBS/src/crtm_v2.3.0/fix

# Set runtime and save directories
export endianness=Big_Endian

# Set variables used in script
#   ncp is cp replacement, currently keep as /bin/cp
ncp=/bin/cp

# Utility to cat NETCDF diag files from all processor
export CATEXEC=$EXECfv3/nc_diag_cat_serial.x

export HYB_ENS=".true."
#export HYB_ENS=".false."

# Get Fv3GDAS Enkf files
# We expect 81 total files to be present (80 enkf + 1 mean)
export nens=81

# Not using FGAT or 4DEnVar, so hardwire nhr_assimilation to 3
export nhr_assimilation=03
##typeset -Z2 nhr_assimilation

UTILfv3=/gpfs/dell6/emc/modeling/noscrub/emc.campara/Shun.Liu/rrfs/ufs-srweather-app/regional_workflow/util/ush
vlddate=$CDATE
python $UTILfv3/getbest_EnKF_FV3GDAS.py -v $vlddate --exact=no --minsize=${nens} -d ${COMINgfs}/enkfgdas -o filelist${nhr_assimilation} --o3fname=gfs_sigf${nhr_assimilation} --gfs_netcdf=yes

#Check to see if ensembles were found 
numfiles=`cat filelist03 | wc -l`

if [ $numfiles -ne 81 ]; then
  echo "Ensembles not found - turning off HYBENS!"
  export HYB_ENS=".false."
fi


echo "HYB_ENS=$HYB_ENS" > $COMOUT/fv3lam.t${CYCrun}z.hybens.${tmmark}

nens=`cat filelist03 | wc -l`
mv filelist03 $COMOUT/${RUN}.t${CYCrun}z.filelist03.${tmmark}
ls ens_* >>filelist03

# Set parameters
export USEGFSO3=.false.
export nhr_assimilation=3
export vs=1.
export fstat=.false.
export i_gsdcldanal_type=0

# Diagnostic files options
export lobsdiag_forenkf=.false.
export netcdf_diag=.true.
export binary_diag=.false.
DIAG_SUFFIX=${DIAG_SUFFIX:-""}
if [ $netcdf_diag = ".true." ] ; then
   DIAG_SUFFIX="${DIAG_SUFFIX}.nc4"
fi

bg_opt=1
nst_gsi=0

# Make gsi namelist
cat << EOF > gsiparm.anl

 &SETUP
   miter=2,niter(1)=50,niter(2)=50,niter_no_qc(1)=20,
   write_diag(1)=.true.,write_diag(2)=.false.,write_diag(3)=.true.,
   gencode=78,qoption=2,
   factqmin=0.0,factqmax=0.0,
   iguess=-1,use_gfs_ozone=${USEGFSO3},
   oneobtest=.false.,retrieval=.false.,
   nhr_assimilation=${nhr_assimilation},l_foto=.false.,
   use_pbl=.false.,gpstop=30.,
   use_gfs_nemsio=.false.,use_gfs_ncio=.true.,
   print_diag_pcg=.true.,
   newpc4pred=.true., adp_anglebc=.true., angord=4,
   passive_bc=.true., use_edges=.false., emiss_bc=.true.,
   diag_precon=.true., step_start=1.e-3, l_reg_update_hydro_delz=.true.,
   netcdf_diag=$netcdf_diag,binary_diag=$binary_diag,
 /
 &GRIDOPTS
   fv3_regional=.true.,grid_ratio_fv3_regional=3.0,nvege_type=20,
 /
 &BKGERR
   hzscl=0.373,0.746,1.50,
   vs=${vs},bw=0.,fstat=${fstat},
 /
 &ANBKGERR
   anisotropic=.false.,
 /
 &JCOPTS
 /
 &STRONGOPTS
   nstrong=0,nvmodes_keep=20,period_max=3.,
   baldiag_full=.true.,baldiag_inc=.true.,
 /
 &OBSQC
   dfact=0.75,dfact1=3.0,noiqc=.false.,c_varqc=0.02,
   vadfile='prepbufr',njqc=.false.,vqc=.true.,
   aircraft_t_bc=.true.,biaspredt=1000.0,upd_aircraft=.true.,cleanup_tail=.true.,
 /
 &OBS_INPUT
   dmesh(1)=120.0,dmesh(2)=60.0,time_window_max=1.5,ext_sonde=.true.,
 /
OBS_INPUT::
!  dfile          dtype       dplat     dsis                 dval    dthin dsfcalc
   prepbufr       ps          null      ps                   1.0     0     0
   prepbufr       t           null      t                    1.0     0     0
   prepbufr       q           null      q                    1.0     0     0
   prepbufr       pw          null      pw                   1.0     0     0
   prepbufr_profl t           null      t                    1.0     0     0
   prepbufr_profl q           null      q                    1.0     0     0
   prepbufr_profl uv          null      uv                   1.0     0     0
   satwndbufr     uv          null      uv                   1.0     0     0
   prepbufr       uv          null      uv                   1.0     0     0
   prepbufr       spd         null      spd                  1.0     0     0
   prepbufr       dw          null      dw                   1.0     0     0
   l2rwbufr       rw          null      l2rw                 1.0     0     0
   prepbufr       sst         null      sst                  1.0     0     0
   gpsrobufr      gps_ref     null      gps                  1.0     0     0
   ssmirrbufr     pcp_ssmi    dmsp      pcp_ssmi             1.0    -1     0
   tmirrbufr      pcp_tmi     trmm      pcp_tmi              1.0    -1     0
   sbuvbufr       sbuv2       n16       sbuv8_n16            0.0     0     0
   sbuvbufr       sbuv2       n17       sbuv8_n17            0.0     0     0
   sbuvbufr       sbuv2       n18       sbuv8_n18            0.0     0     0
   hirs3bufr      hirs3       n16       hirs3_n16            0.0     1     0
   hirs3bufr      hirs3       n17       hirs3_n17            0.0     1     0
   hirs4bufr      hirs4       metop-a   hirs4_metop-a        0.0     2     0
   hirs4bufr      hirs4       n18       hirs4_n18            0.0     1     0
   hirs4bufr      hirs4       n19       hirs4_n19            0.0     2     0
   hirs4bufr      hirs4       metop-b   hirs4_metop-b        0.0     2     0
   gimgrbufr      goes_img    g11       imgr_g11             0.0     1     0
   gimgrbufr      goes_img    g12       imgr_g12             0.0     1     0
   airsbufr       airs        aqua      airs_aqua            0.0     2     0
   amsuabufr      amsua       n15       amsua_n15            0.0     2     0
   amsuabufr      amsua       n18       amsua_n18            0.0     2     0
   amsuabufr      amsua       n19       amsua_n19            0.0     2     0
   amsuabufr      amsua       metop-a   amsua_metop-a        0.0     2     0
   amsuabufr      amsua       metop-b   amsua_metop-b        0.0     2     0
   airsbufr       amsua       aqua      amsua_aqua           0.0     2     0
   amsubbufr      amsub       n17       amsub_n17            0.0     1     0
   mhsbufr        mhs         n18       mhs_n18              0.0     2     0
   mhsbufr        mhs         n19       mhs_n19              0.0     2     0
   mhsbufr        mhs         metop-a   mhs_metop-a          0.0     2     0
   mhsbufr        mhs         metop-b   mhs_metop-b          0.0     2     0
   ssmitbufr      ssmi        f13       ssmi_f13             0.0     2     0
   ssmitbufr      ssmi        f14       ssmi_f14             0.0     2     0
   ssmitbufr      ssmi        f15       ssmi_f15             0.0     2     0
   amsrebufr      amsre_low   aqua      amsre_aqua           0.0     2     0
   amsrebufr      amsre_mid   aqua      amsre_aqua           0.0     2     0
   amsrebufr      amsre_hig   aqua      amsre_aqua           0.0     2     0
   ssmisbufr      ssmis       f16       ssmis_f16            0.0     2     0
   ssmisbufr      ssmis       f17       ssmis_f17            0.0     2     0
   ssmisbufr      ssmis       f18       ssmis_f18            0.0     2     0
   ssmisbufr      ssmis       f19       ssmis_f19            0.0     2     0
   gsnd1bufr      sndrd1      g12       sndrD1_g12           0.0     1     0
   gsnd1bufr      sndrd2      g12       sndrD2_g12           0.0     1     0
   gsnd1bufr      sndrd3      g12       sndrD3_g12           0.0     1     0
   gsnd1bufr      sndrd4      g12       sndrD4_g12           0.0     1     0
   gsnd1bufr      sndrd1      g11       sndrD1_g11           0.0     1     0
   gsnd1bufr      sndrd2      g11       sndrD2_g11           0.0     1     0
   gsnd1bufr      sndrd3      g11       sndrD3_g11           0.0     1     0
   gsnd1bufr      sndrd4      g11       sndrD4_g11           0.0     1     0
   gsnd1bufr      sndrd1      g13       sndrD1_g13           0.0     1     0
   gsnd1bufr      sndrd2      g13       sndrD2_g13           0.0     1     0
   gsnd1bufr      sndrd3      g13       sndrD3_g13           0.0     1     0
   gsnd1bufr      sndrd4      g13       sndrD4_g13           0.0     1     0
   gsnd1bufr      sndrd1      g15       sndrD1_g15           0.0     2     0
   gsnd1bufr      sndrd2      g15       sndrD2_g15           0.0     2     0
   gsnd1bufr      sndrd3      g15       sndrD3_g15           0.0     2     0
   gsnd1bufr      sndrd4      g15       sndrD4_g15           0.0     2     0
   iasibufr       iasi        metop-a   iasi_metop-a         0.0     2     0
   gomebufr       gome        metop-a   gome_metop-a         0.0     2     0
   omibufr        omi         aura      omi_aura             0.0     2     0
   sbuvbufr       sbuv2       n19       sbuv8_n19            0.0     0     0
   tcvitl         tcp         null      tcp                  0.0     0     0
   seviribufr     seviri      m08       seviri_m08           0.0     2     0
   seviribufr     seviri      m09       seviri_m09           0.0     2     0
   seviribufr     seviri      m10       seviri_m10           0.0     2     0
   iasibufr       iasi        metop-b   iasi_metop-b         0.0     2     0
   gomebufr       gome        metop-b   gome_metop-b         0.0     2     0
   atmsbufr       atms        npp       atms_npp             0.0     2     0
   atmsbufr       atms        n20       atms_n20             0.0     2     0
   crisbufr       cris        npp       cris_npp             0.0     2     0
   crisfsbufr     cris-fsr    npp       cris-fsr_npp         0.0     2     0 
   crisfsbufr     cris-fsr    n20       cris-fsr_n20         0.0     2     0 
   abibufr        abi         g16       abi_g16              0.0     2     0
   abibufr        abi         g17       abi_g17              0.0     2     0
   mlsbufr        mls30       aura      mls30_aura           0.0     0     0
   oscatbufr      uv          null      uv                   0.0     0     0
   prepbufr       mta_cld     null      mta_cld              1.0     0     0
   prepbufr       gos_ctp     null      gos_ctp              1.0     0     0
   refInGSI       rad_ref     null      rad_ref              1.0     0     0
   lghtInGSI      lghtn       null      lghtn                1.0     0     0
   larcInGSI      larccld     null      larccld              1.0     0     0
::
 &SUPEROB_RADAR
   del_azimuth=5.,del_elev=.25,del_range=5000.,del_time=.5,elev_angle_max=5.,minnum=50,range_max=100000.,
   l2superob_only=.false.,
 /
 &LAG_DATA
 /
 &HYBRID_ENSEMBLE
   l_hyb_ens=$HYB_ENS,
   n_ens=$nens,
   uv_hyb_ens=.true.,
   beta_s0=0.25,
   s_ens_h=300,
   s_ens_v=5,
   generate_ens=.false.,
   regional_ensemble_option=1,
   fv3sar_bg_opt=${fv3lam_bg_opt},
   aniso_a_en=.false.,
   nlon_ens=0,
   nlat_ens=0,
   jcap_ens=574,
   l_ens_in_diff_time=.true.,
   jcap_ens_test=0,readin_beta=.false.,
   full_ensemble=.true.,pwgtflg=.true.,
   ensemble_path="",
 /
 &RAPIDREFRESH_CLDSURF
   i_gsdcldanal_type=${i_gsdcldanal_type},
   dfi_radar_latent_heat_time_period=20.0,
   l_use_hydroretrieval_all=.false.,
   metar_impact_radius=10.0,
   metar_impact_radius_lowCloud=4.0,
   l_gsd_terrain_match_surfTobs=.false.,
   l_sfcobserror_ramp_t=.false.,
   l_sfcobserror_ramp_q=.false.,
   l_PBL_pseudo_SurfobsT=.false.,
   l_PBL_pseudo_SurfobsQ=.false.,
   l_PBL_pseudo_SurfobsUV=.false.,
   pblH_ration=0.75,
   pps_press_incr=20.0,
   l_gsd_limit_ocean_q=.false.,
   l_pw_hgt_adjust=.false.,
   l_limit_pw_innov=.false.,
   max_innov_pct=0.1,
   l_cleanSnow_WarmTs=.false.,
   r_cleanSnow_WarmTs_threshold=5.0,
   l_conserve_thetaV=.false.,
   i_conserve_thetaV_iternum=3,
   l_cld_bld=.false.,
   cld_bld_hgt=1200.0,
   build_cloud_frac_p=0.50,
   clear_cloud_frac_p=0.1,
   iclean_hydro_withRef=1,
   iclean_hydro_withRef_allcol=0,
 /
 &CHEM
 /
 &SINGLEOB_TEST
 /
 &NST
 /

EOF

echo Shun debug1
anavinfo=$PARMfv3/anavinfo_fv3_65
###anavinfo=$fixgsi/anavinfo_fv3_60
berror=$fixgsi/$endianness/nam_glb_berror.f77.gcv
emiscoef_IRwater=$fixcrtm/Nalli.IRwater.EmisCoeff.bin
emiscoef_IRice=$fixcrtm/NPOESS.IRice.EmisCoeff.bin
emiscoef_IRland=$fixcrtm/NPOESS.IRland.EmisCoeff.bin
emiscoef_IRsnow=$fixcrtm/NPOESS.IRsnow.EmisCoeff.bin
emiscoef_VISice=$fixcrtm/NPOESS.VISice.EmisCoeff.bin
emiscoef_VISland=$fixcrtm/NPOESS.VISland.EmisCoeff.bin
emiscoef_VISsnow=$fixcrtm/NPOESS.VISsnow.EmisCoeff.bin
emiscoef_VISwater=$fixcrtm/NPOESS.VISwater.EmisCoeff.bin
emiscoef_MWwater=$fixcrtm/FASTEM6.MWwater.EmisCoeff.bin
aercoef=$fixcrtm/AerosolCoeff.bin
cldcoef=$fixcrtm/CloudCoeff.bin
#satinfo=$fixgsi/nam_regional_satinfo.txt
satinfo=$PARMfv3/fv3sar_satinfo.txt
scaninfo=$fixgsi/global_scaninfo.txt
satangl=$fixgsi/nam_global_satangbias.txt
atmsbeamdat=$fixgsi/atms_beamwidth.txt
pcpinfo=$fixgsi/nam_global_pcpinfo.txt
ozinfo=$fixgsi/nam_global_ozinfo.txt
errtable=$fixgsi/nam_errtable.r3dv
convinfo=$PARMfv3/rap_nam_regional_convinfo
mesonetuselist=$fixgsi/nam_mesonet_uselist.txt
stnuselist=$fixgsi/nam_mesonet_stnuselist.txt
locinfo=$PARMfv3/nam_hybens_d01_info

# Copy executable and fixed files to $DATA
##$ncp $gsiexec ./gsi.x
$ncp $anavinfo ./anavinfo
$ncp $berror   ./berror_stats
$ncp $emiscoef_IRwater ./Nalli.IRwater.EmisCoeff.bin
$ncp $emiscoef_IRice ./NPOESS.IRice.EmisCoeff.bin
$ncp $emiscoef_IRsnow ./NPOESS.IRsnow.EmisCoeff.bin
$ncp $emiscoef_IRland ./NPOESS.IRland.EmisCoeff.bin
$ncp $emiscoef_VISice ./NPOESS.VISice.EmisCoeff.bin
$ncp $emiscoef_VISland ./NPOESS.VISland.EmisCoeff.bin
$ncp $emiscoef_VISsnow ./NPOESS.VISsnow.EmisCoeff.bin
$ncp $emiscoef_VISwater ./NPOESS.VISwater.EmisCoeff.bin
$ncp $emiscoef_MWwater ./FASTEM6.MWwater.EmisCoeff.bin
$ncp $aercoef  ./AerosolCoeff.bin
$ncp $cldcoef  ./CloudCoeff.bin
$ncp $satangl  ./satbias_angle
$ncp $atmsbeamdat  ./atms_beamwidth.txt
$ncp $satinfo  ./satinfo
$ncp $scaninfo ./scaninfo
$ncp $pcpinfo  ./pcpinfo
$ncp $ozinfo   ./ozinfo
$ncp $convinfo ./convinfo
$ncp $errtable ./errtable
$ncp $mesonetuselist ./mesonetuselist
$ncp $stnuselist ./mesonet_stnuselist
$ncp $fixgsi/prepobs_prep.bufrtable ./prepobs_prep.bufrtable

# Copy CRTM coefficient files based on entries in satinfo file
for file in `awk '{if($1!~"!"){print $1}}' ./satinfo | sort | uniq` ;do
    $ncp $fixcrtm/${file}.SpcCoeff.bin ./
    $ncp $fixcrtm/${file}.TauCoeff.bin ./
done


###export nmmb_nems_obs=${COMINnam}/nam.${PDYrun}
export nmmb_nems_bias=${COMINbias}

# Copy observational data to $tmpdir
# Try RAP first
PDYa=$YYYYMMDD
cya=$HH
export nmmb_nems_obs=${COMINrap}/rap.${PDYa}
$ncp $nmmb_nems_obs/rap.t${cya}z.prepbufr.tm00  ./prepbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.prepbufr.acft_profiles.tm00 prepbufr_profl
$ncp $nmmb_nems_obs/rap.t${cya}z.satwnd.tm00.bufr_d ./satwndbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.1bhrs3.tm00.bufr_d ./hirs3bufr
$ncp $nmmb_nems_obs/rap.t${cya}z.1bhrs4.tm00.bufr_d ./hirs4bufr
$ncp $nmmb_nems_obs/rap.t${cya}z.mtiasi.tm00.bufr_d ./iasibufr
$ncp $nmmb_nems_obs/rap.t${cya}z.1bamua.tm00.bufr_d ./amsuabufr
$ncp $nmmb_nems_obs/rap.t${cya}z.1bamub.tm00.bufr_d ./amsubbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.1bmhs.tm00.bufr_d  ./mhsbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.goesnd.tm00.bufr_d ./gsnd1bufr
$ncp $nmmb_nems_obs/rap.t${cya}z.airsev.tm00.bufr_d ./airsbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.cris.tm00.bufr_d ./crisbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.atms.tm00.bufr_d ./atmsbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.sevcsr.tm00.bufr_d ./seviribufr
$ncp $nmmb_nems_obs/rap.t${cya}z.radwnd.tm00.bufr_d ./radarbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.nexrad.tm00.bufr_d ./l2rwbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.crisf4.tm00.bufr_d ./crisfsbufr
$ncp $nmmb_nems_obs/rap.t${cya}z.gsrcsr.tm00.bufr_d ./abibufr
#new ears data
$ncp $nmmb_nems_obs/rap.t${cya}z.esiasi.tm00.bufr_d ./iasibufrears
$ncp $nmmb_nems_obs/rap.t${cya}z.esamua.tm00.bufr_d ./amsuabufrears
$ncp $nmmb_nems_obs/rap.t${cya}z.esmhs.tm00.bufr_d  ./mhsbufrears
$ncp $nmmb_nems_obs/rap.t${cya}z.esatms.tm00.bufr_d ./atmsbufrears
$ncp $nmmb_nems_obs/rap.t${cya}z.escris.tm00.bufr_d ./crisfsbufrears
#new direct broadcast
$ncp $nmmb_nems_obs/rap.t${cya}z.crsfdb.tm00.bufr_d ./crisfsbufr_db
$ncp $nmmb_nems_obs/rap.t${cya}z.atmsdb.tm00.bufr_d ./atmsbufr_db
$ncp $nmmb_nems_obs/rap.t${cya}z.iasidb.tm00.bufr_d ./iasibufr_db

ls -1 prepbufr
err0=$?

#No RAP obs, get NAM data
if [ $err0 -ne 0 ] ; then
export nmmb_nems_obs=${COMINnam}/nam.${PDYrun}
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.prepbufr.${tmmark}  ./prepbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.prepbufr.acft_profiles.${tmmark} prepbufr_profl
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.satwnd.${tmmark}.bufr_d ./satwndbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.1bhrs3.${tmmark}.bufr_d ./hirs3bufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.1bhrs4.${tmmark}.bufr_d ./hirs4bufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.mtiasi.${tmmark}.bufr_d ./iasibufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.1bamua.${tmmark}.bufr_d ./amsuabufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.esamua.${tmmark}.bufr_d ./amsuabufrears
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.1bamub.${tmmark}.bufr_d ./amsubbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.1bmhs.${tmmark}.bufr_d  ./mhsbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.goesnd.${tmmark}.bufr_d ./gsnd1bufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.airsev.${tmmark}.bufr_d ./airsbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.cris.${tmmark}.bufr_d ./crisbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.atms.${tmmark}.bufr_d ./atmsbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.sevcsr.${tmmark}.bufr_d ./seviribufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.radwnd.${tmmark}.bufr_d ./radarbufr
$ncp $nmmb_nems_obs/nam.t${CYCrun}z.nexrad.${tmmark}.bufr_d ./l2rwbufr
fi

ls -l prepbufr
err10=$?

if [ $err10 -ne 0 ] ; then
  msg="NO PREPBUFR FILE; ABORT GSI JOB"
  err_exit $msg
fi

export GDAS_SATBIAS=NO

if [ $GDAS_SATBIAS = NO ] ; then

# Copy bias correction from prev cycle

GESROOT_HOLD=/gpfs/dell2/emc/modeling/noscrub/emc.campara/nwges/fv3lamdax.hold
cyctm06=06
$ncp $nmmb_nems_bias/fv3lam.${cyctm06}z.satbias.tm01 ./satbias_in
err1=$?
if [ $err1 -ne 0 ] ; then
  cp $GESROOT_HOLD/satbias_in ./satbias_in
fi

$ncp $nmmb_nems_bias/fv3lam.t${cyctm06}z.satbias_pc.tm01 ./satbias_pc
err2=$?
if [ $err2 -ne 0 ] ; then
  cp $GESROOT_HOLD/satbias_pc ./satbias_pc
fi

$ncp $nmmb_nems_bias/fv3lam.t${cyctm06}z.radstat.tm01    ./radstat.gdas
err3=$?
if [ $err3 -ne 0 ] ; then
  cp $GESROOT_HOLD/radstat.nam ./radstat.gdas
fi

else

cp $GESROOT_HOLD/gdas.satbias_out ./satbias_in
cp $GESROOT_HOLD/gdas.satbias_pc ./satbias_pc
cp $GESROOT_HOLD/gdas.radstat_out ./radstat.gdas

fi


#for new type satellite to build the bias coreection file
USE_RADSTAT=NO
USE_RADSTAT=${USE_RADSTAT:-"YES"}
USE_CFP=${USE_CFP:-"NO"}
##############################################################
# If requested, copy and de-tar guess radstat file
if [ $USE_RADSTAT = "YES" ]; then
   if [ $USE_CFP = "YES" ]; then
     rm $DATA/unzip.sh $DATA/mp_unzip.sh
     cat > $DATA/unzip.sh << EOFunzip
#!/bin/sh
   diag_file=\$1
   fname=\$(echo \$diag_file | cut -d'.' -f1)
   fdate=\$(echo \$diag_file | cut -d'.' -f2)
   #$UNCOMPRESS \$diag_file
   gunzip \$diag_file
   fnameges=\$(echo \$fname | sed 's/_ges//g')
   $NMV \$fname.\$fdate \$fnameges
EOFunzip
     chmod 755 $DATA/unzip.sh
   fi
   listdiag=$(tar xvf radstat.gdas | cut -d' ' -f2 | grep _ges)
   for type in $listdiag; do
      diag_file=$(echo $type | cut -d',' -f1)
      if [ $USE_CFP = "YES" ] ; then
         echo "$DATA/unzip.sh $diag_file" | tee -a $DATA/mp_unzip.sh
      else
         fname=$(echo $diag_file | cut -d'.' -f1)
         date=$(echo $diag_file | cut -d'.' -f2)
         #$UNCOMPRESS $diag_file
         gunzip $diag_file
         fnameges=$(echo $fname|sed 's/_ges//g')
         #$NMV $fname.$date $fnameges
         if [ $binary_diag = .true. ] ; then
            mv $fname.$date $fnameges
         elif [ $netcdf_diag = .true. ] ; then
            mv $fname.$date${DIAG_SUFFIX} $fnameges
         fi
      fi
  done
  if [ $USE_CFP = "YES" ] ; then
      chmod 755 $DATA/mp_unzip.sh
      ncmd=$(cat $DATA/mp_unzip.sh | wc -l)
      if [ $ncmd -gt 0 ]; then
         ncmd_max=$((ncmd < npe_node_max ? ncmd : npe_node_max))
         APRUNCFP_UNZIP=$(eval echo $APRUNCFP)
         $APRUNCFP_UNZIP $DATA/mp_unzip.sh
      fi
   fi
fi # if [ $USE_RADSTAT = "YES" ]


# Aircraft bias correction ; get from GDAS/GFS for tm06, cycle through FV3DA
# Try my GDAS dir first

COMINgdas=$COMINgdas/gdas.${PDYa}/${cya}/atmos
COMINgfs=$COMINgfs/gfs.${PDYa}/${cya}/atmos
$ncp $COMINgdas/gdas.t${cya}z.abias_air ./aircftbias_in
#if we don't find aircraft bias file from GDAS dirs, try GFS
if [ -s aircftbias_in ] ; then
  echo "found aircraft bias file from FV3GDAS"
else
  $ncp $COMINgfs/gfs.t${cya}z.abias_air ./aircftbias_in
fi


#if we still don't find aircraft bias file from GDAS/GFS dirs, get best available one
if [ -s aircftbias_in ] ; then
  echo "found aircraft bias file from FV3GDAS"
else
  cp $GESROOT_HOLD/gdas.airbias ./aircftbias_in
fi

cp $COMINrtma/rtma2p5.${PDYa}/rtma2p5.t${cya}z.w_rejectlist ./w_rejectlist
cp $COMINrtma/rtma2p5.${PDYa}/rtma2p5.t${cya}z.t_rejectlist ./t_rejectlist
cp $COMINrtma/rtma2p5.${PDYa}/rtma2p5.t${cya}z.p_rejectlist ./p_rejectlist
cp $COMINrtma/rtma2p5.${PDYa}/rtma2p5.t${cya}z.q_rejectlist ./q_rejectlist



#get coldstart 1st guess

#   This file contains vertical weights for defining hybrid volume hydrostatic pressure interfaces
dothis="false"
if [ $dothis = "true" ]; then
export fv3_case=$GUESSdir
cp $Fix_temp/fv_core.res.nc fv3_akbk
#   This file contains horizontal grid information
cp $fv3_case/user_coupler.res coupler.res
cp $Fix_temp/grid_spec.nc fv3_grid_spec
cp $fv3_case/sfc_data.tile7.nc fv3_sfcdata
cp $fv3_case/gfs_data.tile7.nc .
ln -sf gfs_data.tile7.nc fv3_dynvars
ln -sf gfs_data.tile7.nc fv3_tracer
fi


# Run gsi under Parallel Operating Environment (poe) on NCEP IBM
#export pgm=regional_gsi.x
#. prep_step

startmsg
${APRUN} $EXECfv3/regional_gsi.x < gsiparm.anl > stdout 2> stderr
export err=$?;err_chk


mv fort.201 fit_p1
mv fort.202 fit_w1
mv fort.203 fit_t1
mv fort.204 fit_q1
mv fort.205 fit_pw1
mv fort.207 fit_rad1
mv fort.209 fit_rw1

cat fit_p1 fit_w1 fit_t1 fit_q1 fit_pw1 fit_rad1 fit_rw1 > $COMOUT/fv3lam.t${CYCrun}z.fits.${tmmark}
cat fort.208 fort.210 fort.211 fort.212 fort.213 fort.220 > $COMOUT/fv3lam.t${CYCrun}z.fits2.${tmmark}
echo Shun debug2


if [ ${BKTYPE} -eq 1 ]; then  # cold start, put analysis back to current INPUT 
  cp_vrfy ${analworkdir}/fv3_dynvars                  ${bkpath}/gfs_data.tile7.halo0.nc
  cp_vrfy ${analworkdir}/fv3_sfcdata                  ${bkpath}/sfc_data.tile7.halo0.nc
else                          # cycling
  cp_vrfy ${analworkdir}/fv3_dynvars             ${bkpath}/fv_core.res.tile1.nc
  cp_vrfy ${analworkdir}/fv3_tracer              ${bkpath}/fv_tracer.res.tile1.nc
  cp_vrfy ${analworkdir}/fv3_sfcdata             ${bkpath}/sfc_data.nc
fi

exit

cp satbias_out $GESROOT_HOLD/satbias_in
cp satbias_out $COMOUT/fv3lam.t${CYCrun}z.satbias.${tmmark}
cp satbias_pc.out $GESROOT_HOLD/satbias_pc
cp satbias_pc.out $COMOUT/fv3lam.t${CYCrun}z.satbias_pc.${tmmark}

cp aircftbias_out $COMOUT/fv3lam.t${CYCrun}z.abias_air.${tmmark}
cp aircftbias_out $GESROOT_HOLD/gdas.airbias

cp gfs_data.tile7.nc $ANLdir/.
cp fv3_sfcdata $ANLdir/sfc_data.tile7.nc

##############################################################


#
#-----------------------------------------------------------------------
#
# Print message indicating successful completion of script.
#
#-----------------------------------------------------------------------
#
print_info_msg "
========================================================================
ANALYSIS GSI completed successfully!!!

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

