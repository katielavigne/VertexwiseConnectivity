#/usr/bin/env bash

# Provide analysis stage as argument to the script
STAGE=$1

DSET="Insight"

N_NETS=7
N_SUBS=172
MMS="10 20 40"

#TODO : change this to match the basepath of your current project/system
BP="/project/6008063/katie-gkiar/VertexwiseConnectivity/"
CMD="matlab -nodesktop -nosplash -r "
ALLSTAGE="addpath(\"${BP}data/\",  \"${BP}code/\", \"${BP}code/surfstat/\", \"${BP}code/BrainConnectivityToolbox/2017_01_15_BCT/\");"

mkdir -p ${BP}data/logs/
cd ${BP}data/logs/

SP=${BP}code/slurm_scripts
mkdir -p ${SP}/

if [[ $STAGE == 1 ]]
then
  # Compute residuals
  mkdir -p ${BP}/data/resid/
  for mm in $MMS
  do
    STAGE1="vertexwiseSC(\"${DSET}\", 1, ${mm}, 1, 1);"
    cat << TMP > ${SP}/exec_stage1_${mm}mm.sh
#!/bin/bash
#SBATCH --time 00:05:00
#SBATCH --mem 16G
#SBATCH --account rpp-aevans-ab

cd ${BP}
module load matlab
${CMD} '${ALLSTAGE} ${STAGE1} exit'
TMP
    chmod +x ${SP}/exec_stage1_${mm}mm.sh
    sbatch ${SP}/exec_stage1_${mm}mm.sh
  done
fi

if [[ $STAGE == 2 ]]
then
  # Compute network features
  mkdir -p ${BP}/data/feat/
  for net in `seq 1 ${N_NETS}`
  do
    for sub in `seq 1 ${N_SUBS}`
    do
      for mm in $MMS
      do
        STAGE2="vertexwiseSC(\"${DSET}\", ${sub}, ${mm}, ${net}, 2);"
        cat << TMP > ${SP}/exec_stage2_net${net}_${sub}_${mm}mm.sh
#!/usr/bin/env bash
#SBATCH --time 01:30:00
#SBATCH --mem 30G
#SBATCH --account rpp-aevans-ab

cd ${BP}
module load matlab
${CMD} '${ALLSTAGE} ${STAGE2} exit'
TMP
        chmod +x ${SP}/exec_stage2_net${net}_${sub}_${mm}mm.sh
        sbatch ${SP}/exec_stage2_net${net}_${sub}_${mm}mm.sh
      done
    done
  done
fi

if [[ $STAGE == 3 ]]
then
  # Perform statistical comparison
  mkdir -p ${BP}/data/stats/
  for net in `seq 1 ${N_NETS}`
  do
    for mm in $MMS
    do
      STAGE3="vertexwiseSC(\"${DSET}\", 1, ${mm}, ${net}, 3);"
      cat << TMP > ${SP}/exec_stage3_net${net}_${mm}mm.sh
#!/usr/bin/env bash
#SBATCH --time 00:02:00
#SBATCH --mem 12G
#SBATCH --account rpp-aevans-ab

cd ${BP}
module load matlab
${CMD} '${ALLSTAGE} ${STAGE3} exit'
TMP
        chmod +x ${SP}/exec_stage3_net${net}_${mm}mm.sh
        # sbatch ${SP}/exec_stage3_net${net}_${mm}mm.sh
        echo $net $mm
        ${SP}/exec_stage3_net${net}_${mm}mm.sh
    done
  done
fi
