#/usr/bin/env bash

DSET="Insight"
N_NETS=7
N_SUBS=20
MMS="10 20 40"

#TODO: change this to match the basepath of your current project/system
BP="/project/6008063/katie-gkiar/"
CMD="matlab -nodesktop -nosplash -r "
ALLSTAGE="addpath(\"${BP}data/\",  \"${BP}code/\", \"${BP}code/surfstat/\", \"${BP}code/BrainConnectivityToolbox/2017_01_15_BCT/\");"

mkdir -p ./logs/
cd ./logs/

SP=${BP}slurm_scripts
mkdir -p ${SP}/

# Optional TODO: comment out specific sections, to only run a susbet at a time.
# Compute residuals
for mm in $MMS
do
  STAGE1="vertexwiseSC(\"${DSET}\", 1, ${mm}, 1, \"false\", 1);"
  cat << TMP > ${SP}/exec_stage1_${mm}mm.sh
#!/bin/bash
#SBATCH --time 00:02:00
#SBATCH --mem 12G
#SBATCH --account rpp-aevans-ab

cd ${BP}
module load matlab
${CMD} '${ALLSTAGE} ${STAGE1} exit'
TMP
  chmod +x ${SP}/exec_stage1_${mm}.sh
  sbatch ${SP}/exec_stage1_${mm}.sh
done

# Compute network features
for net in `seq 1 ${N_NETS}`
do
  for sub in `seq 1 ${N_SUBS}`
  do
    for mm in $MMS
    do
      STAGE2="vertexwiseSC(\"${DSET}\", ${sub}, ${mm}, ${net}, \"false\", 2);"
      cat << TMP > ${SP}/exec_stage2_net${net}_${sub}_${mm}mm.sh
#!/usr/bin/env bash
#SBATCH --time 01:30:00
#SBATCH --mem 24G
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

# Perform statistical comparison
for net in `seq 1 ${N_NETS}`
do
  for mm in $MMS
  do
    STAGE3="vertexwiseSC(\"${DSET}\", 1, ${mm}, ${net}, \"false\", 3);"
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
