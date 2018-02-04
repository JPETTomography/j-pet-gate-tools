#!/bin/bash

NR_OF_SPLITS="50"

module load gate/7.2

# GC_DOT_GATE_DIR  : indicates the .Gate directory for splitted mac files
export GC_DOT_GATE_DIR=$PWD
# GC_GATE_EXE_DIR  : indicates the directory with the Gate executable
export GC_GATE_EXE_DIR=$GATE_DIR/bin/

if [[ $PWD/ != $GC_DOT_GATE_DIR ]]; then
	cd $GC_DOT_GATE_DIR
fi

rm -rf .Gate *.submit
# 0. Do not forget to add the /gate/cluster/setTimeSplitHalflife command in your main macro, following the isotope with the biggest activity"

# 1. See the PBS cluster status
qstat -a

# 2. Split the main macro in 2 using the PBS platform and default script
gjs -clusterplatform openPBS -openPBSscript ./array.pbs -numberofsplits $NR_OF_SPLITS main.mac

# 3. See the effect of the gjs program"
#cat main.submit
#cat .Gate/main/main.split

# 4. Launch the jobs using the qsub command
qsub -t 1-$NR_OF_SPLITS array.pbs

# 5. See your jobs in the cluster queue
qstat -t -u $USER
