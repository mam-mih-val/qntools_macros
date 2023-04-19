#!/bin/bash

#SBATCH -p fast
#SBATCH -t 1:00:00
#SBATCH -J QnTools
#SBATCH -o /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/log/%A_%a.log

list_dir=${1}
output_dir=${2}

id=$SLURM_ARRAY_TASK_ID

file_list=$( ls $list_dir | head -n $id | tail -n 1 )

mkdir -p $output_dir
cd $output_dir
mkdir $id
cd $id

source /mnt/pool/nica/7/mam2mih/soft/basov/root-6.24.06/install/bin/thisroot.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/pool/nica/7/mam2mih/soft/basov/QnTools/install/lib:/mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/

echo "/mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correct /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correct.cc $list_dir/$file_list"

# PLAIN
time /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correct /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correct.cc $list_dir/$file_list
# RECENTERING
time /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correct /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correct.cc $list_dir/$file_list
# TWIST AND RESCALING
time /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correct /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correct.cc $list_dir/$file_list

echo "/mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correlate /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correlate.cc correction_out.root"
time /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/build/correlate /mnt/pool/nica/7/mam2mih/soft/basov/qntools_macros/macro/proton_correlate.cc correction_out.root

echo "The End."