#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH -J tensionOfRope
#SBATCH -o %x-%j.log
#SBATCH --mail-user=lichen.simu@outlook.com
#SBATCH --mail-type=all
#SBATCH -p batch
#SBATCH --qos normal
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --mem 80GB
#SBATCH --time=12:00:00

# USELESS TO HAVE LESS THAN 10knodes per CPU!

# onIris=1

### Load latest available ABAQUS
module load cae/ABAQUS
module load lang/Python/3.7.4-GCCcore-8.3.0
module load lang/SciPy-bundle

### Configure environment variables, need to unset SLURM's Global Task ID for ABAQUS's PlatformMPI to work
unset SLURM_GTIDS

### Create ABAQUS environment file for current job, you can set/add your own options (Python syntax)
env_file=abaqus_v6.env

cat << EOF > ${env_file}        # append settings to abaqus env file
#verbose = 3
#ask_delete = OFF
mp_file_system = (SHARED, LOCAL)
EOF

node_list=$(scontrol show hostname ${SLURM_NODELIST} | sort -u)     # current node

mp_host_list="["
for host in ${node_list}; do
    mp_host_list="${mp_host_list}['$host', ${SLURM_CPUS_ON_NODE}],"     # $SLURM_CPUS_ON_NODE is number of CPUs on the node
done

mp_host_list=$(echo ${mp_host_list} | sed -e "s/,$/]/")     # mp_host_list = [['node_name', number of cpus], ['node_name', number of cpus], ...]

echo "mp_host_list=${mp_host_list}"  >> ${env_file}         # append to abaqus env file


initDir=$PWD                                                # record folder of this bash script

# -------------- above are environment settings on cluster -------------

jobName="tension"
initDir=$(pwd)
workDir="$(pwd)/tension"
jobName='tension'
cd ${workDir}
# launch computation
echo ${SLURM_NTASKS}    # print number of CPUs

# using Message Passing Interface (MPI) for parallel computing
abaqus job=${jobName} input=${jobName}.inp cpus=${SLURM_NTASKS} standard_parallel=all mp_mode=mpi interactive;

cd ${initDir}             # go back to folder of sh script

exit 0
