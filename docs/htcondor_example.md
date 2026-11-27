# Running gmx_MMPBSA on HTCondor Clusters

This short guide explains how to run `gmx_MMPBSA` on an HTCondor-managed cluster (e.g., NMRbox), which is a common setup for large-scale MD simulations.

## Sample Submit File (`.sub`)

```condor
universe = vanilla
executable = /usr/software/GROMACS/gpu/bin/gmx_mpi

# Adjust the arguments to your specific GROMACS/MMPBSA setup
arguments = mdrun -s md.tpr -deffnm md -ntmpi 1 -ntomp 4 -nb gpu -pme gpu -bonded gpu -update gpu

should_transfer_files = NO
output = logs/md.out
error = logs/md.err
log = logs/md.log

requirements = (TotalGPUs >= 1) && (OpSys == "LINUX")
request_cpus = 4
request_gpus = 1

queue
Important Notes for OpenCL vs CUDA:

If your GROMACS build uses CUDA, you can safely use -nb gpu -pme gpu -bonded gpu -update gpu.

If your GROMACS build uses OpenCL, the -update gpu and -bonded gpu flags are not supported and will cause the job to crash with Inconsistency in user input.
