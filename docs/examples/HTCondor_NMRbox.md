# Running gmx_MMPBSA on HTCondor (NMRbox)

This example shows how to submit a `gmx_MMPBSA` calculation to the NMRbox HTCondor pool.

## Submit file: `gmx_mmpbsa.sub`

```bash
executable = /usr/local/gromacs/bin/gmx_MMPBSA
arguments  = -O -i mmpbsa.in -cs complex.tpr -ci index.ndx -cg 1 13 -ct complex.xtc -cp topol.top -o FINAL_RESULTS.dat -eo FINAL_RESULTS.csv

log    = logs/mmpbsa_$(PROCID).log
output = logs/mmpbsa_$(PROCID).out
error  = logs/mmpbsa_$(PROCID).err

request_cpus   = 1
request_memory = 4GB
request_disk   = 5GB

+Production = True
should_transfer_files = NO
transfer_executable   = FALSE

queue 4
Notes
request_memory = 4GB is conservative. gmx_MMPBSA can be memory-hungry for per-residue decomposition; monitor the first job and adjust.

+Production = True ensures the job runs on a node with the NMR software stack.

should_transfer_files = NO is safe on NMRbox (shared filesystem).

queue 4 uses $(PROCID) to generate four independent jobs.

Verifying results
bash
condor_q
condor_q -nobatch
tail -f logs/mmpbsa_0.out
Concatenate results:

bash
cat FINAL_RESULTS_*.csv > FINAL_RESULTS_combined.csv
Known limitations
gmx_MMPBSA does not checkpoint per-frame. Split long trajectories into segments.

-eo and -do flags are independent; can be run as separate job arrays.

References
gmx_MMPBSA docs: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/

NMRbox HTCondor tutorial: https://github.com/NMRbox/htcondor-tutorial

MM-PBSA method: Kollman et al., Acc. Chem. Res. 33, 889-897 (2000)
