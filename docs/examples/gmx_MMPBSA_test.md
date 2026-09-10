---
template: main.html
title: gmx_MMPBSA_test
---

## `gmx_MMPBSA_test` command-line
<div class="termy">
    ```bash
    $ gmx_MMPBSA_test -h
    usage: gmx_MMPBSA_test [-h] [-v] 
           [-t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,101,explicit_receptor_waters,gbnsr6,...} ...]]
           [-f FOLDER] [-r] [--examples-dir EXAMPLES_DIR] [--examples-source {clone,local}]
           [--skip-output-check] [-ng] [-n NUM_PROCESSORS] [-j NUM_CONCURRENT]
    
    This program is part of gmx_MMPBSA and will allow you to run various gmx_MMPBSA examples easily.
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    Test options:
      -t [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,101} [{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,101} ...]]
                            The level the test is going to be run at. Multiple systems and analysis can be run at the same 
                            time. Numeric ids, suite ids (`0`/`1`/`2`), legacy `101` (same as `0`), and named aliases
                            such as `explicit_receptor_waters` or `gbnsr6` are supported.
                                  Nr. of Sys  
                            * 0      23     All -- Run all examples (Can take a long time!!!)
                            * 1      12     Minimal -- Does a minimal test with a set of systems and analyzes
                                            that show that gmx_MMPBSA runs correctly. Only exclude 3drism, nmode
                                            protein-ligand MT because take a long time or are redundant
                            * 2       9     Fast -- Only the calculations that take a short time are run (Default)
                            [Systems]:
                                 Slow Frames
                            * 3    . | 10   Protein-Ligand (Single trajectory approximation)
                            * 4    . | 10   Protein-Protein
                            * 5    . | 10   Protein-DNA
                            * 6    x |  4   Protein-Membrane CHARMM-GUI PROA-UQ2
                            * 7    . | 10   Protein-Glycan
                            * 8    x | 10  Metalloprotein-ligand
                            * 9    . | 10  Multicomponent system (Comp_receptor)
                            * 10   x |  4   Protein-Ligand (CHARMM force field)
                            * 11     |      Legacy alias for test 6 (consolidated membrane example)
                            [Analysis]:
                                 Slow Frames
                            * 12   . | 10   Alanine Scanning
                            * 13   . | 10   Stability calculation
                            * 14   . | 10   Decomposition Analysis
                            * 15   . | 10  Interaction Entropy approximation
                            * 16   . | 10   Protein-Ligand (Multiple trajectory approximation)
                            * 17   x | 10  Entropy calculation using Normal Mode approximation
                            * 18   x |  4   Calculations using 3D-RISM approximation
                            * 19          C2 Entropy approximation
                            * 20          LPB Calculation
                            * 21          NLPB Calculation
                            * 22          Protein-Ligand_LPH (CHARMM force field)
                            * 23          QM/MMGBSA Calculation
                            * 24          GBNSR6 Calculation
                            * 25     |  5 AMBER input files
                            * 26     | 10 ST MM/PB(GB)SA with explicit receptor waters
      -f FOLDER, --folder FOLDER
                            Defines the folder to store all data
      -r, --reuse           Defines the existing test forlder will be reuse
      --examples-dir EXAMPLES_DIR
                            Use a local examples directory instead of cloning the repository
      --examples-source {clone,local}
                            Examples source mode. `local` requires `--examples-dir` or `GMXMMPBSA_TEST_EXAMPLES_DIR`
      --skip-output-check   Skip post-run verification of expected output files
      -ng, --nogui          No open gmx_MMPBSA_ana after all calculations finished
      -n NUM_PROCESSORS, --num_processors NUM_PROCESSORS
                            Defines the number of processor cores you want to use with MPI per calculation. If the number 
                            of frames is less than the number of cpus defined, the calculation will be performed with 
                            the number of processors = number of frames.
      -j NUM_CONCURRENT, --num_concurrent NUM_CONCURRENT
                            Defines the number of examples to run concurrently. Each example can use up to
                            --num_processors MPI ranks, so the total rank count can be -j * -n.
    
    
    gmx_MMPBSA is an effort to implement the GB/PB and others calculations in GROMACS. 
    Based on MMPBSA.py (version 16.0) and AmberTools20
    ```
</div>

## Running gmx_MMPBSA_test
gmx_MMPBSA_test can run a predefined set of examples or an individual example.
By default, examples run sequentially. To run multiple examples at the same time, use `-j/--num_concurrent`.
Each example can use up to `-n/--num_processors` MPI ranks.

By default, `gmx_MMPBSA_test` clones the GitHub repository to obtain the `examples/` folder. Developers working
from a local checkout can point directly at that folder with `--examples-dir` (or the
`GMXMMPBSA_TEST_EXAMPLES_DIR` environment variable) to avoid cloning and to test examples that match the
installed code. After each successful run, the tool verifies that expected output files exist unless
`--skip-output-check` is set.

Named selectors such as `-t explicit_receptor_waters` or `-t gbnsr6` are equivalent to their numeric ids.
`-t 101` is a legacy alias for the full `-t 0` suite.
The former CHARMM membrane selector `-t 11` is retained as a legacy alias for the consolidated membrane test 6;
the full suite runs that example only once.

!!! info "Sets in gmx_MMPBSA_test"

    === "Local examples (development)"

            TMP_EXAMPLES=$(mktemp -d)
            cp -a ./examples/. "$TMP_EXAMPLES/"
            gmx_MMPBSA_test -f /tmp/gmx_test --examples-dir "$TMP_EXAMPLES" -t 2 -ng

        This command makes `gmx_MMPBSA_test`:

        * Use the copied temporary examples tree instead of cloning GitHub
        * Run the `Fast` set (`-t 2`) against the checkout you are developing
        * Skip opening `gmx_MMPBSA_ana` at the end (`-ng`)

        `--examples-dir` is the directory whose individual example folders become
        worker directories; it is not redirected by `-f`. Copy the examples to a
        physical temporary tree before local-mode runs so generated files do not
        enter the source checkout. In clone mode, `-f` is the parent directory for
        `gmx_MMPBSA_test/`. That clone is replaced on a normal run when it already
        exists; use `-r/--reuse` to keep and reuse it, and use a new temporary
        parent when an isolated clone is required.

    === "Named selector"

            gmx_MMPBSA_test -f /home/user/Documents -t explicit_receptor_waters

        Equivalent to `-t 26` for the explicit receptor waters example.

    === "Fast"
        
            gmx_MMPBSA_test -f /home/user/Documents -n 10
        
        This command makes `gmx_MMPBSA_test`:
        
        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Use the `Fast` set of examples (`-t 2`, the default)
        * Run nine examples sequentially, using 10 CPUs for each example

    === "Parallel examples"

            gmx_MMPBSA_test -f /home/user/Documents -t 3 5 7 -n 4 -j 2

        This command makes `gmx_MMPBSA_test`:

        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Execute examples `3`, `5`, and `7`
        * Run up to 2 examples at the same time
        * Use up to 4 MPI ranks per example
    
    === "Minimal"
        
            gmx_MMPBSA_test -f /home/user/Documents -n 10 -t 1
        
        This command makes `gmx_MMPBSA_test`:
        
        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Use the `Minimal` set of examples (`-t 1`)
        * Run 12 examples sequentially, using 10 CPUs for each example
    
    === "All"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 0 -n 10
        
        This command makes `gmx_MMPBSA_test`:
        
        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Use the `All` set of examples
        * Run 23 examples sequentially, using 10 CPUs for each example
        
    === "Multiple selection"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 3 5 7
        
        This command makes `gmx_MMPBSA_test`:
        
        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Execute `3` [Protein-Ligand (Single Trajectory approach)], `5` [Protein-DNA], and `7` [Protein-Glycan]
        examples

    === "Single selection"
        
            gmx_MMPBSA_test -f /home/user/Documents -t 3
        
        This command makes `gmx_MMPBSA_test`:
        
        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Run example `3` [Protein-Ligand (Single Trajectory approach)]

    === "Explicit receptor waters"

            gmx_MMPBSA_test -f /home/user/Documents -t 26

        This command makes `gmx_MMPBSA_test`:

        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Run example `26` [ST MM/PB(GB)SA with explicit receptor waters]
        * Run from the `Explicit_receptor_waters` example folder using its local `mmpbsa.in` input

    === "AMBER input files"

            gmx_MMPBSA_test -f /home/user/Documents -t 25

        This command makes `gmx_MMPBSA_test`:

        * Download the gmx_MMPBSA repository to `/home/user/Documents`
        * Run example `25` [AMBER input files]

!!! warning "3D-RISM AmberTools runtime failures"
    Test `18` uses AmberTools 3D-RISM. If this test fails while the other examples pass and the log contains
    `Fortran runtime error: Missing comma between descriptors` from `amber_rism_interface.F90`, the failure is a
    known AmberTools/Fortran runtime compatibility issue. It has been reproduced with conda AmberTools builds linked
    against newer `libgfortran` runtimes. A known working workaround is `gmx_MMPBSA` 1.6.4 with Python 3.9 or 3.10,
    AmberTools 23, and `libgfortran5`/`libgcc-ng` 12.x, or a patched AmberTools build. See the
    [3D-RISM example](3D-RISM/README.md) for details.
