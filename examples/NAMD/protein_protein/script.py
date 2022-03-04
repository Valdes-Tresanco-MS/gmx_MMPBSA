# import ParmEd module
import parmed as pmd

# load psf file
psf = pmd.load_file('step3_input.psf')

#load coordinate file
psf.coordinates = pmd.load_file('step3_input.crd').coordinates

# strip ions and water
psf.strip(':POT, CLA, TIP3, LIT, SOD, RUB, CES, BAR')

# load Charmm Parameter Set. Make sure to include all the necessary force field files in this list
params = pmd.charmm.CharmmParameterSet('toppar/par_all36_carb.prm',
                                       'toppar/par_all36_cgenff.prm',
                                       'toppar/par_all36_lipid.prm',
                                       'toppar/par_all36m_prot.prm',
                                       'toppar/par_all36_na.prm',
                                       'toppar/par_interface.prm',
                                       'toppar/toppar_water_ions.str')
psf.load_parameters(params)

# save GROMACS topology file
psf.save('gromacs.top')
