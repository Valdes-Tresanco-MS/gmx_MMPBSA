"""Convert the bundled CHARMM PSF/CRD system to dry GROMACS inputs."""

import parmed as pmd

# Load the PSF and its matching coordinates.
psf = pmd.load_file('step3_input.psf')
psf.coordinates = pmd.load_file('step3_input.crd').coordinates

# Match the solvent and ion removal used for the XTC generated with cpptraj.
psf.strip(':POT, CLA, TIP3, LIT, SOD, RUB, CES, BAR')

# PDB chain IDs are one character wide, whereas the PSF uses the four-character
# segment IDs PROA and PROB. Map those segments explicitly and give each chain
# positive, sequential residue numbers; PROB begins at residue -3 in the source
# PSF, which is not portable through PDB-based tools.
chain_map = {'PROA': 'A', 'PROB': 'B'}
chain_residue_numbers = {chain: 0 for chain in chain_map.values()}
for residue in psf.residues:
    try:
        residue.chain = chain_map[residue.segid]
    except KeyError as exc:
        raise ValueError(f'Unexpected protein segment {residue.segid!r}') from exc
    chain_residue_numbers[residue.chain] += 1
    residue.number = chain_residue_numbers[residue.chain]

for number, atom in enumerate(psf.atoms, start=1):
    atom.number = number

# Preserve the explicit chain and per-chain residue numbering assigned above.
pmd.formats.PDBFile.write(psf, 'gromacs.pdb', renumber=False)

# Load every CHARMM parameter file required by the PSF.
params = pmd.charmm.CharmmParameterSet('toppar/par_all36_carb.prm',
                                       'toppar/par_all36_cgenff.prm',
                                       'toppar/par_all36_lipid.prm',
                                       'toppar/par_all36m_prot.prm',
                                       'toppar/par_all36_na.prm',
                                       'toppar/par_interface.prm',
                                       'toppar/toppar_water_ions.str')
psf.load_parameters(params)

# Write the GROMACS topology consumed by gmx_MMPBSA.
psf.save('gromacs.top', overwrite=True)
