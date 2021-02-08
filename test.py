import os

current_dir = os.getcwd()

dirs = ['Entropy_calculations/Interaction_Entropy',
        'Protein_ligand/ST',
        'Protein_protein'
        'Stability',
        'Protein_DNA',
        'Protein_membrane',
        'Protein_glycan',
        'Alanine_scanning',
        'Decomposition_analysis',
        'Metalloprotein_peptide',
        'Protein_DNA_RNA_Ion_ligand',
        'Protein_ligand_CHARMMff']

try:
    os.chdir(f'{current_dir}/test_files/Entropy_calculations/Interaction_Entropy')
    print('Testing Interaction Entropy calculations...')
    os.system('python ../../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_ligand/ST')
    print('Testing Protein_ligand/ST BFE calculations...')
    os.system('python ../../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm '
              'ligand.mol2')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_protein')
    print('Testing Protein_protein BFE calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Stability')
    print('Testing Stability calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -s -cs com.tpr -ci index.ndx -cg 19 20 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_DNA')
    print('Testing Protein_DNA BFE calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_membrane')
    print('Testing Membrane Protein PB calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.pdb -ci index.ndx -cg 1 13 -ct com_traj.pdb -lm ligand.mol2')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_glycan')
    print('Testing Protein_glycan BFE calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 26 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Alanine_scanning')
    print('Testing Alanine_scanning calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 12 -ct com_traj.xtc')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Decomposition_analysis')
    print('Testing Decomposition_analysis calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -lm ligand.mol2')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Metalloprotein_peptide')
    print('Testing Metalloprotein_peptide BFE calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 20 21 -ct com_traj.pdb')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_DNA_RNA_Ion_ligand')
    print('Testing Protein_DNA_RNA_Ion_ligand BFE calculations...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 33 14 -ct com_traj.xtc -lm ligand.mol2')
    print('===========================================================================================================')
except:
    pass

try:
    os.chdir(f'{current_dir}/test_files/Protein_ligand_CHARMMff')
    print('Testing Protein_ligand_CHARMMff BFE calculations with CHARMMff files...')
    os.system('python ../../exec.py -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc -cp topol.top')
    print('===========================================================================================================')
except:
    pass


