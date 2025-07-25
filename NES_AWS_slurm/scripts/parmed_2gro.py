import parmed as pmd
import sys

lig = sys.argv[1]
amber = pmd.load_file(f"{lig}.prmtop",xyz=f"{lig}.inpcrd")
amber.save('gromacs.top')
amber.save('gromacs.gro')
amber.save('gromacs.pdb', overwrite=True)
quit()
