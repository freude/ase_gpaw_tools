import os
from ase.lattice.cubic import Diamond
from ase.visualize import view
from ase.io import read, write
from cif_tools import cif2fdf


# cif2fdf('/home/mk/Downloads/POSCAR_cif.cif')

# tc_slab = read('/home/mk/siesta_swarm/tc_molecule/siesta.STRUCT_OUT', index=None, format='struct_out', parallel=True)
# tc_slab1 = read('/home/mk/siesta_swarm/tc_molecule/siesta.STRUCT_NEXT_ITER', index=None, format='struct_out', parallel=True)
#
# tc_slab1.extend(tc_slab)

# read a tetracaene cell and make a slab out of it
tc_slab = read('tetracene.cif', index=None, format='cif', parallel=True)

# cif2fdf('/home/mk/Downloads/POSCAR2.cif')

params = tc_slab.get_cell_lengths_and_angles()
params[3] = 90.0
params[4] = 90.0
params[5] = 90.0
tc_slab.set_cell(params, scale_atoms=True)
tc_slab.wrap()
tc_slab = tc_slab.repeat((7, 1, 2))
tc_slab.translate((-2.488, 0, 0))
tc_slab.rotate(180, 'x', center=(0, 0, 0))
tc_slab.wrap()
# tc.write('tetracene1.cif')
# view(tc_slab)

# read from file a relaxed silicon slab
# silicon_slab = read('/home/mk/siesta_swarm/silicon_slab/siesta.STRUCT_OUT', format='struct_out', parallel=True)
silicon_slab = read('/home/mk/gpaw_swarm/gpaw_comp/relaxed_slab_cons3.gpw')
silicon_slab.set_constraint()
silicon_slab.wrap()
silicon_slab.set_pbc((0, 1, 1))
silicon_slab = silicon_slab.repeat((1, 1, 5))
# view(silicon_slab)

# tc_slab = read('supercell.cif', index=None, format='cif', parallel=True)
# tc_slab = read('/home/mk/tetracene/tetracene.cif', index=None, format='cif', parallel=True)

p_si = silicon_slab.get_positions()
p_tc = tc_slab.get_positions()
shift = 0.7
shift = p_si[:, 0].max() - p_tc[:, 0].min() + shift
tc_slab.translate((shift, 0, 0))

params_tc = tc_slab.get_cell_lengths_and_angles()
params_si = silicon_slab.get_cell_lengths_and_angles()
params_tc[1] = params_si[1]
params_tc[2] = params_si[2]
tc_slab.set_cell(params_tc, scale_atoms=True)


interface = silicon_slab.copy()
interface.extend(tc_slab)
length = -tc_slab.get_positions()[:, 0].max() + silicon_slab.get_positions()[:, 0].min()

interface.center(axis=0, vacuum=length-15)

view(interface)

# define a piece of silicon crystal 1x5x3
# atoms = Diamond(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#                 size=(1, 5, 3),
#                 symbol='Si',
#                 pbc=(1, 1, 0),
#                 latticeconstant=5.4)

# atoms = diamond100('Si', (5, 5, 20), a=5.4, vacuum=10, orthogonal=True)

# define a piece of silicon crystal 1x5x3

# view(atoms)
