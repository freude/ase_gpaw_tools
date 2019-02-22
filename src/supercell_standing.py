from ase.visualize import view
from ase.io import read, write


def supercell_standing(molecules='/home/mk/gpaw_swarm/gpaw_comp1/relaxed_mol.gpw',
                       silicon='/home/mk/gpaw_swarm/gpaw_comp/relaxed_slab_cons3.gpw'):

    # read a tetracaene cell and make a slab out of it
    tc_slab = read(molecules)

    params = tc_slab.get_cell_lengths_and_angles()
    params[3] = 90.0
    params[4] = 90.0
    params[5] = 90.0
    tc_slab.set_cell(params, scale_atoms=True)
    tc_slab.wrap()
    tc_slab = tc_slab.repeat((2, 1, 2))
    tc_slab.translate((-2.488, 0, 0))
    tc_slab.translate((0, 0, 7))
    tc_slab.wrap()

    # read from file a relaxed silicon slab
    silicon_slab = read(silicon)
    silicon_slab.set_constraint()
    cell = silicon_slab.get_cell()
    cell[2, 2] = cell[0, 0]
    cell[0, 0] = cell[1, 1]
    silicon_slab.set_cell(cell)
    silicon_slab.rotate(90, 'y', center=(0, 0, 0))
    silicon_slab.wrap()
    silicon_slab.set_pbc((1, 1, 0))
    silicon_slab = silicon_slab.repeat((3, 1, 1))
    # view(silicon_slab)

    p_si = silicon_slab.get_positions()
    p_tc = tc_slab.get_positions()
    shift = 1.0
    shift = p_si[:, 2].max() - p_tc[:, 2].min() + shift
    tc_slab.translate((0, 0, shift))

    params_tc = tc_slab.get_cell_lengths_and_angles()
    params_si = silicon_slab.get_cell_lengths_and_angles()
    params_tc[1] = params_si[1]
    params_tc[0] = params_si[0]
    tc_slab.set_cell(params_tc, scale_atoms=True)

    interface = silicon_slab.copy()
    interface.extend(tc_slab)
    length = -tc_slab.get_positions()[:, 2].max() + silicon_slab.get_positions()[:, 2].min()

    cell = interface.get_cell()
    cell[2, 2] = length
    interface.set_cell(cell)

    interface.center(axis=2, vacuum=length-10)

    # del interface[interface.positions[:, 2] > 33]

    # view(interface, viewer='vmd')
    view(interface)
    print('hi')


if __name__ == '__main__':

    supercell_standing()
