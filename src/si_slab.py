from __future__ import print_function, division
import numpy as np
import scipy.spatial
from ase.io import read, write
from ase import Atom


def distance_to_surface(coords, surface):
    """
    Computes distance of from the point with coordinates `coords`
    to the surface  defined by the variable `surface`

    :param coords:
    :param surface:
    :return:
    """

    return np.inner(coords - surface, surface) / np.linalg.norm(surface)


def belong_to_surface(coords, surface):
    """
    Check whether the point with coordinates `coords`
    belongs to the surface defined by the variable `surface`

    :param coords:
    :param surface:
    :return:
    """

    dist = distance_to_surface(coords, surface)
    flag = False

    if np.abs(dist) < 0.01:
        flag = True

    return flag


class KdTree(object):
    """
    Wrapper class for `scipy.spatial.cKDTree` class.
    """

    def __init__(self, coords, nn_dist=2.39):
        self._kd_tree = scipy.spatial.cKDTree(coords, leafsize=10)
        self._nn_distance = nn_dist

    def get_neighbours(self, coord):

        ans = self._kd_tree.query(coord,
                                  k=5,
                                  distance_upper_bound=self._nn_distance)

        ans1 = [ans[1][0]]

        for item in zip(ans[0], ans[1]):
            if self._nn_distance * 0.01 < item[0] < self._nn_distance:
                ans1.append(item[1])

        return ans1[1:]


def R_x(degrees):
    """
    Transformation matrix - rotation about the axis x.

    :param degrees:
    :return:
    """

    theta = np.pi / 180 * degrees

    R = np.matrix([[1.0, 0.0, 0],
                   [0.0, np.cos(theta), -np.sin(theta)],
                   [0.0, np.sin(theta), np.cos(theta)]])

    return R


def R_z(degrees):
    """
    Transformation matrix - rotation about the axis x.

    :param degrees:
    :return:
    """

    theta = np.pi / 180 * degrees

    R = np.matrix([[np.cos(theta), -np.sin(theta), 0.0],
                   [np.sin(theta), np.cos(theta), 0.0],
                   [0.0, 0.0, 1.0]])

    return R


def R_matrix(axis, degrees):
    """
    Transformation matrix - rotation about an arbitraty axis
    defined by the variable `axis`.

    :param axis:
    :param degrees:
    :return:
    """

    theta = np.pi / 180 * degrees

    R = np.matrix(np.zeros((3, 3)))

    R[0, 0] = (1 - np.cos(theta)) * axis[0] ** 2 + np.cos(theta)
    R[0, 1] = (1 - np.cos(theta)) * axis[0] * axis[1] - axis[2] * np.sin(theta)
    R[0, 2] = (1 - np.cos(theta)) * axis[0] * axis[2] + axis[1] * np.sin(theta)
    R[1, 0] = (1 - np.cos(theta)) * axis[1] * axis[0] + axis[2] * np.sin(theta)
    R[1, 1] = (1 - np.cos(theta)) * axis[1] ** 2 + np.cos(theta)
    R[1, 2] = (1 - np.cos(theta)) * axis[1] * axis[2] - axis[0] * np.sin(theta)
    R[2, 0] = (1 - np.cos(theta)) * axis[2] * axis[0] - axis[1] * np.sin(theta)
    R[2, 1] = (1 - np.cos(theta)) * axis[2] * axis[1] + axis[0] * np.sin(theta)
    R[2, 2] = (1 - np.cos(theta)) * axis[2] ** 2 + np.cos(theta)

    return R


def add_atoms(coord, coords, length=0.7):
    """
    Loop over all atomic coordinates in the list and
    add to the list an atom with specific requirements
    to its position.

    :param coord:      coordinates of an atom
    :param coords:     coordinates of its neighbours
    :param length:     bond length of passivating atom
    :return:
    """

    a = []

    for item in coords:

        if len(coords) == 2:
            a1 = np.matrix(length * (coord - item))
            a1 = R_z(90) * a1.T
            a.append(coord + np.squeeze(np.array(a1)))

        if len(coords) == 1:
            a0 = np.matrix(length * (coord - item))
            a1 = R_z(90) * a0.T
            a1 = coord + np.squeeze(np.array(a1))
            a.append(a1)

            a1 = R_z(180) * R_z(90) * a0.T
            a1 = coord + np.squeeze(np.array(a1))
            a.append(a1)

    return a


def make_canted(coord, coords, atoms):
    """

    :param coord:      coordinates of an atom
    :param coords:     coordinates of its neighbours
    :param atoms:      coordinates of passivating atoms
    :return:
    """

    atoms_out = []
    origin = np.array([coord[0], coord[1], coords[0][2]])
    axis = np.abs(coords[0] - coord)
    axis[2] = 0.0
    axis /= np.linalg.norm(axis)
    axis = np.squeeze(np.array(R_z(0) * np.matrix(axis).T))

    coord_out = origin + np.squeeze(np.array(R_matrix(axis, 10.0) * np.matrix(coord - origin).T))

    for atom in atoms:
        a1 = origin + np.squeeze(np.array(R_matrix(axis, 10.0) * np.matrix(atom - origin).T))
        atoms_out.append(a1)

    return coord_out, atoms_out


def fold_into(a1, cell):
    for j in range(3):
        if a1[j] < 0:
            a1[j] += cell[j]
        if a1[j] > cell[j]:
            a1[j] -= cell[j]

    return a1


def passivate_surface_ase(atoms, elem, plane):
    atom_list = []
    cell = np.diag(atoms.get_cell())
    coords = atoms.get_scaled_positions()
    coords_cart = atoms.get_positions()

    for j in range(3):
        plane[j] = plane[j](coords[:, j])

    coords_cart = np.array(coords_cart)

    kd_tree = KdTree(coords_cart)

    for j, coord in enumerate(coords_cart):
        if belong_to_surface(coords[j], plane):
            ans = kd_tree.get_neighbours(coord)

            atoms1 = add_atoms(coord, coords_cart[ans])
            # coord, atoms1 = make_canted(coord, coords_cart[ans], atoms1)

            for atom in atoms1:
                a = atom.tolist()
                atoms.append(Atom(elem, a))

    return atoms


def passivate_surface_ase_2x1(atoms, elem, plane):
    atom_list = []
    cell = np.diag(atoms.get_cell())
    coords = atoms.get_scaled_positions()
    coords_cart = atoms.get_positions()

    for j in range(3):
        plane[j] = plane[j](coords[:, j])

    coords_cart = np.array(coords_cart)

    kd_tree = KdTree(coords_cart)

    for j, coord in enumerate(coords_cart):
        if belong_to_surface(coords[j], plane):
            ans = kd_tree.get_neighbours(coord)

            atoms1 = add_atoms(coord, coords_cart[ans])
            coord, atoms1 = make_canted(coord, coords_cart[ans], atoms1)

            for atom in atoms1:
                a = atom.tolist()
                atoms.append(Atom(elem, a))

    return atoms


def make_slab(a_si, width, axis=2, vacuum=10):
    from ase.build import bulk
    from ase.visualize import view

    si = bulk('Si', 'diamond', a_si, cubic=True)
    si = si.repeat((5, 5, width))
    si.center(axis=axis, vacuum=vacuum)

    def zero(x):
        return np.min(x * 0)

    si = passivate_surface_ase(si, 'H', [zero, zero, np.max])
    passivate_surface_ase(si, 'H', [zero, zero, np.min])

    view(si)

    import pickle

    outfile = open('si_slab_' + str(width) + '.pkl', 'wb')
    pickle.dump(si, outfile)
    outfile.close()


def make_slab_111(a_si, width, axis=2, vacuum=10):
    from ase.build import bulk
    from ase.visualize import view
    from ase.build import diamond111

    si = diamond111('Si', size=(2, 2, width), a=a_si)
    si.center(axis=axis, vacuum=vacuum)

    cell = si.get_cell()
    cell[1, 0] = 0.0
    si.set_cell(cell)
    si.wrap()

    bottom_surface = np.min(si.positions[:, 2])
    top_surface = np.max(si.positions[:, 2])

    for atom in si:

        if np.abs(atom.position[2] - bottom_surface) < 0.01:
            atom.symbol = 'H'
            atom.position[2] += 0.5

        if np.abs(atom.position[2] - top_surface) < 0.01:
            atom.symbol = 'H'
            atom.position[2] -= 0.5

    # def zero(x):
    #     return np.min(x * 0)
    #
    # si = passivate_surface_ase(si, 'H', [zero, zero, np.max])
    # passivate_surface_ase(si, 'H', [zero, zero, np.min])

    view(si)

    import pickle

    outfile = open('si_slab111_' + str(width) + '.pkl', 'wb')
    pickle.dump(si, outfile)
    outfile.close()


def extend_by_a_slice(atoms, cut, lattice_vector, axis=2):
    first_cut = cut
    second_cut = cut + lattice_vector

    indices = []

    for atom in atoms:
        if first_cut < atom.position[axis] < second_cut:
            indices.append(atom.index)

    atom_buffer = atoms[indices]

    for atom in atoms:
        if atom.position[axis] > first_cut:
            atom.position[axis] += lattice_vector

    atoms.extend(atom_buffer)
    return atoms


def make_slab_100_1x2(num_periods, vacuum=10):
    atoms = read('/home/mk/ase_gpaw_tools/src/si_slab_100_1x2.xyz', format='xyz')

    num_periods = num_periods - 3
    while num_periods > 0:
        atoms = extend_by_a_slice(atoms, 0.5 * (19.57 - 18.21) + 18.21, 5.473)
        num_periods -= 1
    atoms_extent = np.max(atoms.positions[:, 2]) - np.min(atoms.positions[:, 2])
    cell_z_size = atoms_extent + vacuum
    atoms.set_cell([7.9, 7.9 / 2, cell_z_size], scale_atoms=False)
    atoms.translate((-np.min(atoms.positions[:, 0]), 0, 0.0))
    atoms.translate((0, -np.min(atoms.positions[:, 1]), 0.0))
    atoms.rotate(45, 'z', center=(0, 0, 0))
    atoms.translate((-np.min(atoms.positions[:, 0]) - 0.648, 0, 0.0))
    atoms.translate((0, -np.min(atoms.positions[:, 1]), 0.0))
    atoms.translate((0, 0, -np.min(atoms.positions[:, 2]) - 0.5 * atoms_extent + cell_z_size * 0.5))
    atoms.set_pbc((1, 1, 0))
    del atoms[[atom.index for atom in atoms if atom.position[1] > 2]]

    return atoms


if __name__ == '__main__':

    from ase.visualize import view
    import pickle

    for i in range(4, 20):
        atoms = make_slab_100_1x2(i, vacuum=37)

        write('si_slab_100_1x2_' + str(i) + '.struct', atoms, format='struct')
        si_slab = read('si_slab_100_1x2_' + str(i) + '.struct', format='struct')

        # with open('si_slab_100_1x2_' + str(i) + '.pkl', 'wb') as outfile:
        #     pickle.dump(atoms, outfile)
        #
        # with open('si_slab_100_1x2_' + str(i) + '.pkl', 'rb') as infile:
        #     si_slab = pickle.load(infile)

        view(si_slab)
