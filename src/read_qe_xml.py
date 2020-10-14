import numpy as np
import xml.etree.ElementTree as ET
from ase.io import Trajectory, cif
from ase.io import espresso
from ase.visualize import view
from ase import Atoms, Atom
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
from ase.units import create_units
import matplotlib.pyplot as plt


# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')

atoms = cif.read_cif('/home/mk/tetracene/tetracene.cif', None)

def xml2atoms(xml_elem):

    # symbols = [el.items()[0][1] for el in xml_elem.find('atomic_species').findall('species')]
    #
    # nat = int(xml_elem.find('atomic_structure').items()[0][1])
    # alat = float(xml_elem.find('atomic_structure').items()[1][1])

    atoms = []

    for item in xml_elem.find('atomic_structure').find('atomic_positions').findall('atom'):
        coords = np.array(item.text.split(), dtype=np.float) * units['Bohr']
        atom = Atom(symbol=item.items()[0][1],
                    position=coords)
        atoms.append(atom)

    a1 = np.array(xml_elem.find('atomic_structure').find('cell').find('a1').text.split(), dtype=np.float)
    a2 = np.array(xml_elem.find('atomic_structure').find('cell').find('a2').text.split(), dtype=np.float)
    a3 = np.array(xml_elem.find('atomic_structure').find('cell').find('a3').text.split(), dtype=np.float)

    cell = np.array([a1, a2, a3]) * units['Bohr']

    atoms = Atoms(atoms, cell=cell, pbc=True)

    # ------------------------------------------------------------------------------------------

    # Extract calculation results
    # Energy
    energy = 2*float(xml_elem.find('total_energy').find('etot').text) * units['Ry']

    # Forces
    forces = []
    if xml_elem.find('forces') is not None:
        for item in xml_elem.find('forces').text.strip().split('\n'):
            forces.append(np.fromstring(item.strip(), sep=' ', dtype=np.float))

    forces = np.array(forces) * units['Ry'] / units['Bohr']

    # Stress
    stress = None

    # Magmoms
    magmoms = None

    # Fermi level / highest occupied level
    efermi = None
    if xml_elem.find('band_structure') is not None:
        efermi = 2*float(xml_elem.find('band_structure').find('fermi_energy').text) * units['Ry']

    # # K-points
    ibzkpts = None
    # weights = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print them."
    #
    # for kpts_index in indexes[_PW_KPTS]:
    #     nkpts = int(pwo_lines[kpts_index].split()[4])
    #     kpts_index += 2
    #
    #     if pwo_lines[kpts_index].strip() == kpoints_warning:
    #         continue
    #
    #     # QE prints the k-points in units of 2*pi/alat
    #     # with alat defined as the length of the first
    #     # cell vector
    #     cell = structure.get_cell()
    #     alat = np.linalg.norm(cell[0])
    #     ibzkpts = []
    #     weights = []
    #     for i in range(nkpts):
    #         L = pwo_lines[kpts_index + i].split()
    #         weights.append(float(L[-1]))
    #         coord = np.array([L[-6], L[-5], L[-4].strip('),')],
    #                          dtype=float)
    #         coord *= 2 * np.pi / alat
    #         coord = kpoint_convert(cell, ckpts_kv=coord)
    #         ibzkpts.append(coord)
    #     ibzkpts = np.array(ibzkpts)
    #     weights = np.array(weights)
    #
    # # Bands
    # kpts = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print the bands."
    #
    # for bands_index in indexes[_PW_BANDS] + indexes[_PW_BANDSTRUCTURE]:
    #     if image_index < bands_index < next_index:
    #         bands_index += 2
    #
    #         if pwo_lines[bands_index].strip() == kpoints_warning:
    #             continue
    #
    #         assert ibzkpts is not None
    #         spin, bands, eigenvalues = 0, [], [[], []]
    #
    #         while True:
    #             L = pwo_lines[bands_index].replace('-', ' -').split()
    #             if len(L) == 0:
    #                 if len(bands) > 0:
    #                     eigenvalues[spin].append(bands)
    #                     bands = []
    #             elif L == ['occupation', 'numbers']:
    #                 # Skip the lines with the occupation numbers
    #                 bands_index += len(eigenvalues[spin][0]) // 8 + 1
    #             elif L[0] == 'k' and L[1].startswith('='):
    #                 pass
    #             elif 'SPIN' in L:
    #                 if 'DOWN' in L:
    #                     spin += 1
    #             else:
    #                 try:
    #                     bands.extend(map(float, L))
    #                 except ValueError:
    #                     break
    #             bands_index += 1
    #
    #         if spin == 1:
    #             assert len(eigenvalues[0]) == len(eigenvalues[1])
    #         assert len(eigenvalues[0]) == len(ibzkpts), \
    #             (np.shape(eigenvalues), len(ibzkpts))
    #
    #         kpts = []
    #         for s in range(spin + 1):
    #             for w, k, e in zip(weights, ibzkpts, eigenvalues[s]):
    #                 kpt = SinglePointKPoint(w, s, k, eps_n=e)
    #                 kpts.append(kpt)

    # ------------------------------------------------------------------------------------------

    calc = SinglePointDFTCalculator(atoms, energy=energy,
                                    forces=forces, stress=stress,
                                    magmoms=magmoms, efermi=efermi,
                                    ibzkpts=ibzkpts)
    # calc.kpts = kpts
    atoms.calc = calc

    input_parameters = {}

    input_parameters['ecut'] = None
    if xml_elem.find('basis_set') is not None:
        input_parameters['ecut'] = 2*float(xml_elem.find('basis_set').find('ecutwfc').text)

    input_parameters['input_dft'] = None
    if xml_elem.find('dft') is not None:
        input_parameters['input_dft'] = xml_elem.find('dft').find('functional').text.lower()

    # input_parameters['k_points'] = None
    # if xml_elem.find('band_structure') is not None:
    #     k_points = xml_elem.find('band_structure').find('starting_k_points').find('monkhorst_pack').items()
    #     k_points = [int(item[1]) for item in k_points]
    #     input_parameters['k_points'] = k_points

    return atoms, input_parameters


def read_qe_xml(fileobj, index=-1, results_required=True):
    """Reads Quantum ESPRESSO output files.

    The atomistic configurations as well as results (energy, force, stress,
    magnetic moments) of the calculation are read for all configurations
    within the output file.

    Will probably raise errors for broken or incomplete files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename
    index : slice
        The index of configurations to extract.
    results_required : bool
        If True, atomistic configurations that do not have any
        associated results will not be included. This prevents double
        printed configurations and incomplete calculations from being
        returned as the final configuration with no results data.

    Yields
    ------
    structure : Atoms
        The next structure from the index slice. The Atoms has a
        SinglePointCalculator attached with any results parsed from
        the file.


    """

    root = ET.parse(fileobj).getroot()
    output = root.find('output')
    steps = root.findall('step')

    atoms, input_parameters = xml2atoms(output)

    trajectory = None

    trajectory = Trajectory('t1.traj', 'a')
    atoms_list = []

    for step in steps:
        aaa, _ = xml2atoms(step)
        trajectory.write(aaa)
        atoms_list.append(aaa)

    trajectory.close()

    return atoms, input_parameters, atoms_list


if __name__ == '__main__':

    file_name = '/home/mk/data-file-schema.xml'
    file_name = '/home/mk/tetracene_opt2.xml'
    file_name1 = '/home/mk/tetracene.xml'
    # file_name1 = '/home/mk/tetracene.xml'
    # file_name2 = '/home/mk/tetracene1.xml'

    # atoms, atoms_list1 = read_qe_xml(file_name1)
    # atoms, atoms_list2 = read_qe_xml(file_name2)
    # atoms_list = atoms_list1+atoms_list2
    atoms, ecut, atoms_list = read_qe_xml(file_name)
    atoms1, ecut, atoms_list1 = read_qe_xml(file_name1)

    # print(len(atoms_list))
    # view(atoms_list)

    etot = [item.get_total_energy() for item in atoms_list]
    etot1 = [item.get_total_energy() for item in atoms_list1]

    etot = [np.linalg.norm(item.get_cell()[1]) for item in atoms_list]
    etot1 = [np.linalg.norm(item.get_cell()[1]) for item in atoms_list1]

    plt.plot(etot)
    plt.plot(etot1)
    plt.show()