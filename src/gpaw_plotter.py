import numpy as np
import matplotlib.pyplot as plt
from ase.units import Bohr
from ase.io import read
from gpaw.utilities.ps2ae import PS2AE
from gpaw import GPAW
from invdisttree import Invdisttree


def main():
    # from gpaw import Mixer, restart
    # from gpaw.utilities import h2gpts
    # xc = {'name': 'C09-vdW', 'backend': 'libvdwxc'}
    # interface = read('/home/mk/gpaw_swarm/standing22.gpw')
    #
    # interface.calc = GPAW(mode='lcao',
    # gpts = h2gpts(0.15, interface.get_cell(), idiv=16),
    # basis = 'dzp',
    # kpts = (4, 2, 1),
    # xc = xc,
    # parallel = dict(kpt=4, augment_grids=True),
    # mixer = Mixer(beta=0.001, nmaxold=5, weight=50.0))
    #
    # for j in range(len(interface.calc.setups)):
    #     if interface.positions[j, 2] > 34.5:
    #         interface.calc.setups[j].type = 'ghost'

    # interface = read('/home/mk/gpaw_swarm/gpaw_comp/relaxed_slab_cons3.gpw')
    # calc = GPAW('/home/mk/gpaw_swarm/standing_geometry111/standing111/standing31.gpw', txt=None)
    calc = GPAW('/home/mk/gpaw_swarm/standing22.gpw', txt=None)
    # calc1 = GPAW('/home/mk/gpaw_swarm/si_slab_100_1x2_4.gpw', txt=None)
    # calc = GPAW('/home/mk/gpaw_swarm/si_slab_100_1x2_6.gpw', txt=None)
    calc1 = GPAW('/home/mk/gpaw_swarm/si_slab_100_1x2_8.gpw', txt=None)

    nt = calc.get_pseudo_density()
    nt1 = calc1.get_pseudo_density()

    nt = np.mean(nt, axis=(0, 1))
    nt1 = np.mean(nt1, axis=(0, 1))

    z = calc.density.gd.coords(2)
    z1 = calc1.density.gd.coords(2)

    from scipy.interpolate import interp1d

    my_interpolating_function = interp1d(z, nt, kind='cubic', fill_value="extrapolate")
    z = np.linspace(np.min(z), np.max(z), 5000)
    nt = my_interpolating_function(z)

    my_interpolating_function = interp1d(z1, nt1, kind='cubic', fill_value="extrapolate")
    z1 = np.linspace(np.min(z1), np.max(z1), 5000)
    nt1 = my_interpolating_function(z1)

    from scipy.signal import find_peaks
    peaks, _ = find_peaks(nt, height=0.1 * np.max(nt))
    peaks1, _ = find_peaks(nt1, height=0.1 * np.max(nt1))

    z_peak = z[peaks[0]]
    z_peak1 = z1[peaks1[0]]
    shift = z_peak - z_peak1

    my_interpolating_function = interp1d(z1+shift, nt1, kind='cubic', fill_value="extrapolate")
    nt2 = my_interpolating_function(z)


    # t = PS2AE(calc, h=0.03)


    # Interpolated PS and AE potentials:
    # ps = t.get_electrostatic_potential(ae=False)
    ps = calc.get_electrostatic_potential()
    ps1 = calc1.get_electrostatic_potential()
    # ae = t.get_electrostatic_potential()
    # i = ps.shape[0] // 2
    # x = t.gd.coords(2) * Bohr

    # plt.plot(x, ps[i, i], '--', label=r'$\tilde v$')
    # plt.plot(x, np.mean(ae, axis=(0, 1)), '-', label=r'$v$')

    en = np.linspace(-20, 0)
    num_e = 1000
    coords = calc.atoms.get_positions()
    ind = np.argsort(coords[:, 2])

    num_atoms = len(coords)

    num_points = 80

    z = np.linspace(np.min(coords[:, 2])-10, np.max(coords[:, 2])+10, num_points)
    step = (np.max(coords[:, 2]) - np.min(coords[:, 2])) / num_points
    x = np.linspace(np.min(coords[:, 0]), np.max(coords[:, 0]), (np.max(coords[:, 0]) - np.min(coords[:, 0])) / step)
    y = np.linspace(np.min(coords[:, 1]), np.max(coords[:, 1]), (np.max(coords[:, 1]) - np.min(coords[:, 1])) / step)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # data = np.zeros((num_atoms, num_e))
    data1 = np.zeros((num_atoms, num_e))

    for jj, j in enumerate(ind):
        # energy, pdos = calc.get_orbital_ldos(a=int(j), npts=num_e)
        # data[j, :] = pdos
        energy, pdos = calc.get_lcao_dos(atom_indices=int(j), npts=num_e)
        data1[jj, :] = pdos

    coords = coords[ind, :]
    en_ind = np.arange(len(energy))[500:]
    data = np.zeros((len(x), len(y), len(z), len(en_ind)))

    for jj, j in enumerate(en_ind):
        print(j)
        Fq = Invdisttree(coords, data1[:, j])
        aaa = Fq(np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T, nnear=1, eps=0, p=8).reshape(X.shape)
        data[:, :, :, jj] = aaa

    ef = calc.get_fermi_level()
    plt.contourf(np.arange(num_atoms), energy, np.mean(data, axis=(0, 1)).T, 400)
    plt.plot([0, num_atoms], [ef, ef], 'r-')
    plt.show()

    plt.plot(np.mean(ps, axis=(0, 1)), '-', label=r'$v$')
    plt.show()


def main1():

    ps02 = np.load('/home/mk/nci/ps02.npy')
    nt02 = np.load('/home/mk/nci/nt02.npy')
    calc = GPAW('/home/mk/nci/standing25extended0.gpw', txt=None)
    nt00 = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    ps00 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    calc = GPAW('/home/mk/nci/standing25extended1.gpw', txt=None)
    nt01 = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    ps01 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    calc = GPAW('/home/mk/nci/standing25reduced1.gpw', txt=None)
    nt_1 = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    ps_1 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))

    m1 = np.argmax(ps00[200:-200])
    # m2 = np.argmax(ps[200:-200])
    m2 = np.argmax(ps01[200:-200])
    m3 = np.argmax(ps02[200:-200])

    ps00 = np.pad(ps00, m3-m1, mode='edge')
    # ps = np.pad(ps, m4 - m2, mode='edge')
    ps01 = np.pad(ps01, m3 - m2, mode='edge')

    plt.plot(ps00)
    # plt.plot(ps)
    plt.plot(ps01)
    plt.plot(ps02)
    plt.show()

    # plt.plot(np.mean((nt - nt1 - nt2), axis=(0, 1)), 'k')
    # plt.plot(np.mean((nt - nt3 - nt4), axis=(0, 1)), 'r')
    # plt.fill_between(np.arange(nt.shape[2]), np.mean((nt - nt3 - nt4), axis=(0, 1)),
    #                  np.mean((nt - nt2 - nt1), axis=(0, 1)), color='grey', alpha=0.5)
    # plt.plot(np.mean(nt / 1000, axis=(0, 1)), 'g')
    # plt.show()

    # ps = calc.get_electrostatic_potential()
    ps = np.mean(ps, axis=(0, 1))
    from scipy.interpolate import interp1d
    my_interpolating_function = interp1d(calc.hamiltonian.gd.coords(2), ps[::2], kind='cubic', fill_value="extrapolate")

    en = np.linspace(-20, 0)
    num_e = 1000
    coords = calc.atoms.get_positions()
    ind = np.argsort(coords[:, 2])

    num_atoms = len(coords)

    num_points = 80

    z = np.linspace(np.min(coords[:, 2]) - 10, np.max(coords[:, 2]) + 10, num_points)
    step = (np.max(coords[:, 2]) - np.min(coords[:, 2])) / num_points
    x = np.linspace(np.min(coords[:, 0]), np.max(coords[:, 0]), (np.max(coords[:, 0]) - np.min(coords[:, 0])) / step)
    y = np.linspace(np.min(coords[:, 1]), np.max(coords[:, 1]), (np.max(coords[:, 1]) - np.min(coords[:, 1])) / step)

    ps = my_interpolating_function(z)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # data = np.zeros((num_atoms, num_e))
    data1 = np.zeros((num_atoms, num_e))

    for jj, j in enumerate(ind):
        # energy, pdos = calc.get_orbital_ldos(a=int(j), npts=num_e)
        # data[j, :] = pdos
        energy, pdos = calc.get_lcao_dos(atom_indices=int(j), npts=num_e, width=0.1)
        data1[jj, :] = pdos

    coords = coords[ind, :]
    en_ind = np.arange(len(energy))[500:]
    data = np.zeros((len(x), len(y), len(z), len(en_ind)))

    for jj, j in enumerate(en_ind):
        print(j)
        Fq = Invdisttree(coords, data1[:, j])
        aaa = Fq(np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T, nnear=2, eps=0, p=8).reshape(X.shape)
        data[:, :, :, jj] = aaa

    # ef = calc.get_fermi_level()
    plt.contourf(z, energy[en_ind][220:-120], np.mean(data, axis=(0, 1)).T[220:-120], 400, cmap='terrain')
    plt.xlabel("Distance (angstroms)", fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.xlim([9, 60])
    plt.ylim([-6.1, -2.7])
    # plt.plot([0, num_atoms], [ef, ef], 'r-')
    plt.show()

    # plt.contourf(z, energy[en_ind], np.mean(data, axis=(0, 1)).T, 400, cmap='terrain')
    # plt.xlabel("Distance (angstroms)", fontsize=12)
    # plt.ylabel("Energy (eV)", fontsize=12)
    # # plt.plot([0, num_atoms], [ef, ef], 'r-')
    # plt.show()

    plt.plot(np.mean(ps, axis=(0, 1)), '-', label=r'$v$')
    plt.show()


def main2(files):

    ldos = []

    for file in files:
        calc = GPAW(file, txt=None)

        num_e = 1000
        num_atoms = len(calc.atoms)

        coords = calc.atoms.get_positions()
        ind = np.argsort(coords[:, 2])

        data1 = np.zeros((num_atoms, num_e))

        for jj, j in enumerate(ind):
            energy, pdos = calc.get_lcao_dos(atom_indices=int(j), npts=num_e, width=0.05)
            data1[jj, :] = pdos

        coords = coords[ind, :]
        en_ind = np.arange(len(energy))[500:]

        num_points = 80

        z = np.linspace(np.min(coords[:, 2]) - 10, np.max(coords[:, 2]) + 10, num_points)
        step = (np.max(coords[:, 2]) - np.min(coords[:, 2])) / num_points
        x = np.linspace(np.min(coords[:, 0]), np.max(coords[:, 0]),
                        (np.max(coords[:, 0]) - np.min(coords[:, 0])) / step)
        y = np.linspace(np.min(coords[:, 1]), np.max(coords[:, 1]),
                        (np.max(coords[:, 1]) - np.min(coords[:, 1])) / step)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        data = np.zeros((len(x), len(y), len(z), len(en_ind)))

        for jj, j in enumerate(en_ind):
            print(j)
            Fq = Invdisttree(coords, data1[:, j])
            aaa = Fq(np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T, nnear=1, eps=0, p=8).reshape(X.shape)
            data[:, :, :, jj] = aaa

        ef = calc.get_fermi_level()
        data = np.mean(data, axis=(0, 1)).T

        ldos.append((energy[en_ind], data[:, len(z) // 2]))

        # plt.contourf(z, energy[en_ind], data, 400)
        # plt.plot([0, num_atoms], [ef, ef], 'r-')
        # plt.show()

    return ldos


def main3():

    calc = GPAW('/home/mk/standing40_second.gpw', txt=None)
    nt = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    ps = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))

    plt.plot(ps)
    plt.show()

def main4():

    from ase.io.trajectory import Trajectory
    tr1 = Trajectory('/home/mk/nci/optB88-vdW.traj')
    tr2 = Trajectory('/home/mk/nci/C09-vdW.traj')

    # calc = GPAW('/home/mk/gpaw_swarm/si_slab_100_1x2_8.gpw', txt=None)
    # calc = GPAW('/home/mk/gpaw_swarm/standing22.gpw', txt=None)
    # calc = GPAW('/home/mk/gpaw_swarm/standing_geometry111/standing111/standing31.gpw', txt=None)

    # calc = GPAW('/home/mk/nci/tc_part_ghosts_6x4.gpw', txt=None)
    # nt = calc.get_pseudo_density()
    # hl = calc.get_homo_lumo()
    # fl = calc.get_fermi_level()
    calc = GPAW('/home/mk/nci/tc/tc_part_ghosts.gpw', txt=None)
    ps1 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    nt1 = calc.get_pseudo_density()
    hl1 = calc.get_homo_lumo()
    fl1 = calc.get_fermi_level()
    calc = GPAW('/home/mk/nci/si/si_part_ghosts.gpw', txt=None)
    ps2 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    nt2 = calc.get_pseudo_density()
    calc = GPAW('/home/mk/nci/tc/tc_part.gpw', txt=None)
    nt3 = calc.get_pseudo_density()
    calc = GPAW('/home/mk/nci/si/si_part.gpw', txt=None)
    nt4 = calc.get_pseudo_density()
    calc = GPAW('/home/mk/nci/standing24f.gpw', txt=None)
    nt = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    ps = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))

    x1 = np.linspace(0, calc.atoms.get_cell_lengths_and_angles()[2], len(ps) // 2)
    x2 = np.linspace(0, calc.atoms.get_cell_lengths_and_angles()[2], len(ps))

    plt.figure(figsize=(12, 4))
    plt.plot(x1 / 10 - 3.45, 10 * (nt - np.mean(nt1, axis=(0, 1)) - np.mean(nt2, axis=(0, 1))), 'k')
    plt.plot(x1 / 10 - 3.45, 10 * (nt - np.mean(nt3, axis=(0, 1)) - np.mean(nt4, axis=(0, 1))), 'k--')
    plt.ylabel(r'$\Delta$ n (e/nm$^{-1}$)', fontsize=12)
    plt.xlabel(r'Distance from interface (nm)', fontsize=12)
    plt.legend(["with counterpoise correction", "without counterpoise correction"])
    plt.show()

    # ps = np.mean(ps, axis=(0, 1))
    from scipy.interpolate import interp1d
    my_interpolating_function = interp1d(calc.hamiltonian.gd.coords(2), ps[::2], kind='cubic',
                                         fill_value="extrapolate")

    en = np.linspace(-20, 0)
    num_e = 1000
    coords = calc.atoms.get_positions()
    ind = np.argsort(coords[:, 2])

    num_atoms = len(coords)

    num_points = 80

    z = np.linspace(np.min(coords[:, 2]) - 10, np.max(coords[:, 2]) + 10, num_points)
    step = (np.max(coords[:, 2]) - np.min(coords[:, 2])) / num_points
    x = np.linspace(np.min(coords[:, 0]), np.max(coords[:, 0]),
                    (np.max(coords[:, 0]) - np.min(coords[:, 0])) / step)
    y = np.linspace(np.min(coords[:, 1]), np.max(coords[:, 1]),
                    (np.max(coords[:, 1]) - np.min(coords[:, 1])) / step)

    ps = my_interpolating_function(z)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # data = np.zeros((num_atoms, num_e))
    data1 = np.zeros((num_atoms, num_e))

    for jj, j in enumerate(ind):
        # energy, pdos = calc.get_orbital_ldos(a=int(j), npts=num_e)
        # data[j, :] = pdos
        energy, pdos = calc.get_lcao_dos(atom_indices=int(j), npts=num_e, width=0.1)
        data1[jj, :] = pdos

    coords = coords[ind, :]
    en_ind = np.arange(len(energy))[500:]
    data = np.zeros((len(x), len(y), len(z), len(en_ind)))

    for jj, j in enumerate(en_ind):
        print(j)
        Fq = Invdisttree(coords, data1[:, j])
        aaa = Fq(np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T, nnear=2, eps=0, p=8).reshape(X.shape)
        data[:, :, :, jj] = aaa

    # ef = calc.get_fermi_level()
    plt.contourf(z, energy[en_ind][220:-120], np.mean(data, axis=(0, 1)).T[220:-120], 400, cmap='terrain')
    plt.xlabel("Distance (angstroms)", fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.xlim([9, 60])
    plt.ylim([-6.1, -2.7])
    # plt.plot([0, num_atoms], [ef, ef], 'r-')
    plt.show()


def main5():

    # # load tc slab with ghost atoms
    # calc = GPAW('/home/mk/nci/tc/tc_part_ghosts_4x2.gpw', txt=None)
    # ps1 = calc.get_electrostatic_potential()
    # # ps1 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    # # nt1 = calc.get_pseudo_density()
    # # hl1 = calc.get_homo_lumo()
    # # fl1 = calc.get_fermi_level()
    #
    # # load si slab with ghost atoms
    # calc = GPAW('/home/mk/nci/si/si_part_ghosts.gpw', txt=None)
    # # ps2 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    # # nt2 = calc.get_pseudo_density()
    # ps2 = calc.get_electrostatic_potential()
    #
    # # load tc slab without ghost atoms
    # calc = GPAW('/home/mk/nci/tc/tc_part.gpw', txt=None)
    # # ps3 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    # # nt3 = calc.get_pseudo_density()
    # ps3 = calc.get_electrostatic_potential()
    #
    # # load si slab without ghost atoms
    # calc = GPAW('/home/mk/nci/si/si_part.gpw', txt=None)
    # # ps4 = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    # # nt4 = calc.get_pseudo_density()
    # ps4 = calc.get_electrostatic_potential()

    # load joint slab
    calc = GPAW('/home/mk/nci/standing82_linearsearch.gpw', txt=None)
    # nt = np.mean(calc.get_pseudo_density(), axis=(0, 1))
    # ps = np.mean(calc.get_electrostatic_potential(), axis=(0, 1))
    ps = calc.get_electrostatic_potential()

    # --------------------------- plot density change -----------------------------

    x1 = np.linspace(0, calc.atoms.get_cell_lengths_and_angles()[2], len(ps) // 2)
    x2 = np.linspace(0, calc.atoms.get_cell_lengths_and_angles()[2], len(ps))

    # plt.figure(figsize=(12, 4))
    # # with ghosts atoms
    # plt.plot(x1 / 10 - 3.391, 10 * (nt - np.mean(nt1, axis=(0, 1)) - np.mean(nt2, axis=(0, 1))), 'k')
    # # without ghosts atoms
    # plt.plot(x1 / 10 - 3.391, 10 * (nt - np.mean(nt3, axis=(0, 1)) - np.mean(nt4, axis=(0, 1))), 'k--')
    #
    # plt.ylabel(r'$\Delta$ n (e/nm$^{-1}$)', fontsize=12)
    # plt.xlabel(r'Distance from interface (nm)', fontsize=12)
    # plt.legend(["with counterpoise correction", "without counterpoise correction"])
    # plt.show()

    # --------------------- interpolate and plot potential ------------------------

    ps = np.mean(ps, axis=(0, 1))

    from scipy.interpolate import interp1d
    my_interpolating_function = interp1d(calc.hamiltonian.gd.coords(2), ps[::2],
                                         kind='cubic',
                                         fill_value="extrapolate")

    num_e = 1000
    coords = calc.atoms.get_positions()    # get atomic coordinates
    ind = np.argsort(coords[:, 2])         # get atomic coordinates

    num_atoms = len(coords)                # number of atoms
    num_points = 80                        # number of points along one of dimensions

    z = np.linspace(np.min(coords[:, 2]) - 10, np.max(coords[:, 2]) + 10, num_points)
    step = (np.max(coords[:, 2]) - np.min(coords[:, 2])) / num_points
    x = np.linspace(np.min(coords[:, 0]), np.max(coords[:, 0]),
                    (np.max(coords[:, 0]) - np.min(coords[:, 0])) / step)
    y = np.linspace(np.min(coords[:, 1]), np.max(coords[:, 1]),
                    (np.max(coords[:, 1]) - np.min(coords[:, 1])) / step)

    ps = my_interpolating_function(z)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    data1 = np.zeros((num_atoms, num_e))

    for jj, j in enumerate(ind):
        # energy, pdos = calc.get_orbital_ldos(a=int(j), npts=num_e)
        # data[j, :] = pdos
        energy, pdos = calc.get_lcao_dos(atom_indices=int(j), npts=num_e, width=0.136)
        data1[jj, :] = pdos

    coords = coords[ind, :]
    en_ind = np.arange(len(energy))[500:]
    data = np.zeros((len(x), len(y), len(z), len(en_ind)))

    for jj, j in enumerate(en_ind):
        print(j)
        Fq = Invdisttree(coords, data1[:, j])
        aaa = Fq(np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T, nnear=4, eps=0, p=2).reshape(X.shape)
        data[:, :, :, jj] = aaa

    # ef = calc.get_fermi_level()
    plt.contourf(z, energy[en_ind][220:-120], np.mean(data, axis=(0, 1)).T[220:-120], 400, cmap='terrain')
    plt.xlabel("Distance (angstroms)", fontsize=12)
    plt.ylabel("Energy (eV)", fontsize=12)
    plt.xlim([9, 60])
    plt.ylim([-6.1, -2.7])
    plt.show()


if __name__ == "__main__":
    # from ase.visualize import view
    # # calc = GPAW('/home/mk/gpaw_swarm/si_part_ghosts.gpw')
    # calc = GPAW('/home/mk/gpaw_swarm/standing22.gpw')
    # view(calc.atoms)

    # files = ['/home/mk/gpaw_swarm/si_slab_100_1x2_8.gpw',
    #          '/home/mk/gpaw_swarm/si_slab_100_1x2_7.gpw',
    #          '/home/mk/gpaw_swarm/si_slab_100_1x2_6.gpw',
    #          '/home/mk/gpaw_swarm/si_slab_100_1x2_4.gpw',
    #          '/home/mk/gpaw_swarm/si_slab_100_1x2_3.gpw']
    #
    # # files = ['/home/mk/gpaw_swarm/si_slab_100_1x2_8.gpw']
    #
    # ldos = main2(files)
    #
    # for item in ldos:
    #     # plt.plot(item[0], item[1])
    #     plt.fill_between(item[0], item[1], alpha=0.0)
    #     plt.ylabel("LDOS (a.u.)")
    #     plt.xlabel("Energy (eV)")
    #
    # plt.show()

    main5()
