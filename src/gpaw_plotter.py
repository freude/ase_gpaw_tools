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

    # calc = GPAW('/home/mk/gpaw_swarm/si_slab_100_1x2_8.gpw', txt=None)
    # calc = GPAW('/home/mk/gpaw_swarm/standing22.gpw', txt=None)
    # calc = GPAW('/home/mk/gpaw_swarm/standing_geometry111/standing111/standing31.gpw', txt=None)
    calc = GPAW('/home/mk/standing24f.gpw', txt=None)

    nt = calc.get_pseudo_density()
    ps = calc.get_electrostatic_potential()
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

    main1()
