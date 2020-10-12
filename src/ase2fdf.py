from ase.units import Ry, eV, Ha
from ase.io.trajectory import Trajectory
from ase.visualize import view
from ase.calculators.siesta import Siesta


traj = Trajectory('/home/mk/nci/C09-vdW.traj')
atoms = traj[-1]

siesta = Siesta(
    mesh_cutoff=450 * Ry,
    basis_set='DZP',
    xc="BH",
    energy_shift=(10 * 10**-3) * eV,
    pseudo_path = '.',
    fdf_arguments={
        'SCFMustConverge': False,
        'COOP.Write': True,
        'WriteDenchar': True,
        'PAO.BasisType': 'split',
        "PAO.SoftDefault": True,
        'DM.MixingWeight': 0.01,
        'MaxSCFIterations': 300,
        'DM.NumberPulay': 4,
        'XML.Write': True,
        'DM.UseSaveDM': True})

print(siesta.prefix)
print(siesta.getpath(ext='fdf'))
siesta.write_input(atoms, properties=['energy'])
