from ase.io import read
from ase.visualize import view


si = read('si_slab_100_1x2.cif', format='cif')
view(si)