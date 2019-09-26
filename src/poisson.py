import numpy as np
from gpaw.poisson import FastPoissonSolver, FDPoissonSolver, create_poisson_solver, BasePoissonSolver
from gpaw.grid_descriptor import GridDescriptor


class PoissonSolver(BasePoissonSolver):
    def __init__(self, solver_name='fast', max_num_iter=10, **kwargs):
        BasePoissonSolver.__init__(self, **kwargs)
        self.max_num_iter = max_num_iter
        self.solver_name = solver_name
        self.solver = create_poisson_solver(name='fast', **kwargs)
        self.delta = []

    def set_grid_descriptor(self, gd):
        self.solver.set_grid_descriptor(gd)

    def solve(self, pot, dens0, permitivity):

        self.delta = []

        conv_criteria = np.max(dens0) * 0.001
        alpha = 0.01
        dens = dens0

        for j in range(self.max_num_iter):
            print(j)
            self.solver.solve(pot, dens)

            aaa, _, _ = np.gradient(eps, dx)
            bbb, _, _ = np.gradient(pot, dx)

            delta = aaa * bbb / eps

            self.delta.append(np.max(delta))
            dens = dens0 - delta * alpha

            if j > 1:
                if np.abs(self.delta[-1] - self.delta[-2]) < conv_criteria:
                    return

        raise RuntimeError('PoissonSolver have not converged.')

    def todict(self):
        self.solver.todict()

    def get_description(self):
        return self.__class__.__name__

    def estimate_memory(self, mem):
        self.solver.estimate_memory(mem)


def laplacian(input_mat, dx):

    top = input_mat[0:-2, 1:-1, 1:-1]
    left = input_mat[1:-1, 0:-2, 1:-1]
    bottom = input_mat[2:, 1:-1, 1:-1]
    right = input_mat[1:-1, 2:, 1:-1]
    center = input_mat[1:-1, 1:-1, 1:-1]
    back = input_mat[1:-1, 1:-1, 0:-2]
    forward = input_mat[1:-1, 1:-1, 2:]
    return np.pad(top + left + bottom + right + back + forward - 6 * center,
                  ((1, 1), (1, 1), (1, 1)),
                  'constant',
                  constant_values=((0, 0), (0, 0), (0, 0))) / dx**2


def gf(x, y, z, x0, y0, z0, sigma):
    return np.exp(-((x-x0)**2+(y-y0)**2+(z-z0)**2) / (2*sigma**2))

num_points = 128
x = np.linspace(0, 100, num_points - 1)
dx = x[3] - x[2]
X, Y, Z = np.meshgrid(x, x, x)
gd = GridDescriptor([num_points, num_points, num_points], cell_cv=[100, 100, 100], pbc_c=False)

dens0 = 10*gf(X, Y, Z, 50, 50, 50, 3)
eps = 1.0*(1.0 - 0.5*np.tanh(3000.0*(Y - 70))-0.5) + 22.0*(0.5*np.tanh(3000.0*(Y - 70))+0.5)
dens0 = dens0 / eps
pot = np.zeros(dens0.shape)


ps = PoissonSolver()
ps.set_grid_descriptor(gd)
ps.solve(pot, dens0, eps)


import matplotlib.pyplot as plt
# plt.plot(np.array(delta))
plt.contour(pot[:, :, 64], 50)
plt.show()