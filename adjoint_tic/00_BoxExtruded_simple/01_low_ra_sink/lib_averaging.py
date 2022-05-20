# noqa: F405
from firedrake import *  # noqa: F403
import numpy as np
import scipy.linalg


class LayerAveraging:
    """Manager to compute a vertical profile of horizontal layer averages"""
    def __init__(self, mesh, r1d, cartesian=True, quad_degree=None):
        """
        :arg mesh: source mesh
        :arg r1d:  array of either y coordinates, or radii self.r, at which depths
                   the averages should be calculated
        :kwarg cartesian: if cartesian, the depths are y values, if not they
                          are radii
        :kwarg quad_degree: quadrature degree for the integrals. If the depths
                            do not exactly match layers in the source mesh, you
                            may have to yank up this number to avoid a noisy
                            result."""
        self.mesh = mesh
        x, y = SpatialCoordinate(mesh)
        if cartesian:
            self.r = y
        else:
            self.r = sqrt(x**2 + y**2)
        if quad_degree is None:
            self.dx = dx
        else:
            self.dx = dx(degree=quad_degree)

        self.r1d = r1d
        # symmetric tri-diagonal mass matrix stored as 2 rows
        # first row contains main diagonal
        # second row contains off diagonal (ignoring last entry)
        self.mass = np.zeros((2, len(r1d)))
        self.rhs = np.zeros(len(r1d))
        self._assemble_mass()

    def _assemble_mass(self):
        """Assembles the mass matrix for the layer average projection"""
        # compute main diagonal of mass matrix
        rc = Constant(self.r1d[0])
        rn = Constant(self.r1d[1])
        rp = Constant(0.)
        # radial P1 hat function between rp<self.r<rn with max. at rc
        Phi = Max(Min((self.r-rp)/(rc-rp), (rn-self.r)/(rn-rc)), 0)
        for i, rin in enumerate(self.r1d[1:]):
            rn.assign(rin)
            self.mass[0, i] = assemble(Phi**2*self.dx)
            rp.assign(rc)
            rc.assign(rn)
        # at this point rp and rc=rn are now the last two vertical nodes
        # final basis function increases from 0 to 1 between rp and rn,
        # and is constant 1 for self.r>rn
        Phi = Max(Min(1, (self.r-rp)/(rn-rp)), 0)
        self.mass[0, -1] = assemble(Phi**2*self.dx)

        # now compute off-diagonal
        rp = Constant(self.r1d[0])
        rn = Constant(self.r1d[1])
        # overlapping product between two basis functions between rp<self.r<rn
        overlap = Max((rn-self.r)/(rn-rp), 0) * Max((self.r-rp)/(rn-rp), 0) * self.dx
        for i, rin in enumerate(self.r1d[1:]):
            rn.assign(rin)
            self.mass[1, i] = assemble(overlap)
            rp.assign(rn)

        # sum of mass matrix should equal domain volume
        # (need to count off-diagonal twice!)
        np.testing.assert_almost_equal(self.mass.sum() + self.mass[1, :].sum(),
                                       assemble(Constant(1, domain=self.mesh)*self.dx))

    def _assemble_rhs(self, T):
        """Assembles the rhs for the layer averaging projection into self.rhs"""
        rc = Constant(self.r1d[0])
        rn = Constant(self.r1d[1])
        rp = Constant(0.)
        # radial P1 hat function between rp<self.r<rn with max. at rc
        Phi = Max(Min((self.r-rp)/(rc-rp), (rn-self.r)/(rn-rc)), 0)
        for i, rin in enumerate(self.r1d[1:]):
            rn.assign(rin)
            self.rhs[i] = assemble(Phi*T*self.dx)
            rp.assign(rc)
            rc.assign(rn)
        # at this point rp and rc=rn are now the last two vertical nodes
        # final basis function increases from 0 to 1 between rp and rn, and is constant 1 for self.r>rn
        Phi = Max(Min(1, (self.r-rp)/(rn-rp)), 0)
        self.rhs[-1] = assemble(Phi*T*self.dx)

    def get_layer_average(self, T):
        """Compute the layer average of Function T at the specified depths and return the array of averages"""
        self._assemble_rhs(T)
        return scipy.linalg.solveh_banded(self.mass, self.rhs, lower=True)

    def extrapolate_layer_average(self, u, avg, DirBCs=None):
        """Given an array of layer averages (such as returned by get_layer_average(), extrapolate these to Function u"""
        # Make sure the field is empty
        u.assign(0.0)

        rc = Constant(self.r1d[0])
        rn = Constant(self.r1d[1])
        rp = Constant(0.)
        # radial P1 hat function between rp<self.r<rn with max. at rc
        Phi = Max(Min((self.r-rp)/(rc-rp), (rn-self.r)/(rn-rc)), 0)
        value = Constant(0.)
        for a, rin in zip(avg[:-1], self.r1d[1:]):
            value.assign(a)
            rn.assign(rin)
            u.interpolate(u + value*Phi)
            rp.assign(rc)
            rc.assign(rn)
        # at this point rp and rc=rn are now the last two vertical nodes
        # final basis function increases from 0 to 1 between rp and rn, and is constant 1 for self.r>rn
        Phi = Max(Min(1, (self.r-rp)/(rn-rp)), 0)
        value.assign(avg[-1])
        u.interpolate(u + value*Phi)

        # todo: Look at why at the boundaries we have inaccurate values
        if DirBCs:
            for bc in DirBCs:
                bc.apply(u)
