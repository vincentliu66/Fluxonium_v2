"""
Fluxonium circuit object

Author: Long Nguyen, lbnguyen@lbl.gov
Date: July 1, 2020
"""

import numpy as np
import qutip as qt
import sys
sys.dont_write_bytecode = True

class Fluxonium_qubit(object):
    """A class representing fluxonium qubit"""

    def __init__(self, E_J, E_C, E_L, phi_ext,
                 nhilbert=25, nlevels=20):
        """The Hamiltonian is H = 4E_C*n^2 +0.5E_L*phi^2 -E_J*cos(phi+phi_ext)"""
        """Energies have unit of GHz"""
        self.E_J = E_J #Josephson energy
        self.E_C = E_C #Charging energy
        self.E_L = E_L #Inductive energy
        self.phi_ext = phi_ext #External normalized flux
        self.nhilbert = nhilbert #Size of matrix
        self.nlevels = nlevels #Number of energy levels simulated, nlevels < nhilbert
        self.type = 'qubit'

    def __str__(self):
        text = ('A fluxonium qubit with E_L = {} '.format(self.E_L) + 'GHz'
             + ', E_C = {} '.format(self.E_C) + 'GHz'
             + ', and E_J = {} '.format(self.E_J) + 'GHz'
             + '. The external flux is {} (Phi_o).'.format(
             self.phi_ext/np.pi/2.0))
        return text

    @property
    def E_J(self):
        return self._E_J

    @property
    def E_C(self):
        return self._E_C

    @property
    def E_L(self):
        return self._E_L

    @property
    def phi_ext(self):
        return self._phi_ext

    @property
    def nhilbert(self):
        return self._nhilbert

    @property
    def nlevels(self):
        return self._nlevels

    @E_J.setter
    def E_J(self, value):
        if value <= 0:
            raise Exception('Josephson energy must be positive.')
        else:
            self._E_J = value

    @E_C.setter
    def E_C(self, value):
        if value <= 0:
            raise Exception('Charging energy must be positive.')
        else:
            self._E_C = value

    @E_L.setter
    def E_L(self, value):
        if value <= 0:
            raise Exception('Inductive energy must be positive.')
        else:
            self._E_L = value

    @phi_ext.setter
    def phi_ext(self, value):
        self._phi_ext = value

    @nhilbert.setter
    def nhilbert(self, value):
        if value <= 0:
            raise Exception('Matrix dimension must be positive.')
        else:
            self._nhilbert = value

    @nlevels.setter
    def nlevels(self, value):
        if value <= 0:
            raise Exception('Number of levels must be positive.')
        else:
            self._nlevels = value

    def _phi_lc(self):
        """Phase operator in the LC basis."""
        return (8 * self.E_C / self.E_L) ** (0.25) * qt.position(self.nhilbert)


    def _n_lc(self):
        """Charge operator in the LC basis."""
        return (self.E_L / (8 * self.E_C)) ** (0.25) * qt.momentum(self.nhilbert)


    def _hamiltonian_lc(self):
        """Qubit Hamiltonian."""
        E_C = self.E_C
        E_L = self.E_L
        E_J = self.E_J
        phi = self._phi_lc()
        n = self._n_lc()
        delta_phi = phi + self.phi_ext
        return 4 * E_C * n ** 2 + 0.5 * E_L * phi ** 2 - E_J * delta_phi.cosm()


    def _eigenspectrum_lc(self):
        """Eigenenergies and eigenstates in the LC basis.
        Hold values previously calculated, unless cache is reset"""
        if self._eigvals is None or self._eigvecs is None:
            H_lc = self._hamiltonian_lc()
            self._eigvals, self._eigvecs = H_lc.eigenstates()
        return self._eigvals, self._eigvecs

    def level(self, level_index):
        """Energy of a single level of the qubit.

        Parameters
        ----------
        level_ind : int
            The qubit level starting from zero.

        Returns
        -------
        float
            Energy of the level.
        """
        if level_index < 0 or level_index >= self.nhilbert:
            raise Exception('The level is out of bounds')
        return self._eigenspectrum_lc()[level_index]

    def freq(self, level1, level2):
        """Transition energy/frequency between two levels of the qubit.

        Parameters
        ----------
        level1, level2 : int
            The qubit levels.

        Returns
        -------
        float
            Transition energy/frequency between `level1` and `level2` defined
            as the difference of energies. Positive if `level1` < `level2`.
        """
        return self.level(level2) - self.level(level1)

    def H(self, nlev=None):
        """Qubit Hamiltonian in its eigenbasis.

        Parameters
        ----------
        nlev : int, optional
            The number of qubit eigenstates if different from `self.nlev`.

        Returns
        -------
        :class:`qutip.Qobj`
            The Hamiltonian operator.
        """
        return qt.Qobj(np.diag(self.levels(nlev=nlev)))

    def one(self, nlev=None):
        """Identity operator in the qubit eigenbasis.

        Parameters
        ----------
        nlev : int, optional
            The number of qubit eigenstates if different from `self.nlev`.

        Returns
        -------
        :class:`qutip.Qobj`
            The identity operator.
        """
        if nlev is None:
            nlev = self.nlev
        if nlev < 1 or nlev > self.nhilbert:
            raise Exception('`nlev` is out of bounds.')
        return qt.qeye(nlev)

    def phi(self, nlev=None):
        """Generalized-flux operator in the qubit eigenbasis.

        Parameters
        ----------
            The number of qubit eigenstates if different from `self.nlev`.

        Returns
        -------
        :class:`qutip.Qobj`
            The flux operator.
        """
        if nlev is None:
            nlev = self.nlev
        if nlev < 1 or nlev > self.nhilbert:
            raise Exception('`nlev` is out of bounds.')
        _, evecs = self._eigenspectrum_lc()
        phi_op = np.zeros((nlev, nlev), dtype=complex)
        for index1 in range(nlev):
            for index2 in range(nlev):
                phi_op[ind1, ind2] = self._phi_lc().matrix_element(
                    evecs[index1].dag(), evecs[index2])
        return qt.Qobj(phi_op)

    def n(self, nlev=None):
        """Charge operator in the qubit eigenbasis.

        Parameters
        ----------
            The number of qubit eigenstates if different from `self.nlev`.

        Returns
        -------
        :class:`qutip.Qobj`
            The charge operator.
        """
        if nlev is None:
            nlev = self.nlev
        if nlev < 1 or nlev > self.nhilbert:
            raise Exception('`nlev` is out of bounds.')
        _, evecs = self._eigenspectrum_lc()
        n_op = np.zeros((nlev, nlev), dtype=complex)
        for ind1 in range(nlev):
            for ind2 in range(nlev):
                n_op[ind1, ind2] = self._n_lc().matrix_element(
                    evecs[ind1].dag(), evecs[ind2])
        return qt.Qobj(n_op)

    def phi_ij(self, level1, level2):
        """The flux matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            The qubit levels.

        Returns
        -------
        complex
            The matrix element of the flux operator.
        """
        if (level1 < 0 or level1 > self.nlev_lc
                or level2 < 0 or level2 > self.nlev_lc):
            raise Exception('Level index is out of bounds.')
        _, evecs = self._eigenspectrum_lc(eigvecs_flag=True)
        return self._phi_lc().matrix_element(
            evecs[level1].dag(), evecs[level2])

    def n_ij(self, level1, level2):
        """The charge matrix element between two eigenstates.

        Parameters
        ----------
        level1, level2 : int
            The qubit levels.

        Returns
        -------
        complex
            The matrix element of the charge operator.
        """
        if (level1 < 0 or level1 > self.nlev_lc
                or level2 < 0 or level2 > self.nlev_lc):
            raise Exception('Level index is out of bounds.')
        _, evecs = self._eigenspectrum_lc(eigvecs_flag=True)
        return self._n_lc().matrix_element(evecs[level1].dag(), evecs[level2])

    def potential(self, phi_points=None):
        """Returns the fluxonium potential as a function of flux.

        Parameters
        ----------
        phi_points : float or ndarray, optional
            Value(s) of the flux variable.

        Returns
        -------
        phi_points : ndarray, if not given as a parameter
            Values of the flux variable.
        V_phi : float or ndarray
            Potential energy at the flux values from `phi_points`.

        Examples
        --------
        >>> plt.plot(phi_values, qubit.potential(phi_values))
        >>> plt.plot(qubit.potential()[0], qubit.potential()[1])
        """
        if phi_points is None:
            phi_points = np.linspace(-np.pi, np.pi, 201)
            phi_flag = False
        else:
            phi_flag = True
        if not (isinstance(phi_points, np.ndarray)
                or isinstance(phi_points, float)
                or isinstance(phi_points, int)):
            raise Exception('Inappropriate type of `phi`')
        V_phi = (0.5 * self.E_L * phi_points ** 2
                 - self.E_J * np.cos(phi_points - self.phi_ext))
        if phi_flag:
            return V_phi
        else:
            return phi_points, V_phi

    def wavefunc(self, level_ind, phi_points=None):
        """Returns the flux-dependent wavefunction of an eigenstate.

        Parameters
        ----------
        level_ind : int
            The level of the qubit.
        phi_points : float or ndarray, optional
            Value(s) of the flux variable.

        Returns
        -------
        phi_points : ndarray, if not given as a parameter
            Values of the flux variable.
        wfunc : float or ndarray
            The wavefunction at the flux value from `phi`.

        Examples
        --------
        >>> plt.plot(phi_values, qubit.wavefunc(1, phi_values))
        >>> plt.plot(qubit.wavefunc(2)[0], qubit.wavefunc(2)[1])
        """
        if level_ind < 0 or level_ind > self.nlev_lc:
            raise Exception('The level is out of bounds.')
        if phi_points is None:
            phi_points = np.linspace(-np.pi, np.pi, 201)
            phi_flag = False
        else:
            phi_flag = True
        if not (isinstance(phi_points, np.ndarray)
                or isinstance(phi_points, float)
                or isinstance(phi_points, int)):
            raise Exception('Inappropriate type of `phi`')

        def ho_wf(phi, l, E_C, E_L):
            # Calcultes wave functions of the harmonic oscillator.
            ratio = (8.0 * E_C / E_L) ** (0.25)
            coeff = (2.0 ** l * np.math.factorial(l)
                     * np.sqrt(np.pi) * ratio) ** (-0.5)
            return coeff * np.exp(-0.5 * (phi / ratio) ** 2) * hpoly(l, phi / ratio)

        wfunc = 0
        _, evecs = self._eigenspectrum_lc(eigvecs_flag=True)
        for lvl_idx in range(self._nlev_lc):
            coeff = np.real(evecs[level_ind].full()[lvl_idx, 0])
            wfunc = wfunc + coeff * ho_wf(
                phi_points, lvl_idx, self.E_C, self.E_L)
        if phi_flag:
            return wfunc
        else:
            return phi_points, wfunc