"""
Class to compute pseudo analytical solution for suspension development

from Hjelmfelt and Lenau
Journal of the Hydraulic Division  (1970)
"""


import numpy as np
import matplotlib.pyplot as plt
from math import factorial
from scipy.special import gamma, digamma
from scipy.optimize import minimize


class suspensionDevSol:

    def __init__(self, RouseNumber, Aref, nRoots=5, nSeries=15):
        """
        RouseNumber: float, Rouse number
        Aref: float, dimensionless reference level, 0 < Aref < 1
        nRoots: int, number of roots alphaK computed for solution
        nSeries: int, number of element in infinite sum, default is 15
        """
        # dimension less reference level
        self.Aref = Aref
        # compute Rouse number
        self.Ro = RouseNumber
        # number of alphaK roots, P(A,alphaK)=0
        self.nRoots = nRoots
        # number of term to compute to approximate infinite sum
        self.nSeries = nSeries
        self._checkParameters()
        # roots alphaK for solution, Sturm Liouville
        self.alphasKroots = np.zeros(self.nRoots)

# - - - - - PRIVATE MEMBER FUNCTIONS - - - - - #
    def _checkParameters(self):
        """
        verify parameters of suspesionDevSol instance
        raise error in case of a wrong input
        """
        # check type of entries
        if not isinstance(self.Aref, float):
            raise TypeError("wrong type for Aref, must be a float")
        if not isinstance(self.Ro, float):
            raise TypeError("wrong type for RouseNumber, must be a float")
        if not isinstance(self.nRoots, int):
            raise TypeError("wrong type for nRoots, must be an int")
        if not isinstance(self.nSeries, int):
            raise TypeError("wrong type for nSeries, must be an int")
        # check value of Aref and RouseNumber
        if self.Aref <= 0 or self.Aref >= 1:
            raise ValueError("Aref value is in interval ]0, 1[")
        if self.Ro <= 0:
            raise ValueError("Rouse number value has to be greater than 0")
        return

    def _coefHyper(self, a, n):
        """
        compute falling factorial of a, n
        a: float
        n: int
        ----------------
        (a)n = 1, if n=0
        (a)n+1 = (a+n) * (a)n
        """
        coef = 1. + np.zeros_like(a)
        for i in range(n):
            coef *= a + i
        return coef

    def _hyperGeoFunc(self, a, b, c, Z):
        """compute the hyper geometric function"""
        Fout = 0.
        for i in range(self.nSeries):
            ai = self._coefHyper(a, i)
            bi = self._coefHyper(b, i)
            ci = self._coefHyper(c, i)
            ifac = factorial(i)
            fi = ((ai * bi) / ci) * (Z**i / ifac)
            Fout += fi
        return Fout

    def _F1(self, Z, alpha):
        """
        Z: float
        alpha: float
        eq 62
        """
        Fout = self._hyperGeoFunc(
            0.5+self.Ro+alpha, 0.5+self.Ro-alpha, 1+self.Ro, Z)
        return Fout

    def _dF1dAlpha(self, alpha):
        """
        partial derivative of F1 as regard to alpha
        evaluated in A
        eq 73
        """
        dFout = 0.
        for i in range(1, self.nSeries+1):
            dfi = self._coefHyper(0.5+self.Ro+alpha, i)
            dfi *= self._coefHyper(
                0.5+self.Ro-alpha, i) / self._coefHyper(1+self.Ro, i)
            dfi *= self.Aref**i / factorial(i)
            sumOverj = 0.
            for j in range(1, i+1):
                sumOverj += 1/(
                    (self.Ro+alpha+j-0.5)*(self.Ro+-alpha+j-0.5))
            dfi *= sumOverj
            dFout += dfi
        return -2 * alpha * dFout

    def _dF1dZ(self, alpha):
        """
        partial derivative of F1 with respect to alpha
        evaluated in A
        eq 70
        """
        dFout = (0.5+self.Ro+alpha)*(0.5+self.Ro-alpha)/(1+self.Ro)
        for i in range(1, self.nSeries+1):
            dfi = self._coefHyper(0.5+self.Ro+alpha, i)
            dfi *= self._coefHyper(
                0.5+self.Ro-alpha, i) / self._coefHyper(1+self.Ro, i)
            dfi *= self.Aref**(i-1) / factorial(i-1)
            dFout += dfi
        return dFout

    def _F2(self, Z, alpha):
        """
        Z: float
        alpha: float
        eq 63
        """
        Fout = self._hyperGeoFunc(0.5+alpha, 0.5-alpha, 1-self.Ro, Z)
        return Fout

    def _dF2dAlpha(self, alpha):
        """
        partial derivative of F2 with respect to alpha
        evaluated in A
        eq 74
        """
        dFout = 0.
        for i in range(1, self.nSeries+1):
            dfi = self._coefHyper(0.5+alpha, i)
            dfi *= self._coefHyper(0.5-alpha, i)

            dfi *= (self.Aref**i / factorial(i))
            dfi /= self._coefHyper(1-self.Ro, i)

            sumOverj = 0.
            for j in range(1, i+1):
                sumOverj += 1/((alpha+j-0.5)*(-alpha+j-0.5))
            dfi *= sumOverj
            dFout += dfi
        return -2*alpha * dFout

    def _dF2dZ(self, alpha):
        """
        partial derivative of F2 with respect to Z
        evaluated in A
        eq 71
        """
        dFout = 0.
        for i in range(1, self.nSeries+1):
            dfi = self._coefHyper(
                0.5+alpha, i) * self._coefHyper(0.5-alpha, i)
            dfi *= (self.Aref**(i-1) / factorial(i-1)) / self._coefHyper(
                1-self.Ro, i)
            dFout += dfi
        return dFout

    def _PZalpha(self, Z, alpha):
        """
        P function as a function of Z and alpha
        eq 21
        """
        Pout = (1-Z)**self.Ro
        Pout *= self._hyperGeoFunc(
            0.5+self.Ro+alpha, 0.5+self.Ro-alpha, 1+self.Ro, 1-Z)
        return Pout

    def _PAalpha(self, alpha):
        """
        P function evaluated in A (reference level)
        eq 68
        """
        Pout = (gamma(1+self.Ro) * gamma(-self.Ro)/np.pi)
        Pout *= np.sin(np.pi*(0.5+alpha)) * self._F1(self.Aref, alpha)
        Pout += (self.Aref**(-self.Ro)
                 * (gamma(1+self.Ro) * gamma(self.Ro)
                    / (np.pi * gamma(0.5+self.Ro+alpha)))
                 * np.sin(np.pi*(0.5+alpha-self.Ro)) * gamma(
                     0.5+alpha-self.Ro)*self._F2(self.Aref, alpha))
        Pout *= (1-self.Aref)**self.Ro
        return Pout

    def _dPAdAlpha(self, alpha):
        """
        partial derivative of P with respect to alpha
        evaluated in A (dimensionless reference level)
        eq 72
        """
        dP1 = np.pi * np.cos(np.pi*(0.5+alpha))*self._F1(self.Aref, alpha)
        dP1 += np.sin(np.pi*(0.5+alpha))*self._dF1dAlpha(alpha)
        dP1 *= (gamma(1+self.Ro) * gamma(-self.Ro) / np.pi)
        dP2 = (digamma(0.5-self.Ro+alpha)
               - digamma(0.5+self.Ro+alpha))*self._F2(self.Aref, alpha)
        dP2 += self._dF2dAlpha(alpha)
        dP2 *= np.sin(np.pi*(0.5-self.Ro+alpha))
        dP2 += np.pi*np.cos(
            np.pi*(0.5-self.Ro+alpha))*self._F2(self.Aref, alpha)
        dP2 *= (gamma(0.5-self.Ro+alpha)/gamma(0.5+self.Ro+alpha))
        dP2 *= self.Aref**(-self.Ro) * (gamma(1+self.Ro)*gamma(self.Ro)/np.pi)
        dPout = (1-self.Aref)**self.Ro * (dP1 + dP2)
        return dPout

    def _dPAdZ(self, alpha):
        """
        partial derivative of P with respect to Z
        evaluated in A (dimensionless reference level)
        eq 69
        """
        dP1 = gamma(self.Ro+1)*gamma(-self.Ro)/np.pi
        dP1 *= np.sin(np.pi*(0.5+alpha))*self._dF1dZ(alpha)
        dP2 = gamma(1+self.Ro)*gamma(self.Ro)/(np.pi*gamma(0.5+self.Ro+alpha))
        dP2 *= np.sin(np.pi*(0.5-self.Ro+alpha))*gamma(0.5-self.Ro+alpha)
        dP2 *= (-self.Ro * self.Aref**(-self.Ro-1)
                * self._F2(self.Aref, alpha)
                + self.Aref**(-self.Ro)*self._dF2dZ(alpha))

        dPout = dP1 + dP2
        dPout *= (1-self.Aref)**(self.Ro)
        dPout += -self.Ro*self._PAalpha(alpha)/(1-self.Aref)
        return dPout

    def _integralCav(self):
        """
        compute integral from A to 1 of: Y -> ((1-Y) / Y))^Ro
        necessary to compute suspension averaged over water depth
        """
        # evaluate integral from 0 to 1
        int_01 = np.pi*self.Ro / np.sin(np.pi*self.Ro)
        # evaluate integral from 0 to 1, binomial expansion
        int_0A = 0.
        for i in range(self.nSeries):
            binomCoef = self._coefHyper(self.Ro, i) / factorial(i)
            int_0Ai = (binomCoef * (-1)**i
                       * self.Aref**(1+i-self.Ro) / (1+i-self.Ro))

            int_0A += int_0Ai
        return int_01 - int_0A

    def _getAlphasKroots(self):
        """
        Obtain alphaK coefficients minimizing cost functions
        in sequence, find 1st root and move to the next one
        """
        if np.all(self.alphasKroots != 0):
            return np.copy(self.alphasKroots)
        print("find root of P(A, alpha)")

        def costFunction(alpha):
            return np.abs(self._PAalpha(alpha))
        # initial values of alpha coef
        self.alphasKroots[0] = self.Ro + 0.6
        print(f"start point for alpha0 : {self.alphasKroots[0]}")
        self.alphasKroots[0] = minimize(
            costFunction, self.alphasKroots[0]).x[0]
        for i in range(1, self.nRoots):
            self.alphasKroots[i] = self.alphasKroots[i-1] + 1.1
            self.alphasKroots[i] = minimize(
                costFunction, self.alphasKroots[i]).x[0]
        return np.copy(self.alphasKroots)

# - - - - - PUBLIC MEMBER FUNCTIONS - - - - - #
    def get_solution(self, X, Z):
        """
        return dimensionless concentration field
        suspension development as a function of X, Z
        X: dimensionless horizontal distance
        Z: dimensionless vertical distance
        eq 29
        """
        alphasK = self._getAlphasKroots()
        print(f"roots of P(A, alpha), alphasK: {alphasK}")
        C = (self.Aref/(1-self.Aref))**self.Ro * ((1-Z)/Z)**self.Ro
        for alphaK in alphasK:
            C += 2 * (
                alphaK*self._PZalpha(Z, alphaK)
                / ((alphaK**2 - 0.25) * self._dPAdAlpha(alphaK))
            ) * np.exp(-X*(alphaK**2 - 0.25))
        return C

    def get_Caveraged(self, X):
        """
        Return depth integrated average suspended concentration
        eq 31
        """
        Cav = np.zeros_like(X)
        alphasK = self._getAlphasKroots()
        # integral from A to 1 of, Z -> ((1-Z)/Z)^Ro
        int_A1 = self._integralCav()
        print(r"int_A^1 ((1-Y)/Y)^Ro = " + f"{int_A1}")
        for k in range(self.nRoots):
            alphak = alphasK[k]
            coefi = alphak * self._dPAdZ(alphak)
            coefi /= (alphak**2 - 0.25)*(self._dPAdAlpha(alphak))
            coefi *= np.exp(-X*(alphak**2 - 0.25))
            print(f"coefi: max={np.max(coefi)}; min = {np.min(coefi)}")
            Cav += 2*self.Aref * coefi
        Cav += int_A1 * self.Aref**self.Ro / (1-self.Aref)**(1+self.Ro)
        return Cav

    def getDimCoords(self, X, Z, Hwater, Umean, ustar, Sc=1.):
        """
        return dimensioned coordinates x and z
        X: array-like, dimensionless horizontal coordinate
        Z: array like, dimensionless vertical coordinate
        Hwater: float, water depth, in meters
        Umean, float, mean flow velocity
        ustar: float, friction velocity, m/s
        Sc: float, Schmidt number, defailt is 1
        """
        kappa = 0.41  # von Karman constant
        x = X * (Hwater * Umean) / (kappa * Sc * ustar)
        z = Z * Hwater
        return x, z
