"""
Numeric strategy for calculation. The exp, sin and cos functions are calculated using numpy,
numpy with truncation of the arguments, or arbitrary precision using mpmath.
"""

import numpy

class CalculationStrategy(object):
    """Abstract strategy for calculation. Can be plain python or arbitrary precision like mpmath."""
    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def sin(self, power):
        """Sin to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")

    def cos(self, power):
        """Cos to the power.

        Parameters
        ----------
        power :
            The power to raise to.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Raises
        ------
        Exception
            Must override this method.

        """
        raise Exception("Must override this method.")


class CalculationStrategyMPMath(CalculationStrategy):
    """Use mpmath for calculation."""

    def __init__(self):
        """
        Constructor.
        """
        import mpmath
        self.mpmath_sin = numpy.vectorize(mpmath.sin)
        self.mpmath_cos = numpy.vectorize(mpmath.cos)
        self.mpmath_exp = numpy.vectorize(mpmath.exp)
        self.mpmath_mpc = mpmath.mpc
        # Use 32 digits in mpmath calculations.
        mpmath.mp.dps = 32

    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value : float, complex or numpy array
            Initial value of the variable.

        Returns
        -------
        instance of CalculationStrategyMPMath
            variable.

        """

        if not(isinstance(initial_value, numpy.ndarray)):
            initial_value = numpy.array(initial_value)

        if initial_value.size == 1:
            mpc = self.mpmath_mpc(complex(initial_value.real) + 1j * complex(initial_value.imag))
        else:
            mpc = self.mpmath_mpc(complex(1) + 1j * complex(0)) * initial_value

        return mpc

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power : float
            The power to raise to.

        Returns
        -------
        mpmath variable
            Exponential.

        """
        return self.mpmath_exp(power)

    def sin(self, power):
        """Sin function.

        Parameters
        ----------
        power : float or numpy array
            The arg of sin.

        Returns
        -------
        mpmath variable
            Sin.

        """
        return self.mpmath_sin(power)

    def cos(self, power):
        """Cos function.

        Parameters
        ----------
        power : float or numpy array
            The arg of cos.

        Returns
        -------
        mpmath variable
            Cos.

        """
        return self.mpmath_cos(power)

    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            variable to convert.

        Returns
        -------
        numpy array
            Native python complex variable.

        """
        return numpy.array(variable, dtype=complex)


class CalculationStrategyNumpy(CalculationStrategy):
    """Use plain python for calculation."""
    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Returns
        -------
        instance of CalculationStrategyNumpy
            variable.

        """
        return initial_value + 0j # complex(initial_value)

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power : float
            The power to raise to.

        Returns
        -------
        numpy array
            Exponential.

        """
        try:
            ans =  numpy.exp(power)
        except:
            ans = float("Inf")
        return ans

    def sin(self, power):
        """Sin function.

        Parameters
        ----------
        power :
            The sin argument.

        Returns
        -------
        numpy array
            Sin.

        """
        return numpy.sin(power)

    def cos(self, power):
        """Cos function.

        Parameters
        ----------
        power :
            The coa argument.

        Returns
        -------
        numpy array
            Cos.

        """
        return numpy.cos(power)


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Returns
        -------
        numpy array (complex)
            Native python complex variable.

        """
        return complex(variable)

class CalculationStrategyNumpyTruncated(CalculationStrategy):
    """Use plain python for calculation."""
    def __init__(self, limit=1000):
        self.limit = limit

    def createVariable(self, initial_value):
        """Factory method for calculation variable.

        Parameters
        ----------
        initial_value :
            Initial value of the variable.

        Returns
        -------
        instance of CalculationStrategyNumpy
            variable.

        """
        return initial_value + 0j # complex(initial_value)

    def exponentiate(self, power):
        """Exponentiates to the power.

        Parameters
        ----------
        power : float
            The power to raise to.

        Returns
        -------
        numpy array
            Exponential.

        """
        if power.size == 1:
            power1 = numpy.array([power], dtype=numpy.complex128)
        else:
            power1 = numpy.array(power, dtype=numpy.complex128)

        if numpy.any(power1.real > self.limit):
            ii = numpy.where(power1.real > self.limit)
            power1[ii] = self.limit + power1.imag[ii] * 1j

        try:
            ans =  numpy.exp(power1)
        except:
            ans = float("Inf")

        return ans

    def sin(self, power):
        """Sin function.

        Parameters
        ----------
        power :
            The sin argument.

        Returns
        -------
        numpy array
            Sin.

        """
        if power.size == 1:
            power1 = numpy.array([power])
        else:
            power1 = numpy.array(power)

        if numpy.any(power1.imag > self.limit):
            ii = numpy.where(power1.imag > self.limit)
            power1[ii] = power1.real[ii] + self.limit * 1j
        if numpy.any(power1.imag < -self.limit):
            ii = numpy.where(power1.imag < self.limit)
            power1[ii] = power1.real[ii] - self.limit * 1j
        ans = numpy.sin(power1)

        return ans

    def cos(self, power):
        """Cos function.

        Parameters
        ----------
        power :
            The coa argument.

        Returns
        -------
        numpy array
            Cos.

        """
        if power.size == 1:
            power1 = numpy.array([power])
        else:
            power1 = numpy.array(power)

        if numpy.any(power1.imag > self.limit):
            ii = numpy.where(power1.imag > self.limit)
            power1[ii] = power1.real[ii] + self.limit * 1j
        if numpy.any(power1.imag < -self.limit):
            ii = numpy.where(power1.imag < self.limit)
            power1[ii] = power1.real[ii] - self.limit * 1j
        ans = numpy.cos(power1)

        return ans


    def toComplex(self, variable):
        """Converts calculation variable to native python complex.

        Parameters
        ----------
        variable :
            Calculation variable to convert.

        Returns
        -------
        numpy array (complex)
            Native python complex variable.

        """
        return complex(variable)


if __name__ == "__main__":
    a = CalculationStrategyMPMath()

    a0 = a.createVariable(numpy.array(0))
    print(type(a))
    print("cos(0)", a.cos(a0))

    api = a.createVariable(numpy.array(numpy.pi))
    print("cos(pi)", a.cos(api))

    api = a.createVariable(numpy.array([numpy.pi] * 10))
    print("cos(pi)", a.cos(api))
    # print("exp(pi)", a.exponentiate(api))

