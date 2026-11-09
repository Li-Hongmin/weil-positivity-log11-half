"""Exact algebraic checks of the two boundary-response identities."""
import sympy as sp


def check():
    t, x, y = sp.symbols("t x y", positive=True)
    # t=exp(-c*a), x=exp(c*physical_x), y=exp(c*physical_y), x<=y.
    green_even = (x + 1 / x) * (1 / (t * y) + t * y) / (1 / t - t)
    green_odd = (x - 1 / x) * (1 / (t * y) + t * y) / (1 / t + t)
    convolution_even = (x + 1 / x) / y
    convolution_odd = (x - 1 / x) / y
    response_even = t**2 * (x + 1 / x) * (y + 1 / y) / (1 - t**2)
    response_odd = t**2 * (x - 1 / x) * (y - 1 / y) / (1 + t**2)
    assert sp.cancel(green_even - convolution_even - response_even) == 0
    assert sp.cancel(green_odd - convolution_odd - response_odd) == 0
    c, frequency = sp.symbols("c frequency", positive=True)
    assert sp.cancel(2 / c - 2 * c / (c**2 + frequency**2)
                     - 2 * frequency**2 / (c * (c**2 + frequency**2))) == 0
    print("even and odd Green boundary-response identities: PASS")
    print("digamma summand rational identity: PASS")


if __name__ == "__main__":
    check()
