import math
from typing import Union
from decimal import Decimal


# This is a sample Python script.

def rationalize(number: Union[Decimal, str, float]):
    def gcd_decimal(a: Decimal, b: Decimal):
        return math.gcd(int(a), int(b))

    if isinstance(number, float):
        decimal_places = len(str(number) % 1) - 2
        denominator = int(10 ** decimal_places)
        nominator = int(number * denominator)
    else:
        decimal_places = len(str(Decimal(number) % 1)) - 2
        denominator = int(10 ** decimal_places)
        nominator = int(Decimal(number) * denominator)

    gcd = gcd_decimal(nominator, denominator)

    while gcd > 1:
        nominator = nominator / gcd
        denominator = denominator / gcd
        gcd = gcd_decimal(nominator, denominator)

    return nominator, denominator


def bezout1(a: int, b: int):
    r = b
    x = a  # becomes gcd(a, b)
    s = 0
    y = 1  # the coefficient of a
    t = 1
    z = 0  # the coefficient of b
    while r:
        q = x // r
        x, r = r, x % r
        y, s = s, y - q * s
        z, t = t, z - q * t
    return y % (b / x), z % (-a / x)  # modulus in this way so that y is positive and z is negative

def bezout2(a: int, b: int):

    class Step():
        def __init__(self, quotient, remainder, bezout_s:int, bezout_t:int):
            self.quotient = quotient
            self.remainder = remainder
            self.bezout_s = bezout_s
            self.bezout_t = bezout_t

    algorithm = [Step(None, a, 1, 0), Step(None, b, 0, 1)]

    while algorithm[-2].remainder % algorithm[-1].remainder > 0:
        quotient = algorithm[-2].remainder // algorithm[-1].remainder
        algorithm.append(Step(quotient=quotient,
                              remainder=algorithm[-2].remainder % algorithm[-1].remainder,
                              bezout_s=algorithm[-2].bezout_s - quotient * algorithm[-1].bezout_s,
                              bezout_t=algorithm[-2].bezout_t - quotient * algorithm[-1].bezout_t
                              )
                         )
    assert algorithm[-1].remainder == math.gcd(a, b)
    return algorithm[-1].bezout_s, algorithm[-1].bezout_t, algorithm[-1].remainder


def padic(nominator: Union[Decimal, int], denominator: Union[Decimal, int], prime: Union[int, Decimal] = 5):
    bezout = bezout2(int(denominator), int(prime))

    return ''


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # p = padic(Decimal(1), Decimal(3))
    # print(p)

    print(bezout2(3, 5))