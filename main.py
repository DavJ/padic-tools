import math
from typing import Union
from decimal import Decimal

# This is a sample Python script.

def rationalize(number: Union[Decimal, str]):

    def gcd_decimal(a: Decimal, b: Decimal):
        return math.gcd(int(a), int(b))

    decimal_places = len(str(Decimal(number) % 1)) - 2
    denominator = Decimal(10 ** decimal_places)
    nominator = Decimal(number) * denominator

    gcd = gcd_decimal(nominator, denominator)

    while gcd > 1:
        nominator = nominator / gcd
        denominator = denominator / gcd
        gcd = gcd_decimal(nominator, denominator)

    return Decimal(nominator), Decimal(denominator)


def padic(number: Union[Decimal, str], base=5):
    pass


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    padic('123.4')

    print(rationalize('123.4'))

