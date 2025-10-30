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


def padic_valuation(n: int, p: int) -> int:
    """
    Compute the p-adic valuation of n.
    Returns the largest power of p that divides n.
    """
    if n == 0:
        return float('inf')
    
    valuation = 0
    n = abs(n)
    while n % p == 0:
        n //= p
        valuation += 1
    return valuation


def padic_norm(n: int, p: int) -> float:
    """
    Compute the p-adic norm of n.
    Returns p^(-v_p(n)) where v_p(n) is the p-adic valuation.
    """
    if n == 0:
        return 0.0
    
    val = padic_valuation(n, p)
    return p ** (-val)


def padic_expansion(numerator: int, denominator: int, prime: int = 5, terms: int = 10) -> list:
    """
    Compute the p-adic expansion of numerator/denominator.
    Returns a list of coefficients [a_0, a_1, a_2, ...] where
    numerator/denominator = a_0 + a_1*p + a_2*p^2 + ...
    """
    # Normalize the fraction
    gcd = math.gcd(abs(numerator), abs(denominator))
    numerator //= gcd
    denominator //= gcd
    
    # Check if denominator is coprime to p
    if math.gcd(denominator, prime) != 1:
        raise ValueError(f"Denominator must be coprime to prime {prime}")
    
    # Find the inverse of denominator mod p^terms using Hensel's lemma
    # Start with inverse mod p and lift it
    s, t, g = bezout2(denominator, prime)
    if g != 1:
        raise ValueError(f"Denominator not invertible mod {prime}")
    
    inv = s % prime
    
    # Hensel lifting: lift inverse from mod p to mod p^terms
    for k in range(1, terms):
        # inv_new = inv * (2 - denominator * inv) mod p^(k+1)
        # But we use a simpler method: inv_new = inv + p^k * t where
        # t is chosen so that denominator * inv_new ≡ 1 (mod p^(k+1))
        mod_pk = prime ** k
        mod_pk1 = prime ** (k + 1)
        
        # Find correction term
        residual = (1 - denominator * inv) // mod_pk
        s_k, t_k, g_k = bezout2(denominator, prime)
        correction = (s_k * residual) % prime
        inv = (inv + correction * mod_pk) % mod_pk1
    
    # Now compute the expansion
    mod = prime ** terms
    inv_full = inv
    
    # Recompute inverse for full precision
    s, t, g = bezout2(denominator, mod)
    if g != 1:
        raise ValueError(f"Denominator not invertible mod {prime}^{terms}")
    inv_full = s % mod
    
    # Compute numerator * inverse mod p^terms
    result = (numerator * inv_full) % mod
    
    # Extract coefficients
    coefficients = []
    for k in range(terms):
        coeff = result % prime
        coefficients.append(coeff)
        result = result // prime
    
    return coefficients


def padic(nominator: Union[Decimal, int], denominator: Union[Decimal, int], prime: Union[int, Decimal] = 5, terms: int = 10) -> str:
    """
    Compute p-adic expansion and return as a formatted string.
    """
    try:
        expansion = padic_expansion(int(nominator), int(denominator), int(prime), terms)
        
        # Format the output
        result_parts = []
        for i, coeff in enumerate(expansion):
            if coeff != 0:
                if i == 0:
                    result_parts.append(str(coeff))
                elif i == 1:
                    result_parts.append(f"{coeff}*{prime}")
                else:
                    result_parts.append(f"{coeff}*{prime}^{i}")
        
        if not result_parts:
            return "0"
        
        return " + ".join(result_parts)
    except ValueError as e:
        return f"Error: {str(e)}"


class PadicNumber:
    """
    Represents a p-adic number with its expansion coefficients.
    """
    def __init__(self, numerator: int, denominator: int, prime: int = 5, terms: int = 10):
        """
        Initialize a p-adic number from a rational number.
        """
        self.prime = prime
        self.terms = terms
        
        # Normalize the fraction
        gcd = math.gcd(abs(numerator), abs(denominator))
        self.numerator = numerator // gcd
        self.denominator = denominator // gcd
        
        # Compute the expansion
        self.coefficients = padic_expansion(self.numerator, self.denominator, self.prime, self.terms)
    
    def __str__(self):
        """String representation of the p-adic number."""
        result_parts = []
        for i, coeff in enumerate(self.coefficients):
            if coeff != 0:
                if i == 0:
                    result_parts.append(str(coeff))
                elif i == 1:
                    result_parts.append(f"{coeff}*{self.prime}")
                else:
                    result_parts.append(f"{coeff}*{self.prime}^{i}")
        
        if not result_parts:
            return "0"
        
        return " + ".join(result_parts)
    
    def __repr__(self):
        return f"PadicNumber({self.numerator}/{self.denominator}, p={self.prime})"
    
    def valuation(self):
        """Return the p-adic valuation."""
        return padic_valuation(self.numerator, self.prime) - padic_valuation(self.denominator, self.prime)
    
    def norm(self):
        """Return the p-adic norm."""
        val = self.valuation()
        if val == float('inf'):
            return 0.0
        return self.prime ** (-val)
    
    def __add__(self, other):
        """Add two p-adic numbers."""
        if not isinstance(other, PadicNumber):
            raise TypeError("Can only add PadicNumber to PadicNumber")
        if self.prime != other.prime:
            raise ValueError("Can only add p-adic numbers with the same prime")
        
        # Add as fractions
        num = self.numerator * other.denominator + other.numerator * self.denominator
        den = self.denominator * other.denominator
        
        return PadicNumber(num, den, self.prime, max(self.terms, other.terms))
    
    def __sub__(self, other):
        """Subtract two p-adic numbers."""
        if not isinstance(other, PadicNumber):
            raise TypeError("Can only subtract PadicNumber from PadicNumber")
        if self.prime != other.prime:
            raise ValueError("Can only subtract p-adic numbers with the same prime")
        
        # Subtract as fractions
        num = self.numerator * other.denominator - other.numerator * self.denominator
        den = self.denominator * other.denominator
        
        return PadicNumber(num, den, self.prime, max(self.terms, other.terms))
    
    def __mul__(self, other):
        """Multiply two p-adic numbers."""
        if not isinstance(other, PadicNumber):
            raise TypeError("Can only multiply PadicNumber with PadicNumber")
        if self.prime != other.prime:
            raise ValueError("Can only multiply p-adic numbers with the same prime")
        
        # Multiply as fractions
        num = self.numerator * other.numerator
        den = self.denominator * other.denominator
        
        return PadicNumber(num, den, self.prime, max(self.terms, other.terms))
    
    def __truediv__(self, other):
        """Divide two p-adic numbers."""
        if not isinstance(other, PadicNumber):
            raise TypeError("Can only divide PadicNumber by PadicNumber")
        if self.prime != other.prime:
            raise ValueError("Can only divide p-adic numbers with the same prime")
        if other.numerator == 0:
            raise ZeroDivisionError("Cannot divide by zero")
        
        # Divide as fractions
        num = self.numerator * other.denominator
        den = self.denominator * other.numerator
        
        return PadicNumber(num, den, self.prime, max(self.terms, other.terms))
    
    def __neg__(self):
        """Negate a p-adic number."""
        return PadicNumber(-self.numerator, self.denominator, self.prime, self.terms)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # Example usage
    print("=== P-adic Arithmetic Tools ===\n")
    
    # Test p-adic valuation
    print("P-adic valuation examples:")
    print(f"v_5(25) = {padic_valuation(25, 5)}")  # Should be 2
    print(f"v_5(15) = {padic_valuation(15, 5)}")  # Should be 1
    print(f"v_5(7) = {padic_valuation(7, 5)}")    # Should be 0
    print()
    
    # Test p-adic norm
    print("P-adic norm examples:")
    print(f"||25||_5 = {padic_norm(25, 5)}")  # Should be 1/25
    print(f"||15||_5 = {padic_norm(15, 5)}")  # Should be 1/5
    print(f"||7||_5 = {padic_norm(7, 5)}")    # Should be 1
    print()
    
    # Test p-adic expansion
    print("P-adic expansion examples:")
    print(f"1/3 in 5-adic: {padic(1, 3, 5, 10)}")
    print(f"1/2 in 5-adic: {padic(1, 2, 5, 10)}")
    print()
    
    # Test PadicNumber class
    print("PadicNumber class examples:")
    p1 = PadicNumber(1, 3, 5, 10)
    print(f"p1 = {p1}")
    print(f"repr: {repr(p1)}")
    print(f"valuation: {p1.valuation()}")
    print(f"norm: {p1.norm()}")
    print()
    
    p2 = PadicNumber(1, 2, 5, 10)
    print(f"p2 = {p2}")
    print(f"repr: {repr(p2)}")
    print()
    
    # Test arithmetic operations
    print("Arithmetic operations:")
    p3 = p1 + p2  # 1/3 + 1/2 = 5/6
    print(f"p1 + p2 = {repr(p3)}")
    print(f"  = {p3}")
    print()
    
    p4 = p1 - p2  # 1/3 - 1/2 = -1/6
    print(f"p1 - p2 = {repr(p4)}")
    print(f"  = {p4}")
    print()
    
    p5 = p1 * p2  # 1/3 * 1/2 = 1/6
    print(f"p1 * p2 = {repr(p5)}")
    print(f"  = {p5}")
    print()
    
    p6 = p1 / p2  # 1/3 / 1/2 = 2/3
    print(f"p1 / p2 = {repr(p6)}")
    print(f"  = {p6}")
    print()
    
    p7 = -p1  # -1/3
    print(f"-p1 = {repr(p7)}")
    print(f"  = {p7}")
    print()