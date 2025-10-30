# padic-tools
some tools for padic number system calculations

## Features

This library provides a comprehensive set of tools for working with p-adic numbers:

### Core Functions

- **`padic_valuation(n, p)`**: Compute the p-adic valuation of an integer n (the largest power of p that divides n)
- **`padic_norm(n, p)`**: Compute the p-adic norm of an integer n
- **`padic_expansion(numerator, denominator, prime, terms)`**: Compute the p-adic expansion of a rational number
- **`padic(numerator, denominator, prime, terms)`**: Format and display a p-adic expansion as a string

### PadicNumber Class

The `PadicNumber` class represents p-adic numbers and supports:
- Initialization from rational numbers (numerator/denominator)
- Pretty printing of p-adic expansions
- Arithmetic operations: addition (+), subtraction (-), multiplication (*), division (/), negation (-)
- Valuation and norm calculations

## Usage Examples

```python
from main import *

# Basic p-adic functions
print(padic_valuation(25, 5))  # Returns 2 (since 25 = 5^2)
print(padic_norm(25, 5))       # Returns 0.04 (= 5^-2)

# P-adic expansions
print(padic(1, 3, 5, 10))  # 1/3 in 5-adic representation

# Using the PadicNumber class
p1 = PadicNumber(1, 3, 5, 10)  # 1/3 in 5-adic
p2 = PadicNumber(1, 2, 5, 10)  # 1/2 in 5-adic

# Arithmetic operations
p3 = p1 + p2  # 1/3 + 1/2 = 5/6
p4 = p1 * p2  # 1/3 * 1/2 = 1/6
p5 = p1 / p2  # 1/3 / 1/2 = 2/3
p6 = -p1      # -1/3

print(p3)  # Displays the 5-adic expansion
```

## Mathematical Background

P-adic numbers are an alternative way to complete the rational numbers, different from the real numbers. They provide a notion of "closeness" based on divisibility by a prime p rather than the usual distance. This makes them useful in number theory and various areas of mathematics.

## Running the Examples

Run the included examples:
```bash
python3 main.py
```

