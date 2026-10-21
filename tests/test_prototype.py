"""Verify the source-tree compatibility helpers against exact integer oracles."""
import importlib.util
import math
from pathlib import Path


def test_prototype_helpers():
    spec=importlib.util.spec_from_file_location("prototype",Path(__file__).resolve().parents[1]/"main.py")
    prototype=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prototype)
    for a in range(-9,10):
        for b in range(-9,10):
            s,t,g=prototype.bezout2(a,b)
            assert s*a+t*b==g==math.gcd(a,b)
            assert prototype.bezout1(a,b)==(s,t)
    assert prototype.rationalize("0.125")== (1,8)
    assert prototype.padic(1,3).contains("1/3")
