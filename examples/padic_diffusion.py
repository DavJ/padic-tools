"""Reproduce exact spectra and numerical conservative diffusion checks.

Run python examples/padic_diffusion.py after installing .[complex].
All three lattices contain the same initial density 9*1_(9*Z3).
The b=2 and b=3 grids share the same ball; a=-1 changes the boundary.
"""
from fractions import Fraction
import json

from padic_tools import BallVladimirov, PadicTestFunction, Qp


def run_demo():
    import mpmath
    ctx = mpmath.mp.clone(); ctx.dps = 70
    K = Qp(3)
    initial = PadicTestFunction.ball(K,exponent=2)*9
    times = [Fraction(0),Fraction(1,100),Fraction(1,10),Fraction(1),Fraction(5)]
    grids,states = [],{}
    for a,b in [(0,2),(0,3),(-1,2)]:
        f = initial.on_lattice(support_exponent=a,constancy_exponent=b)
        operator = BallVladimirov(K,support_exponent=a,constancy_exponent=b)
        rows = []
        for t in times:
            state = operator.evolve(f,t,dps=60)
            states[a,b,t] = state
            rows.append({
                'time': str(t), 'mass_error': ctx.nstr(abs(state.integral()-1),12),
                'minimum_density': ctx.nstr(min(v.real for v in state.values),20),
                'density_at_zero': ctx.nstr(state(0).real,20),
                'norm_squared': ctx.nstr(state.norm_squared(),20),
                'maximum_imaginary_part': ctx.nstr(max(abs(v.imag) for v in state.values),12),
            })
        grids.append({'support_exponent':a,'constancy_exponent':b,
            'cosets':len(f.values),'initial_mass':str(f.integral().as_fraction()),
            'escape_rate':str(operator.escape_rate),
            'spectrum':sorted(set(str(v) for v in operator.eigenvalues)),
            'exact_zero_row_sums':all(sum(row) == 0 for row in operator.matrix()),'samples':rows})
    coarse,fine = states[0,2,Fraction(1,10)],states[0,3,Fraction(1,10)]
    error = max(abs(v-coarse.values[j % 9]) for j,v in enumerate(fine.values))
    return {'prime':3,'order':1,'working_dps':60,
        'initial_density':'9 * indicator(9*Z3)',
        'boundary':'regional conservative operator on p**a * Zp',
        'refinement_error_at_time_1_over_10':ctx.nstr(error,12),'grids':grids}


if __name__ == '__main__':
    print(json.dumps(run_demo(),indent=2))
