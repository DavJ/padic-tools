from fractions import Fraction as F
import random

import pytest

from padic_tools import (
    BallVladimirov, CyclotomicValue as CV, DomainError, KozyrevBasis,
    PadicTestFunction as TF, PrecisionError, Qp, additive_character,
    haar_volume, kozyrev_wavelet, valuation, vladimirov,
)


def pp(p,n):
    return F(p**n) if n >= 0 else F(1,p**(-n))


@pytest.mark.parametrize('p,scale,center', [(2,-1,F(1,4)),(2,2,1),
    (3,0,F(1,3)),(3,1,0),(5,-1,0),(7,0,0)])
def test_wavelet_definition_norm_mean_and_digit_contract(p,scale,center):
    K = Qp(p,12)
    for j in range(1,p):
        wave = kozyrev_wavelet(K,scale=scale,center=center,frequency=j)
        assert wave.integral() == 0
        assert wave.norm_squared() == pp(p,-scale)
        for k in range(p):
            x = center+pp(p,scale)*k
            assert wave(x) == additive_character(F(j*k,p),K)
        assert wave(center+pp(p,scale-1)) == 0
        assert kozyrev_wavelet(K,scale=scale,center=K(center,absprec=scale+1),frequency=j) == wave
        with pytest.raises(PrecisionError):
            kozyrev_wavelet(K,scale=scale,center=K(center,absprec=scale),frequency=j)


@pytest.mark.parametrize('p,a,b', [(2,-1,2),(3,0,2),(3,-1,1),(5,0,1),(2,2,2)])
def test_basis_gram_matrix_reconstruction_and_weighted_parseval(p,a,b):
    basis = KozyrevBasis(p,support_exponent=a,constancy_exponent=b)
    n = p**(b-a)
    assert len(basis.modes) == len(basis.labels) == n
    for i,left in enumerate(basis.modes):
        for j,right in enumerate(basis.modes):
            assert left.inner_product(right) == (basis.norms_squared[i] if i == j else 0)
    rng = random.Random(171+p+b)
    values = [F(rng.randint(-5,7),rng.randint(1,5))+CV.root(p,1,j) for j in range(n)]
    f = TF(p,values,support_exponent=a,constancy_exponent=b)
    coefficients = basis.decompose(f)
    assert basis.reconstruct(coefficients) == f
    assert basis.energy(coefficients) == f.norm_squared()
    assert coefficients[0] == f.integral()/haar_volume(p,a)


@pytest.mark.parametrize('p,a,b,order', [(2,0,3,1),(3,-1,1,2),(3,1,3,1),(5,0,1,3)])
def test_regional_kernel_matrix_spectrum_wavelets_and_positive_form(p,a,b,order):
    op = BallVladimirov(p,order,support_exponent=a,constancy_exponent=b)
    basis = KozyrevBasis(p,support_exponent=a,constancy_exponent=b)
    matrix = op.matrix()
    for i,row in enumerate(matrix):
        assert sum(row) == 0 and row[i] >= 0
        for j,v in enumerate(row):
            assert v == matrix[j][i]
            if i != j: assert v < 0
    for label,mode in zip(basis.labels,basis.modes):
        rate = 0 if label is None else pp(p,order*(label[0]+1))-op.escape_rate
        expected = mode*rate
        applied = op.apply(mode)
        assert applied == expected
        assert applied.values == tuple(sum(matrix[j][k]*mode.values[k]
            for k in range(len(mode.values))) for j in range(len(mode.values)))
        assert rate >= 0
    f = basis.template._replace_values([F(j*j-3*j+1,2) for j in range(len(basis.modes))])
    applied = op.apply(f)
    assert applied.integral() == 0
    weight = haar_volume(p,b)
    energy = sum(-matrix[j][k]*weight*(f.values[j]-f.values[k])**2/2
        for j in range(len(f.values)) for k in range(len(f.values)))
    assert f.inner_product(applied) == energy
    assert energy.as_fraction() >= 0


@pytest.mark.parametrize('p,a,b,order', [(2,-1,2,1),(3,0,2,2),(5,0,1,1)])
def test_whole_qp_zero_mean_derivative_matches_independent_fourier_multiplier(p,a,b,order):
    basis = KozyrevBasis(p,support_exponent=a,constancy_exponent=b)
    for label,wave in zip(basis.labels[1:],basis.modes[1:]):
        image = vladimirov(wave,order)
        assert image.tail_coefficient == 0
        assert image.as_test_function() == wave*pp(p,order*(label[0]+1))
    f = basis.reconstruct([0]+[F(j+1,j+2) for j in range(len(basis.modes)-1)])
    hat = f.fourier()
    assert hat.values[0] == 0
    multiplied = hat._replace_values([0]+[
        value*pp(p,order*(b-valuation(k,p))) for k,value in enumerate(hat.values[1:],1)])
    image = vladimirov(f,order)
    assert image.as_test_function() == multiplied.fourier(inverse=True)
    assert image.integral() == 0


@pytest.mark.parametrize('p,a,b,order', [(2,0,2,1),(3,-1,1,2),(5,1,2,1)])
def test_nonzero_mass_whole_qp_tail_and_full_integral_norm(p,a,b,order):
    f = TF(p,[j+1 for j in range(p**(b-a))],support_exponent=a,constancy_exponent=b)
    image = vladimirov(f,order)
    C = (p**order-1)/(1-pp(p,-1-order))
    for n in range(a-4,a):
        x = pp(p,n)
        # Outside B the ultrametric kernel has constant distance to every y in B.
        assert image(x) == -C*f.integral()*pp(p,(1+order)*n)
        assert image(Qp(p)(x,absprec=n+1)) == image(x)
    assert image.integral() == 0
    with pytest.raises(DomainError): image.as_test_function()
    # An independent Parseval integral of |xi|^(2*order)*|fhat|^2,
    # integrating the whole zero-frequency ball analytically.
    hat = f.fourier()
    norm = sum(value*value.conjugate()*pp(p,2*order*(b-valuation(k,p)))*pp(p,a)
               for k,value in enumerate(hat.values[1:],1))
    zero_cell = (1-F(1,p))*pp(p,a*(2*order+1))/(1-pp(p,-2*order-1))
    norm += hat.values[0]*hat.values[0].conjugate()*zero_cell
    assert image.norm_squared() == norm
    larger = f.on_lattice(support_exponent=a-1)
    larger_image = vladimirov(larger,order)
    assert larger_image.integral() == 0
    assert larger_image.norm_squared() == image.norm_squared()
    for n in range(-8,5): assert larger_image(pp(p,n)) == image(pp(p,n))


@pytest.mark.parametrize('p,a,b', [(2,-1,2),(3,0,2),(5,0,1)])
def test_diffusion_mass_positivity_energy_stationarity_and_refinement(p,a,b):
    mp = pytest.importorskip('mpmath')
    ctx = mp.mp.clone(); ctx.dps = 70
    before = mp.mp.dps
    op = BallVladimirov(p,support_exponent=a,constancy_exponent=b)
    f = op.template._replace_values([p**b if b >= 0 else pp(p,b)]+[0]*(len(op.template.values)-1))
    energy = f.norm_squared().complex(70).real
    assert f.integral() == 1
    for time in [0,F(1,100),F(1,10),1,10]:
        state = op.evolve(f,time,dps=65)
        assert abs(state.integral()-1) < ctx.mpf('1e-60')
        assert min(value.real for value in state.values) > -ctx.mpf('1e-60')
        assert max(abs(value.imag) for value in state.values) < ctx.mpf('1e-60')
        assert state.norm_squared() <= energy+ctx.mpf('1e-60')
        energy = state.norm_squared()
        assert state(pp(p,a-1)) == 0
    # Matrix exponential is an independent algorithm for the Markov solution.
    matrix = ctx.matrix([[ctx.mpf(v.numerator)/v.denominator for v in row] for row in op.matrix()])
    oracle = ctx.expm(-matrix/10)*ctx.matrix([ctx.mpf(value.as_fraction().numerator)/value.as_fraction().denominator
                                            for value in f.values])
    state = op.evolve(f,F(1,10),dps=65)
    assert max(abs(value-oracle[j]) for j,value in enumerate(state.values)) < ctx.mpf('1e-60')
    fine = f.on_lattice(constancy_exponent=b+1)
    fine_op = BallVladimirov(p,support_exponent=a,constancy_exponent=b+1)
    fine_state = fine_op.evolve(fine,F(1,10),dps=65)
    assert max(abs(fine_state.values[j]-state.values[j % len(f.values)])
               for j in range(len(fine.values))) < ctx.mpf('1e-60')
    constant = op.template._replace_values([1]*len(f.values))
    assert op.apply(constant) == constant*0
    assert all(abs(v-1) < ctx.mpf('1e-60') for v in op.evolve(constant,100,dps=65).values)
    assert mp.mp.dps == before


def test_lattice_refinement_and_explicit_domains():
    f = TF(3,[1,2,4])
    enlarged = f.on_lattice(support_exponent=-1,constancy_exponent=2)
    assert f.integral() == enlarged.integral()
    assert f.norm_squared() == enlarged.norm_squared()
    for x in [F(1,3),0,1,2,3,4]: assert f(x) == enlarged(x)
    with pytest.raises(DomainError): f.on_lattice(constancy_exponent=0)
    with pytest.raises(PrecisionError): f.on_lattice(constancy_exponent=8,max_size=10)
    with pytest.raises(TypeError): f.on_lattice(support_exponent=True)
    with pytest.raises(DomainError): kozyrev_wavelet(3,scale=1,constancy_exponent=1)
    with pytest.raises(ValueError): kozyrev_wavelet(3,frequency=0)
    with pytest.raises(ValueError): kozyrev_wavelet(3,frequency=True)
    with pytest.raises(PrecisionError): KozyrevBasis(3,constancy_exponent=4,max_size=10)
    basis = KozyrevBasis(3)
    with pytest.raises(ValueError): basis.reconstruct([1,2])
    with pytest.raises(DomainError): basis.decompose(TF.ball(3))
    for order in [0,-1,True,F(1,2)]:
        with pytest.raises(ValueError): BallVladimirov(3,order)
    with pytest.raises(PrecisionError): BallVladimirov(3,constancy_exponent=5,max_size=20)
    with pytest.raises(TypeError): vladimirov(1)
    op = BallVladimirov(3)
    with pytest.raises(DomainError): op.evolve(f,-1)
    with pytest.raises(TypeError): op.evolve(f,0.1)
    with pytest.raises(ValueError): op.evolve(f,1,dps=True)


def test_reproducible_diffusion_example():
    mp = pytest.importorskip('mpmath')
    import json
    from pathlib import Path
    import subprocess
    import sys
    script = Path(__file__).resolve().parents[1]/'examples/padic_diffusion.py'
    result = subprocess.run([sys.executable,str(script)],capture_output=True,text=True,
                            check=True,timeout=60)
    report = json.loads(result.stdout)
    assert mp.mpf(report['refinement_error_at_time_1_over_10']) < mp.mpf('1e-50')
    for grid in report['grids']:
        assert grid['initial_mass'] == '1' and grid['exact_zero_row_sums']
        previous = mp.inf
        for sample in grid['samples']:
            assert mp.mpf(sample['mass_error']) < mp.mpf('1e-50')
            assert mp.mpf(sample['minimum_density']) > -mp.mpf('1e-50')
            energy = mp.mpf(sample['norm_squared'])
            assert energy <= previous
            previous = energy
