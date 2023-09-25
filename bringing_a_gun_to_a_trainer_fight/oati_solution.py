from fractions import Fraction
from math import sqrt, floor, ceil
from collections import Counter


def solveLatticeIntersectionProblem(L, A1, A2, base, D):
    lx, ly = L
    a1x, a1y = A1
    a2x, a2y = A2
    bx, by = base

    d1x, d1y = a1x - b1x, a1y - b1y
    d2x, d2y = a2x - b2x, a2y - b2y

    try:
        k, n = solveVectorCongruence((d1x, d1y), (-d2x, -d21y), L)
        ct = Counter()
        for m in range(k, D, n):
            zx = (d1x - m * d2x) // lx
            zy = (d1y - m * d2y) // ly

    except ValueError:
        print("Lattices are incompatible")
        raise


def solveLatticeProblem(base, anchor, diag_basis, dist):
    Bx, By = -base[0] + anchor[0], -base[1] + anchor[1]
    lx, ly = diag_basis
    e = Fraction(lx, ly)
    qx, qy = Fraction(Bx, lx), Fraction(By, ly)
    a, b = Fraction(dist, lx), Fraction(dist, ly)

    mmin, mmax = ceil(qx - a), floor(qx + a)

    total = 0
    for m in range(mmin, mmax + 1):
        s = e * sqrt(a**2 - (m - qx)**2)
        nmin, nmax = ceil(qy - b * s), floor(qy + b * s)
        total += nmax - nmin + 1

    return total


def solveVisibleLatticeProblem(base, anchor, diag_basis, dist):
    total = solveLatticeProblem(base, anchor, diag_basis, dist)
    sieve = [0 for _ in range(dist - 1)]
    for k in range(2, dist + 1):
        if sieve[k - 2] == 1:
            continue
        d_anchor = ((k - 1) * diag_basis[0] * base[0],
                    (k - 1) * diag_basis[1] * base[1])
        redundant = solveLatticeProblem(
            base, (anchor[0] + d_anchor[0], anchor[1] + d_anchor[1]),
            (k * diag_basis[0], k * diag_basis[1]), dist)
        multiplicity = 1 - sieve[k - 2]
        total -= multiplicity * redundant
        for j in range(k - 2, len(sieve), k):
            sieve[j] += multiplicity
    return total


def solution(dims, pos, tpos, dist):
    basis = (2 * dims[0], 2 * dims[1])
    numTargetBearings = 0
    numSelfOwns = 0
    for target in tpos, (-tpos[0], tpos[1]), (-tpos[0], -tpos[1]), (tpos[0],
                                                                    -tpos[1]):
        print("----------------------------------------")
        print(f"Finding target bearings for target point {target}...")
        target_bearings = solveVisibleLatticeProblem(pos, target, basis, dist)
        numTargetBearings += target_bearings
        print(f"\t{target_bearings} bearings found")
        target_disp = (target[0] - pos[0], target[1] - pos[1])

        for image in pos, (-pos[0], pos[1]), (-pos[0], -pos[1]), (pos[0],
                                                                  -pos[1]):
            print(f"Searching for self-owns with image point {image}...")
            image_disp = (image[0] - pos[0], image[1] - pos[1])
            try:
                k, n = solveVectorCongruence(target_disp,
                                             (-image_disp[0], -image_disp[1]),
                                             basis)
                self_owns = solveVisibleLatticeProblem(pos, image, basis,
                                                       Fraction(dist, k))
                print(f"\tSelf-owns found for image {image}: {self_owns}")
                numSelfOwns += self_owns
            except ValueError:
                print(
                    f"\tLattice congruence problem is invalid.  No self-owns")
        print("----------------------------------------")
        print()
    net = numTargetBearings - numSelfOwns
    print(
        f"numTargetBearings: {numTargetBearings}\nnumSelfOwns: {numSelfOwns}\nnet: {net}"
    )
    return net


def solveVectorCongruence(u, v, diag_basis):
    """Find the minimal k and modulus n such that diag_basis.xhat | (u + k'*v).xhat and diag_basis.yhat | (u + k'*v).yhat for all k' == k (mod n).

    Raise a ValueError if this is not possible."""
    ux, uy = u
    vx, vy = v
    lx, ly = diag_basis

    try:
        kx, nx = solveScalarCongruence(ux, vx, lx)
        ky, ny = solveScalarCongruence(uy, vy, ly)
        d, r, s = euclideanAlgorithm(nx, ny)

        if kx % d != ky % d:
            raise ValueError("2D congruence has no solution")

        Q = (ky - kx) // d
        assert Q * r * nx + kx == -Q * s * ny + \
            ky, "Inconsistent results obtained for vector congruence"
        lcm = (nx * ny // d)
        return (Q * r * nx + kx) % lcm, lcm
    except ValueError:
        raise


def solveScalarCongruence(x, y, z):
    """Find the minimal k that solves z | (x + k*y), along with the modulus n such that k' == k (mod n) guarantees that x + k*y is a solution.

    Raise a ValueError if this is not possible."""
    d, a, b = euclideanAlgorithm(y, z)
    if x % d != 0:
        raise ValueError("Congruence has no solution")
    x //= d
    y //= d
    z //= d
    yinv = a % z
    return (-a * x) % z, z


def euclideanAlgorithm(a, b):
    """Computes r and s such that r*a + s*b = gcd(a, b), returning (gcd(a, b), r, s)"""
    if a == b:
        return a, 1, 0
    if a > b:
        d, r, s = euclideanAlgorithm(b, a)
        return d, s, r
    ca, cb = [1, 0], [0, 1]
    while a % b != 0:
        d = a // b
        a -= d * b
        ca[0] -= d * cb[0]
        ca[1] -= d * cb[1]
        if a < b:
            a, b, ca, cb = b, a, cb, ca
    return b, cb[0], cb[1]
