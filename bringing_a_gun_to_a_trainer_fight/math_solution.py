"""
After failing with brute-force and Farey-star-based solutions, I decided to finally sit down and solve the problem "the old-fashioned way."

Fortunately, it's not actually too difficult, and the solution cuts the runtime from O(dist**2) to O(dist).
"""

from math import ceil, floor, sqrt, gcd
from .vector import IntegerVector, BilinearForm
from itertools import combinations
from fractions import Fraction


def lcm(a, b):
    """Compute the LCM of two integers"""
    return a * b // gcd(a, b)


def numLatticePointsInEllipse(basis, anchor, center, axes):
    Lx, Ly = basis
    lx, ly = anchor
    x0, y0 = center
    a, b = axes
    total = 0
    for m in range(ceil(x0 - a), floor(x0 + a) + 1):
        if (m - lx) % Lx != 0:
            continue
        s = sqrt(1 - Fraction(m - x0, a)**2)
        total += (floor(y0 + b * s) - ly) // Ly - (ceil(y0 - b * s) -
                                                   ly) // Ly + 1
    return total


def helper(base, L, H, dist):
    """Compute the number of integer lattice points satisfying abs(base + [[2*L, 0], [0, 2*H]]*[m, n]) <= dist"""
    print(
        f"Calculating lattice points in critical ellipse for given parameters:"
    )
    print(f"\tbase = {base}\n\tL, H = {L, H}\n\tdist = {dist}")

    # (delta_x, delta_y) := center of critical ellipse
    # k0 := min(n in integers such that n*delta[i] in integers for i = 0, 1)
    delta_x, delta_y = Fraction(base[0], 2 * L), Fraction(base[1], 2 * H)
    k0 = lcm(delta_x.denominator, delta_y.denominator)

    # Axes of ellipse
    a, b = Fraction(dist, 2 * L), Fraction(dist, 2 * H)

    print(
        f"\tCalculated parameters:\n\t\tdelta = {(delta_x, delta_y)}\n\t\t(a, b) = {(a, b)}\n\t\tk0 = {k0}"
    )

    total = 0

    # Loop over all possible values of m
    for m in range(ceil(-delta_x - a), floor(-delta_x + a) + 1):
        s = sqrt(1 - Fraction(m + delta_x, a)**2)
        c, d = ceil(-delta_y - b * s), floor(-delta_y + b * s)
        total += d - c + 1
        print(
            f"\tLattice points with m = {m}:\n\t\trange = {(c, d)}\n\t\tnum = {d - c + 1}"
        )
        print(
            f"\t\tBearings: {[(base[0] + 2*L*m, base[1] + 2*H*n) for n in range(c, d + 1)]}"
        )

    return total


def find_num_target_points(pos, tpos, L, H, dist):
    total = 0
    for target in tpos, (tpos[0], -tpos[1]), (-tpos[0], tpos[1]), (-tpos[0],
                                                                   -tpos[1]):
        base = (target[0] - pos[0], target[1] - pos[1])
        total += helper(base, L, H, dist)
    return total


def solution(dims, pos, tpos, dist):
    L, H = dims
    me = find_num_target_points(pos, pos, L, H, dist)
    him = find_num_target_points(pos, tpos, L, H, dist)
    print("", f"Number of target points:\n\tme --> {me}\n\thim --> {him}")
    return him - me
