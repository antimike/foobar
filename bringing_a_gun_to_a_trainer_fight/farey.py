def farey(N: int):
    """
    Generate the Nth Farey sequence.
    """
    a, b, c, d = 0, 1, 1, N
    yield (0, 1)
    while c <= N:
        yield (c, d)
        k = (N + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b


def solution(dims, pos, tpos, dist):
    a, b, c, d = 0, 1, 1, dist - 1
    while c < dist:

        k = (dist + b - 1) // d
        a, b, c, d = c, d, k * c - a, k * d - b


def coprime_integers(N):
    """
    Generate all pairs of coprime integers (a, b) with a < b <= N.

    Approach 1:
    gcd(a, b) == 1 iff there are integers r, s such that a*r + b*s == 1.
    Note that this implies that r, the inverse of a mod b, is also coprime to b; thus these numbers can be computed in pairs.

    Approach 2:
    1. Precompute all primes in the range [1, N].
        * For this problem, N <= 10**4, so numPrimes ~ 250
    2. Expand the list of primes to a list of prime powers in the range [1, N].
        * This might be too much to ask.  Perhaps it's enough to just store the maximum exponent for each prime...
    3 ...

    Approach 3:
    Use the Farey Sequence!
    --> [(a, b) for (a, b) in farey(N) if a**2 + b**2 <= D2]

    Approach 4:
    Algorithm to generate all coprime (a, b) as nodes in two distinct ternary trees rooted at (1, 2) and (1, 3)
    """
    ...
