"""
Key insight: symmetry

"Self-owns" are the same as "backward redundant" lattice points
"""

from fractions import Fraction
from math import sqrt, gcd, floor, ceil
from itertools import chain


def get_bearing(p):
    return p // gcd(p.x, p.y)


def solution(dims, pos, tpos, dist):
    L = ReflectionLattice(*dims)
    me, trainer = Point(*pos), Point(*tpos)
    bearings = {get_bearing(trainer - me)}
    hits = 1
    for image, src in sorted(chain(
        map(lambda q: (q, trainer), L.reflections(trainer, me, dist)),
        map(lambda q: (q, me), L.reflections(me, me, dist))),
            key=lambda x: (x[0] - me) * (x[0] - me)):
        b = get_bearing(image - me)
        if b not in bearings:
            bearings.add(b)
            if src != me:
                hits += 1
    return hits


class Point:
    """
    Adds arithmetic structure to pairs (a, b) via dunder methods
    """

    def __init__(self, x, y):
        self.x, self.y = x, y

    def __neg__(self):
        return self.__class__(-self.x, -self.y)

    def __add__(self, other):
        return self.__class__(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return self.__class__(self.x - other.x, self.y - other.y)

    def __rmul__(self, num):
        return self.__class__(num * self.x, num * self.y)

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            return self.x * other.x + self.y * other.y
        else:
            return self.__class__(other * self.x, other * self.y)

    def __truediv__(self, other):
        if isinstance(other, self.__class__):
            return NotImplemented
        else:
            return self.__class__(self.x / other, self.y / other)

    def __floordiv__(self, other):
        if isinstance(other, self.__class__):
            return NotImplemented
        else:
            return self.__class__(self.x // other, self.y // other)

    def __mod__(self, other):
        if isinstance(other, self.__class__):
            return self.__class__(self.x % other.x, self.y % other.y)
        else:
            return self.__class__(self.x % other, self.y % other)

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y)

    def __repr__(self):
        return "<" + ",".join([repr(self.x), repr(self.y)]) + ">"

    @classmethod
    @property
    def origin(cls):
        return cls(0, 0)


class ReflectionLattice:
    """
    The lattice of reflections in the sides of a rectangle aligned with the coordinate axes.
    """

    def __init__(self, lx, ly):
        self._dims = Point(lx, ly)
        self._basis = 2 * Point(Point(lx, 0), Point(0, ly))
        self._inverse_basis = Point(Point(Fraction(1, 2 * lx), 0),
                                    Point(0, Fraction(1, 2 * ly)))
        self._eccentricity = Fraction(lx, ly)

    def primary_reflections(self, p):
        """Generates the four reflections of p through the sides of the lattice's defining rectangle."""
        # yield Point(p.x, 2 * self._dims.y - p.y)
        yield Point(p.x, -p.y)
        yield Point(-p.x, p.y)
        yield Point(-p.x, -p.y)
        # yield Point(2 * self._dims.x - p.x, p.y)

    def reflections(self, source, observer, maxDist):
        """Generates all image-points of source within a given distance of observer."""
        for q in self.primary_reflections(source):
            for z in map(lambda p: p + q, self.in_disk(observer - q, maxDist)):
                yield z

    def in_disk(self, center, radius):
        xmin, xmax = center.x - radius, center.x + radius
        xmin += (-xmin) % (2 * self._dims.x)
        xmax -= xmax % (2 * self._dims.x)
        for x in range(xmin, xmax + 1, 2 * self._dims.x):
            dy = floor(sqrt(radius**2 - (x - center.x)**2))
            ymin, ymax = center.y - dy, center.y + dy
            ymin += (-ymin) % (2 * self._dims.y)
            ymax -= ymax % (2 * self._dims.y)
            for y in range(ymin, ymax + 1, 2 * self._dims.y):
                yield Point(x, y)

    def numVisibleInDisk(self, center, radius):
        z0 = center * self._inverse_basis
        d1, d2 = z0.x.denominator, z0.y.denominator
        lcm = d1 * d2 // gcd(d1, d2)

        numVisible = self.numInDisk(center, radius)
        numMasked = 0

        sieve = [0] * (radius + 1)
        for n in range(1, radius // lcm):
            k = 1 + n * lcm
            if sieve[k - 1] != 1:
                redundant = self.numInDisk(center, radius // k)
                multiplicity = 1 - sieve[k - 1]
                correction = multiplicity * redundant
                numVisible -= correction
                numMasked += correction
                for j in range(k - 1, len(sieve), k):
                    sieve[j] += multiplicity
        return numVisible, numMasked

    def numInDisk(self, center, radius):
        """
        Computes the number of integer lattice points lying within the given radius of the point center.

        Note that the "anchor point" of the lattice (i.e., its "origin") can be taken to be the same as center without loss of generality.
        This is not true for the equivalent problem on visilibity classes, however.
        """
        # center of ellipse
        ecenter = center * self._inverse_basis

        # axes of ellipse
        axes = radius * self._inverse_basis

        mmin, mmax = ceil(ecenter.x - axes.x.x), floor(ecenter.x + axes.x.x)
        print(f"Center of ellipse: {ecenter}\nAxes: {axes}")
        print(
            f"Checking x in range {2*self._dims.x*mmin + center.x, 2*self._dims.x*mmax + center.x}..."
        )

        total = 0
        for m in range(mmin, mmax + 1):
            s = (axes.y.y / axes.x.x) * sqrt(axes.x.x**2 - (m - ecenter.x)**2)
            nmin, nmax = ceil(ecenter.y - s), floor(ecenter.y + s)
            print(
                f"Num points found with x = {2*self._dims.x*m + center.x}: {nmax - nmin + 1}",
            )
            print(f"\tyrange = {nmin, nmax}")
            total += nmax - nmin + 1

        return total

    def _find_linear_combo(self, u, v):
        """Finds the minimal k and modulus n such that u + k'*v is in self for all k' == k (mod n).

        Raises a ValueError if this is not possible."""
        try:
            kx, nx = self.solveScalarCongruence(u.x, v.x, self._basis.x.x)
            ky, ny = self.solveScalarCongruence(u.y, v.y, self._basis.y.y)
            d, r, s = self.euclideanAlgorithm(nx, ny)

            if kx % d != ky % d:
                raise ValueError("2D congruence has no solution")

            Q = (ky - kx) // d
            assert Q * r * nx + kx == -Q * s * ny + \
                ky, "Inconsistent results obtained for vector congruence"
            lcm = (nx * ny // d)
            return (Q * r * nx + kx) % lcm, lcm
        except ValueError:
            raise

    @staticmethod
    def solveScalarCongruence(x, y, z):
        """Find the minimal k that solves z | (x + k*y), along with the modulus n such that k' == k (mod n) guarantees that x + k*y is a solution.

        Raise a ValueError if this is not possible."""
        d, a, b = ReflectionLattice.euclideanAlgorithm(y, z)
        if x % d != 0:
            raise ValueError("Congruence has no solution")
        x //= d
        y //= d
        z //= d
        yinv = a % z
        return (-a * x) % z, z

    @staticmethod
    def euclideanAlgorithm(a, b):
        """Computes r and s such that r*a + s*b = gcd(a, b), returning (gcd(a, b), r, s)"""
        if a == b:
            return a, 1, 0
        if a > b:
            d, r, s = ReflectionLattice.euclideanAlgorithm(b, a)
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


if __name__ == "__main__":
    t1 = [300, 275], [150, 150], [185, 100], 500
    # [1, 0], [1, 2], [1, -2], [3, 2], [3, -2], [-3, 2], [-3, -2]
    t2 = [3, 2], [1, 1], [2, 1], 4
    print(solution(*t1))
    print(solution(*t2))
