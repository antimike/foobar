"""
Key insight: symmetry

"Self-owns" are the same as "backward redundant" lattice points
"""

from fractions import Fraction
from math import sqrt, gcd, floor
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


if __name__ == "__main__":
    t1 = [300, 275], [150, 150], [185, 100], 500
    # [1, 0], [1, 2], [1, -2], [3, 2], [3, -2], [-3, 2], [-3, -2]
    t2 = [3, 2], [1, 1], [2, 1], 4
    print(solution(*t1))
    print(solution(*t2))
