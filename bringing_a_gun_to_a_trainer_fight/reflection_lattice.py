from math import sqrt, floor
from itertools import chain


def gcd(a, b):
    # Python 2.7 math module doesn't implement this :/
    while b != 0:
        a, b = b, a % b
    return abs(a)


def get_bearing(p):
    if p == Point(0, 0):
        return p
    return p // gcd(p.x, p.y)


def solution(dims, pos, tpos, dist):
    L = ReflectionLattice(*dims)
    me, trainer = Point(*pos), Point(*tpos)
    bearings = set()
    hits = set()
    for image, src in sorted(chain(
        map(lambda q: (q, trainer), L.reflections(trainer, me, dist)),
        map(lambda q: (q, me), L.reflections(me, me, dist))),
            key=lambda x: (x[0] - me) * (x[0] - me)):
        b = get_bearing(image - me)
        if b not in bearings:
            bearings.add(b)
            if src != me:
                hits.add(b)
    return len(hits)


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


class ReflectionLattice:
    """
    The lattice of reflections in the sides of a rectangle aligned with the coordinate axes.
    """

    def __init__(self, lx, ly):
        self._dims = Point(lx, ly)

    def primary_reflections(self, p):
        """Generates representative members of the four equivalence classes of images of p under reflections through the sides of the rectangle."""
        yield p
        yield Point(p.x, -p.y)
        yield Point(-p.x, p.y)
        yield Point(-p.x, -p.y)

    def reflections(self, source, observer, maxDist):
        """Generates all image-points of source within a given distance of observer."""
        for q in self.primary_reflections(source):
            for z in map(lambda p: p + q, self.in_disk(observer - q, maxDist)):
                yield z

    def in_disk(self, center, radius):
        """Generates all lattice points within a distance radius of center."""
        xmin, xmax = center.x - radius, center.x + radius
        xmin += (-xmin) % (2 * self._dims.x)
        xmax -= xmax % (2 * self._dims.x)
        for x in range(xmin, xmax + 1, 2 * self._dims.x):
            dy = int(floor(sqrt(radius**2 - (x - center.x)**2)))
            ymin, ymax = center.y - dy, center.y + dy
            ymin += (-ymin) % (2 * self._dims.y)
            ymax -= ymax % (2 * self._dims.y)
            for y in range(ymin, ymax + 1, 2 * self._dims.y):
                yield Point(x, y)


if __name__ == "__main__":
    testcases = {
        ((3, 2), (1, 1), (2, 1), 4): 7,
        ((300, 275), (150, 150), (185, 100), 500): 9,
        ((2, 5), (1, 2), (1, 4), 11): 27,
        ((23, 10), (6, 4), (3, 2), 23): 8,
        ((1250, 1250), (1000, 1000), (500, 400), 10000): 196,
        ((10, 10), (4, 4), (3, 3), 5000): 739323,
        ((3, 2), (1, 1), (2, 1), 7): 19,
        ((2, 3), (1, 1), (1, 2), 4): 7,
        ((3, 4), (1, 2), (2, 1), 7): 10,
        ((4, 4), (2, 2), (3, 1), 6): 7,
        ((300, 275), (150, 150), (180, 100), 500): 9,
        ((3, 4), (1, 1), (2, 2), 500): 54243
    }
    for t, ans in testcases.items():
        s = solution(*t)
        if s != ans:
            print("FAILED test:")
            print("\tsolution(%s) == %s, not %s" %
                  (", ".join([str(x) for x in t]), ans, s))
