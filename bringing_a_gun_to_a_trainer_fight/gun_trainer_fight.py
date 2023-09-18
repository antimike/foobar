from itertools import product
from fractions import Fraction


class ReflectionPoints:

    def __init__(self, x, y, length, height):
        self._x, self._y = x, y
        self._l, self._h = length, height

    def __contains__(self, pt):
        x, y = pt
        if x % (2 * self._l) in (self._x % (2 * self._l),
                                 -self._x % (2 * self._l)):
            if y % (2 * self._h) in (self._y % (2 * self._h),
                                     -self._y % (2 * self._h)):
                return True
        return False


class Displacement:

    def __init__(self, dx, dy):
        self._dx, self._dy = dx, dy
        self.bearing = self.getBearing(dx, dy)
        self.d2 = dx**2 + dy**2

    def displace(self, pos):
        x, y = pos
        return (x + self._dx, y + self._dy)

    @staticmethod
    def getBearing(dx, dy):
        m = -1 if dx < 0 or dx == 0 and dy < 0 else 1
        if dy == 0:
            return (m, 0)
        else:
            f = Fraction(dx, dy)
            return (m * f.numerator, m * f.denominator)


class Bearing:

    def __init__(self, x, y):
        if y == 0:
            self._ratio = float("inf")
        else:
            self._ratio = Fraction(x, y)


def solution(dimensions, your_position, trainer_position, distance):
    """
    Correct, but TLE.
    """
    s_lattice = ReflectionPoints(*your_position, *dimensions)
    t_lattice = ReflectionPoints(*trainer_position, *dimensions)

    targets = {}
    valid = set()
    D2 = distance**2

    for dx in range(distance + 1):
        for dy in range(distance + 1):
            if dx == dy == 0:
                continue
            if dx**2 + dy**2 > D2:
                break
            for mx, my in product((1, -1), (1, -1)):
                d = Displacement(mx * dx, my * dy)
                p = d.displace(your_position)
                targets.setdefault(d.bearing, [float("inf"), float("inf")])

                if p in t_lattice:
                    targets[d.bearing][0] = min(targets[d.bearing][0], d.d2)
                if p in s_lattice:
                    targets[d.bearing][1] = min(targets[d.bearing][1], d.d2)

                t_dist, s_dist = targets[d.bearing]
                if s_dist > t_dist:
                    valid.add(d.bearing)
                else:
                    valid.discard(d.bearing)

    return len(valid)
