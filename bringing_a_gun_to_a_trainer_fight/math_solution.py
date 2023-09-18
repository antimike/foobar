"""
After failing with brute-force and Farey-star-based solutions, I decided to finally sit down and solve the problem "the old-fashioned way."

Fortunately, it's not actually too difficult, and the solution cuts the runtime from O(dist**2) to O(dist).
"""

from .vector import IntegerVector, BilinearForm
from itertools import combinations


class Range:
    from bisect import bisect_left

    def __init__(self, include=[], exclude=[]):
        self._intervals = [-float("inf"), float("inf")]
        for interval in include:
            self.add_interval(*interval)
        for interval in exclude:
            self.remove_interval(*interval)

    def _get_endpoint_index(self, num):
        """Helper function wrapping bisect_left"""
        return self.bisect_left(self._intervals, num)

    def add_interval(self, a, b):
        i, j = self._get_endpoint_index(a), self._get_endpoint_index(b)
        if i == j:
            if i % 2 == 0:
                self._intervals.insert(i, b)
                self._intervals.insert(i, a)
        else:
            if i % 2 == 0:
                self._intervals[i] = a
            if j % 2 == 0:
                # we can assume j > 0 here due to the case handled above (i == j == 0)
                self._intervals[j - 1] = b
            self._intervals = self._intervals[:i + 1 - (i % 2)] + \
                self._intervals[j - 1 + (j % 2):]

    def remove_interval(self, a, b):
        i, j = self._get_endpoint_index(a), self._get_endpoint_index(b)
        if i == j:
            if i % 2 == 1:
                self._intervals.insert(i, b)
                self._intervals.insert(i, a)
        else:
            if i % 2 == 1:
                self._intervals[i] = a
            if j % 2 == 1:
                self._intervals[j - 1] = b
            self._intervals = self._intervals[:i + (i % 2)] + \
                self._intervals[j - (j % 2):]

    def size(self):
        if len(self._intervals) > 0:
            if self._intervals[0] == -float(
                    "inf") or self._intervals[-1] == float("inf"):
                return float("inf")
            else:
                return sum([
                    x * (-1)**(j + 1) for j, x in enumerate(self._intervals)
                ]) + len(self._intervals)
        else:
            return 0

    def clear(self):
        self._intervals = []

    def copy(self):
        ret = self.__class__()
        ret._intervals = list(self._intervals)
        return ret

    def union(self, other):
        ret = self.copy()
        for a, b in other.intervals:
            ret.add_interval(a, b)
        return ret

    def difference(self, other):
        ret = self.copy()
        for a, b in other.intervals:
            ret.remove_interval(a, b)
        return ret

    def intersect(self, other):
        ret = self.copy()
        ret.clear()
        for a, b in other.intervals:
            i, j = self._get_endpoint_index(a), self._get_endpoint_index(b)
            if i == j:
                if i % 2 == 1:
                    ret.add_interval(a, b)
            else:
                if i % 2 == 1:
                    ret.add_interval(a, self._intervals[i])
                if j % 2 == 1:
                    ret.add_interval(self._intervals[j - 1], b)

    @property
    def intervals(self):
        for i in range(0, len(self._intervals) - 1, 2):
            yield (self._intervals[i], self._intervals[i + 1])

    def __add__(self, other):
        return self.union(other)

    def __sub__(self, other):
        return self.difference(other)

    def __contains__(self, x):
        return self._get_endpoint_index(x) % 2 == 1


def solution(dims, pos, tpos, dist):
    me, target = IntegerVector(*pos), IntegerVector(*tpos)
    lattice = IntegerVector.Lattice(2 * dims[0] * IntegerVector.xhat,
                                    2 * dims[1] * IntegerVector.yhat)
    collision_vectors = []
    for p, q in combinations(
            [v - u for u in me.reflections() for v in target.reflections()], 2):
        if p - q in lattice:
            collision_vectors.append(lattice.coefficients(p - q))
