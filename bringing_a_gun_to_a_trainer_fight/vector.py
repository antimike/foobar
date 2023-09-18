from .farey import farey


class BilinearForm:

    def __init__(self, col1, col2):
        self._cols = (Vector(*col1), Vector(*col2))
        self._rows = (Vector(col1[0], col2[0]), Vector(col1[1], col2[1]))

    def det(self):
        return self[0, 0] * self[1, 1] - self[0, 1] * self[1, 0]

    @property
    def transpose(self):
        return self.__class__(*self._rows)

    def __mul__(self, other):
        if isinstance(other, Vector):
            return other[0] * self._cols[0] + other[1] * self._cols[1]
        elif isinstance(other, self.__class__):
            return self.__class__(self * other._cols[0], self * other._cols[1])
        elif isinstance(other, int) or isinstance(other, float):
            return self.__class__(*[other * col for col in self._cols])

    def __getitem__(self, item):
        try:
            m, n = iter(item)
            return self._rows[m][n]
        except TypeError:
            return self._rows[item]


class Vector:
    from math import sqrt

    class Lattice:

        def __init__(self, u, v):
            self._form = BilinearForm(u, v)
            self._det = self._form.det()

        def __contains__(self, vec):
            m1, m2 = BilinearForm(vec, self._form._cols[1]), BilinearForm(
                self._form._cols[0], vec)
            return m1.det() % self._det == 0 and m2.det() % self._det == 0

    def __init__(self, x, y):
        self._tuple = (x, y)

    @property
    def x(self):
        return self._tuple[0]

    @property
    def y(self):
        return self._tuple[1]

    def __hash__(self):
        return hash(self._tuple)

    def __eq__(self, other):
        return self._tuple == other._tuple

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, idx):
        return self._tuple[idx]

    def __add__(self, other):
        return self.__class__(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return self.__class__(self.x - other.x, self.y - other.y)

    def __mul__(self, other):
        return self.x * other.x + self.y * other.y

    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self.__class__(other * self.x, other * self.y)
        else:
            return NotImplemented

    def __neg__(self):
        return self.__class__(-self.x, -self.y)

    def __iter__(self):
        """Allows unpacking"""
        return iter(self._tuple)

    def __repr__(self):
        return f"<{self.x}, {self.y}>"

    def __abs__(self):
        return self.sqrt(self * self)

    def reflect(self, other):
        return 2 * (self * other) / (other * other) * other - self

    def reflections(self):
        return set(
            (self, self.reflect(self.xhat), self.reflect(self.yhat), -self))

    @classmethod
    @property
    def xhat(cls):
        return cls(1, 0)

    @classmethod
    @property
    def yhat(cls):
        return cls(0, 1)


class IntegerVector(Vector):

    def __init__(self, x, y):
        super().__init__(int(x), int(y))

    @classmethod
    def FareyStar(cls, maxDenom):
        for a, b in farey(maxDenom):
            v = cls(a, b)
            yield from v.reflections().union(
                v.reflect(cls.xhat + cls.yhat).reflections())

    @classmethod
    def printFarey(cls,
                   N,
                   blank=" ",
                   star="*",
                   spacing_horiz=4,
                   spacing_vert=2):
        S = set(cls.FareyStar(N))
        print(S)
        lines = []
        for y in range(N, -N - 1, -1):
            line = []
            for x in range(-N, N + 1):
                if (x, y) == (0, 0):
                    line.append("O")
                elif Vector(x, y) in S:
                    line.append(star)
                else:
                    line.append(blank)
            lines.append(line)
        print(("\n" * spacing_vert).join([(blank * spacing_horiz).join(l)
                                          for l in lines]))
