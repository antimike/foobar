from math import gcd


def bearing(a, b):
    if a == b == 0:
        return (0, 0)
    d = gcd(a, b)
    return (a // d, b // d)


def get_all_bearings(dims, pos, tpos, dist):
    hits = {}
    for x in range(pos[0] - dist, pos[0] + dist + 1):
        dx = x - pos[0]
        dy_max = floor(sqrt(dist**2 - dx**2))
        ymin, ymax = pos[1] - dy_max, pos[1] + dy_max
        for y in range(ymin, ymax + 1):
            dy = y - pos[1]
            D2 = dx**2 + dy**2
            b = bearing(dx, dy)
            for tx, ty in tpos, [-tpos[0],
                                 tpos[1]], [tpos[0],
                                            -tpos[1]], [-tpos[0], -tpos[1]]:
                if (x - tx) % (2 * dims[0]) == 0 and (y - ty) % (2 *
                                                                 dims[1]) == 0:
                    hits.setdefault(b, [float("inf"), float("inf")])
                    hits[b][0] = min(hits[b][0], D2)
            for mx, my in pos, [-pos[0],
                                pos[1]], [pos[0], -pos[1]], [-pos[0], -pos[1]]:
                if b in hits and (x - mx) % (2 * dims[0]) == 0 and (y - my) % (
                        2 * dims[1]) == 0:
                    hits[b][1] = min(hits[b][1], D2)
    return [b for b, x in hits.items() if x[0] < x[1]]
