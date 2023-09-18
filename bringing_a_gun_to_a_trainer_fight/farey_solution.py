from .vector import IntegerVector


def solution(dims, pos, tpos, dist):
    """
    Both incorrect and TLE :/

    Not every valid target point is on the Farey star.  In addition, the size of the star grows quadratically, so this isn't much better than the "naive" solution.
    """
    D2 = dist**2
    me, target = IntegerVector(*pos), IntegerVector(*tpos)
    lattice = IntegerVector.Lattice(2 * dims[0] * IntegerVector.xhat,
                                    2 * dims[1] * IntegerVector.yhat)
    ans = 0
    for bearing in filter(lambda v: v * v <= D2,
                          IntegerVector.FareyStar(dist)):
        if any(map(lambda v: me + bearing - v in lattice,
                   target.reflections())):
            print(f"Found: {bearing}")
            ans += 1
    return ans


if __name__ == "__main__":
    t1 = [300, 275], [150, 150], [185, 100], 500
    t2 = [3, 2], [1, 1], [2, 1], 4
    print(solution(*t1))
    print(solution(*t2))
