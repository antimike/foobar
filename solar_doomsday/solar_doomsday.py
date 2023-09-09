from math import sqrt, floor


def solution(area):
    ans = []
    while area > 0:
        # Note that Python 2.7 doesn't have math.isqrt
        nearestSquare = int(floor(sqrt(area)))**2
        area -= nearestSquare
        ans.append(nearestSquare)
    return ans
