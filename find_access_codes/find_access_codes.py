from collections import Counter


def solution(l):
    """
    Passes all but the last testcase (TLE?)
    O(N**3)
    """
    N = len(l)
    ct = 0
    nums = sorted([(x, i) for i, x in enumerate(l)])
    for pos, pair in enumerate(nums):
        n, i = pair
        for L in range(pos + 1, N):
            a, j = nums[L]
            if j <= i or a % n != 0:
                continue
            for R in range(N - 1, L, -1):
                b, k = nums[R]
                if k <= j or b % a != 0:
                    continue
                ct += 1
    return ct


def solution(l):
    """
    A stripped-down version of the above.  Still TLEs on the last testcase.
    """
    N = len(l)
    ct = 0
    for i, x in enumerate(l):
        for j in range(i + 1, N):
            if l[j] % l[i] != 0:
                continue
            for k in range(j + 1, N):
                if l[k] % l[j] == 0:
                    ct += 1
    return ct


def solution(l):
    """
    Asymptotically superior to the above, at O(N**2).  The key is to iterate backwards through the list.
    """
    ct = Counter()
    total = 0
    for i, x in enumerate(l[-1::-1], 1):
        for j, y in enumerate(l[-i - 1::-1], i + 1):
            if x % y == 0:
                ct[j] += 1
                total += ct[i]
    return total
