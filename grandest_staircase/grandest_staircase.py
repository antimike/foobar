def numStaircases(N, k, memo):
    """
    Number of staircases with N bricks and k steps.
    """
    if (N, k) in memo:
        return memo[N, k]
    elif N < k * (k + 1) // 2 or k < 2:
        return 0
    elif N == k * (k + 1) // 2:
        return 1
    elif k == 2:
        return 1 + (N - 3) // 2

    memo[N, k] = 0
    for m in range(N - k, k * (k - 1) // 2 - 1, -k):
        memo[N, k] += numStaircases(m, k - 1, memo)

    return memo[N, k]


def solution(n):
    memo = {}
    total = 0
    k = 2
    while k * (k + 1) // 2 <= n:
        total += numStaircases(n, k, memo)
        k += 1
    return total
