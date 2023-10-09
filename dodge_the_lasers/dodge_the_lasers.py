from math import sqrt, floor
from os import path


def convergent_of_sqrt2(num_decimal_places):
    """
    Computes a rational approximation of sqrt(2) to a specified accuracy.
    Args:
        num_decimal_places (int): The number of decimal places to which the
            return value is guaranteed to be accurate.

    Returns:
        An ordered pair (A, B) such that gcd(A, B) == 1 and
            abs(sqrt(2) - (A/B)) <= 10**-num_decimal_places.

    """
    a, b, c, d = 1, 2, 1, 1
    num_it = 0
    while num_it < 1000:
        num_it += 1
        a, b, c, d = a**2 + b * c, a * b + b * d, a * c + c * d, b * c + d**2
        A, B = a + b, c + d
        if 5 * 10**num_decimal_places * abs(A**2 - 2 * B**2) < 14 * B**2:
            return A, B


def brute_force(num_str):
    N = int(num_str)
    return int(sum([floor(j * sqrt(2)) for j in range(N + 1)]))


def solution(num_str):
    A, B = convergent_of_sqrt2(300)
    N = int(num_str)
    S = 0
    sgn = 1
    while N > 10**4:
        M = N * A // B
        N = M - N
        S += sgn * (M * (M + 1) // 2 - N * (N + 1))
        sgn *= -1
    return str(S + sgn * brute_force(N))


if __name__ == "__main__":
    with open(path.join(path.dirname(__file__), "testcases.txt"),
              "r") as tests:
        test_num = 0
        while tests.readline():
            test_num += 1
            num_str = tests.readline().strip()
            ans = tests.readline().strip()
            print("Running testcase %s:\n\tnum_str = %s\n\tans = %s" %
                  (test_num, num_str, ans))
            computed = solution(num_str)
            if computed != ans:
                print("\tFAILED:\n\t\treturned %s" % computed)
