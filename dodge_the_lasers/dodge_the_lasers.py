"""
The solution is short but somewhat opaque, so I've included an explanation
on the off-chance that anyone actually reads this.

Notation:
    For irrational a, let
        S(a) := {floor(k*a) for k a positive integer}.
    For integer N, let
        Q(N) := sum([floor(sqrt(2)*j) for j in range(N + 1)])
    (i.e., the quantity we seek to calculate).

It's a well-known fact from elementary mathematics (which, iirc, I first
encountered in Arthur Engel's "Problem-Solving Strategies") that, for
irrational a and b,
    (1/a) + (1/b) == 1 iff
    1. S(a) and S(b) do not intersect, and
    2. the union of S(a) and S(b) contains all positive integers.

Clearly we are interested in S(sqrt(2)).  By the above result, S(sqrt(2))
contains all positive integers except those appearing in the complimentary
sequence S(x), where
    x = 1/(1 - 1/sqrt(2)) = 2 + sqrt(2).
But since floor(N + x) == N + floor(x) for any integer N, we have
    S(x) == {2*j + floor(j*sqrt(2)) for j a positive integer}
            == {2*k + t_k for k >= 1},
where t_k denotes the kth element of the sequence S(sqrt(2)).

In particular, if we let
    M := floor(N * sqrt(2))
    K := max pos. integer such that 2*K + t_K <= M
    T(n) := sum of the first n terms of S(2 + sqrt(2))
then by the above, we have
    Q(N) + T(K) == M*(M + 1) // 2
        (since all n from 1 to M appear exactly once in the LHS)
and
    T(K) == K*(K + 1) + Q(K)
        (by the identity given 12 lines above).
Thus,
    Q(N) == M*(M + 1) // 2 - K*(K + 1) - Q(K).
This recurrence allows a recursive calculation of Q(N).

It remains to calculate K.  To that end, note that
    k * (2 + sqrt(2)) <= N * sqrt(2)
<=> k <= N * (sqrt(2) - 1);
thus, we have
    K == floor(N * (sqrt(2) - 1)) == M - N.

Assuming that we can calculate M in constant time, this implies a total
runtime of O(lg(N)), since at each application of the recurrence we
replace N with N * (sqrt(2) - 1).  Fortunately, Python 2.5+ includes
builtin arbitrary-precision integer arithmetic, so we don't have to
worry about implementing bigint multiplication; however, we do need
access to a rational approximation of sqrt(2) accurate to 10**-100.
Empirically, I found that an appropriate convergent of sqrt(2) takes on
the order of microseconds to calculate, while the main algorithm runs in
time on the order of milliseconds per 100-digit input, so no significant
performance benefit results from precomputing a good approximation and
copy-pasting it into the solution.
"""

from os import path


def convergent_of_sqrt2(num_decimal_places):
    """
    Computes a rational approximation of sqrt(2) to a specified accuracy.

    The convergents of sqrt(2) are given by a_n / b_n, where
        [a_n, b_n]^T := [[1, 2], [1, 1]]**n * [1, 1]^T
    (where "X^T" refers to the transpose of X, and M**n to the normal matrix
    power of M).  By a well-known result in the theory of continued fractions,
    the error R_n := abs((a_n/b_n) - sqrt(2)) is bounded by
        R_n < 1/(b_n * b_{n+1}) = 1/(b_n * (a_n + b_n))
    Thus, in order to guarantee the requested accuracy, it is sufficient to
    ensure that
        R_n < 10**-num_decimal_places
    <=> 10**num_decimal_places < b_n * (a_n + b_n).
    The strategy used is simply to successively square the matrix
    [[1, 2], [1, 1]], obtaining a_{2**k} and b_{2**k} as the sums of the first
    and second row (respectively), and stopping when the above error condition
    is met.  This will obviously "overshoot" in the majority of cases,
    yielding return values significantly longer than strictly necessary.
    However, I found empirically that the performance hit associated with the
    larger convergents returned by this algorithm was insignificant compared
    to the penalty incurred by locating the "optimal" convergent (i.e., the
    first one satisfying the given error condition).

    Args:
        num_decimal_places (int): The number of decimal places to which the
            return value is guaranteed to be accurate.

    Returns:
        An ordered pair (A, B) such that gcd(A, B) == 1 and
            abs(sqrt(2) - (A/B)) <= 10**-num_decimal_places.
    """
    a, b, c, d = 1, 2, 1, 1
    while 10**num_decimal_places >= (c + d) * (a + b + c + d):
        a, b, c, d = a**2 + b * c, a * b + b * d, a * c + c * d, b * c + d**2
    return a + b, c + d


def solution(num_str):
    # We need sqrt(2) accurate to at least 100 digits
    A, B = convergent_of_sqrt2(100)
    N = int(num_str)
    S = 0
    sgn = 1
    # I tried short-circuiting and returning S + sgn * brute_force_solution(N)
    # when N dips below some cutoff, but found the performance benefit to be
    # negligible.
    while N > 0:
        M = N * A // B
        N = M - N
        S += sgn * (M * (M + 1) // 2 - N * (N + 1))
        sgn *= -1
    return str(S)


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
