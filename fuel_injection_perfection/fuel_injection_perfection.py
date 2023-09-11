def solution(n):
    n = bin(int(n))[2:]
    pos = n.rfind('1')
    ans = len(n) - pos - 1
    while pos > 0:
        z = n.rfind('0', 0, pos)
        if z == -1:
            if pos == 1:
                return ans + 2
            else:
                return ans + pos + 2
        elif pos - z == 1:
            pos = n.rfind('1', 0, z)
            ans += 2 + z - pos
        else:
            ans += 1 + pos - z
            pos = z
    return ans
