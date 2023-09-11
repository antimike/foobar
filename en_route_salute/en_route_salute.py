def solution(s):
    numRight = 0
    ct = 0
    for c in s:
        if c == '>':
            numRight += 1
        elif c == '<':
            ct += 2 * numRight
    return ct
