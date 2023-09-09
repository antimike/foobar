from collections import deque
from itertools import product, chain


def knightlyNeighbors(pos):
    for dx, dy in chain(product((1, -1), (2, -2)), product((2, -2), (1, -1))):
        if 0 <= (pos % 8) + dx < 8 and 0 <= pos + 8 * dy < 64:
            yield pos + dx + 8 * dy


def solution(src, dest):
    seen = {src: 0}
    queue = deque()
    loc = src
    while dest not in seen:
        for c in knightlyNeighbors(loc):
            if c not in seen:
                queue.append(c)
                seen[c] = seen[loc] + 1
        loc = queue.popleft()
    return seen[dest]
