def bellman_ford(adj, start):
    """
    Implementation of the Bellman-Ford single-source shortest path algorithm.

    For this problem, we don't need to recover the actual shortest paths, but
    we do need to account for the possibility of a negative cycle.
    See https://cp-algorithms.com/graph/bellman_ford.html.

    Args:
      adj: adjacency matrix of weighted digraph
      start: source node

    Returns:
      dist: list of minimal path-weights from source to each vertex in adj
      pred: list of each node's predecessor on an optimal path from source
      negative_cycle: list of nodes in a negative cycle, if any is found
    """
    N = len(adj)
    dist = [float("inf")] * N
    pred = [-1] * N
    dist[start] = 0
    for phase in range(N):
        last = None
        for a, b, c in [(i, j, adj[i][j]) for i in range(N) for j in range(N)
                        if i != j]:
            if dist[a] + c < dist[b]:
                dist[b] = dist[a] + c
                pred[b] = a
                last = b
    negative_cycle = []
    if last is not None:
        # negative cycle
        for i in range(N):
            last = pred[last]
        curr = pred[last]
        negative_cycle.append(last)
        while curr != last:
            negative_cycle.insert(0, curr)
            curr = pred[curr]
    return dist, pred, negative_cycle


def solution(times, time_limit):
    # Approach: run Bellman-Ford starting at each vertex to get the "actual"
    # edge weights, i.e., the minimum amount of time required to travel
    # between vertices i and j (along a path of any length) in the case of an
    # edge (i, j).  Then, run a DFS on the updated digraph (a sort of
    # "weighted transitive closure") starting at the source node 0.
    N = len(times)
    START, EXIT = 0, N - 1
    min_path_lengths = []

    for i in range(N):
        dist, pred, neg = bellman_ford(times, i)
        if len(neg) > 0:
            # If there's a negative cycle, then all bunnies are accessible
            return list(range(N - 2))
        min_path_lengths.append(dist)

    saved = set()

    def helper(path, time):
        if len(path) > N:
            return
        curr = path[-1]
        if curr == EXIT and time >= 0:
            # bunnies are 0-indexed
            bunnies = [x - 1 for x in filter(lambda x: START < x < EXIT, path)]
            saved.add(tuple(sorted(bunnies)))
        for node in range(N):
            if node not in path and min_path_lengths[curr][node] + \
                    min_path_lengths[node][EXIT] <= time:
                path.append(node)
                helper(path, time - min_path_lengths[curr][node])
                path.pop()

    helper([START], time_limit)
    num_saved = max(len(l) for l in saved)
    return list(
        sorted([l for l in saved if len(l) == num_saved], reverse=True).pop())

    # Let M(i1, ..., ik) := minimal path-weight among paths connecting i1, ..., ik **beginning with i1 and ending with ik**
    # M(i, j, k) >= max(M(i, j), M(i, k))
    # However, M(i, j, k) >= M(j, k) is **not** true, since M(i, j) could be negative


if __name__ == "__main__":
    import os
    with open(os.path.join(os.path.dirname(__file__), "testcases.txt"),
              "r") as testfile:
        while testfile.readline():
            dim = int(testfile.readline().strip())
            times = []
            for row in range(dim):
                times.append(
                    [int(x) for x in testfile.readline().strip().split()])
            limit = int(testfile.readline().strip())
            ans = [int(x) for x in testfile.readline().strip().split()]

            print("===================================")
            print("Running testcase:\ntimes = %s\nlimit = %s" % (times, limit))
            print("")
            try:
                computed = solution(times, limit)
                assert ans == computed, "Expected: %s\tReceived: %s" % (
                    ans, computed)
            except Exception as e:
                print("FAILED: %s" % e)
            print("===================================")
