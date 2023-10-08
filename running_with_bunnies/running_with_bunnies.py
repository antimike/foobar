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
        list of minimal path-weights from source to each vertex in adj

    Raises:
        ValueError: if a negative cycle is found
    """
    N = len(adj)
    dist = [float("inf")] * N
    dist[start] = 0
    for phase in range(N):
        last = None
        for a, b, c in [(i, j, adj[i][j]) for i in range(N) for j in range(N)
                        if i != j]:
            if dist[a] + c < dist[b]:
                dist[b] = dist[a] + c
                last = b
    if last is not None:
        # negative cycle
        raise ValueError("Digraph contains a negative cycle")
    return dist


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
        try:
            min_path_lengths.append(bellman_ford(times, i))
        except ValueError:
            # If there's a negative cycle, then all bunnies are accessible
            return list(range(N - 2))

    def best_solution(s1, s2):
        if len(s1) != len(s2):
            return max(s1, s2, key=len)
        else:
            return min(s1, s2)

    def dfs(path, time_remaining, can_save):
        if len(path) <= N:
            curr = path[-1]
            if curr == EXIT:
                if time_remaining >= 0 and len(path) >= len(can_save) + 2:
                    # bunnies are 0-indexed
                    bunnies = sorted([x - 1 for x in path if START < x < EXIT])
                    can_save = best_solution(can_save, bunnies)
            else:
                for node in range(N):
                    if node not in path and min_path_lengths[curr][node] + \
                            min_path_lengths[node][EXIT] <= time_remaining:
                        path.append(node)
                        can_save = best_solution(
                            dfs(path,
                                time_remaining - min_path_lengths[curr][node],
                                can_save), can_save)
                        path.pop()
        return can_save

    return dfs([START], time_limit, [])


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
