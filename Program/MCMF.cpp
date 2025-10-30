//
// Created by celio on 25/10/2025.
//

#include "MCMF.h"

#include <queue>

#include "Searcher.h"

#define INITIAL_COST 0
#define CAPACITY 1
#define NULL_ID (-1)

#define INF 99999999

#define vector std::vector

MCMF::MCMF(int verticesSize, int deliveriesSize)
{
    size = verticesSize;
    graph = vector<Edge> (verticesSize);
}

Searcher GetSearcher(TSol &s, long nodesSize)
{
    vector<int> nodeSols = vector<int>(nodesSize);

    for (int i = 0; i < nodesSize; i++)
    {
        nodeSols[i] = s.vec[i].sol;
    }

    return {nodesSize, nodeSols};
}

void MCMF::build(TSol &s, const vector<TDelivery> &del, const int Q)
{
    for (int i = 0; i < size - 1; i++)
    {
        int v1 = s.vec[i].sol, v2 = s.vec[i + 1].sol;
        graph.emplace_back(v1, v2, Q, INITIAL_COST);
    }

    for (int i = 0; i < del.size(); i++)
    {
        TDelivery delivery = del[i];

        Searcher searcher = GetSearcher(s, s.vec.size() - del.size());
        int origin = searcher.search(delivery.origin);
        int destination = searcher.search(delivery.destination);

        // printf("I: %d, Origin: %d, Destination: %d\n", i, origin, destination);

        if (origin < destination || destination == 0 && origin <= size - 1)
        {
            graph.emplace_back(origin, destination, CAPACITY, -delivery.value);
        }
    }
}

int MCMF::solve(vector<int> &edges) {
    return minCostFlow(INF, 0, edges.size());
}

vector<vector<int>> adj, cost, cap;

void MCMF::shortestPaths(int n, int v0, vector<int>& d, vector<int>& p) {
    d.assign(n, INF);
    d[v0] = 0;
    vector<bool> inq(n, false);
    std::queue<int> q;
    q.push(v0);
    p.assign(n, -1);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        inq[u] = false;
        for (int v : adj[u]) {
            if (cap[u][v] > 0 && d[v] > d[u] + cost[u][v]) {
                d[v] = d[u] + cost[u][v];
                p[v] = u;
                if (!inq[v]) {
                    inq[v] = true;
                    q.push(v);
                }
            }
        }
    }
}


int MCMF::minCostFlow(const int K, const int s, const int t) {
    adj.assign(size, vector<int>());
    cost.assign(size, vector<int>(size, 0));
    cap.assign(size, vector<int>(size, 0));
    for (Edge e : graph) {
        adj[e.from].push_back(e.to);
        adj[e.to].push_back(e.from);
        cost[e.from][e.to] = e.cost;
        cost[e.to][e.from] = -e.cost;
        cap[e.from][e.to] = e.capacity;
    }

    int flow = 0;
    int cost = 0;
    vector<int> d, p;
    while (flow < K) {
        shortestPaths(size, s, d, p);
        if (d[t] == INF)
            break;

        int f = K - flow;
        int cur = t;
        while (cur != s) {
            f = std::min(f, cap[p[cur]][cur]);
            cur = p[cur];
        }

        flow += f;
        cost += f * d[t];
        cur = t;
        while (cur != s) {
            cap[p[cur]][cur] -= f;
            cap[cur][p[cur]] += f;
            cur = p[cur];
        }
    }

    return flow < K ? -1 : cost;
}
