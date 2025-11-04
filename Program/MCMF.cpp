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

MCMF::MCMF(int nodes, int deliveries) : nodes(nodes), deliveries(deliveries) {}

Searcher GSearcher(TSol &s, long nodesSize)
{
    std::vector<int> nodeSols = std::vector<int>(nodesSize);

    for (int i = 0; i < nodesSize; i++)
    {
        nodeSols[i] = s.vec[i].sol;
    }

    return {nodesSize, nodeSols};
}

void MCMF::build(TSol &s, const std::vector<TDelivery> &del, const int Q)
{
    for (int i = 0; i <= nodes + 1; i++) {
        network.add_node();
    }

    for (int i = 0; i <= nodes; i++) {
        network.add_arc(NULL_ID, i, i + 1, 0, Q, 0);
    }

    Searcher searcher = GSearcher(s, nodes);

    for (int i = 0; i < deliveries; i++)
    {
        TDelivery delivery = del[i];

        int origin = searcher.search(delivery.origin);
        int destination = searcher.search(delivery.destination);

        if (destination == 0) {
            destination = nodes;
        }

        if (origin < destination)
        {
            network.add_arc(i, origin, destination, 0, CAPACITY, -delivery.value);
        }
    }
}

int MCMF::solve(int del[])
{
    int result = network.min_cost_max_flow(0, nodes + 1);

    for (int i = 0; i < network.arcs_.size(); i++) {
        auto arc = network.arcs_[i];
        if (arc.id_ != NULL_ID && arc.flow_ == CAPACITY) {
            del[arc.id_] = 1;
        }
    }

    return result;

}