//
// Created by celio on 25/10/2025.
//

#ifndef BRKGA_QL_V2_0_NODE_H
#define BRKGA_QL_V2_0_NODE_H
#include <memory>


class Edge {
    public:
        int from;
        int to;
        int capacity;
        int cost;
        Edge();
        Edge(const int from, const int to, const int q, const int cost): from(from), to(to), capacity(q), cost(cost) {}
};


#endif //BRKGA_QL_V2_0_NODE_H