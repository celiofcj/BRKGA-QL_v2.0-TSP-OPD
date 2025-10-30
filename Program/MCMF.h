//
// Created by celio on 25/10/2025.
//

#ifndef BRKGA_QL_V2_0_MCMF_H
#define BRKGA_QL_V2_0_MCMF_H
#include <vector>

#include "Edge.h"
#include "Problem.h"


class MCMF {
    private:
        int size;
        std::vector<Edge> graph;
        static void shortestPaths(int n, int v0, std::vector<int>& d, std::vector<int>& p);
        int minCostFlow(int K, int s, int t);
    public:
        explicit MCMF(int verticesSize, int deliveriesSize);
        void build(TSol &s, const std::vector<TDelivery> &deliveries, int Q);
        int solve(std::vector<int> &edges);
};




#endif //BRKGA_QL_V2_0_MCMF_H