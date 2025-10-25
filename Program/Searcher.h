//
// Created by celio on 25/10/2025.
//

#ifndef BRKGA_QL_V2_0_SEARCHER_H
#define BRKGA_QL_V2_0_SEARCHER_H
#include <utility>
#include <vector>


class Searcher {
    private:
        std::vector<int> cache;
        int i;
        long n;
        std::vector <int> vec;
    public:
        Searcher(long n, std::vector<int> vec) : cache(n, -1), i(0), n(n), vec(std::move(vec)) {};
        int search(int id);
};


#endif //BRKGA_QL_V2_0_SEARCHER_H