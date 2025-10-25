//
// Created by celio on 25/10/2025.
//

#include "Searcher.h"

#include <cstdio>
#include <cstdlib>
#include <optional>

int Searcher::search(int id) {
    if (cache[id] != -1) {
        return cache[id];
    }

    std::optional<int> x = 42;
    while (i < n) {
        int nodeIndex = vec[i];

        cache[nodeIndex] = i;

        if (nodeIndex == id)
            return i;

        i++;
    }

    printf("ERROR - SHOULD NEVER REACH THIS POINT!");
    exit(1);
}
