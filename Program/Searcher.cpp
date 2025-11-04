//
// Created by celio on 25/10/2025.
//

#include "Searcher.h"

#include <cstdio>
#include <cstdlib>

int Searcher::search(int id) {
    if (cache[id] != -1) {
        return cache[id];
    }

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
