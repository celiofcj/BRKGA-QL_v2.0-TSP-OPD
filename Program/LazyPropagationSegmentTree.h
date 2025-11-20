//
// Created by celio on 11/10/2025.
//

#ifndef BRKGA_QL_V2_0_LAZYPROPAGATIONSEGMENTTREE_H
#define BRKGA_QL_V2_0_LAZYPROPAGATIONSEGMENTTREE_H
#include <vector>


class LazyPropagationSegmentTree {
    private:
        std::vector<int> t;
        std::vector<int> lazy;
        int n;
        int treeSize;

        void build(int a[], int index, int start, int end);
        int query(int v, int start, int end, int l, int r);
        void update(int v, int start, int end, int l, int r, int addend);
        void push(int index);


    public:
        explicit LazyPropagationSegmentTree(int size);
        void build(int a[]);
        void update(int l, int r, int addend);
        int query(int l, int r);
};


#endif //BRKGA_QL_V2_0_LAZYPROPAGATIONSEGMENTTREE_H