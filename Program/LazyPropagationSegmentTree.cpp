#include "LazyPropagationSegmentTree.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <vector>

LazyPropagationSegmentTree::LazyPropagationSegmentTree(int size) {
    n = size;
    treeSize = n * 4 + 1;
    t = std::vector(treeSize, 0);
    lazy = std::vector(treeSize, 0);
}

void LazyPropagationSegmentTree::build(int a[]) {
    build(a, 1, 0, n);
}

void LazyPropagationSegmentTree::build(int a[], int index, int start, int end) {
    if (start == end) {
        t[index] = a[start];
    } else {
        int middle = (start + end) / 2;
        build(a, index*2, start, middle);
        build(a, index*2+1, middle+1, end);
        t[index] = std::max(t[index*2], t[index*2 + 1]);
    }
}

void LazyPropagationSegmentTree::push(int index)
{
    if (index * 2 < treeSize) {

        t[index*2] += lazy[index];
        lazy[index*2] += lazy[index];
    }
    if (index * 2 + 1 < treeSize) {

        t[index*2+1] += lazy[index];
        lazy[index*2+1] += lazy[index];
    }
    lazy[index] = 0;
}

void LazyPropagationSegmentTree::update(int l, int r, int addend) {
    update(1, 0, n, l, r, addend);
}


void LazyPropagationSegmentTree::update(int index, int start, int end, int left, int right, int addend) {
    if (right < start || left > end)
        return;

    if (left <= start && right >= end) {
        t[index] += addend;
        lazy[index] += addend;
    }
    else {
        push(index);
        int middle = (start + end) / 2;
        update(index*2, start, middle, left, right, addend);
        update(index*2+1, middle+1, end, left, right, addend);

        t[index] = std::max(t[index*2], t[index*2+1]);
    }
}

int LazyPropagationSegmentTree::query(int l, int r) {
    return query(1, 0, n, l, r);
}

int LazyPropagationSegmentTree::query(int index, int start, int end, int left, int right) {
    assert (left <= right);

    push(index);

    if (right < start || left > end)
        return 0;

    if (left <= start && right >= end) {
        return t[index];
    }

    int middle = (start + end) / 2;

    return std::max(query(index*2, start, middle, left, right),
        query(index*2+1, middle+1, end, left, right));
}