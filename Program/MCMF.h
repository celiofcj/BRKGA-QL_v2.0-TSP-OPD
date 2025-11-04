//
// Created by celio on 25/10/2025.
//

#ifndef BRKGA_QL_V2_0_MCMF_H
#define BRKGA_QL_V2_0_MCMF_H
#include <queue>
#include <vector>
#include<vector>
#include<deque>
#include<queue>
#include<limits>
#include<tuple>
#include<algorithm>
#include<iostream>

#include "Problem.h"


class MCMF {
    struct SuccessiveShortestPathFlowNetwork {
        struct Arc {
            int id_;
            int start_;
            int end_;
            int flow_;
            int capacity_;
            int cost_;

            int get_dest(int from) {
                // Return the destination node if we were to traverse the arc from node `from`
                if(from == start_) {
                    return end_;
                }
                return start_;
            }

            void add_flow(int from, int to_add) {  // Adds flow from originating vertex `from`
                if(from == start_) {
                    flow_ += to_add;
                } else {
                    flow_ -= to_add;
                }
            }

            int get_capacity(int from) {
                // Gets the capacity of the edge if the originating vertex is `from`
                if(from == start_) {
                    return capacity_ - flow_;
                }
                return flow_;
            }

            int get_cost_from(int from) {
                if(from == start_) {
                    return cost_;
                }
                return -cost_;
            }
        };

        struct Node {
            int index_;
            std::vector<Arc*> connected_arcs_;
        };

        std::vector<Node> nodes_;
        std::deque<Arc> arcs_;

        void add_node() {
            nodes_.push_back({int(nodes_.size()), {}});
        }

        int add_arc(int id, int start, int end, int flow, int capacity, int cost) {
            arcs_.push_back({id, start, end, flow, capacity, cost});
            nodes_[start].connected_arcs_.push_back(&arcs_.back());
            nodes_[end].connected_arcs_.push_back(&arcs_.back());
            return (int)arcs_.size() - 1;
        }

        // Successive shortest paths min-cost max-flow algorithm
        // If there is an negative cost cycle initially, then it goes into infinite loop
        int min_cost_max_flow(int source_i, int sink_i) {
            int result = 0;

            // First calculate the potentials with Bellman–Ford derivative
            // It starts from a single vertex and is optimized to only operate on active vertices on each layer
            // Thus, it works more like a BFS derivative
            std::vector<int> potentials(nodes_.size(), std::numeric_limits<int>::max());
            {

                std::deque<std::pair<int, int> > front;
                front.push_back({0, source_i});

                while(front.size() > 0) {
                    int potential;
                    int cur_i;
                    std::tie(potential, cur_i) = front.front();
                    front.pop_front();

                    if(potential >= potentials[cur_i]) {
                        continue;
                    }
                    potentials[cur_i] = potential;

                    for(Arc* arc : nodes_[cur_i].connected_arcs_) if(arc->get_capacity(cur_i) > 0) {
                        // Traverse the arc if there is some remaining capacity
                        front.push_back({potential + arc->get_cost_from(cur_i), arc->get_dest(cur_i)});
                    }
                }
            }

            // Next loop Dijkstra to saturate flow. Once we subtract the difference in potential, every arc will have a
            // non-negative cost in both directions, so using Dijkstra is safe

            while(1) {
                std::priority_queue<std::tuple<int, int, Arc*> > frontier;
                std::vector<bool> explr(nodes_.size(), false);
                std::vector<int> cost_to_node(nodes_.size(), -1);
                std::vector<Arc*> arc_used(nodes_.size(), NULL);
                frontier.push({0, source_i, NULL});

                while(frontier.size() > 0) {
                    int path_cost;
                    int cur_i;
                    Arc* cur_arc_used;
                    std::tie(path_cost, cur_i, cur_arc_used) = frontier.top();
                    path_cost = -path_cost;
                    frontier.pop();

                    if(!explr[cur_i]) {
                        explr[cur_i] = true;
                        arc_used[cur_i] = cur_arc_used;
                        cost_to_node[cur_i] = path_cost;

                        for(Arc* arc : nodes_[cur_i].connected_arcs_) if(arc->get_capacity(cur_i) > 0) {
                            int next_i = arc->get_dest(cur_i);
                            // As priority_queue is a max-heap, we use the negative of the path cost for convenience
                            // We subtract the difference of potentials from the arc cost to ensure all arcs have positive
                            // cost
                            frontier.push({
                                -path_cost - (arc->get_cost_from(cur_i) - potentials[next_i] + potentials[cur_i]),
                                next_i,
                                arc
                            });
                        }
                    }
                }

                if(arc_used[sink_i] == NULL) {
                    return result;  // We didn't find a path, so return
                }
                std::vector<Arc*> arcs;
                int flow_pushed = std::numeric_limits<int>::max();
                {
                    // Here we counstruct the path of arcs from source to sink
                    int cur_i = sink_i;
                    while(cur_i != source_i) {
                        Arc* arc = arc_used[cur_i];
                        cur_i = arc->get_dest(cur_i);
                        flow_pushed = std::min(flow_pushed, arc->get_capacity(cur_i));
                        arcs.push_back(arc);
                    }

                    // Next we push flow back across all the arcs
                    for(auto arc_it = arcs.rbegin(); arc_it != arcs.rend(); arc_it++) {
                        Arc* arc = *arc_it;
                        arc->add_flow(cur_i, flow_pushed);
                        result += arc->get_cost_from(cur_i) * flow_pushed;
                        cur_i = arc->get_dest(cur_i);
                    }
                }

                // Finally, update the potentials so all edge-traversal costs remain non-negative
                for(int i=0; i<(int)nodes_.size(); i++) if(cost_to_node[i] != -1) {
                    potentials[i] += cost_to_node[i];
                }
            }
        }
    };
    SuccessiveShortestPathFlowNetwork network;
    int nodes;
    int deliveries;
    public:
        explicit MCMF(int nodes, int deliveries);
        void build(TSol &s, const std::vector<TDelivery> &deliveries, int Q);
        int solve(int del[]);
};

#endif //BRKGA_QL_V2_0_MCMF_H