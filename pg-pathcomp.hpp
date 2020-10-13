#pragma once

#include <iostream>
#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/path_position_handle_graph.hpp"

// compute pairwise path jaccard distance and print to table in os
std::pair<std::vector<std::string>, std::vector<std::vector<float>>> all_paths_jaccard(const handlegraph::PathPositionHandleGraph& graph);
