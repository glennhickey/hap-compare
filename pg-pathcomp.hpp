#pragma once

#include <iostream>
#include <map>
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/path_position_handle_graph.hpp"
#include "vg/vg.pb.h"
#include "handlegraph/util.hpp"

using namespace std;
using namespace vg;

// compute pairwise path jaccard distance and print to table in os
pair<vector<string>, vector<vector<float>>> all_paths_jaccard(const handlegraph::PathPositionHandleGraph& graph);

// compute wiggle with the number of sv's vs reference path
vector<float> sv_pileup(const handlegraph::PathPositionHandleGraph& graph,
                        const vector<Snarl>& snarls,
                        const string& ref_path,
                        size_t window_size,
                        size_t min_sv_len);

// get the sv coverage from an alt traversal
void process_alt(const handlegraph::PathPositionHandleGraph& graph,
                 const set<handlegraph::handle_t>& ref_handles,
                 const SnarlTraversal& ref_traversal,
                 handlegraph::path_handle_t ref_path_handle,
                 size_t ref_pos,
                 const pair<handlegraph::step_handle_t, handlegraph::step_handle_t>& alt_traversal,
                 size_t min_sv_len,
                 vector<float>& ref_pileup);

// find traversals through a snarl
pair<vector<SnarlTraversal>, vector<pair<handlegraph::step_handle_t, handlegraph::step_handle_t> > > find_path_traversals(const handlegraph::PathPositionHandleGraph& graph, const Snarl& site);

// get the traversal lengths
vector<size_t> get_trav_lengths(const handlegraph::HandleGraph& graph, const vector<SnarlTraversal>& travs);
