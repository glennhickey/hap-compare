#include <unordered_map>
#include <unordered_set>
#include "pg-pathcomp.hpp"
#include "vg.pb.h"
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/util.hpp"
#include "vg/io/json2pb.h"

//#define debug
//#define debug1

using namespace std;
using namespace handlegraph;
using namespace vg;

pair<vector<string>, vector<vector<float>>> all_paths_jaccard(const PathPositionHandleGraph& graph) {

    // cannonical ordering
    vector<path_handle_t> path_handles;
    
    graph.for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });

    if (path_handles.empty()) {
        cerr << "Warning: no paths found" << endl;
    }

    // the output
    vector<vector<float>> out_matrix(path_handles.size());
    vector<string> path_names(path_handles.size());

    // todo: matrix is symmetric -- only need to iterate half!
#pragma omp parallel for
    for (size_t i = 0; i < path_handles.size(); ++i) {

        // row coverage
        unordered_map<path_handle_t, size_t> row_coverage;

        size_t path_len = graph.get_path_length(path_handles[i]);

        // traverse the given path
        graph.for_each_step_in_path(path_handles[i], [&] (step_handle_t step_handle) {
            handle_t handle = graph.get_handle_of_step(step_handle);
            size_t handle_len = graph.get_length(handle);
            vector<step_handle_t> steps = graph.steps_of_handle(handle);
            // add up the coverage of all other paths that share the node
            for (size_t j = 0; j < steps.size(); ++j) {
                path_handle_t step_path = graph.get_path_handle_of_step(steps[j]);
                row_coverage[step_path] += handle_len;
            }
            
            });

        path_names[i] = graph.get_path_name(path_handles[i]);
        out_matrix[i].resize(path_handles.size());
        // compute the overlap with all the other paths
        for (size_t j = 0; j < path_handles.size(); ++j) {
            double jaccard = 2. * (double)row_coverage[path_handles[j]] / (double)(graph.get_path_length(path_handles[j]) + path_len);
            out_matrix[i][j] = jaccard;
        }
    }

    return make_pair(path_names, out_matrix);
}

vector<float> sv_pileup(const PathPositionHandleGraph& graph,
                        const vector<Snarl>& snarls,
                        const string& ref_path,
                        size_t window_size,
                        size_t min_sv_len) {

    path_handle_t ref_path_handle = graph.get_path_handle(ref_path);
    
    // output sv pileup
    vector<float> pileup(graph.get_path_length(ref_path_handle), 0.0);
    
    for (const Snarl& snarl : snarls) {

        // find our snarl endpoints in the grpah
        handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
        handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());

        // find the snarl endpoints in the path
        step_handle_t start_step;
        step_handle_t end_step;

        size_t start_count = 0;
        graph.for_each_step_on_handle(start_handle, [&](step_handle_t step_handle) {
                if (graph.get_path_handle_of_step(step_handle) == ref_path_handle) {
                    start_step = step_handle;
                    ++start_count;
                }
            });

        size_t end_count = 0;
        graph.for_each_step_on_handle(end_handle, [&](step_handle_t step_handle) {
                if (graph.get_path_handle_of_step(step_handle) == ref_path_handle) {
                    end_step = step_handle;
                    ++end_count;
                }
            });

        if (start_count > 1 || end_count > 1) {
            cerr << "Error (svpileup) : cycle not supported in reference path" << endl;
            exit(1);
        }

        if (start_count != 1 || end_count != 1) {
            continue;
        }

        size_t start_pos = graph.get_position_of_step(start_step);
        size_t end_pos = graph.get_position_of_step(end_step);

        // orient snarl on our path
        if (start_pos > end_pos) {
            swap(start_handle, end_handle);
            start_handle = graph.flip(start_handle);
            end_handle = graph.flip(end_handle);
            start_pos = graph.get_position_of_step(start_step);
            end_pos = graph.get_position_of_step(end_step);
        }

        pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs = find_path_traversals(graph, snarl);

        vector<size_t> trav_lengths = get_trav_lengths(graph, path_travs.first);

        // find the reference traversal and index it
        int ref_trav_idx = -1;
        set<handle_t> ref_handles;
        for (int i = 0; i < path_travs.second.size(); ++i) {
            if (graph.get_path_handle_of_step(path_travs.second[i].first) == ref_path_handle) {
                assert(ref_trav_idx < 0);
                ref_trav_idx = i;
                for (step_handle_t ref_step = path_travs.second[i].first; ref_step != graph.get_next_step(path_travs.second[i].second);
                     ref_step = graph.get_next_step(ref_step)) {
                    ref_handles.insert(graph.get_handle_of_step(ref_step));
                }
            }
        }

#ifdef debug
        cerr << "found " << trav_lengths.size() << " traversals with refidx " << ref_trav_idx << " for snarl " << pb2json(snarl) << endl;

        cerr << "the alt paths with different lengths are ";
        for (int i = 0; i < trav_lengths.size(); ++i) {
            if (trav_lengths[i] != trav_lengths[ref_trav_idx]) {
                cerr << "i=" << i << ": " << trav_lengths[i] << " vs ref=" << trav_lengths[ref_trav_idx] << endl;
            }
        }
#endif

        size_t ref_pos = start_pos + graph.get_length(start_handle);

        // process each alt traversal
        for (int i = 0; i < path_travs.second.size(); ++i) {
            // todo: can remove this check to make more general
            if (i != ref_trav_idx && trav_lengths[i] != trav_lengths[ref_trav_idx]) {                
                process_alt(graph, ref_handles, path_travs.first[ref_trav_idx], ref_path_handle, ref_pos, path_travs.second[i], min_sv_len, pileup);
            }
        }

    }
        
    return pileup;

}

void process_alt(const PathPositionHandleGraph& graph,
                 const set<handle_t>& ref_handles,
                 const SnarlTraversal& ref_traversal,
                 path_handle_t ref_path_handle,
                 size_t ref_pos,
                 const pair<step_handle_t, step_handle_t>& alt_traversal,
                 size_t min_sv_len,
                 vector<float>& ref_pileup) {

    bool on_ref = true;
    bool prev_ref = false;
    step_handle_t prev_step;
    step_handle_t bubble_start;

    // hacky nested bubbles
    vector<pair<step_handle_t, step_handle_t>> bubbles;

#ifdef debug
    cerr << "process alt on snarl with ref nodes ";
    for (auto rh : ref_handles) {
        cerr << graph.get_id(rh) << ", ";
    }
    cerr << endl;
#endif
    
    for (step_handle_t alt_step = alt_traversal.first; alt_step != graph.get_next_step(alt_traversal.second);
         alt_step = graph.get_next_step(alt_step)) {

        handle_t handle = graph.get_handle_of_step(alt_step);

        if (on_ref) {
            if (!ref_handles.count(handle)) {
#ifdef debug
                cerr << "leaving ref " << graph.get_id(handle) << endl;
#endif
                // starting new bubble
                bubble_start = prev_step;
                on_ref = false;
            } else if (prev_ref) {
                // ref-ref = deletion
                bubbles.push_back(make_pair(prev_step, alt_step));                    
            }
          prev_ref = true;
        } else {
            if (ref_handles.count(handle)) {
#ifdef debug
                cerr << "back on ref " << graph.get_id(handle) << endl;
#endif
                // finishing bubble
                bubbles.push_back(make_pair(bubble_start, alt_step));
                on_ref = true;
            }
            prev_ref = false;
        }
            
        prev_step = alt_step;
    }

    for (auto& bubble : bubbles) {
        size_t alt_length = 0;
        for (step_handle_t alt_step = bubble.first; alt_step != graph.get_next_step(bubble.second);
             alt_step = graph.get_next_step(alt_step)) {
            alt_length += graph.get_length(graph.get_handle_of_step(alt_step));
        }
        size_t ref_length;
        vector<step_handle_t> start_steps = graph.steps_of_handle(graph.get_handle_of_step(bubble.first));
        size_t ref_start = 0;
        for (step_handle_t step : start_steps) {
            if (graph.get_path_handle_of_step(step) == ref_path_handle) {
                ref_start = graph.get_position_of_step(step) + graph.get_length(graph.get_handle_of_step(step)) - 1;
                break;
            }
        }
        vector<step_handle_t> end_steps = graph.steps_of_handle(graph.get_handle_of_step(bubble.second));
        size_t ref_end = 0;
        for (step_handle_t step : end_steps) {
            if (graph.get_path_handle_of_step(step) == ref_path_handle) {
                ref_end = graph.get_position_of_step(step);
                break;
            }
        }
        ref_length = ref_end - ref_start + 1;
#ifdef debug
        cerr << "found bubble " << ref_start << " " << ref_end << " with alt length " << alt_length << endl;
#endif

        if (ref_length >= alt_length + min_sv_len || alt_length >= ref_length + min_sv_len) {
            for (size_t i = ref_start; i < ref_start + ref_length; ++i) {
                ref_pileup[i] += 1;
            }
        }
    }
}


vector<size_t> get_trav_lengths(const HandleGraph& graph, const vector<SnarlTraversal>& travs) {

    vector<size_t> sizes;

    for (const auto& trav : travs) {
        size_t trav_len = 0;
        for (size_t i = 0; i < trav.visit_size(); ++i) {
            trav_len += graph.get_length(graph.get_handle(trav.visit(i).node_id()));
        }
        sizes.push_back(trav_len);
    }

    return sizes;    
}
    
///////////////////////////////////////////////////////////////////////////////////
// copied from https://github.com/vgteam/vg/blob/master/src/traversal_finder.cpp //
///////////////////////////////////////////////////////////////////////////////////
static pair<unordered_set<nid_t>, set<edge_t> > deep_contents(const Snarl* snarl, const HandleGraph& graph,
                                                                        bool include_boundary_nodes) {
        
    pair<unordered_set<nid_t>, set<edge_t> > to_return;
        
    unordered_set<nid_t> already_stacked;
        
    // initialize stack for DFS traversal of site
    vector<handle_t> stack;

    handle_t start_node = graph.get_handle(snarl->start().node_id());
    handle_t end_node = graph.get_handle(snarl->end().node_id());
        
    // mark the boundary nodes as already stacked so that paths will terminate on them
    already_stacked.insert(graph.get_id(start_node));
    already_stacked.insert(graph.get_id(end_node));
        
    // add boundary nodes as directed
    if (include_boundary_nodes) {
        to_return.first.insert(graph.get_id(start_node));
        to_return.first.insert(graph.get_id(end_node));
    }

    // stack up the nodes one edge inside the snarl from the start
    graph.follow_edges(start_node, snarl->start().backward(), [&](const handle_t& node) {            

            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->start().backward()) {
                to_return.second.insert(graph.edge_handle(node, start_node));
            } else {
                to_return.second.insert(graph.edge_handle(start_node, node));
            }
        });
      
    // stack up the nodes one edge inside the snarl from the end
    graph.follow_edges(end_node, !snarl->end().backward(), [&](const handle_t& node) {
            
            if (!already_stacked.count(graph.get_id(node))) {
                stack.push_back(node);
                already_stacked.insert(graph.get_id(node));
            }
            if (snarl->end().backward()) {
                to_return.second.insert(graph.edge_handle(end_node, node));
            } else {
                to_return.second.insert(graph.edge_handle(node, end_node));
            }
        });
        
    // traverse the snarl with DFS, skipping over any child snarls
    // do not pay attention to valid walks since we also want to discover any tips
    while (stack.size()) {
            
        // pop the top node off the stack
        handle_t node = stack.back();
        stack.pop_back();
            
        // record that this node is in the snarl
        to_return.first.insert(graph.get_id(node));

        graph.follow_edges(node, false, [&] (const handle_t& next_node) {
            edge_t edge = graph.edge_handle(node, next_node);
            to_return.second.insert(edge);
            if (!already_stacked.count(graph.get_id(next_node))) {
                stack.push_back(next_node);
                already_stacked.insert(graph.get_id(next_node));
            }
            });
        
        graph.follow_edges(node, true, [&] (const handle_t& prev_node) {
            edge_t edge = graph.edge_handle(prev_node, node);
            to_return.second.insert(edge);
            if (!already_stacked.count(graph.get_id(prev_node))) {
                stack.push_back(prev_node);
                already_stacked.insert(graph.get_id(prev_node));
            }
            });
    }
        
    return to_return;
}

pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > find_path_traversals(const PathPositionHandleGraph& graph, const Snarl& site) {

    handle_t start_handle = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t end_handle = graph.get_handle(site.end().node_id(), site.end().backward());
    
    vector<step_handle_t> start_steps = graph.steps_of_handle(start_handle);
    vector<step_handle_t> end_steps = graph.steps_of_handle(end_handle);

    pair<unordered_set<nid_t>, set<edge_t> > snarl_contents = deep_contents(&site, graph, true);
    
    // use this to skip paths that don't reach the end node
    unordered_set<path_handle_t> end_path_handles;
    for (const step_handle_t& step : end_steps) {
        end_path_handles.insert(graph.get_path_handle_of_step(step));
    }

#ifdef debug1
    cerr << "Finding traversals of " << pb2json(site) << " using PathTraversalFinder" << endl
         << " - there are " << start_steps.size() << " start_steps, " << end_steps.size() << " end_steps"
         << " and " << end_path_handles.size() << " end_path_handles" << endl;
#endif

    vector<SnarlTraversal> out_travs;
    vector<pair<step_handle_t, step_handle_t> > out_steps;

    unordered_set<path_handle_t> paths;
    
    for (const step_handle_t& start_step : start_steps) {
        path_handle_t start_path_handle = graph.get_path_handle_of_step(start_step);
        // only crawl paths that have a chance of reaching the end
        if ((paths.empty() || paths.count(start_path_handle)) && end_path_handles.count(start_path_handle)) {

            handle_t end_check = end_handle;

#ifdef debug1
            cerr << " - considering path " << graph.get_path_name(start_path_handle) << endl;
#endif
            // try to make a traversal by walking forward
            SnarlTraversal trav;
            bool can_continue = true;
            step_handle_t step = start_step;
            while (can_continue) {
                handle_t handle = graph.get_handle_of_step(step);
                Visit* start_visit = trav.add_visit();
                start_visit->set_node_id(graph.get_id(handle));
                start_visit->set_backward(graph.get_is_reverse(handle));

                can_continue = false;
                if (graph.has_next_step(step) && handle != end_handle) {
                    step_handle_t next_step = graph.get_next_step(step);
                    handle_t next_handle = graph.get_handle_of_step(next_step);
                    if (snarl_contents.first.count(graph.get_id(next_handle)) &&
                        snarl_contents.second.count(graph.edge_handle(handle, next_handle))) {
                        step = next_step;
                        can_continue = true;
                    } 
                }
            }

            if (graph.get_handle_of_step(step) != end_check) {
#ifdef debug1
                cerr << "     - failed to find forward traversal of path " << graph.get_path_name(start_path_handle) << endl;
#endif
                // try to make a traversal by walking backward
                end_check = graph.flip(end_handle);
                
                trav.Clear();
                can_continue = true;
                step = start_step;
                while (can_continue) {
                    handle_t handle = graph.flip(graph.get_handle_of_step(step));
                    
                    Visit* start_visit = trav.add_visit();
                    start_visit->set_node_id(graph.get_id(handle));
                    start_visit->set_backward(graph.get_is_reverse(handle));

                    can_continue = false;
                    if (graph.has_previous_step(step) && handle != end_handle) {
                        step_handle_t prev_step = graph.get_previous_step(step);
                        handle_t prev_handle = graph.flip(graph.get_handle_of_step(prev_step));

                        if (snarl_contents.first.count(graph.get_id(prev_handle)) &&
                            snarl_contents.second.count(graph.edge_handle(handle, prev_handle))) {
                            step = prev_step;
                            can_continue = true;
                        } 
                    }
                }
            }
            if (graph.get_handle_of_step(step) == end_check) {
                out_travs.push_back(trav);
                out_steps.push_back(make_pair(start_step, step));
            } 
        }
    }
    
    return make_pair(out_travs, out_steps);
}

void print_path_coverage(const handlegraph::PathPositionHandleGraph& graph, const string& name_prefix) {
    cout << "path-name" << "\t" << "pct-coverage" << "\t" << "max-gap" << "\t" << "avg-gap" << endl;
    cout << "----------" << "\t" << "------------" << "\t" << "-------" << "\t" << "-------" << endl;

    graph.for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph.get_path_name(path_handle);
            if (path_name.substr(0, name_prefix.length()) == name_prefix) {
                int64_t last_covered = -1;
                int64_t total_covered = 0;
                int64_t total_length = 0;
                vector<int64_t> gaps;
                graph.for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph.get_handle_of_step(step_handle);
                        int64_t length = graph.get_length(handle);
                        bool covered = false;
                        graph.for_each_step_on_handle(handle, [&](step_handle_t step_handle_2) {
                                path_handle_t path_handle_2 = graph.get_path_handle_of_step(step_handle_2);
                                string path_name_2 = graph.get_path_name(path_handle_2);
                                if (path_name_2.substr(0, name_prefix.length()) != name_prefix) {
                                    covered = true;
                                }
                                return !covered;
                            });
                        total_length += length;
                        if (covered) {
                            total_covered += length;
                            if (total_length - last_covered > 1) {
                                gaps.push_back(total_length - last_covered - 1);
                            }
                            last_covered = length - 1;
                        }
                    });

                if (total_length - last_covered > 1) {
                    gaps.push_back(total_length - last_covered - 1);
                }

                int64_t max_gap = 0;
                int64_t total_gap = 0;
                for (int i = 0; i < gaps.size(); ++i) {
                    max_gap = std::max(gaps[i], max_gap);
                    total_gap += gaps[i];
                }

                cout << path_name << "\t"
                    << ((float)total_covered / total_length) << "\t"
                    << max_gap << "\t"
                    << (total_gap / gaps.size()) << endl;
            }   
        });
        
}
