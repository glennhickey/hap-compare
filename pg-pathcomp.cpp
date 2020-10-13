#include <unordered_map>
#include "pg-pathcomp.hpp"

using namespace std;
using namespace handlegraph;

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
