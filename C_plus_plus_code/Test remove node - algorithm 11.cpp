
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath>
#include <set>

using namespace std;

const double INF = numeric_limits<double>::infinity();

// Alias for adjacency matrix
using Matrix = vector<vector<double>>;

// Floyd-Warshall Algorithm
vector<vector<double>> cal_APSP_floyd_warshall(const vector<vector<double>> &distance_matrix) {
    int N = distance_matrix.size();
    vector<vector<double>> shortest_path_matrix = distance_matrix;

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                shortest_path_matrix[i][j] = min(shortest_path_matrix[i][j], shortest_path_matrix[i][k] + shortest_path_matrix[k][j]);
            }
        }
    }
    return shortest_path_matrix;
}

// Dijkstra's Algorithm (One-to-All)
vector<double> dijkstra_one_to_all(const vector<vector<double>> &distance_matrix, int src, const vector<int> &need_up_i, int remove_node) {
    int N = distance_matrix.size();
    vector<double> dist(N, INF);
    vector<bool> visited(N, false);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    dist[src] = 0;
    pq.emplace(0, src);

    int kk = need_up_i.size();
    int counter = 0;
    visited[remove_node] = true;

    while (!pq.empty()) {
        auto [current_dist, u] = pq.top();
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        if (find(need_up_i.begin(), need_up_i.end(), u) != need_up_i.end()) {
            counter++;
            if (counter == kk) return dist;
        }

        for (int v = 0; v < N; ++v) {
            if (!visited[v] && distance_matrix[u][v] > 0) {
                double new_dist = current_dist + distance_matrix[u][v];
                if (new_dist < dist[v]) {
                    dist[v] = new_dist;
                    pq.emplace(new_dist, v);
                }
            }
        }
    }

    return dist;
}


vector<int> cal_candidate_node_list(int i, int j, 
        const vector<vector<double>> &X_APSP_matrix, 
        const vector<vector<double>> &X_distance_matrix) {

    int N = X_distance_matrix.size();
    vector<int> remaining_list;

    for (int k = 0; k < N; ++k) {
        if (k != i && k != j) {
            remaining_list.push_back(k);
        }
    }

    vector<int> candidate_node_list;
    candidate_node_list.push_back(i);

    for (int t : remaining_list) {
        if (X_APSP_matrix[i][j] >= X_APSP_matrix[i][t] + X_APSP_matrix[t][j]) {
            candidate_node_list.push_back(t);
        }
    }

    candidate_node_list.push_back(j);
    return candidate_node_list;
}

vector<vector<double>> output_small_matrix(const vector<vector<double>>& X_distance_matrix, const vector<int>& candidate_node_list) {
    int size = candidate_node_list.size();
    vector<vector<double>> small_matrix(size, vector<double>(size, 0));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            small_matrix[i][j] = X_distance_matrix[candidate_node_list[i]][candidate_node_list[j]];
        }
    }
    return small_matrix;
}

typedef pair<double, int> Pair;
vector<int> dijkstra_show_path(const vector<vector<double>>& adj_matrix, int start, int end) {
    int n = adj_matrix.size();
    vector<double> distances(n, INF);
    vector<int> previous_nodes(n, -1);
    distances[start] = 0;
    priority_queue<Pair, vector<Pair>, greater<>> priority_queue;
    priority_queue.emplace(0, start);

    while (!priority_queue.empty()) {
        auto [current_distance, current_node] = priority_queue.top();
        priority_queue.pop();

        if (current_node == end) {
            break;
        }

        if (current_distance > distances[current_node])
            continue;

        for (int neighbor = 0; neighbor < n; ++neighbor) {
            double weight = adj_matrix[current_node][neighbor];
            if (weight > 0) {
                double distance = current_distance + weight;
                if (distance < distances[neighbor]) {
                    distances[neighbor] = distance;
                    previous_nodes[neighbor] = current_node;
                    priority_queue.emplace(distance, neighbor);
                }
            }
        }
    }

    vector<int> path;
    for (int current = end; current != -1; current = previous_nodes[current]) {
        path.push_back(current);
    }
    reverse(path.begin(), path.end());

    return (path.front() == start) ? path : vector<int>();
}



vector<int> dijkstra_cal_previous_nodes(const vector<vector<double>>& adj_matrix, int start) {
    int n = adj_matrix.size();
    vector<double> distances(n, INF);
    vector<int> previous_nodes(n, -1);
    distances[start] = 0;
    priority_queue<Pair, vector<Pair>, greater<>> priority_queue;
    priority_queue.emplace(0, start);

    while (!priority_queue.empty()) {
        auto [current_distance, current_node] = priority_queue.top();
        priority_queue.pop();

        if (current_distance > distances[current_node])
            continue;

        for (int neighbor = 0; neighbor < n; ++neighbor) {
            double weight = adj_matrix[current_node][neighbor];
            if (weight > 0) {
                double distance = current_distance + weight;
                if (distance < distances[neighbor]) {
                    distances[neighbor] = distance;
                    previous_nodes[neighbor] = current_node;
                    priority_queue.emplace(distance, neighbor);
                }
            }
        }
    }
    return previous_nodes;
}


vector<int> translate_path(const vector<int> &path, const vector<int> &candidate_node_list) {
    vector<int> new_path;
    for (int i : path) {
        new_path.push_back(candidate_node_list[i]);
    }
    return new_path;
}

bool check_if_temp_path_already_in(const vector<vector<int>> &all_paths_list, const vector<int> &temp_path) {
    set<int> aa(temp_path.begin(), temp_path.end());
    for (const auto &bb : all_paths_list) {
        set<int> bb_set(bb.begin(), bb.end());
        if (aa == bb_set) {
            return true;
        }
    }
    return false;
}

vector<int> dijkstra_show_path_from_previous_nodes(const vector<int> &previous_nodes, int start, int end) {
    vector<int> path;
    int current = end;
    while (current != -1) {
        path.push_back(current);
        current = previous_nodes[current];
    }
    reverse(path.begin(), path.end());
    return (path.front() == start) ? path : vector<int>();
}

vector<vector<int>> cal_all_paths_warm_start(int i, int j, 
        const vector<vector<double>> &X_APSP_matrix, 
        const vector<vector<double>> &X_distance_matrix) {

    vector<int> candidate_node_list = cal_candidate_node_list(i, j, X_APSP_matrix, X_distance_matrix);
    vector<vector<double>> small_matrix = output_small_matrix(X_distance_matrix, candidate_node_list);
    
    vector<vector<int>> all_paths_list;
    all_paths_list.push_back(dijkstra_show_path(small_matrix, 0, candidate_node_list.size() - 1));

    int K = candidate_node_list.size();
    vector<int> temp_list;
    
    for (int q = 0; q < K; ++q) {
        if (find(all_paths_list[0].begin(), all_paths_list[0].end(), q) == all_paths_list[0].end()) {
            temp_list.push_back(q);
        }
    }

    vector<int> previous_nodes = dijkstra_cal_previous_nodes(small_matrix, 0);

    for (int m : temp_list) {
        vector<int> temp_path1 = dijkstra_show_path_from_previous_nodes(previous_nodes, 0, m);
        vector<int> temp_path2 = dijkstra_show_path(small_matrix, m, K - 1);

        vector<int> temp_path = temp_path1;
        temp_path.pop_back(); // Remove the last element to avoid duplication
        temp_path.insert(temp_path.end(), temp_path2.begin(), temp_path2.end());

        if (!check_if_temp_path_already_in(all_paths_list, temp_path)) {
            all_paths_list.push_back(temp_path);
        }
    }

    for (auto &path : all_paths_list) {
        path = translate_path(path, candidate_node_list);
    }

    return all_paths_list;
}

vector<int> cal_key_node_list_one(int i, int j, const vector<vector<double>>& X_APSP_matrix, const vector<vector<double>>&  X_distance_matrix) {
    vector<int> key_node_list;
    
    if (i == j) {
        return key_node_list;
    } else {
        vector<vector<int>> all_paths_list = cal_all_paths_warm_start(i, j, X_APSP_matrix, X_distance_matrix);
        set<int> key_node_set(all_paths_list[0].begin(), all_paths_list[0].end());
        
        for (size_t m = 1; m < all_paths_list.size(); ++m) {
            set<int> path_set(all_paths_list[m].begin(), all_paths_list[m].end());
            set<int> intersection;
            
            set_intersection(key_node_set.begin(), key_node_set.end(), path_set.begin(), path_set.end(),
                             inserter(intersection, intersection.begin()));
            key_node_set = intersection;
        }
        
        key_node_list.assign(key_node_set.begin(), key_node_set.end());
        key_node_list.erase(remove(key_node_list.begin(), key_node_list.end(), i), key_node_list.end());
        key_node_list.erase(remove(key_node_list.begin(), key_node_list.end(), j), key_node_list.end());
    }
    
    return key_node_list;
}

// Need Update List Calculation
pair<vector<vector<int>>, int> cal_need_update_list_mix(int remove_node, const vector<vector<double>> &X_APSP_matrix, const vector<vector<double>> &X_distance_matrix) {
    int N = X_APSP_matrix.size();
    vector<int> remaining_list;
    for (int i = 0; i < N; ++i) {
        if (i != remove_node) remaining_list.push_back(i);
    }

    vector<vector<int>> need_update_list(N);
    vector<int> how_many_counter(N, 0);

    for (int i : remaining_list) {
        for (int j : remaining_list) {
            if (i != j && X_APSP_matrix[i][j] >= X_APSP_matrix[i][remove_node] + X_APSP_matrix[remove_node][j]) {
                    vector<int> key_node_list = cal_key_node_list_one(i, j, X_APSP_matrix, X_distance_matrix);
                    if (find(key_node_list.begin(), key_node_list.end(), remove_node) != key_node_list.end()) {
                        need_update_list[i].push_back(j);
                        how_many_counter[i]++;
                    }
            }
        }
    }

    int how_many = count_if(how_many_counter.begin(), how_many_counter.end(), [](int x) { return x > 0; });
    return {need_update_list, how_many};
}

// Generate New Distance Matrix
vector<vector<double>> generate_new_distance_matrix(const vector<vector<double>> &X_distance_matrix, int remove_node) {
    int N = X_distance_matrix.size();
    vector<vector<double>> X_distance_matrix_new = X_distance_matrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == remove_node || j == remove_node) {
                X_distance_matrix_new[i][j] = INF;
            }
        }
    }

    X_distance_matrix_new[remove_node][remove_node] = 0;
    return X_distance_matrix_new;
}


// Calculate New APSP Matrix using Dijkstra
vector<vector<double>> cal_new_APSP_matrix_warm(const vector<vector<double>> &X_distance_matrix,
                                                    const vector<vector<double>> &X_APSP_matrix,
                                                    int remove_node) {
    int N = X_distance_matrix.size();
    vector<vector<double>> X_APSP_matrix_new_warm = X_APSP_matrix;
    auto [need_update_list, _] = cal_need_update_list_mix(remove_node, X_APSP_matrix, X_distance_matrix);

    auto X_distance_matrix_new = generate_new_distance_matrix(X_distance_matrix, remove_node);

    for (int i = 0; i < N; ++i) {
        if (i == remove_node) continue;
        if (need_update_list[i].size() > 0) {
            auto temp = dijkstra_one_to_all(X_distance_matrix_new, i, need_update_list[i], remove_node);
            for (int j : need_update_list[i]) {
                X_APSP_matrix_new_warm[i][j] = temp[j];
            }
        } 
    }


    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == remove_node || j == remove_node) {
                X_APSP_matrix_new_warm[i][j] = INF;
            }
        }
    }

    X_APSP_matrix_new_warm[remove_node][remove_node] = 0;
    return X_APSP_matrix_new_warm;
}

// Load Distance Matrix from CSV
Matrix load_distance_matrix(const string& filename) {
    ifstream file(filename);
    string line;
    Matrix matrix;

    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> row;
        string cell;
        while (getline(ss, cell, ',')) {
            row.push_back(stod(cell));
        }
        matrix.push_back(row);
    }
    return matrix;
}

// Main Function
int main() {


    int remove_node = 3;
    string file_path = "./X_100_distance_matrix.csv";

    vector<vector<double>> X_distance_matrix = load_distance_matrix(file_path);
    vector<vector<double>> X_APSP_matrix = cal_APSP_floyd_warshall(X_distance_matrix);

    vector<vector<double>> new_X_distance_matrix = generate_new_distance_matrix(X_distance_matrix, remove_node);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<double>> X_APSP_matrix_new_warm = cal_new_APSP_matrix_warm(X_distance_matrix, X_APSP_matrix, remove_node);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> time_used1 = end - start;

    cout << "time_used1 warm: " << time_used1.count() << endl;

    start = chrono::high_resolution_clock::now();
    vector<vector<double>> X_APSP_matrix_new_floyd_warshall = cal_APSP_floyd_warshall(new_X_distance_matrix);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> time_used2 = end - start;
 
    cout << "time_used2 floyd: " << time_used2.count() << endl;

    cout << "Time Ratio: " << round(time_used1.count() / time_used2.count() * 100) / 100 << endl;

    for (int i = 0; i < 10; ++i) {
        cout << X_APSP_matrix_new_warm[10][i] << " ";
    }
    cout << endl;

    cout << "Matrix size: " << X_APSP_matrix_new_warm.size() << "x" << X_APSP_matrix_new_warm[0].size() << endl;

    // Compare results
    bool matrices_equal = true;
    for (size_t i = 0; i < X_APSP_matrix_new_warm.size(); ++i) {
        for (size_t j = 0; j < X_APSP_matrix_new_warm[i].size(); ++j) {
            if (fabs(X_APSP_matrix_new_warm[i][j] - X_APSP_matrix_new_floyd_warshall[i][j]) > 1e-6) {
                matrices_equal = false;
                break;
            }
        }
    }

    cout << "Matrices are " << (matrices_equal ? "equal" : "not equal") << endl;
 
}
