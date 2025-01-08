#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath> 


using namespace std;



// Alias for adjacency matrix
using Matrix = vector<vector<double>>;

// Constants
const double INF = numeric_limits<double>::infinity();

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

// Need Update List Calculation
pair<vector<vector<int>>, int> cal_need_update_list(int remove_node, const vector<vector<double>> &X_APSP_matrix) {
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
                need_update_list[i].push_back(j);
                how_many_counter[i]++;
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

// Calculate New APSP Matrix using Dijkstra
vector<vector<double>> cal_new_APSP_matrix_warm(const vector<vector<double>> &X_distance_matrix,
                                                    const vector<vector<double>> &X_APSP_matrix,
                                                    int remove_node) {
    int N = X_distance_matrix.size();
    vector<vector<double>> X_APSP_matrix_new_warm = X_APSP_matrix;
    auto [need_update_list, _] = cal_need_update_list(remove_node, X_APSP_matrix);

    auto X_distance_matrix_new = generate_new_distance_matrix(X_distance_matrix, remove_node);

    for (int i = 0; i < N; ++i) {
        if (i == remove_node) continue;
        if (need_update_list[i].size() > 0) {
            auto temp = dijkstra_one_to_all(X_distance_matrix_new, i, need_update_list[i], remove_node);
            for (int j : need_update_list[i]) {
                X_APSP_matrix_new_warm[i][j] = temp[j];
            }
        } 
        // else {
        //     for (int j : need_update_list[i]) {
        //         X_APSP_matrix_new_warm[i][j] = dijkstra_one_to_one(X_distance_matrix_new, i, j);
        //     }
        // }

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

// Helper function to get the minimum of two values
template <typename T>
T get_min(T a, T b) {
    return (a < b) ? a : b;
}

// Helper function to get the minimum element in a vector
template <typename T>
T get_min_element(const vector<T>& vec) {
    return *min_element(vec.begin(), vec.end());
}

double shortest_path_n_to_r_new(const vector<vector<double>>& distance_matrix,
                                const vector<vector<double>>& shortest_path_matrix,
                                const vector<int>& remaining_list, int r, int remove_node) {
    vector<double> max_jump_list;
    for (int t : remaining_list) {
        double m_jump = distance_matrix[remove_node][t] + shortest_path_matrix[t][r];
        max_jump_list.push_back(m_jump);
    }
    return get_min_element(max_jump_list);
}

double shortest_path_r_to_n_new(const vector<vector<double>>& distance_matrix,
                                const vector<vector<double>>& shortest_path_matrix,
                                const vector<int>& remaining_list, int r, int remove_node) {
    vector<double> max_jump_list;
    for (int t : remaining_list) {
        double m_jump = shortest_path_matrix[r][t] + distance_matrix[t][remove_node];
        max_jump_list.push_back(m_jump);
    }
    return get_min_element(max_jump_list);
}

double update_shortest_path_ij_new(const vector<vector<double>>& distance_matrix, 
                                    const vector<vector<double>>& shortest_path_matrix,
                                    int remove_node, int i, int j) {
    double m1 = shortest_path_matrix[i][j];
    double m2 = shortest_path_matrix[i][remove_node] + shortest_path_matrix[remove_node][j];
    return get_min(m1, m2);
}

void cal_n_shortest_path_new(const vector<vector<double>>& distance_matrix, vector<vector<double>>& shortest_path_matrix,
                              int remove_node, const vector<int>& remaining_list) {
    for (int r : remaining_list) {
        shortest_path_matrix[remove_node][r] = shortest_path_n_to_r_new(distance_matrix, shortest_path_matrix, remaining_list, r, remove_node);
        shortest_path_matrix[r][remove_node] = shortest_path_r_to_n_new(distance_matrix, shortest_path_matrix, remaining_list, r, remove_node);
    }

    for (int i : remaining_list) {
        for (int j : remaining_list) {
            if (i < j) {
                shortest_path_matrix[i][j] = update_shortest_path_ij_new(distance_matrix, shortest_path_matrix, remove_node, i, j);
                shortest_path_matrix[j][i] = update_shortest_path_ij_new(distance_matrix, shortest_path_matrix, remove_node, j, i);
            }
        }
    }
}

int cal_which_to_remove(const vector<vector<double>>& X_APSP_matrix, const vector<int>& edge_nodes) {

    auto [need_update_list0, how_many0] = cal_need_update_list(edge_nodes[0], X_APSP_matrix);
    auto [need_update_list1, how_many1] = cal_need_update_list(edge_nodes[1], X_APSP_matrix);

    return (how_many0 < how_many1) ? edge_nodes[0] : edge_nodes[1];
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

vector<vector<double>> cal_APSP_matrix_after_graph_update(const vector<vector<double>>& X_distance_matrix,
                                                           const vector<vector<double>>& X_APSP_matrix,
                                                           const vector<int>& edge_nodes,
                                                           const vector<vector<double>>& X_distance_matrix_new) {
    int N = X_distance_matrix.size();
    
    int remove_node = cal_which_to_remove(X_APSP_matrix, edge_nodes);

    vector<vector<double>> X_APSP_matrix_updated = cal_new_APSP_matrix_warm(X_distance_matrix, X_APSP_matrix, remove_node);  // You need to implement this

    vector<int> remaining_list;
    for (int i = 0; i < N; i++) {
        if (i != remove_node) {
            remaining_list.push_back(i);
        }
    }

    cal_n_shortest_path_new(X_distance_matrix_new, X_APSP_matrix_updated, remove_node, remaining_list);

    return X_APSP_matrix_updated;
}



// Entry point
int main() {

    string file_path = "./X_100_distance_matrix.csv";

    vector<vector<double>> X_distance_matrix = load_distance_matrix(file_path);
    vector<vector<double>> X_APSP_matrix = cal_APSP_floyd_warshall(X_distance_matrix);

    vector<int> edge_nodes = {3, 9};
    double edge_weight = 50;

    vector<vector<double>> X_distance_matrix_new = X_distance_matrix;

    X_distance_matrix_new[edge_nodes[0]][edge_nodes[1]]= edge_weight;

    auto start = chrono::high_resolution_clock::now();
    Matrix X_APSP_matrix_new_warm = cal_APSP_matrix_after_graph_update(X_distance_matrix, X_APSP_matrix, edge_nodes, X_distance_matrix_new);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> time_used1_warm = end - start;

    start = chrono::high_resolution_clock::now();
    Matrix X_APSP_matrix_new_floyd_warshall = cal_APSP_floyd_warshall(X_distance_matrix_new);
    end = chrono::high_resolution_clock::now();

    chrono::duration<double> time_used2_floyd = end - start;

    cout << "Time Used (warm): " << time_used1_warm.count() << "s\n";
    cout << "Time Used (Floyd-Warshall): " << time_used2_floyd.count() << "s\n";
 

    cout << "Time Ratio: " << round(time_used1_warm.count() / time_used2_floyd.count() * 100) / 100 << endl;

    // for (int i = 0; i < 10; ++i) {
    //     cout << X_APSP_matrix_new_warm[0][i] << " ";
    // }
    // cout << endl;

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

    return 0;
}
