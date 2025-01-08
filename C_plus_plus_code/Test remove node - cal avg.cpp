
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <chrono>
#include <fstream>
#include <sstream>
#include <cmath>

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

// Function to create a distance matrix
vector<vector<double>> create_distance_matrix(int N) {
    // Initialize random seed
    srand(time(nullptr));

    // Create a distance matrix filled with zeros
    vector<vector<double>> dist_matrix(N, vector<double>(N, 0.0));

    // Fill the matrix with random weights
    for (int u = 0; u < N; ++u) {
        for (int v = 0; v < N; ++v) {
            if (u != v) {
                int w = rand() % 1000 + 1; // Random weight between 1 and 1000
                dist_matrix[u][v] = static_cast<double>(w);
            }
        }
    }

    return dist_matrix;
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
 
    vector<int> N_list = {100, 200, 300};

    for (int N : N_list){

    cout << "N: " << N << endl;

    vector<double> uuu;

    for (int j = 0; j < 20; ++j) {

    int remove_node = 3;
 
    // Create the initial distance matrix
    vector<vector<double>> X_distance_matrix = create_distance_matrix(N);
 
    vector<vector<double>> X_APSP_matrix = cal_APSP_floyd_warshall(X_distance_matrix);

    vector<vector<double>> new_X_distance_matrix = generate_new_distance_matrix(X_distance_matrix, remove_node);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<double>> X_APSP_matrix_new_warm = cal_new_APSP_matrix_warm(X_distance_matrix, X_APSP_matrix, remove_node);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> time_used1 = end - start;

    start = chrono::high_resolution_clock::now();
    vector<vector<double>> X_APSP_matrix_new_floyd_warshall = cal_APSP_floyd_warshall(new_X_distance_matrix);
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> time_used2 = end - start;
 
    // Calculate ratio and add to the list
    double kkk = round(time_used1.count() / time_used2.count() * 100) / 100;
    uuu.push_back(kkk);
    }

    // Calculate and print the mean
    double sum = 0;
    for (double value : uuu) {
        sum += value;
    }
    double mean = sum / uuu.size();
    cout << "Mean: " << mean << endl;

    // Print all values
    cout << "Values: ";
    for (double value : uuu) {
        cout << value << ", ";
    }
    cout << endl;
    }
    return 0;   
 
}
