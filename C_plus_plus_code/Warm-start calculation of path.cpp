#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;

// Alias for matrix
using Matrix = vector<vector<double>>;
using Path = vector<int>;

// Floyd-Warshall Algorithm
Matrix floyd_warshall(const Matrix& distance_matrix) {
    int N = distance_matrix.size();
    Matrix shortest_path_matrix = distance_matrix;

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                shortest_path_matrix[i][j] = min(shortest_path_matrix[i][j], 
                                                shortest_path_matrix[i][k] + shortest_path_matrix[k][j]);
            }
        }
    }
    return shortest_path_matrix;
}

// Dijkstra's Algorithm to find the shortest path
Path dijkstra_show_path(const Matrix& adj_matrix, int start, int end) {
    int n = adj_matrix.size();
    vector<double> distances(n, numeric_limits<double>::infinity());
    vector<int> previous_nodes(n, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    distances[start] = 0;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [current_distance, current_node] = pq.top();
        pq.pop();

        if (current_node == end) break;
        if (current_distance > distances[current_node]) continue;

        for (int neighbor = 0; neighbor < n; ++neighbor) {
            double weight = adj_matrix[current_node][neighbor];
            if (weight > 0) {
                double distance = current_distance + weight;
                if (distance < distances[neighbor]) {
                    distances[neighbor] = distance;
                    previous_nodes[neighbor] = current_node;
                    pq.push({distance, neighbor});
                }
            }
        }
    }

    Path path;
    for (int current = end; current != -1; current = previous_nodes[current]) {
        path.push_back(current);
    }
    reverse(path.begin(), path.end());

    return (path.front() == start) ? path : Path();
}

// Translate path
Path translate_path(const Path& path, const vector<int>& candidate_node_list) {
    Path new_path;
    for (int index : path) {
        new_path.push_back(candidate_node_list[index]);
    }
    return new_path;
}

// Create small matrix based on candidate nodes
Matrix output_small_matrix(const Matrix& dist_matrix, const vector<int>& candidate_node_list) {
    int size = candidate_node_list.size();
    Matrix small_matrix(size, vector<double>(size, 0.0));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            small_matrix[i][j] = dist_matrix[candidate_node_list[i]][candidate_node_list[j]];
        }
    }

    return small_matrix;
}

// Calculate candidate node list
vector<int> cal_candidate_node_list(int i, int j, const Matrix& shortest_path_matrix, const Matrix& dist_matrix) {
    int N = dist_matrix.size();
    vector<int> remaining_list;
    for (int k = 0; k < N; ++k) {
        if (k != i && k != j) {
            remaining_list.push_back(k);
        }
    }

    vector<int> candidate_node_list = {i};
    for (int t : remaining_list) {
        if (shortest_path_matrix[i][j] >= shortest_path_matrix[i][t] + shortest_path_matrix[t][j]) {
            candidate_node_list.push_back(t);
        }
    }
    candidate_node_list.push_back(j);

    return candidate_node_list;
}

// Smart path calculation
Path cal_path_smart(int i, int j, const Matrix& shortest_path_matrix, const Matrix& dist_matrix) {
    vector<int> candidate_node_list = cal_candidate_node_list(i, j, shortest_path_matrix, dist_matrix);
    Matrix small_matrix = output_small_matrix(dist_matrix, candidate_node_list);
    Path path = dijkstra_show_path(small_matrix, 0, candidate_node_list.size() - 1);
    return translate_path(path, candidate_node_list);
}

// Load distance matrix from file
Matrix load_distance_matrix(const string& filename) {
    ifstream file(filename);
    string line;
    Matrix matrix;

    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> row;
        double value;
        while (ss >> value) {
            row.push_back(value);
            if (ss.peek() == ',') ss.ignore();
        }
        matrix.push_back(row);
    }

    return matrix;
}

int main() {
    string filepath = "./X_100_distance_matrix.csv";
    Matrix dist_matrix = load_distance_matrix(filepath);

    Matrix shortest_path_matrix = floyd_warshall(dist_matrix);
    int N = dist_matrix.size();

    cout << "N: " << N << endl;

    int i = 2, j = 5;
    Path path1 = dijkstra_show_path(dist_matrix, i, j);
    Path path2 = cal_path_smart(i, j, shortest_path_matrix, dist_matrix);

    cout << "Path1: ";
    for (int node : path1) cout << node << " ";
    cout << endl;

    cout << "Path2: ";
    for (int node : path2) cout << node << " ";
    cout << endl;

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == 2) {
                cal_path_smart(i, j, shortest_path_matrix, dist_matrix);
            }
        }
    }
    auto end = chrono::high_resolution_clock::now();
    auto time_used1_warm = chrono::duration<double>(end - start).count();
    cout << "time_used1_warm: " << time_used1_warm << " seconds" << endl;

    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == 2) {
                dijkstra_show_path(dist_matrix, i, j);
            }
        }
    }
    end = chrono::high_resolution_clock::now();
    auto time_used2_cold = chrono::duration<double>(end - start).count();
    cout << "time_used2_cold: " << time_used2_cold << " seconds" << endl;

    cout << "Ratio: " << time_used1_warm / time_used2_cold << endl;

    return 0;
}
