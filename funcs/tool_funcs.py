import numpy as np
import random
import time
import heapq

# create distance matrix of a directed graph

def create_distance_matrix(N):
    edges = []
    for u in range(N):  
        for v in range(N): 
            if u != v:
                w = random.randint(1, 1000)            
                edges.append((u, v, w))    
    
    dist_matrix = np.zeros((N, N))
    for u, v, w in edges:
        dist_matrix[u, v] = w
        
    
    return dist_matrix.astype('float64')

def floyd_warshall(distance_matrix):
    N = len(distance_matrix)
    shortest_path_matrix = distance_matrix.copy()
   
    for k in range(N):
        for i in range(N):
            for j in range(N):
                shortest_path_matrix[i][j] = min(shortest_path_matrix[i][j], shortest_path_matrix[i][k] + shortest_path_matrix[k][j])
    return shortest_path_matrix

def floyd_warshall_minimax_path(distance_matrix):
    N = len(distance_matrix)
    shortest_path_matrix = distance_matrix.copy()
   
    for k in range(N):
        for i in range(N):
            for j in range(N):
                shortest_path_matrix[i][j] = min(shortest_path_matrix[i][j], max(shortest_path_matrix[i][k], shortest_path_matrix[k][j]))
    return shortest_path_matrix

def floyd_warshall_widest_path(distance_matrix):
    N = len(distance_matrix)
    shortest_path_matrix = distance_matrix.copy()
   
    for k in range(N):
        for i in range(N):
            for j in range(N):
                shortest_path_matrix[i][j] = max(shortest_path_matrix[i][j], min(shortest_path_matrix[i][k], shortest_path_matrix[k][j]))
    return shortest_path_matrix





