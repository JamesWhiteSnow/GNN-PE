#pragma once


#include "graph/graph.h"
class GraphOperations {
public:
    static void getKCore(const Static_Graph* graph, int* core_table);
    static void match_bfs(int* col_ptrs, int* col_ids, int* match, int* row_match, int* visited,
        int* queue, int* previous, int n, int m);
    static void bfsTraversal(const Static_Graph* graph, VertexID root_vertex, TreeNode*& tree, VertexID*& bfs_order);
    static void dfsTraversal(TreeNode* tree, VertexID root_vertex, ui node_num, VertexID*& dfs_order);
private:
    static void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m);
    static void dfs(TreeNode* tree, VertexID cur_vertex, VertexID* dfs_order, ui& count);
};

