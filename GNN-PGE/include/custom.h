#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./blockfile/blk_file.h"
#include "./blockfile/cache.h"
#include "./linlist/linlist.h"
#include "./rtree/rtree_cmd.h"
#include "rand.h"
#include "cdf.h"
#include "./graph/graph.h"

#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <set>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <chrono>
#include <limits>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iomanip>

#define NANOSECTOSEC(elapsed_time) ((elapsed_time) / (double)1000000000)
#define NANOSECTOMSEC(elapsed_time) ((elapsed_time) / (double)1000000)

using namespace std;
using namespace chrono;

const double epsilon = 1e-6;

ui MAX_LIMIT = UINT_MAX;

ui vde_dim = 2;			  
ui path_length = 1 + 1;  
ui pde_dim = vde_dim * (path_length);
ui partition_num = 5;

void dfs(ui node, ui depth, vector<ui> &path, const Static_Graph *static_data_graph, vector<vector<ui>> &all_paths)
{
	if (depth == path_length)
	{
		all_paths.push_back(path);
		return;
	}

	ui count = 0;
	const ui *neighbors = static_data_graph->getVertexNeighbors(node, count);
	for (ui i = 0; i < count; i++)
	{
		if (find(path.begin(), path.end(), neighbors[i]) == path.end())
		{
			path.push_back(neighbors[i]);
			dfs(neighbors[i], depth + 1, path, static_data_graph, all_paths);
			path.pop_back();
		}
	}
}

class Vertex
{
public:
	ui vid;
	ui label;
	ui degree;

	vector<double> x;
	vector<double> nx;
	vector<double> vde;

	vector<double> path_group;
	vector<double> path_label_group;

	double key;
};

class Auxiliary_Index
{
public:
	double key;
	ui degree;
	vector<double> label_mbr;

	Auxiliary_Index()
	{
		this->key = 0;
		this->degree = 0;
		this->label_mbr.resize(2 * pde_dim, 0);
	}
};

class Query_Plan
{
public:
	double key;
	vector<Vertex> query_vertices;

	Query_Plan() {}

	Query_Plan(vector<Vertex> query_vertices)
	{
		this->query_vertices = query_vertices;
		calculate_key();
	}

	void calculate_key()
	{
		key = -1e20;
		for (ui i = 0; i < query_vertices.size(); i++)
		{
			if (query_vertices[i].key > key)
			{
				key = query_vertices[i].key;
			}
		}
	}
};

class Partition
{
public:
	vector<Vertex> vertices;
	ui node_num;
	vector<Auxiliary_Index> auxiliary_index;

	RTree *rtree_index;

	Partition(const vector<Vertex> &data_vertices, const string partition_path, const vector<ui> &partition_vertices, double &index_time)
	{
		for (ui node : partition_vertices)
		{
			this->vertices.push_back(data_vertices[node]);
		}

		string index_path = partition_path + "index.dat";
		char *index_file_name = new char[index_path.length() + 1];
		strcpy(index_file_name, index_path.c_str());

		fstream in(index_file_name, ios_base::in);
		bool generate_tree = true;
		if (!in)
		{
			in.close();
			generate_tree = true;
		}
		else
		{
			generate_tree = false;
			in.close();
		}

		if (generate_tree)
		{
			Entry *node;
			rtree_index = new RTree(index_file_name, 4096, NULL, pde_dim);
			for (ui i = 0; i < this->vertices.size(); i++)
			{
				node = new Entry(pde_dim, NULL);
				node->son = i;
				for (ui j = 0; j < pde_dim; j++)
				{
					node->bounces[2 * j] = vertices[i].path_group[2 * j];
					node->bounces[2 * j + 1] = vertices[i].path_group[2 * j + 1];
				}
				double ind_st = clock();

				rtree_index->insert(node);

				double ind_ed = clock();
				index_time += (ind_ed - ind_st) / CLOCKS_PER_SEC;
			}
			delete rtree_index;
		}
		rtree_index = new RTree(index_file_name, NULL);

		rtree_index->load_root();
		this->node_num = rtree_index->root_ptr->get_num_of_node();
		delete[] index_file_name;

		auxiliary_index.resize(node_num);
		build_auxiliary_index(rtree_index->root_ptr, rtree_index->root);
	}

	void build_auxiliary_index(RTNode *rtn, ui current_ID)
	{
		if (rtn->level == 0)
		{
			for (ui i = 0; i < rtn->num_entries; i++)
			{
				ui data_ID = rtn->entries[i].son;

				if (i == 0)
				{
					auxiliary_index[current_ID].degree = vertices[data_ID].degree;
				}
				else
				{
					if (auxiliary_index[current_ID].degree < vertices[data_ID].degree)
					{
						auxiliary_index[current_ID].degree = vertices[data_ID].degree;
					}
				}

				for (ui j = 0; j < pde_dim; j++)
				{
					if (i == 0)
					{
						auxiliary_index[current_ID].label_mbr[2 * j] = vertices[data_ID].path_label_group[2 * j];
						auxiliary_index[current_ID].label_mbr[2 * j + 1] = vertices[data_ID].path_label_group[2 * j + 1];
					}
					else
					{
						if (auxiliary_index[current_ID].label_mbr[2 * j] > vertices[data_ID].path_label_group[2 * j])
						{
							auxiliary_index[current_ID].label_mbr[2 * j] = vertices[data_ID].path_label_group[2 * j];
						}
						if (auxiliary_index[current_ID].label_mbr[2 * j + 1] < vertices[data_ID].path_label_group[2 * j + 1])
						{
							auxiliary_index[current_ID].label_mbr[2 * j + 1] = vertices[data_ID].path_label_group[2 * j + 1];
						}
					}
				}
			}
		}
		else
		{
			for (ui i = 0; i < rtn->num_entries; i++)
			{
				RTNode *child = new RTNode(rtn->my_tree, rtn->entries[i].son);
				build_auxiliary_index(child, rtn->entries[i].son);

				auxiliary_index[rtn->entries[i].son].key = 0;
				for (ui j = 0; j < pde_dim; j++)
				{
					auxiliary_index[rtn->entries[i].son].key -= rtn->entries[i].bounces[2 * j + 1];
				}

				ui node_ID = rtn->entries[i].son;

				for (ui j = 0; j < path_length; j++)
				{
					if (i == 0)
					{
						auxiliary_index[current_ID].degree = auxiliary_index[node_ID].degree;
					}
					else
					{
						if (auxiliary_index[current_ID].degree < auxiliary_index[node_ID].degree)
						{
							auxiliary_index[current_ID].degree = auxiliary_index[node_ID].degree;
						}
					}
				}

				for (ui j = 0; j < pde_dim; j++)
				{
					if (i == 0)
					{
						auxiliary_index[current_ID].label_mbr[2 * j] = auxiliary_index[node_ID].label_mbr[2 * j];
						auxiliary_index[current_ID].label_mbr[2 * j + 1] = auxiliary_index[node_ID].label_mbr[2 * j + 1];
					}
					else
					{
						if (auxiliary_index[current_ID].label_mbr[2 * j] > auxiliary_index[node_ID].label_mbr[2 * j])
						{
							auxiliary_index[current_ID].label_mbr[2 * j] = auxiliary_index[node_ID].label_mbr[2 * j];
						}
						if (auxiliary_index[current_ID].label_mbr[2 * j + 1] < auxiliary_index[node_ID].label_mbr[2 * j + 1])
						{
							auxiliary_index[current_ID].label_mbr[2 * j + 1] = auxiliary_index[node_ID].label_mbr[2 * j + 1];
						}
					}
				}
				delete child;
			}
		}
	}

	double query(vector<set<ui>> &candidates, Query_Plan Q)
	{
		Heap *hp = new Heap();
		hp->init(rtree_index->dimension);

		HeapEntry *he = new HeapEntry();
		rtree_index->load_root();
		he->son1 = rtree_index->root;
		he->key = -1e20;
		he->level = 1;
		hp->insert(he);
		delete he;

		vector<Query_Plan> Q_map(node_num);
		Q_map[rtree_index->root] = Q;

		RTNode *child;

		ui traversal_index_nodes = 0;

		auto start = chrono::high_resolution_clock::now();
		while (hp->used > 0)
		{
			he = new HeapEntry();
			hp->remove(he);
			ui son = he->son1;
			ui level = he->level;
			double key = he->key;
			delete he;

			if (Q_map[son].key < key)
			{
				break;
			}

			if (level == 0)
			{
				child = new RTNode(rtree_index, son);
				for (ui i = 0; i < child->num_entries; i++)
				{
					ui vertex_ID = child->entries[i].son;
					for (ui j = 0; j < Q_map[son].query_vertices.size(); j++)
					{
						if (Q_map[son].query_vertices[j].degree <= vertices[vertex_ID].degree && Q_map[son].query_vertices[j].label == vertices[vertex_ID].label)
						{
							ui k = 0;
							k = vde_dim;
							for (; k < vde_dim; k++)
							{
								if (Q_map[son].query_vertices[j].vde[k] > vertices[vertex_ID].vde[k] && abs(Q_map[son].query_vertices[j].vde[k] - vertices[vertex_ID].vde[k]) > epsilon)
								{
									break;
								}
							}
							if (k == vde_dim)
							{
								for (k = 0; k < pde_dim; k++)
								{
									if (vertices[vertex_ID].path_label_group[2 * k + 1] < Q_map[son].query_vertices[j].path_label_group[2 * k] || vertices[vertex_ID].path_label_group[2 * k] > Q_map[son].query_vertices[j].path_label_group[2 * k + 1])
									{
										break;
									}
								}
								if (k == pde_dim)
								{
									for (k = 0; k < pde_dim; k++)
									{
										if (vertices[vertex_ID].path_group[2 * k + 1] < Q_map[son].query_vertices[j].path_group[2 * k])
										{
											break;
										}
									}
									if (k == pde_dim)
									{
										candidates[Q_map[son].query_vertices[j].vid].insert(vertices[vertex_ID].vid);
									}
								}
							}
						}
					}
				}
				delete child;
			}
			else
			{
				child = new RTNode(rtree_index, son);
				for (ui i = 0; i < child->num_entries; i++)
				{
					ui node_ID = child->entries[i].son;
					for (ui j = 0; j < Q_map[son].query_vertices.size(); j++)
					{
						ui k = 0;
						for (; k < pde_dim; k++)
						{
							if (Q_map[son].query_vertices[j].path_label_group[2 * k] > auxiliary_index[node_ID].label_mbr[2 * k + 1] || Q_map[son].query_vertices[j].path_label_group[2 * k + 1] < auxiliary_index[node_ID].label_mbr[2 * k])
							{
								break;
							}
						}
						if (k == pde_dim)
						{
							k = 0;
							for (; k < pde_dim; k++)
							{
								if (Q_map[son].query_vertices[j].path_group[2 * k] > child->entries[i].bounces[2 * k + 1] && abs(Q_map[son].query_vertices[j].path_group[2 * k] - child->entries[i].bounces[2 * k + 1]) > epsilon)
								{
									break;
								}
							}
							if (k == pde_dim)
							{
								Q_map[node_ID].query_vertices.push_back(Q_map[son].query_vertices[j]);
							}
						}
					}
					if (Q_map[node_ID].query_vertices.size() != 0)
					{
						Q_map[node_ID].calculate_key();
						he = new HeapEntry();
						he->son1 = child->entries[i].son;
						he->level = child->level - 1;
						he->key = auxiliary_index[he->son1].key;
						hp->insert(he);
						delete he;
						traversal_index_nodes++;
					}
				}
				delete child;
			}
		}
		delete hp;
		auto end = chrono::high_resolution_clock::now();
		return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
	}
};

vector<double> gen_vde_x(const ui &vertex_label)
{
	ui seed = vertex_label;
	mt19937 psr_gen(seed);
	uniform_real_distribution<double> dis(0.0, 1.0);

	vector<double> embedding_x(vde_dim);
	for (ui i = 0; i < vde_dim; i++)
	{
		embedding_x[i] = dis(psr_gen);
	}

	double sum = std::accumulate(embedding_x.begin(), embedding_x.end(), 0.0);
	for (ui i = 0; i < vde_dim; i++)
	{
		embedding_x[i] = embedding_x[i] / sum;
	}

	return embedding_x;
}

vector<Vertex> gen_vde(const Static_Graph *graph)
{
	vector<Vertex> vertices(graph->getVerticesCount());
	for (ui i = 0; i < graph->getVerticesCount(); i++)
	{
		vertices[i].vid = i;
		vertices[i].label = graph->getVertexLabel(i);
		vertices[i].degree = graph->getVertexDegree(i);
		vertices[i].x = gen_vde_x(vertices[i].label);
	}

	for (ui i = 0; i < graph->getVerticesCount(); i++)
	{
		ui cnt = 0;
		const ui *neighbors = graph->getVertexNeighbors(i, cnt);
		vertices[i].nx = vector<double>(vde_dim, 0);
		for (ui j = 0; j < cnt; j++)
		{
			for (ui k = 0; k < vde_dim; k++)
			{
				vertices[i].nx[k] += vertices[neighbors[j]].x[k];
			}
		}

		vertices[i].vde = vector<double>(vde_dim, 0);
		for (ui j = 0; j < vde_dim; j++)
		{
			vertices[i].vde[j] = vertices[i].x[j] + vertices[i].nx[j];
		}
	}

	return vertices;
}

VertexID selectGQLStartVertex(const Static_Graph *query_graph, ui *candidates_count)
{
	ui start_vertex = 0;

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID cur_vertex = i;

		if (candidates_count[cur_vertex] < candidates_count[start_vertex])
		{
			start_vertex = cur_vertex;
		}
		else if (candidates_count[cur_vertex] == candidates_count[start_vertex] && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex))
		{
			start_vertex = cur_vertex;
		}
	}

	return start_vertex;
}

void updateValidVertices(const Static_Graph *query_graph, VertexID query_vertex, std::vector<bool> &visited,
						 std::vector<bool> &adjacent)
{
	visited[query_vertex] = true;
	ui nbr_cnt;
	const ui *nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i)
	{
		ui nbr = nbrs[i];
		adjacent[nbr] = true;
	}
}

double generateGQLQueryPlan(const Static_Graph *data_graph, const Static_Graph *query_graph, ui *candidates_count, ui *&order, ui *&pivot)
{
	std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
	std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
	order = new ui[query_graph->getVerticesCount()];
	pivot = new ui[query_graph->getVerticesCount()];

	auto start = chrono::high_resolution_clock::now();
	VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
	order[0] = start_vertex;
	updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID next_vertex;
		ui min_value = data_graph->getVerticesCount() + 1;
		for (ui j = 0; j < query_graph->getVerticesCount(); ++j)
		{
			VertexID cur_vertex = j;

			if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex])
			{
				if (candidates_count[cur_vertex] < min_value)
				{
					min_value = candidates_count[cur_vertex];
					next_vertex = cur_vertex;
				}
				else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex))
				{
					next_vertex = cur_vertex;
				}
			}
		}
		updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
		order[i] = next_vertex;
	}

	for (ui i = 1; i < query_graph->getVerticesCount(); ++i)
	{
		VertexID u = order[i];
		for (ui j = 0; j < i; ++j)
		{
			VertexID cur_vertex = order[j];
			if (query_graph->checkEdgeExistence(u, cur_vertex))
			{
				pivot[i] = cur_vertex;
				break;
			}
		}
	}
	auto end = chrono::high_resolution_clock::now();
	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

void generateBN(const Static_Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count)
{
	ui query_vertices_num = query_graph->getVerticesCount();
	bn_count = new ui[query_vertices_num];
	std::fill(bn_count, bn_count + query_vertices_num, 0);
	bn = new ui *[query_vertices_num];
	for (ui i = 0; i < query_vertices_num; ++i)
	{
		bn[i] = new ui[query_vertices_num];
	}

	std::vector<bool> visited_vertices(query_vertices_num, false);
	visited_vertices[order[0]] = true;
	for (ui i = 1; i < query_vertices_num; ++i)
	{
		VertexID vertex = order[i];

		ui nbrs_cnt;
		const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
		for (ui j = 0; j < nbrs_cnt; ++j)
		{
			VertexID nbr = nbrs[j];

			if (visited_vertices[nbr] && nbr != pivot[i])
			{
				bn[i][bn_count[i]++] = nbr;
			}
		}

		visited_vertices[vertex] = true;
	}
}

void generateValidCandidates(const Static_Graph *query_graph, const Static_Graph *data_graph, ui depth, vector<ui> &embedding,
							 ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot)
{
	VertexID u = order[depth];
	LabelID u_label = query_graph->getVertexLabel(u);
	ui u_degree = query_graph->getVertexDegree(u);

	idx_count[depth] = 0;

	VertexID p = embedding[pivot[depth]];
	ui nbr_cnt;
	const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

	for (ui i = 0; i < nbr_cnt; ++i)
	{
		VertexID v = nbrs[i];

		if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
			u_degree <= data_graph->getVertexDegree(v))
		{
			bool valid = true;

			for (ui j = 0; j < bn_cnt[depth]; ++j)
			{
				VertexID u_nbr = bn[depth][j];
				VertexID u_nbr_v = embedding[u_nbr];

				if (!data_graph->checkEdgeExistence(v, u_nbr_v))
				{
					valid = false;
					break;
				}
			}

			if (valid)
			{
				valid_candidate[depth][idx_count[depth]++] = v;
			}
		}
	}
}

double exploreQuickSIStyle(const Static_Graph *data_graph, const Static_Graph *query_graph, ui **candidates,
						   ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &embedding_cnt)
{
	int cur_depth = 0;
	int max_depth = query_graph->getVerticesCount();
	VertexID start_vertex = order[0];

	ui **bn;
	ui *bn_count;

	ui *idx;
	ui *idx_count;
	VertexID **valid_candidate;
	bool *visited_vertices;

	idx = new ui[max_depth];
	idx_count = new ui[max_depth];
	vector<ui> embedding(max_depth);
	visited_vertices = new bool[data_graph->getVerticesCount()];
	std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
	valid_candidate = new ui *[max_depth];

	ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
	for (ui i = 0; i < max_depth; ++i)
	{
		valid_candidate[i] = new VertexID[max_candidate_count];
	}

	idx[cur_depth] = 0;
	idx_count[cur_depth] = candidates_count[start_vertex];
	std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
			  valid_candidate[cur_depth]);

	auto start = chrono::high_resolution_clock::now();

	generateBN(query_graph, order, pivot, bn, bn_count);

	while (true)
	{
		while (idx[cur_depth] < idx_count[cur_depth])
		{
			VertexID u = order[cur_depth];
			VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
			embedding[u] = v;
			visited_vertices[v] = true;
			idx[cur_depth] += 1;

			if (cur_depth == max_depth - 1)
			{
				embedding_cnt += 1;

				visited_vertices[v] = false;
				if (embedding_cnt >= output_limit_num)
				{
					goto EXIT;
				}
			}
			else
			{
				cur_depth += 1;
				idx[cur_depth] = 0;
				generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
										visited_vertices, bn, bn_count, order, pivot);
			}
		}

		cur_depth -= 1;
		if (cur_depth < 0)
			break;
		else
			visited_vertices[embedding[order[cur_depth]]] = false;
	}

EXIT:
	auto end = chrono::high_resolution_clock::now();
	delete[] bn_count;
	delete[] idx;
	delete[] idx_count;
	delete[] visited_vertices;
	for (ui i = 0; i < max_depth; ++i)
	{
		delete[] bn[i];
		delete[] valid_candidate[i];
	}

	delete[] bn;
	delete[] valid_candidate;

	return double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
}

double refinement(const Static_Graph *static_data_graph, const Static_Graph *static_query_graph, const vector<set<ui>> &candidate_set, ui &answer_num)
{
	size_t call_count = 0;
	size_t output_limit = MAX_LIMIT;

	ui **candidates = new ui *[static_query_graph->getVerticesCount()];
	ui *candidates_count = new ui[static_query_graph->getVerticesCount()];
	ui *matching_order = NULL;
	ui *pivots = NULL;

	for (ui i = 0; i < static_query_graph->getVerticesCount(); i++)
	{
		candidates[i] = new ui[candidate_set[i].size()];
		candidates_count[i] = candidate_set[i].size();
		ui count = 0;
		for (set<ui>::iterator iter = candidate_set[i].begin(); iter != candidate_set[i].end(); iter++)
		{
			candidates[i][count] = *iter;
			count++;
		}
	}

	vector<vector<ui>> answers;

	double query_plan_time = generateGQLQueryPlan(static_data_graph, static_query_graph, candidates_count, matching_order, pivots);

	size_t embedding_cnt = 0;

	double refinement_time = exploreQuickSIStyle(static_data_graph, static_query_graph, candidates, candidates_count, matching_order, pivots, output_limit, embedding_cnt);

	answer_num = embedding_cnt;
	
	for (ui i = 0; i < static_query_graph->getVerticesCount(); i++)
	{
		delete[] candidates[i];
	}
	delete[] candidates;
	delete[] candidates_count;
	delete[] matching_order;
	delete[] pivots;

	return query_plan_time + refinement_time;
}