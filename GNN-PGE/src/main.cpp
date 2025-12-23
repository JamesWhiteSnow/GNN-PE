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
#include "custom.h"
#include "CLI11.hpp"

#define NOMINMAX
#undef min
#undef max

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
#include <omp.h>
#include <chrono>
#include <cstdio>
#include <numeric>

using namespace std;
using namespace chrono;

int main(int argc, char *argv[])
{
	string data_graph_path = "../Test/data_graph.graph";
	string mode = "offline";
	string query_graph_path = "../Test/query_graph.graph";
	string answer_limit = "MAX";
	string dataset_path = "../Test/";

	CLI::App app{"App description"};
	app.add_option("-f,--file", dataset_path, "dataset path");
	app.add_option("-d,--data", data_graph_path, "data graph path");
	app.add_option("-q,--query", query_graph_path, "query graph path");
	app.add_option("-m,--mode", mode, "offline or online mode");
	app.add_option("-p,--partition", partition_num, "partition number");
	app.add_option("-l,--length", path_length, "path length");
	app.add_option("-e,--embedding", vde_dim, "embedding dimension");
	app.add_option("-n,--answers", answer_limit, "max answer number");

	CLI11_PARSE(app, argc, argv);

	pde_dim = vde_dim * (path_length);

	string partitions_path = dataset_path + "gnn-pge/partitions/";
	if (answer_limit == "MAX")
	{
		MAX_LIMIT = UINT_MAX;
	}
	else
	{
		MAX_LIMIT = stoi(answer_limit);
	}

	Static_Graph *static_data_graph = new Static_Graph(true);
	static_data_graph->loadGraphFromFile(data_graph_path);

	vector<Vertex> data_vertices;
	vector<vector<ui>> partition_vertices(partition_num);

	vector<ui> membership(static_data_graph->getVerticesCount());
	vector<ui> sorted_nodes(static_data_graph->getVerticesCount());

	ifstream fin(dataset_path + "gnn-pge/membership.txt");
	for (ui i = 0; i < static_data_graph->getVerticesCount(); i++)
	{
		fin >> sorted_nodes[i] >> membership[sorted_nodes[i]];
	}
	fin.close();

	for (ui node : sorted_nodes)
	{
		partition_vertices[membership[node]].push_back(node);
	}

	if (mode == "offline")
	{
		data_vertices = gen_vde(static_data_graph);

		ui total_paths_num = 0;

		for (ui node : sorted_nodes)
		{
			vector<ui> path = {node};
			vector<vector<ui>> all_paths;
			dfs(node, 1, path, static_data_graph, all_paths);

			total_paths_num += all_paths.size();

			if (all_paths.size() == 0)
			{
				for (ui i = 0; i < vde_dim; i++)
				{
					data_vertices[node].path_group.push_back(data_vertices[node].vde[i]);
					data_vertices[node].path_group.push_back(data_vertices[node].vde[i]);
					data_vertices[node].path_label_group.push_back(data_vertices[node].x[i]);
					data_vertices[node].path_label_group.push_back(data_vertices[node].x[i]);
				}
				for (ui i = vde_dim; i < pde_dim; i++)
				{
					data_vertices[node].path_group.push_back(0);
					data_vertices[node].path_group.push_back(0);
					data_vertices[node].path_label_group.push_back(0);
					data_vertices[node].path_label_group.push_back(0);
				}
				continue;
			}

			vector<vector<double>> path_embeddings;
			vector<vector<double>> path_label_embeddings;
			for (ui i = 0; i < all_paths.size(); i++)
			{
				vector<double> path_embedding;
				vector<double> path_label_embedding;
				for (ui j = 0; j < path_length; j++)
				{
					for (ui k = 0; k < vde_dim; k++)
					{
						path_embedding.push_back(data_vertices[all_paths[i][j]].vde[k]);
						path_label_embedding.push_back(data_vertices[all_paths[i][j]].x[k]);
					}
				}
				path_embeddings.push_back(path_embedding);
				path_label_embeddings.push_back(path_label_embedding);
			}

			for (ui i = 0; i < pde_dim; i++)
			{
				data_vertices[node].path_group.push_back(path_embeddings[0][i]);
				data_vertices[node].path_group.push_back(path_embeddings[0][i]);
			}
			for (ui i = 0; i < pde_dim; i++)
			{
				data_vertices[node].path_label_group.push_back(path_label_embeddings[0][i]);
				data_vertices[node].path_label_group.push_back(path_label_embeddings[0][i]);
			}
			for (ui i = 1; i < path_embeddings.size(); i++)
			{
				for (ui j = 0; j < pde_dim; j++)
				{
					if (data_vertices[node].path_group[2 * j] > path_embeddings[i][j])
					{
						data_vertices[node].path_group[2 * j] = path_embeddings[i][j];
					}
					if (data_vertices[node].path_group[2 * j + 1] < path_embeddings[i][j])
					{
						data_vertices[node].path_group[2 * j + 1] = path_embeddings[i][j];
					}
				}
				for (ui j = 0; j < pde_dim; j++)
				{
					if (data_vertices[node].path_label_group[2 * j] > path_label_embeddings[i][j])
					{
						data_vertices[node].path_label_group[2 * j] = path_label_embeddings[i][j];
					}
					if (data_vertices[node].path_label_group[2 * j + 1] < path_label_embeddings[i][j])
					{
						data_vertices[node].path_label_group[2 * j + 1] = path_label_embeddings[i][j];
					}
				}
			}
		}

		ofstream fout(dataset_path + "gnn-pge/data_vertices.bin", std::ios::binary);
		ui count = data_vertices.size();
		fout.write(reinterpret_cast<const char *>(&count), sizeof(ui));
		for (ui i = 0; i < data_vertices.size(); i++)
		{
			fout.write(reinterpret_cast<const char *>(&data_vertices[i].vid), sizeof(ui));
			fout.write(reinterpret_cast<const char *>(&data_vertices[i].label), sizeof(ui));
			fout.write(reinterpret_cast<const char *>(&data_vertices[i].degree), sizeof(ui));
			fout.write(reinterpret_cast<const char *>(&data_vertices[i].key), sizeof(double));
			fout.write(reinterpret_cast<const char *>(data_vertices[i].x.data()), vde_dim * sizeof(double));
			fout.write(reinterpret_cast<const char *>(data_vertices[i].nx.data()), vde_dim * sizeof(double));
			fout.write(reinterpret_cast<const char *>(data_vertices[i].vde.data()), vde_dim * sizeof(double));
			fout.write(reinterpret_cast<const char *>(data_vertices[i].path_group.data()), pde_dim * 2 * sizeof(double));
			fout.write(reinterpret_cast<const char *>(data_vertices[i].path_label_group.data()), pde_dim * 2 * sizeof(double));
		}
		fout.close();
	}

	if (mode == "online")
	{
		fin.open(dataset_path + "gnn-pge/data_vertices.bin", std::ios::binary);

		ui count;
		fin.read(reinterpret_cast<char *>(&count), sizeof(ui));

		for (ui i = 0; i < count; i++)
		{
			Vertex data_vertex;

			fin.read(reinterpret_cast<char *>(&data_vertex.vid), sizeof(ui));
			fin.read(reinterpret_cast<char *>(&data_vertex.label), sizeof(ui));
			fin.read(reinterpret_cast<char *>(&data_vertex.degree), sizeof(ui));
			fin.read(reinterpret_cast<char *>(&data_vertex.key), sizeof(double)); 

			std::vector<double> x(vde_dim);
			fin.read(reinterpret_cast<char *>(x.data()), vde_dim * sizeof(double));
			data_vertex.x = x;

			std::vector<double> nx(vde_dim);
			fin.read(reinterpret_cast<char *>(nx.data()), vde_dim * sizeof(double));
			data_vertex.nx = nx;

			std::vector<double> vde(vde_dim);
			fin.read(reinterpret_cast<char *>(vde.data()), vde_dim * sizeof(double));
			data_vertex.vde = vde;

			std::vector<double> path_group(pde_dim * 2);
			fin.read(reinterpret_cast<char *>(path_group.data()), pde_dim * 2 * sizeof(double));
			data_vertex.path_group = path_group;

			std::vector<double> path_label_group(pde_dim * 2);
			fin.read(reinterpret_cast<char *>(path_label_group.data()), pde_dim * 2 * sizeof(double));
			data_vertex.path_label_group = path_label_group;

			data_vertices.push_back(data_vertex);
		}
		fin.close();
	}

	double index_build_time = 0;
	vector<Partition> partitions;
	for (ui i = 0; i < partition_num; i++)
	{
		string partition_path = partitions_path + "partition-" + to_string(i) + "/";
		Partition partition(data_vertices, partition_path, partition_vertices[i], index_build_time);
		partitions.push_back(partition);
	}

	if (mode == "online")
	{
		string query_path = query_graph_path;
		Static_Graph *static_query_graph = new Static_Graph(true);
		static_query_graph->loadGraphFromFile(query_path);

		vector<Vertex> query_vertices = gen_vde(static_query_graph);

		ui total_paths_num = 0;

		auto start = chrono::high_resolution_clock::now();
		for (ui node = 0; node < query_vertices.size(); node++)
		{
			vector<ui> path = {node};
			vector<vector<ui>> all_paths;
			dfs(node, 1, path, static_query_graph, all_paths);

			total_paths_num += all_paths.size();

			vector<vector<double>> path_embeddings;
			vector<vector<double>> path_label_embeddings;
			for (ui i = 0; i < all_paths.size(); i++)
			{
				vector<double> path_embedding;
				vector<double> path_label_embedding;
				for (ui j = 0; j < path_length; j++)
				{
					for (ui k = 0; k < vde_dim; k++)
					{
						path_embedding.push_back(query_vertices[all_paths[i][j]].vde[k]);
						path_label_embedding.push_back(query_vertices[all_paths[i][j]].x[k]);
					}
				}
				path_embeddings.push_back(path_embedding);
				path_label_embeddings.push_back(path_label_embedding);
			}

			if (path_embeddings.size() == 0)
			{
				continue;
			}

			for (ui i = 0; i < pde_dim; i++)
			{
				query_vertices[node].path_group.push_back(path_embeddings[0][i]);
				query_vertices[node].path_group.push_back(path_embeddings[0][i]);
			}
			for (ui i = 0; i < pde_dim; i++)
			{
				query_vertices[node].path_label_group.push_back(path_label_embeddings[0][i]);
				query_vertices[node].path_label_group.push_back(path_label_embeddings[0][i]);
			}
			for (ui i = 1; i < path_embeddings.size(); i++)
			{
				for (ui j = 0; j < pde_dim; j++)
				{
					if (query_vertices[node].path_group[2 * j] > path_embeddings[i][j])
					{
						query_vertices[node].path_group[2 * j] = path_embeddings[i][j];
					}
					if (query_vertices[node].path_group[2 * j + 1] < path_embeddings[i][j])
					{
						query_vertices[node].path_group[2 * j + 1] = path_embeddings[i][j];
					}
				}
				for (ui j = 0; j < pde_dim; j++)
				{
					if (query_vertices[node].path_label_group[2 * j] > path_label_embeddings[i][j])
					{
						query_vertices[node].path_label_group[2 * j] = path_label_embeddings[i][j];
					}
					if (query_vertices[node].path_label_group[2 * j + 1] < path_label_embeddings[i][j])
					{
						query_vertices[node].path_label_group[2 * j + 1] = path_label_embeddings[i][j];
					}
				}
			}

			query_vertices[node].key = 0;
			for (ui i = 0; i < pde_dim; i++)
			{
				query_vertices[node].key -= query_vertices[node].path_group[2 * i];
			}
		}
		auto end = chrono::high_resolution_clock::now();

		double query_plan_time = double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());

		Query_Plan Q(query_vertices);

		vector<vector<set<ui>>> candidate_sets(partition_num, vector<set<ui>>(static_query_graph->getVerticesCount()));

		vector<set<ui>> candidate_set(static_query_graph->getVerticesCount());

		vector<double> search_time(partition_num, 0);
#pragma omp parallel for
		for (ui pid = 0; pid < partition_num; pid++)
		{
			search_time[pid] = partitions[pid].query(candidate_sets[pid], Q);
		}

		for (ui pid = 0; pid < partition_num; pid++)
		{
			for (ui i = 0; i < static_query_graph->getVerticesCount(); i++)
			{
				candidate_set[i].insert(candidate_sets[pid][i].begin(), candidate_sets[pid][i].end());
			}
		}

		double max_searh_time = *max_element(search_time.begin(), search_time.end());

		ui answer_num = 0;
		double refine_time = refinement(static_data_graph, static_query_graph, candidate_set, answer_num);

		cout << "Answer Num: " << answer_num << " Query Time (ms): " << NANOSECTOMSEC(max_searh_time + refine_time + query_plan_time) << endl;
	}

	return 0;
}