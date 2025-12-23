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
	string dataset_path = "../Test/";
	string data_graph_name = "../Test/data_graph.graph";
	string query_graph_name = "../Test/query_graph.graph";
	string mode = "offline";
	string answer_limit = "MAX";

	CLI::App app{"App description"};
	app.add_option("-f,--file", dataset_path, "dataset path")->required(false);
	app.add_option("-d,--data", data_graph_name, "data graph path")->required(false);
	app.add_option("-q,--query", query_graph_name, "query graph path")->required(false);
	app.add_option("-m,--mode", mode, "offline or online mode")->required(false);
	app.add_option("-p,--partition", partition_num, "partition number")->required(false);
	app.add_option("-l,--length", path_length, "path length")->required(false);
	app.add_option("-e,--embedding", vde_dim, "embedding dimension")->required(false);
	app.add_option("-n,--answers", answer_limit, "max answer number")->required(false);

	CLI11_PARSE(app, argc, argv);

	path_length += 1;
	pde_dim = vde_dim * (path_length);

	string partitions_path = dataset_path + "gnn-pe/partitions/";
	if (answer_limit == "MAX")
	{
		MAX_LIMIT = UINT_MAX;
	}
	else
	{
		MAX_LIMIT = stoi(answer_limit);
	}

	Static_Graph *static_data_graph = new Static_Graph(true);
	static_data_graph->loadGraphFromFile(data_graph_name);
	static_data_graph->printGraphMetaData();

	if (mode == "offline")
	{
		vector<ui> membership(static_data_graph->getVerticesCount());
		vector<ui> sorted_nodes(static_data_graph->getVerticesCount());

		ifstream fin(dataset_path + "gnn-pe/membership.txt");
		for (ui i = 0; i < static_data_graph->getVerticesCount(); i++)
		{
			fin >> sorted_nodes[i] >> membership[sorted_nodes[i]];
		}
		fin.close();

		vector<vector<ui>> partitions_paths(partition_num);

		vector<vector<ui>> all_paths;
		unordered_set<vector<ui>, VectorHash> all_paths_set;

		for (ui node : sorted_nodes)
		{
			vector<ui> path = {node};
			dfs(node, path_length - 2, path, static_data_graph, all_paths, all_paths_set, partitions_paths[membership[node]]);
		}

		for (ui i = 0; i < partition_num; i++)
		{
			string partition_path = partitions_path + "partition-" + to_string(i) + "/";
			ofstream fout(partition_path + "partition_paths.txt");
			fout << partitions_paths[i].size() << endl;
			for (ui j = 0; j < partitions_paths[i].size(); j++)
			{
				fout << partitions_paths[i][j] << endl;
			}
			fout.close();
		}

		ofstream fout(dataset_path + "/gnn-pe/all_paths.txt");
		fout << all_paths.size() << endl;
		for (ui i = 0; i < all_paths.size(); i++)
		{
			for (ui j = 0; j < path_length; j++)
			{
				fout << all_paths[i][j] << " ";
			}
			fout << endl;
		}
	}

	if (mode == "online")
	{
		vector<Vertex> data_vertices = gen_vde(static_data_graph);

		vector<Path> data_paths = gen_pde(data_vertices, dataset_path + "/gnn-pe/all_paths.txt");

		vector<Partition> partitions;
		for (ui i = 0; i < partition_num; i++)
		{
			string partition_path = partitions_path + "partition-" + to_string(i) + "/";
			Partition partition(data_paths, partition_path);
			partitions.push_back(partition);
		}

		Static_Graph *static_query_graph = new Static_Graph(true);
		static_query_graph->loadGraphFromFile(query_graph_name);

		vector<vector<ui>> all_paths;
		unordered_set<vector<ui>, VectorHash> all_paths_set;

		for (ui node = 0; node < static_query_graph->getVerticesCount(); node++)
		{
			vector<ui> path = {node};
			dfs_query(node, path_length - 2, path, static_query_graph, all_paths, all_paths_set);
		}

		auto start = chrono::high_resolution_clock::now();
		vector<Vertex> query_vertices = gen_vde(static_query_graph);
		vector<Query_Path> query_paths = gen_query_pde(query_vertices, all_paths);
		auto end = chrono::high_resolution_clock::now();
		double query_plan_time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
		Query_Plan Q(query_paths);

		vector<vector<set<ui>>> candidate_sets(partition_num, vector<set<ui>>(static_query_graph->getVerticesCount()));

		vector<set<ui>> candidate_set(static_query_graph->getVerticesCount());

		vector<double> search_time(partition_num, 0);
#pragma omp parallel for
		for (ui pid = 0; pid < partition_num; pid++)
		{
			search_time[pid] = partitions[pid].query(static_query_graph->getVerticesCount(), candidate_sets[pid], Q);
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

		cout << "Answer Number: " << answer_num << " Query Time (ms): " << NANOSECTOMSEC(query_plan_time + max_searh_time + refine_time) << endl;

		delete static_query_graph;
	}

	return 0;
}