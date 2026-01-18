# GNN-based Path Embeddings for Efficient and Exact Subgraph Matching

This repository contains the implementation for the papers "Efficient Exact Subgraph Matching via GNN-based Path Dominance Embedding" and "GNN-based Path Embeddings for Efficient and Exact Subgraph Matching". These papers present two GNN-based method for solving the exact subgraph matching query problem: given a data graph $G$ and a query graph $q$, the subgraph matching query aims to enumerate all subgraphs $g\in G$ that are isomorphic to $q$.

For details of our methods, please see our papers at [VLDB 2024](https://vldb.org/pvldb/volumes/17/paper/Efficient%20Exact%20Subgraph%20Matching%20via%20GNN-based%20Path%20Dominance%20Embedding) and [arxiv 2025](https://arxiv.org/abs/2309.15641).

Welcome to cite our paper!

```
@inproceedings{ye2024efficient,
  title={Efficient Exact Subgraph Matching via GNN-based Path Dominance Embedding},
  author={Ye, Yutong and Lian, Xiang and Chen, Mingsong},
  booktitle={Proceedings of the International Conference on Very Large Data Bases (PVLDB)},
  pages={1628--1641},
  year={2024}
}
```

```
@misc{ye2025efficientexactsubgraphmatching,
      title={Efficient Exact Subgraph Matching via GNN-based Path Dominance Embedding (Technical Report)}, 
      author={Yutong Ye and Xiang Lian and Mingsong Chen},
      year={2025},
      eprint={2309.15641},
      archivePrefix={arXiv},
      primaryClass={cs.DB},
      url={https://arxiv.org/abs/2309.15641}, 
}
```

## Workflow

To run a subgraph matching query, follow these three steps:

1. Data set Preparation: Prepare data set by Python
2. Offline Pre-computation: Generate vertex/path embeddings and construct index
3. Online Query Process: Traverse index and refine the candidates

## Getting Started

### Dependencies
1. The codes require the following dependences:

* A modern C++ compiler compliant with the C++17 standard (gcc/g++ >= 12.2)
* CMake (>= 3.28)

2. Under the root directory, execute the following conda commands to configure the Python environment.

```
conda create --name <new_environment_name> --file requirements.txt
conda activate <new_environment_name>
```

### GNN-PE

1. Turn into the directory of GNN-PE, execute the following command to prepare the data set.

```
python gnnpe.py --f ../Test/ --d ../Test/data_graph.gpickle.gz
```

|Parameter|Value|Description|
|:----:|:----:|:----:|
|-f|../Test/|Dataset Path|
|-d|../Test/data_graph.gpickle.gz|Data Graph Path|
|-p|5|Partition Number|
|-l|2|Path Length|


2. Build the project.

```
mkdir build
cd build
cmake ..
make
```

3. Return to the root directory of GNN-PE, and execute the following command to run a quick start example.

```
./build/src/main -f ../Test/ -d ../Test/data_graph.graph -m offline

./build/src/main -f ../Test/ -d ../Test/data_graph.graph -q ../Test/query_graph.graph -m online
```

|Parameter|Value|Description|
|:----:|:----:|:----:|
|-f|../Test/|Dataset Path|
|-d|../Test/data_graph.gpickle.gz|Data Graph Path|
|-q|../Test/query_graph.graph|Query Graph Path|
|-m|offline/online|Execution Mode|
|-p|5|Partition Number|
|-l|2|Path Length|
|-e|2|Vertex Embedding Dimension|
|-n|MAX|Max Answer Limit|

### GNN-PGE

1. Trun into the directory of GNN-PE, execute the following command to prepare the data set.

```
python gnnpge.py --f ../Test/ --d ../Test/data_graph.gpickle.gz
```

|Parameter|Value|Description|
|:----:|:----:|:----:|
|-f|../Test/|Dataset Path|
|-d|../Test/data_graph.gpickle.gz|Data Graph Path|
|-p|5|Partition Number|


2. Build the project.

```
mkdir build
cd build
cmake ..
make
```

3. Return to the root directory of GNN-PE execute the following command to run a quick start example.

```
./build/src/main -f ../Test/ -d ../Test/data_graph.graph -m offline

./build/src/main -f ../Test/ -d ../Test/data_graph.graph -q ../Test/query_graph.graph -m online
```

|Parameter|Value|Description|
|:----:|:----:|:----:|
|-f|../Test/|Dataset Path|
|-d|../Test/data_graph.gpickle.gz|Data Graph Path|
|-q|../Test/query_graph.graph|Query Graph Path|
|-m|offline/online|Execution Mode|
|-p|5|Partition Number|
|-l|2|Path Length|
|-e|2|Vertex Embedding Dimension|
|-n|MAX|Max Answer Limit|

