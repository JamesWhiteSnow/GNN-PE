import numpy as np
import networkx as nx
import os
import time
import random
import pymetis
import pickle
import argparse
import shutil

os.environ["CUDA_VISIBLE_DEVICES"] = "0"
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

seed=2022
random.seed(seed)
np.random.seed(seed)


def delete_folder(folder_path):
    try:
        shutil.rmtree(folder_path)
        print('Successfully Delete!')
    except FileNotFoundError:
        print('File not Found!')
    except PermissionError:
        print('Permission Error!')
    except Exception as e:
        print('Errors during Deletion!')

def create_folder(filename):
    filename = filename.strip()
    filename = filename.rstrip("\\")
    isExists = os.path.exists(filename)

    if not isExists:
        os.makedirs(filename)
        print(filename+" Create Successful")
        return  True
    else:
        print(filename+" already exists")
        return False


parser=argparse.ArgumentParser()
parser.add_argument('--f',type=str,default='../Test/',help='file path')
parser.add_argument('--d',type=str,default='../Test/data_graph.gpickle.gz',help='data graph')
parser.add_argument('--p',type=int,default=5,help='partition number')
parser.add_argument('--l',type=int,default=2,help='path length')
args=parser.parse_args()

dataset_path=str(args.f)
data_graph_path=str(args.d)
partition_num=args.p
path_length=args.l

with open(data_graph_path,'rb') as f:
    data_graph=pickle.load(f)

#创建分区文件路径
delete_folder(dataset_path+'gnn-pe')
create_folder(dataset_path+'gnn-pe/')
create_folder(dataset_path+'gnn-pe/partitions')
for i in range(partition_num):
    create_folder(dataset_path+'gnn-pe/partitions/'+'partition-{0}'.format(i))

adjacency_list = []
for node in data_graph.nodes():
    adjacency_list.append(np.array(list(data_graph.neighbors(node))))
n_cuts, membership = pymetis.part_graph(partition_num, adjacency=adjacency_list,recursive=True)

nodes_degree=dict(data_graph.degree())
nodes_degree_sorted=sorted(nodes_degree.items(),key=lambda x:x[1])

with open(dataset_path+'gnn-pe/membership.txt','w') as f:
    for node,degree in nodes_degree_sorted:
        f.write(str(node)+' '+str(membership[node])+'\n')