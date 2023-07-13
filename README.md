# GNN-PE

In this project, we open-source the source code of our GNN-PE approach, including both offline and online processes (i.e., offline.tar.gz and online.tar.gz).

The real-world and default synthetic datasets used in our paper are stored in the Datasets directory. As some synthetic datasets are large, we do not upload them. You can easily generate them by following the instruction in our paper.

Due to space limitations, we will introduce how to reproduce the results of our experiments over the Yeast dataset.

## Offline Process
1. Download the offline.tar.gz, execute the following commands to decompress the source code and go to the root directory of the file.

```
tar -xzvf offline.tar.gz
cd offline
```

2. Under the root directory of the offline directory, execute the following conda commands to configure the python environment.

```

```




## Online Process
1. Download the online.tar.gz, execute the following commands to decompress the source code and go to the root directory of the file.

```
tar -xzvf online.tar.gz
cd online
```

3. Under the root directory of the online directory, execute the following commands to compile the source code.

```
mkdir build
cd build
cmake ..
make
```

3. Execute the following command to run the experiment over the Yeast dataset (the necessary pre-computed data produced by the offline process is already stored in the Datasets directory).

```
./build/src/main
```




