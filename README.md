# GNN-PE

In this project, we open-source the source code of our GNN-PE approach, including both offline and online processes (i.e., offline.tar.gz and online.tar.gz).

Next, we will introduce how to reproduce the results of our experiments over the Yeast dataset.

## Offline Process


## Online Process
1. Download the online.tar.gz, execute the following commands to decompress the source code and go to the root directory of the file

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




