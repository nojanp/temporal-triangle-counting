# temporal-triangle-counting

Instructions:

1) Build the code by running "make" in the main directory

2) Run the code (ettc: exact temporal triangle counting) using the following command:
```
./ettc  <path to input> <delta_{1,3}> <\delta_{1,2}> <\delta_{2,3}> <path to output>
```
The first line in the input file should contain the number of vertices and edges. Starting from the second line, each line is a directed temporal edge consisting of three numbers, starting point, end point, and its timestamp.

Example:
```
./ettc ./example_graph.txt 100 80 80 out.txt
```
