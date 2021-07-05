#ifndef ETTC_GRAPH_H_
#define ETTC_GRAPH_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

//type aliases
typedef int64_t VertexEdgeId;  // type for ID of vertices or edges
typedef int64_t TemporalTime;  // type for temporal time


class Edge
{
  public:
    VertexEdgeId src_;
    VertexEdgeId dst_;

    bool operator<(const Edge& edge) const  // comparison: first sources and then destinations
    {
        if (src_ == edge.src_)
            return dst_ < edge.dst_;
        return src_ < edge.src_;
    }
};

class TemporalEdge: public Edge 
{
  public:
    TemporalTime time_;

    bool operator<(const TemporalEdge& edge) const // comparison: first sources then destinations and then timestamps
    {
        if (src_ == edge.src_)
        {
            if(dst_ == edge.dst_)
                return time_ < edge.time_;
            return dst_ < edge.dst_;
        }
        return src_ < edge.src_;
    }
};


class CSRGraph
{
  public:
    VertexEdgeId num_vertices_; 
    VertexEdgeId num_edges_;      
    VertexEdgeId *offsets_;  // for each vertex, determines the starting point of the sequence of neighbors in the nbrs_ array
    VertexEdgeId *nbrs_;     // list of neighbors
    VertexEdgeId *degen_order_; // place of each vertex in degeneracy ordering
    VertexEdgeId *sort_by_degen_; // vertices sorted by degeneracy ordering

    VertexEdgeId *temporal_start_pos_; 
    
    // constructor
    CSRGraph(){};
    CSRGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges);


    // copy constructor

    // destructor
    // ~CSRGraph();

    VertexEdgeId maxDegree();    // finds vertex with max degree (static)
    void sortNbrs();             // sort the neighbor sequnce of each vertex in nbrs_ arrray by destination ID
    void findDegenOrdering();    // find degeneracy ordering and stores it in degen_order_
    void relabelByDegenOrder();  
    VertexEdgeId getEdgeIndx(VertexEdgeId u, VertexEdgeId v);  // find the index of the edge in the nbrs_ array using a binary search on the sequence of neighbors of u

    
    void printGraph();
    
    void deleteGraph();
};


class CSRTemporalNbr  // a neighbor in a temporal graph: destination and timestamp of the edge
{
  public:
    VertexEdgeId dst_;
    TemporalTime time_;

    // // constructors
    // CSRTemporalNbr(){};
    // CSRTemporalNbr(VertexEdgeId dst, TemporalTime time);

    bool operator<(const CSRTemporalNbr& edge) const // comparison: first destinations and then timestamps
    {
        if (dst_ == edge.dst_)
            return time_ < edge.time_;
        return dst_ < edge.dst_;
    }
};

class CSRTemporalGraph
{
  public:
    VertexEdgeId num_vertices_;
    VertexEdgeId num_edges_; 
    VertexEdgeId *offsets_;    // for each vertex, determines the starting point of the sequence of neighbors in the temporal_nbrs_ array    
    CSRTemporalNbr *temporal_nbrs_;  // list of temporal neighbors
    TemporalTime *times_;  // list of timestamps (in same order with the temporal_nbrs_)
    
    // constructor
    CSRTemporalGraph();
    CSRTemporalGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges);


    // copy constructor

    // destructor
    // ~CSRTemporalGraph();

    VertexEdgeId numberOfStaticDirectedEdges();  // number of directed static edges ignoring the timestamps
    void sortNbrs();  // sort the neighbor sequnce of each vertex in temporal_nbrs_ arrray by destination ID and then by timestamp
    void relabelByDegenOrder(VertexEdgeId *degen_order, VertexEdgeId *sort_by_degen); 
    CSRGraph extractMultGraph(); // multGraph is a temporary static grraph (DAG). For each edge in multGraph, temporal_start_pos_ indicates
                                 // the starting point of sequence of temporal edges corresponding to this static edge

    void reloadTimes(); // writes the timestamp of all temporal edges into the times_ array again (used after sorting the vertices)

    // lower bounds and upper bounds on the times_ array
    VertexEdgeId edgeTimeMinLimitSearch(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value);
    VertexEdgeId edgeTimeMaxLimitSearch(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime max_value);

    VertexEdgeId edgeTimeMinLimitCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value);
    VertexEdgeId edgeTimeMaxLimitCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime max_value);
    VertexEdgeId edgeTimeIntervalCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value, TemporalTime max_value);

    void printTimeSpan();
    void printGraph();
    
    void deleteGraph();

};

class CSRDAG
{
  public:
    CSRGraph out_edge_dag_;  // DAG stored as a CSR graph. For each vertex we store out-neighbors as neighbors
    CSRGraph in_edge_dag_;   // DAG stored as a CSR graph. For each vertex we store in-neighbors as neighbors

    // constructor
    CSRDAG(){};
    CSRDAG(CSRGraph &csr_graph);

};

class Graph
{
  public:
    VertexEdgeId num_vertices_; 
    VertexEdgeId num_edges_; 
    Edge *edges_; 

    // Constructors
    Graph(){}
    Graph(VertexEdgeId num_vertices, VertexEdgeId num_edges);

    //Copy constructor

    //Destructor
    // ~Graph(); // change to delete graph : when we return a graph then the destructor deletes the array and then returns a copy of the pointer to the array

    void sortEdges(); // sort the list of edges by source and then destination
    CSRGraph convertToCSR();

    void printGraph();

    void deleteGraph();
};


class TemporalGraph
{
  public:
    VertexEdgeId num_vertices_;
    VertexEdgeId num_edges_; 
    TemporalEdge *temporal_edges_; 

    // Constructor
    TemporalGraph(){}
    TemporalGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges); 

    //Copy constructor

    //Destructor
    // ~TemporalGraph(); // change to delete graph : when we return a graph then the destructor deletes the array and then returns a copy of the pointer to the array
    
    void sortEdges(); // sort the list of temporal edges first by sources then destinations and then timestamps
    VertexEdgeId numberOfStaticDirectedEdges();
    Graph ExtractStaticGraph();
    CSRTemporalGraph convertToCSR();
    void printGraph();

    void deleteGraph();
};


TemporalGraph loadTemporalGraph(const char *path);


#endif


