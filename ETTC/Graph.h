#ifndef ETTC_GRAPH_H_
#define ETTC_GRAPH_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

//type aliases
typedef int64_t VertexEdgeId;
typedef int64_t TemporalTime;


class Edge
{
  public:
    VertexEdgeId src_;
    VertexEdgeId dst_;

    bool operator<(const Edge& edge) const
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

    bool operator<(const TemporalEdge& edge) const
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
    VertexEdgeId *offsets_;         
    VertexEdgeId *nbrs_;     
    VertexEdgeId *degen_order_; 
    VertexEdgeId *sort_by_degen_;

    VertexEdgeId *temporal_start_pos_; 
    
    // constructor
    CSRGraph(){};
    CSRGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges);


    // copy constructor

    // destructor
    // ~CSRGraph();

    VertexEdgeId maxDegree();
    void sortNbrs();
    void findDegenOrdering();
    void relabelByDegenOrder();
    VertexEdgeId getEdgeIndx(VertexEdgeId u, VertexEdgeId v);

    
    void printGraph();
    
    void deleteGraph();
};

class CSRTemporalNbr
{
  public:
    VertexEdgeId dst_;
    TemporalTime time_;

    // // constructors
    // CSRTemporalNbr(){};
    // CSRTemporalNbr(VertexEdgeId dst, TemporalTime time);

    bool operator<(const CSRTemporalNbr& edge) const
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
    VertexEdgeId *offsets_;        
    CSRTemporalNbr *temporal_nbrs_;  
    TemporalTime *times_; 
    
    // constructor
    CSRTemporalGraph();
    CSRTemporalGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges);


    // copy constructor

    // destructor
    // ~CSRTemporalGraph();

    VertexEdgeId numberOfStaticDirectedEdges();
    void sortNbrs();
    void relabelByDegenOrder(VertexEdgeId *degen_order, VertexEdgeId *sort_by_degen); 
    CSRGraph extractMultGraph();

    void reloadTimes();

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
    CSRGraph out_edge_dag_;
    CSRGraph in_edge_dag_;

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

    void sortEdges();
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
    
    void sortEdges();
    VertexEdgeId numberOfStaticDirectedEdges();
    Graph ExtractStaticGraph();
    CSRTemporalGraph convertToCSR();
    void printGraph();

    void deleteGraph();
};


TemporalGraph loadTemporalGraph(const char *path);


#endif


