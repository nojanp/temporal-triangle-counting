#ifndef T_TRIANGLE_H_
#define T_TRIANGLE_H_


#include "Graph.h"
#include <map>

//type aliases
typedef int64_t Count;

class MotifCounter
{
  public:

    VertexEdgeId *edge_count; 
    VertexEdgeId *edge_count_cum;
    TemporalTime delta_, delta1_, delta2_;
    VertexEdgeId highest_mult_;

    
    Count static_triangles_count_ = 0;

    Count motif_type_counts_[8] = {0};
    Count motif_counts_[6][8] = {0};
    Count triangle_motif_counts_[6][8] = {0};
    
    // void populateEdgeCountIndexMap(VertexEdgeId highest_mult);
    void populateEdgeCount(CSRTemporalGraph &t_graph,
                          VertexEdgeId search_for_start_pos,
                          VertexEdgeId search_for_end_pos,
                          VertexEdgeId search_from_start_pos,
                          VertexEdgeId search_from_end_pos,
                          VertexEdgeId minus_inf,
                          VertexEdgeId minus_delta,
                          VertexEdgeId plus_delta,
                          VertexEdgeId plus_inf,
                          VertexEdgeId minus_delta1,
                          VertexEdgeId plus_delta1,
                          VertexEdgeId minus_delta2,
                          VertexEdgeId plus_delta2);

    void countTemporalTriangle(CSRGraph &static_graph, CSRTemporalGraph &temporal_graph, TemporalTime delta, TemporalTime delta1, TemporalTime delta2);
    
    Count countCaseA(CSRTemporalGraph &t_graph,
                    VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                    VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                    VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                    VertexEdgeId t1_on_t3_plus_delta,
                    VertexEdgeId t1_on_t3_plus_inf,
                    VertexEdgeId t2_on_t3_plus_inf,
                    VertexEdgeId t2_on_t3_plus_delta2);
    
    Count countCaseB(CSRTemporalGraph &t_graph,
                    VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                    VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                    VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                    VertexEdgeId t1_on_t2_plus_delta1,
                    VertexEdgeId t1_on_t2_minus_inf,
                    VertexEdgeId t3_on_t2_minus_delta2,
                    VertexEdgeId t3_on_t2_minus_inf);
    
    Count countCaseC(CSRTemporalGraph &t_graph,
                    VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                    VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                    VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                    VertexEdgeId t3_on_t1_minus_delta,
                    VertexEdgeId t3_on_t1_minus_inf,
                    VertexEdgeId t2_on_t1_minus_inf,
                    VertexEdgeId t2_on_t1_minus_delta1);
    
    void printCounts();
    void printCountsFile(const char *path);
    void freeMemory();
};  




#endif

