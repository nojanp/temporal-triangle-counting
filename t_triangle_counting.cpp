#include "ETTC/t_triangle_counting.h"

#include <algorithm>

using namespace std;


void MotifCounter::populateEdgeCount(CSRTemporalGraph &t_graph,
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
                          VertexEdgeId plus_delta2)
{

    for (VertexEdgeId t_pos = search_for_start_pos; t_pos < search_for_end_pos; t_pos++)
    {   
        edge_count[minus_inf + t_pos - search_for_start_pos] = 
        edge_count_cum[minus_inf + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeMaxLimitCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos]);

        edge_count[minus_delta + t_pos - search_for_start_pos] = 
        edge_count_cum[minus_delta + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos] - delta_, t_graph.times_[t_pos]);

        edge_count[plus_delta + t_pos - search_for_start_pos] =
        edge_count_cum[plus_delta + t_pos - search_for_start_pos] =
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos], t_graph.times_[t_pos] + delta_);

        edge_count[plus_inf + t_pos - search_for_start_pos] = 
        edge_count_cum[plus_inf + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeMinLimitCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos]);

        edge_count[minus_delta1 + t_pos - search_for_start_pos] = 
        edge_count_cum[minus_delta1 + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos] - delta1_, t_graph.times_[t_pos]);

        edge_count[plus_delta1 + t_pos - search_for_start_pos] = 
        edge_count_cum[plus_delta1 + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos], t_graph.times_[t_pos] + delta1_);

        edge_count[minus_delta2 + t_pos - search_for_start_pos] = 
        edge_count_cum[minus_delta2 + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos] - delta2_, t_graph.times_[t_pos]);

        edge_count[plus_delta2 + t_pos - search_for_start_pos] = 
        edge_count_cum[plus_delta2 + t_pos - search_for_start_pos] = 
            t_graph.edgeTimeIntervalCount(search_from_start_pos, search_from_end_pos, t_graph.times_[t_pos], t_graph.times_[t_pos] + delta2_);
    }

    // cumulative
    for (VertexEdgeId t_pos = search_for_start_pos + 1; t_pos < search_for_end_pos; t_pos++)
    {
        edge_count_cum[minus_inf + t_pos - search_for_start_pos] +=
             edge_count_cum[minus_inf + t_pos - search_for_start_pos - 1];

        edge_count_cum[minus_delta + t_pos - search_for_start_pos] += 
            edge_count_cum[minus_delta + t_pos - search_for_start_pos - 1];

        edge_count_cum[plus_delta + t_pos - search_for_start_pos] +=
            edge_count_cum[plus_delta + t_pos - search_for_start_pos - 1];

        edge_count_cum[plus_inf + t_pos - search_for_start_pos] += 
            edge_count_cum[plus_inf + t_pos - search_for_start_pos - 1];


        edge_count_cum[minus_delta1 + t_pos - search_for_start_pos] +=
            edge_count_cum[minus_delta1 + t_pos - search_for_start_pos - 1];

        edge_count_cum[plus_delta1 + t_pos - search_for_start_pos] +=
            edge_count_cum[plus_delta1 + t_pos - search_for_start_pos - 1];

        edge_count_cum[minus_delta2 + t_pos - search_for_start_pos] +=
            edge_count_cum[minus_delta2 + t_pos - search_for_start_pos - 1];

        edge_count_cum[plus_delta2 + t_pos - search_for_start_pos] +=
            edge_count_cum[plus_delta2 + t_pos - search_for_start_pos - 1];

    }

}



Count MotifCounter::countCaseA(CSRTemporalGraph &t_graph,
                            VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                            VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                            VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                            VertexEdgeId t1_on_t3_plus_delta,
                            VertexEdgeId t1_on_t3_plus_inf,
                            VertexEdgeId t2_on_t3_plus_inf,
                            VertexEdgeId t2_on_t3_plus_delta2)
{
    
    Count sum = 0, sum_1 = 0, sum_2 = 0;


    for (VertexEdgeId t2_pos = t2_start_pos; t2_pos < t2_end_pos; t2_pos++)
    {
        // sum_2
        VertexEdgeId t2_on_t1 = t_graph.edgeTimeIntervalCount(t1_start_pos, t1_end_pos, t_graph.times_[t2_pos] + delta2_ - delta_ + 1, t_graph.times_[t2_pos]);
        sum_2 += t2_on_t1 * edge_count[t2_on_t3_plus_delta2  + t2_pos - t2_start_pos];

        // sum_1

        VertexEdgeId start = t_graph.edgeTimeMinLimitSearch(t1_start_pos, t1_end_pos, t_graph.times_[t2_pos] - delta1_);
        VertexEdgeId finish = t_graph.edgeTimeMaxLimitSearch(t1_start_pos, t1_end_pos, t_graph.times_[t2_pos] + delta2_ - delta_);

        if(start == -1)
        {
            continue;
        }
        if(finish == -1)
        {
            continue;
        }
        start -= t1_start_pos;
        finish -= t1_start_pos;
 


        Count cur_sum = edge_count_cum[t1_on_t3_plus_delta + finish - 1] -
                        edge_count_cum[t1_on_t3_plus_delta + start] +
                        edge_count[t1_on_t3_plus_delta + start];

        cur_sum -= edge_count_cum[t1_on_t3_plus_inf + finish - 1] -
            edge_count_cum[t1_on_t3_plus_inf + start] +
            edge_count[t1_on_t3_plus_inf + start];

        cur_sum += (finish - start) * edge_count[t2_on_t3_plus_inf + t2_pos - t2_start_pos];

        sum_1 += cur_sum;

    }


    sum = sum_1 + sum_2;
    return sum;

}

Count MotifCounter::countCaseB(CSRTemporalGraph &t_graph,
                            VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                            VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                            VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                            VertexEdgeId t1_on_t2_plus_delta1,
                            VertexEdgeId t1_on_t2_minus_inf,
                            VertexEdgeId t3_on_t2_minus_delta2,
                            VertexEdgeId t3_on_t2_minus_inf)
{

    Count sum = 0, sum_1 = 0, sum_2 = 0, sum_3 = 0;


    for (VertexEdgeId t1_pos = t1_start_pos; t1_pos < t1_end_pos; t1_pos++)
    {
        //sum_1

        VertexEdgeId start = t_graph.edgeTimeMinLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t1_pos]);
        VertexEdgeId finish = t_graph.edgeTimeMaxLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t1_pos] + min(delta1_, delta2_));

        if(start == -1)
        {
            continue;
        }
        if(finish == -1)
        {
            continue;
        }
        start -= t3_start_pos;
        finish -= t3_start_pos;

        
        Count cur_sum = edge_count_cum[t3_on_t2_minus_inf + finish - 1] -
                        edge_count_cum[t3_on_t2_minus_inf + start] +
                        edge_count[t3_on_t2_minus_inf + start];

        cur_sum -= (finish - start) * edge_count[t1_on_t2_minus_inf + t1_pos - t1_start_pos];

        cur_sum += (finish - start) * t_graph.edgeTimeIntervalCount(t2_start_pos, t2_end_pos, t_graph.times_[t1_pos], t_graph.times_[t1_pos]);
        
        sum_1 += cur_sum;


        //sum_2

        start = finish;
        finish = t_graph.edgeTimeMinLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t1_pos] + max(delta1_, delta2_));
        if(start == -1)
        {
            continue;
        }
        if(finish == -1)
        {
            continue;
        }

        finish -= t3_start_pos;

        if(delta1_ <= delta2_)
        {
            sum_2 += (finish - start) * edge_count[t1_on_t2_plus_delta1 + t1_pos - t1_start_pos];
        }
        else
        {
            sum_2 += edge_count_cum[t3_on_t2_minus_delta2 + finish - 1] -
                        edge_count_cum[t3_on_t2_minus_delta2 + start] +
                        edge_count[t3_on_t2_minus_delta2 + start];
        }


        start = finish;
        finish = t_graph.edgeTimeMaxLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t1_pos] + delta_);
 
        if(start == -1)
        {
            continue;
        }
        if(finish == -1)
        {
            continue;
        }

        finish -= t3_start_pos;

        cur_sum = edge_count_cum[t3_on_t2_minus_delta2 + finish - 1] -
                        edge_count_cum[t3_on_t2_minus_delta2 + start] +
                        edge_count[t3_on_t2_minus_delta2 + start];

        cur_sum -= edge_count_cum[t3_on_t2_minus_inf + finish - 1] -
                        edge_count_cum[t3_on_t2_minus_inf + start] +
                        edge_count[t3_on_t2_minus_inf + start];

        cur_sum += (finish - start) * t_graph.edgeTimeMaxLimitCount(t2_start_pos, t2_end_pos, t_graph.times_[t1_pos] + delta1_);
        sum_3 += cur_sum;
    }

    sum= sum_1 + sum_2 + sum_3;
    return sum;

}

Count MotifCounter::countCaseC(CSRTemporalGraph &t_graph,
                            VertexEdgeId t1_start_pos, VertexEdgeId t1_end_pos,
                            VertexEdgeId t2_start_pos, VertexEdgeId t2_end_pos,
                            VertexEdgeId t3_start_pos, VertexEdgeId t3_end_pos,
                            VertexEdgeId t3_on_t1_minus_delta,
                            VertexEdgeId t3_on_t1_minus_inf,
                            VertexEdgeId t2_on_t1_minus_inf,
                            VertexEdgeId t2_on_t1_minus_delta1)
{
    
    Count sum = 0, sum_1 = 0, sum_2 = 0;

    for (VertexEdgeId t2_pos = t2_start_pos; t2_pos < t2_end_pos; t2_pos++)
    {

        // sum_2
        VertexEdgeId t2_on_t3 = t_graph.edgeTimeIntervalCount(t3_start_pos, t3_end_pos, t_graph.times_[t2_pos], t_graph.times_[t2_pos] + delta_ - delta1_-1);
        sum_2 += t2_on_t3 * edge_count[t2_on_t1_minus_delta1  + t2_pos - t2_start_pos];


        //sum_1
        VertexEdgeId start = t_graph.edgeTimeMinLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t2_pos] - delta1_ + delta_);
        VertexEdgeId finish = t_graph.edgeTimeMaxLimitSearch(t3_start_pos, t3_end_pos, t_graph.times_[t2_pos] + delta2_);
 
        if(start == -1)
        {
            continue;
        }
        if(finish == -1)
        {
            continue;
        }
        start -= t3_start_pos;
        finish -= t3_start_pos;

        Count cur_sum = edge_count_cum[t3_on_t1_minus_delta + finish - 1] -
                        edge_count_cum[t3_on_t1_minus_delta + start] +
                        edge_count[t3_on_t1_minus_delta + start];
        
        cur_sum -= edge_count_cum[t3_on_t1_minus_inf + finish - 1] -
            edge_count_cum[t3_on_t1_minus_inf + start] +
            edge_count[t3_on_t1_minus_inf + start];
        
        cur_sum += (finish - start) * edge_count[t2_on_t1_minus_inf + t2_pos - t2_start_pos];
        
        sum_1 += cur_sum;

    }

    sum = sum_1 + sum_2;
    return sum;

}


void MotifCounter::countTemporalTriangle(CSRGraph &s_dag, CSRTemporalGraph &t_graph, TemporalTime delta, TemporalTime delta1, TemporalTime delta2)
{
    CSRGraph mult_graph = t_graph.extractMultGraph();
    VertexEdgeId highest_mult = 0;
    for (VertexEdgeId i = 0; i < mult_graph.num_edges_; i++)
        highest_mult = max(highest_mult, mult_graph.temporal_start_pos_[i+1] - mult_graph.temporal_start_pos_[i]);
    highest_mult_ = highest_mult;


    // reload edge times from temporal_Nbrs_ to times_
    t_graph.reloadTimes();


    edge_count = new VertexEdgeId[64 * highest_mult];
    edge_count_cum = new VertexEdgeId[64 * highest_mult];
    delta_ = delta;
    delta1_ = delta1;
    delta2_ = delta2;

    VertexEdgeId index_array[64];
    for(int i=0; i < 64; i++)
        index_array[i] = i * highest_mult;


    for (VertexEdgeId i=0; i < s_dag.num_vertices_; i++)
    {
        VertexEdgeId u = i;
        for (VertexEdgeId j = s_dag.offsets_[i]; j < s_dag.offsets_[i+1]; j++)
        {
            VertexEdgeId v = s_dag.nbrs_[j];
            for (VertexEdgeId k = j+1; k < s_dag.offsets_[i+1]; k++)
            {
                VertexEdgeId w = s_dag.nbrs_[k];
                if(w == v)
                    continue;

                VertexEdgeId w_index = s_dag.getEdgeIndx(v,w);
                if(w_index == -1)
                    continue;
                
                static_triangles_count_++;
                              

                VertexEdgeId mult_graph_indx_u_v = mult_graph.getEdgeIndx(u, v);
                VertexEdgeId mult_graph_indx_u_w = mult_graph.getEdgeIndx(u, w);
                VertexEdgeId mult_graph_indx_v_u = mult_graph.getEdgeIndx(v, u);
                VertexEdgeId mult_graph_indx_v_w = mult_graph.getEdgeIndx(v, w);
                VertexEdgeId mult_graph_indx_w_u = mult_graph.getEdgeIndx(w, u);
                VertexEdgeId mult_graph_indx_w_v = mult_graph.getEdgeIndx(w, v);


                VertexEdgeId start_pos_u_v = mult_graph_indx_u_v == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_u_v];
                VertexEdgeId start_pos_u_w = mult_graph_indx_u_w == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_u_w];
                VertexEdgeId start_pos_v_u = mult_graph_indx_v_u == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_v_u];
                VertexEdgeId start_pos_v_w = mult_graph_indx_v_w == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_v_w];
                VertexEdgeId start_pos_w_u = mult_graph_indx_w_u == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_w_u];
                VertexEdgeId start_pos_w_v = mult_graph_indx_w_v == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_w_v];


                // one past the actual end index
                VertexEdgeId end_pos_u_v = mult_graph_indx_u_v == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_u_v+1];
                VertexEdgeId end_pos_u_w = mult_graph_indx_u_w == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_u_w+1];
                VertexEdgeId end_pos_v_u = mult_graph_indx_v_u == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_v_u+1];
                VertexEdgeId end_pos_v_w = mult_graph_indx_v_w == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_v_w+1];
                VertexEdgeId end_pos_w_u = mult_graph_indx_w_u == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_w_u+1];
                VertexEdgeId end_pos_w_v = mult_graph_indx_w_v == -1 ? -1 : mult_graph.temporal_start_pos_[mult_graph_indx_w_v+1];

                // the actual temporal multiplicity of each static edge
                VertexEdgeId mult_u_v = start_pos_u_v == -1 ? 0 : end_pos_u_v - start_pos_u_v;
                VertexEdgeId mult_u_w = start_pos_u_w == -1 ? 0 : end_pos_u_w - start_pos_u_w;
                VertexEdgeId mult_v_u = start_pos_v_u == -1 ? 0 : end_pos_v_u - start_pos_v_u;
                VertexEdgeId mult_v_w = start_pos_v_w == -1 ? 0 : end_pos_v_w - start_pos_v_w;
                VertexEdgeId mult_w_u = start_pos_w_u == -1 ? 0 : end_pos_w_u - start_pos_w_u;
                VertexEdgeId mult_w_v = start_pos_w_v == -1 ? 0 : end_pos_w_v - start_pos_w_v;
                

                // we are looking at a triangle u,v,w. Directions in the degeneracy DAG are (u,v), (u,w), and (v,w)
                
                for(int p = 0; p < 6; p++)
                    for(int q = 0; q < 8; q++)
                        triangle_motif_counts_[p][q] = 0;
 
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                // Here we populate the auxilary arrays for edge interval counting
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                // for each temporal edge with time t in (u,v) from (v,w)
                populateEdgeCount(t_graph, start_pos_u_v, end_pos_u_v, start_pos_v_w, end_pos_v_w,
                                index_array[0], index_array[1], index_array[2], index_array[3], index_array[32], index_array[40], index_array[48], index_array[56]);

                // for each temporal edge with time t in (u,v) from (w,v)
                populateEdgeCount(t_graph, start_pos_u_v, end_pos_u_v, start_pos_w_v, end_pos_w_v,
                                index_array[4], index_array[5], index_array[6], index_array[7], index_array[33], index_array[41], index_array[49], index_array[57]);

                // for each temporal edge with time t in (v,u) from (v,w)
                populateEdgeCount(t_graph, start_pos_v_u, end_pos_v_u, start_pos_v_w, end_pos_v_w,
                                index_array[8], index_array[9], index_array[10], index_array[11], index_array[34], index_array[42], index_array[50], index_array[58]);

                // for each temporal edge with time t in (v,u) from (w,v)
                populateEdgeCount(t_graph, start_pos_v_u, end_pos_v_u, start_pos_w_v, end_pos_w_v,
                                index_array[12], index_array[13], index_array[14], index_array[15], index_array[35], index_array[43], index_array[51], index_array[59]);

                // for each temporal edge with time t in (u,w) from (v,w)
                populateEdgeCount(t_graph, start_pos_u_w, end_pos_u_w, start_pos_v_w, end_pos_v_w,
                                index_array[16], index_array[17], index_array[18], index_array[19], index_array[36], index_array[44], index_array[52], index_array[60]);

                // for each temporal edge with time t in (u,w) from (w,v)
                populateEdgeCount(t_graph, start_pos_u_w, end_pos_u_w, start_pos_w_v, end_pos_w_v,
                                index_array[20], index_array[21], index_array[22], index_array[23], index_array[37], index_array[45], index_array[53], index_array[61]);

                // for each temporal edge with time t in (w,u) from (v,w)
                populateEdgeCount(t_graph, start_pos_w_u, end_pos_w_u, start_pos_v_w, end_pos_v_w,
                                index_array[24], index_array[25], index_array[26], index_array[27], index_array[38], index_array[46], index_array[54], index_array[62]);

                // for each temporal edge with time t in (w,u) from (w,v)
                populateEdgeCount(t_graph, start_pos_w_u, end_pos_w_u, start_pos_w_v, end_pos_w_v,
                                index_array[28], index_array[29], index_array[30], index_array[31], index_array[39], index_array[47], index_array[55], index_array[63]);


                // %%%%%%%%%%%%%%%%%%%%%%%%%                
                //Case A1 (\pi_1)
                // %%%%%%%%%%%%%%%%%%%%%%%%%

                /////////////////////////////////////////////////////////////////////
                //A1 Dir D1: M7
                triangle_motif_counts_[0][0] = countCaseA(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_w, end_pos_v_w,
                            index_array[2], index_array[3], index_array[19], index_array[60]);
                // cout << "case A1 D1 : M7 addition: " << triangle_motif_counts_[0][0] << endl;
                motif_counts_[0][0] += triangle_motif_counts_[0][0];

                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D2: M6
                triangle_motif_counts_[0][1] = countCaseA(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_w_v, end_pos_w_v,
                            index_array[6], index_array[7], index_array[31], index_array[63]);
                // cout << "case A1 D2 : M6 addition: " << triangle_motif_counts_[0][1] << endl;
                motif_counts_[0][1] += triangle_motif_counts_[0][1];

                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D3: M5
                triangle_motif_counts_[0][2] = countCaseA(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_w_v, end_pos_w_v,
                            index_array[6], index_array[7], index_array[23], index_array[61]);
                // cout << "case A1 D3 : M5 addition: " << triangle_motif_counts_[0][2] << endl;
                motif_counts_[0][2] += triangle_motif_counts_[0][2];
             
                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D4: M8
                triangle_motif_counts_[0][3] = countCaseA(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_w, end_pos_v_w,
                            index_array[2], index_array[3], index_array[27], index_array[62]);
                // cout << "case A1 D4 : M8 addition: " << triangle_motif_counts_[0][3] << endl;
                motif_counts_[0][3] += triangle_motif_counts_[0][3];
             
                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D5: M3
                triangle_motif_counts_[0][4] = countCaseA(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_w, end_pos_v_w,
                            index_array[10], index_array[11], index_array[19], index_array[60]);
                // cout << "case A1 D5 : M3 addition: " << triangle_motif_counts_[0][4] << endl;
                motif_counts_[0][4] += triangle_motif_counts_[0][4];
             
                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D6: M2
                triangle_motif_counts_[0][5] = countCaseA(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_w_v, end_pos_w_v,
                            index_array[14], index_array[15], index_array[31], index_array[63]);
                // cout << "case A1 D6 : M2 addition: " << triangle_motif_counts_[0][5] << endl;
                motif_counts_[0][5] += triangle_motif_counts_[0][5];
             
                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D7: M4
                triangle_motif_counts_[0][6] = countCaseA(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_w_v, end_pos_w_v,
                            index_array[14], index_array[15], index_array[23], index_array[61]);
                // cout << "case A1 D7 : M4 addition: " << triangle_motif_counts_[0][6] << endl;
                motif_counts_[0][6] += triangle_motif_counts_[0][6];
             
                // /////////////////////////////////////////////////////////////////////
                // //A1 Dir D8: M1
                triangle_motif_counts_[0][7] = countCaseA(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_w, end_pos_v_w,
                            index_array[10], index_array[11], index_array[27], index_array[62]);
                // cout << "case A1 D8 : M1 addition: " << triangle_motif_counts_[0][7] << endl;
                motif_counts_[0][7] += triangle_motif_counts_[0][7];

                // %%%%%%%%%%%%%%%%%%%%%%%%%
                //Case A2 (\pi_2)
                // %%%%%%%%%%%%%%%%%%%%%%%%%                

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D1: M5
                triangle_motif_counts_[1][0] = countCaseA(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_v_w, end_pos_v_w,
                            index_array[18], index_array[19], index_array[3], index_array[56]);
                // cout << "case A2 D1 : M5 addition: " << triangle_motif_counts_[1][0] << endl;
                motif_counts_[1][0] += triangle_motif_counts_[1][0];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D2: M3
                triangle_motif_counts_[1][1] = countCaseA(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_v, end_pos_w_v,
                            index_array[30], index_array[31], index_array[7], index_array[57]);
                // cout << "case A2 D2 : M3 addition: " << triangle_motif_counts_[1][1] << endl;
                motif_counts_[1][1] += triangle_motif_counts_[1][1];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D3: M7
                triangle_motif_counts_[1][2] = countCaseA(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_v, end_pos_w_v,
                            index_array[22], index_array[23], index_array[7], index_array[57]);
                // cout << "case A2 D3 : M7 addition: " << triangle_motif_counts_[1][2] << endl;
                motif_counts_[1][2] += triangle_motif_counts_[1][2];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D4: M4
                triangle_motif_counts_[1][3] = countCaseA(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_v_w, end_pos_v_w,
                            index_array[26], index_array[27], index_array[3], index_array[56]);
                // cout << "case A2 D4 : M4 addition: " << triangle_motif_counts_[1][3] << endl;
                motif_counts_[1][3] += triangle_motif_counts_[1][3];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D5: M6
                triangle_motif_counts_[1][4] = countCaseA(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_v_w, end_pos_v_w,
                            index_array[18], index_array[19], index_array[11], index_array[58]);
                // cout << "case A2 D5 : M6 addition: " << triangle_motif_counts_[1][4] << endl;
                motif_counts_[1][4] += triangle_motif_counts_[1][4];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D6: M1
                triangle_motif_counts_[1][5] = countCaseA(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_v, end_pos_w_v,
                            index_array[30], index_array[31], index_array[15], index_array[59]);
                // cout << "case A2 D6 : M1 addition: " << triangle_motif_counts_[1][5] << endl;
                motif_counts_[1][5] += triangle_motif_counts_[1][5];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D7: M8
                triangle_motif_counts_[1][6] = countCaseA(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_v, end_pos_w_v,
                            index_array[22], index_array[23], index_array[15], index_array[59]);
                // cout << "case A2 D7 : M8 addition: " << triangle_motif_counts_[1][6] << endl;
                motif_counts_[1][6] += triangle_motif_counts_[1][6];

                /////////////////////////////////////////////////////////////////////
                //A2 Dir D8: M2
                triangle_motif_counts_[1][7] = countCaseA(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_v_w, end_pos_v_w,
                            index_array[26], index_array[27], index_array[11], index_array[58]);
                // cout << "case A2 D8 : M2 addition: " << triangle_motif_counts_[1][7] << endl;
                motif_counts_[1][7] += triangle_motif_counts_[1][7];
                

                // %%%%%%%%%%%%%%%%%%%%%%%%%
                //Case B1 (\pi_3)
                // %%%%%%%%%%%%%%%%%%%%%%%%%                

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D1: M3
                triangle_motif_counts_[2][0] = countCaseB(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_w, end_pos_u_w,
                            index_array[40], index_array[0], index_array[52], index_array[16]);
                // cout << "case B1 D1 : M3 addition: " << triangle_motif_counts_[2][0] << endl;
                motif_counts_[2][0] += triangle_motif_counts_[2][0];


                /////////////////////////////////////////////////////////////////////
                //B1 Dir D2: M2
                triangle_motif_counts_[2][1] = countCaseB(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_w_u, end_pos_w_u,
                            index_array[41], index_array[4], index_array[55], index_array[28]);
                // cout << "case B1 D2 : M2 addition: " << triangle_motif_counts_[2][1] << endl;
                motif_counts_[2][1] += triangle_motif_counts_[2][1];


                /////////////////////////////////////////////////////////////////////
                //B1 Dir D3: M1
                triangle_motif_counts_[2][2] = countCaseB(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_w, end_pos_u_w,
                            index_array[41], index_array[4], index_array[53], index_array[20]);
                // cout << "case B1 D3 : M1 addition: " << triangle_motif_counts_[2][2] << endl;
                motif_counts_[2][2] += triangle_motif_counts_[2][2];
             

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D4: M4
                triangle_motif_counts_[2][3] = countCaseB(t_graph,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_w_u, end_pos_w_u,
                            index_array[40], index_array[0], index_array[54], index_array[24]);
                // cout << "case B1 D4 : M4 addition: " << triangle_motif_counts_[2][3] << endl;
                motif_counts_[2][3] += triangle_motif_counts_[2][3];

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D5: M7
                triangle_motif_counts_[2][4] = countCaseB(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_w, end_pos_u_w,
                            index_array[42], index_array[8], index_array[52], index_array[16]);
                // cout << "case B1 D5 : M7 addition: " << triangle_motif_counts_[2][4] << endl;
                motif_counts_[2][4] += triangle_motif_counts_[2][4];
             

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D6: M6
                triangle_motif_counts_[2][5] = countCaseB(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_w_u, end_pos_w_u,
                            index_array[43], index_array[12], index_array[55], index_array[28]);
                // cout << "case B1 D6 : M6 addition: " << triangle_motif_counts_[2][5] << endl;
                motif_counts_[2][5] += triangle_motif_counts_[2][5];

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D7: M8
                triangle_motif_counts_[2][6] = countCaseB(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_w, end_pos_u_w,
                            index_array[43], index_array[12], index_array[53], index_array[20]);
                // cout << "case B1 D7 : M8 addition: " << triangle_motif_counts_[2][6] << endl;
                motif_counts_[2][6] += triangle_motif_counts_[2][6];
             

                /////////////////////////////////////////////////////////////////////
                //B1 Dir D8: M5
                triangle_motif_counts_[2][7] = countCaseB(t_graph,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_w_u, end_pos_w_u,
                            index_array[42], index_array[8], index_array[54], index_array[24]);
                // cout << "case B1 D8 : M5 addition: " << triangle_motif_counts_[2][7] << endl;
                motif_counts_[2][7] += triangle_motif_counts_[2][7];


                // %%%%%%%%%%%%%%%%%%%%%%%%%
                //Case B2 (\pi_4)
                // %%%%%%%%%%%%%%%%%%%%%%%%%                
                
                /////////////////////////////////////////////////////////////////////
                //B2 Dir D1: M1
                triangle_motif_counts_[3][0] = countCaseB(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_v, end_pos_u_v,
                            index_array[44], index_array[16], index_array[48], index_array[0]);
                // cout << "case B2 D1 : M1 addition: " << triangle_motif_counts_[3][0] << endl;
                motif_counts_[3][0] += triangle_motif_counts_[3][0];


                /////////////////////////////////////////////////////////////////////
                //B2 Dir D2: M7
                triangle_motif_counts_[3][1] = countCaseB(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_v, end_pos_u_v,
                            index_array[47], index_array[28], index_array[49], index_array[4]);
                // cout << "case B2 D2 : M7 addition: " << triangle_motif_counts_[3][1] << endl;
                motif_counts_[3][1] += triangle_motif_counts_[3][1];

                /////////////////////////////////////////////////////////////////////
                //B2 Dir D3: M3
                triangle_motif_counts_[3][2] = countCaseB(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_v, end_pos_u_v,
                            index_array[45], index_array[20], index_array[49], index_array[4]);
                // cout << "case B2 D3 : M3 addition: " << triangle_motif_counts_[3][2] << endl;
                motif_counts_[3][2] += triangle_motif_counts_[3][2];

                /////////////////////////////////////////////////////////////////////
                //B2 Dir D4: M8
                triangle_motif_counts_[3][3] = countCaseB(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_v, end_pos_u_v,
                            index_array[46], index_array[24], index_array[48], index_array[0]);
                // cout << "case B2 D4 : M8 addition: " << triangle_motif_counts_[3][3] << endl;
                motif_counts_[3][3] += triangle_motif_counts_[3][3];

                /////////////////////////////////////////////////////////////////////
                //B2 Dir D5: M2
                triangle_motif_counts_[3][4] = countCaseB(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_v_u, end_pos_v_u,
                            index_array[44], index_array[16], index_array[50], index_array[8]);
                // cout << "case B2 D5 : M2 addition: " << triangle_motif_counts_[3][4] << endl;
                motif_counts_[3][4] += triangle_motif_counts_[3][4];

                /////////////////////////////////////////////////////////////////////
                //B2 Dir D6: M5
                triangle_motif_counts_[3][5] = countCaseB(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_v_u, end_pos_v_u,
                            index_array[47], index_array[28], index_array[51], index_array[12]);
                // cout << "case B2 D6 : M5 addition: " << triangle_motif_counts_[3][5] << endl;
                motif_counts_[3][5] += triangle_motif_counts_[3][5];


                /////////////////////////////////////////////////////////////////////
                //B2 Dir D7: M4
                triangle_motif_counts_[3][6] = countCaseB(t_graph,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_v_u, end_pos_v_u,
                            index_array[45], index_array[20], index_array[51], index_array[12]);
                // cout << "case B2 D7 : M4 addition: " << triangle_motif_counts_[3][6] << endl;
                motif_counts_[3][6] += triangle_motif_counts_[3][6];


                /////////////////////////////////////////////////////////////////////
                //B2 Dir D8: M6
                triangle_motif_counts_[3][7] = countCaseB(t_graph,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_v_u, end_pos_v_u,
                            index_array[46], index_array[24], index_array[50], index_array[8]);
                // cout << "case B2 D8 : M6 addition: " << triangle_motif_counts_[3][7] << endl;
                motif_counts_[3][7] += triangle_motif_counts_[3][7];


                // %%%%%%%%%%%%%%%%%%%%%%%%%
                //Case C1 (\pi_5)
                // %%%%%%%%%%%%%%%%%%%%%%%%%                

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D1: M6
                triangle_motif_counts_[4][0] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_u_w, end_pos_u_w,
                            index_array[17], index_array[16], index_array[0], index_array[32]);
                // cout << "case C1 D1 : M6 addition: " << triangle_motif_counts_[4][0] << endl;
                motif_counts_[4][0] += triangle_motif_counts_[4][0];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D2: M1
                triangle_motif_counts_[4][1] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_u, end_pos_w_u,
                            index_array[29], index_array[28], index_array[4], index_array[33]);
                // cout << "case C1 D2 : M1 addition: " << triangle_motif_counts_[4][1] << endl;
                motif_counts_[4][1] += triangle_motif_counts_[4][1];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D3: M2
                triangle_motif_counts_[4][2] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_u_w, end_pos_u_w,
                            index_array[21], index_array[20], index_array[4], index_array[33]);
                // cout << "case C1 D3 : M2 addition: " << triangle_motif_counts_[4][2] << endl;
                motif_counts_[4][2] += triangle_motif_counts_[4][2];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D4: M8
                triangle_motif_counts_[4][3] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_v, end_pos_u_v,
                            start_pos_w_u, end_pos_w_u,
                            index_array[25], index_array[24], index_array[0], index_array[32]);
                // cout << "case C1 D4 : M8 addition: " << triangle_motif_counts_[4][3] << endl;
                motif_counts_[4][3] += triangle_motif_counts_[4][3];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D5: M5
                triangle_motif_counts_[4][4] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_u_w, end_pos_u_w,
                            index_array[17], index_array[16], index_array[8], index_array[34]);
                // cout << "case C1 D5 : M5 addition: " << triangle_motif_counts_[4][4] << endl;
                motif_counts_[4][4] += triangle_motif_counts_[4][4];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D6: M3
                triangle_motif_counts_[4][5] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_u, end_pos_w_u,
                            index_array[29], index_array[28], index_array[12], index_array[35]);
                // cout << "case C1 D6 : M3 addition: " << triangle_motif_counts_[4][5] << endl;
                motif_counts_[4][5] += triangle_motif_counts_[4][5];
               

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D7: M4
                triangle_motif_counts_[4][6] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_u_w, end_pos_u_w,
                            index_array[21], index_array[20], index_array[12], index_array[35]);
                // cout << "case C1 D7 : M4 addition: " << triangle_motif_counts_[4][6] << endl;
                motif_counts_[4][6] += triangle_motif_counts_[4][6];

                /////////////////////////////////////////////////////////////////////
                //C1 Dir D8: M7
                triangle_motif_counts_[4][7] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_v_u, end_pos_v_u,
                            start_pos_w_u, end_pos_w_u,
                            index_array[25], index_array[24], index_array[8], index_array[34]);
                // cout << "case C1 D8 : M7 addition: " << triangle_motif_counts_[4][7] << endl;
                motif_counts_[4][7] += triangle_motif_counts_[4][7];



                // %%%%%%%%%%%%%%%%%%%%%%%%%
                //Case C2 (\pi_6)
                // %%%%%%%%%%%%%%%%%%%%%%%%%                

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D1: M2
                triangle_motif_counts_[5][0] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_u_v, end_pos_u_v,
                            index_array[1], index_array[0], index_array[16], index_array[36]);
                // cout << "case C2 D1 : M2 addition: " << triangle_motif_counts_[5][0] << endl;
                motif_counts_[5][0] += triangle_motif_counts_[5][0];

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D2: M5
                triangle_motif_counts_[5][1] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_u_v, end_pos_u_v,
                            index_array[5], index_array[4], index_array[28], index_array[39]);
                // cout << "case C2 D2 : M5 addition: " << triangle_motif_counts_[5][1] << endl;
                motif_counts_[5][1] += triangle_motif_counts_[5][1];

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D3: M6
                triangle_motif_counts_[5][2] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_u_v, end_pos_u_v,
                            index_array[5], index_array[4], index_array[20], index_array[37]);
                // cout << "case C2 D3 : M6 addition: " << triangle_motif_counts_[5][2] << endl;
                motif_counts_[5][2] += triangle_motif_counts_[5][2];

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D4: M4
                triangle_motif_counts_[5][3] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_u_v, end_pos_u_v,
                            index_array[1], index_array[0], index_array[24], index_array[38]);
                // cout << "case C2 D4 : M4 addition: " << triangle_motif_counts_[5][3] << endl;
                motif_counts_[5][3] += triangle_motif_counts_[5][3];


                /////////////////////////////////////////////////////////////////////
                //C2 Dir D5: M1
                triangle_motif_counts_[5][4] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_u, end_pos_v_u,
                            index_array[9], index_array[8], index_array[16], index_array[36]);
                // cout << "case C2 D5 : M1 addition: " << triangle_motif_counts_[5][4] << endl;
                motif_counts_[5][4] += triangle_motif_counts_[5][4];


                /////////////////////////////////////////////////////////////////////
                //C2 Dir D6: M7
                triangle_motif_counts_[5][5] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_u, end_pos_v_u,
                            index_array[13], index_array[12], index_array[28], index_array[39]);
                // cout << "case C2 D6 : M7 addition: " << triangle_motif_counts_[5][5] << endl;
                motif_counts_[5][5] += triangle_motif_counts_[5][5];

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D7: M8
                triangle_motif_counts_[5][6] = countCaseC(t_graph,
                            start_pos_w_v, end_pos_w_v,
                            start_pos_u_w, end_pos_u_w,
                            start_pos_v_u, end_pos_v_u,
                            index_array[13], index_array[12], index_array[20], index_array[37]);
                // cout << "case C2 D7 : M8 addition: " << triangle_motif_counts_[5][6] << endl;
                motif_counts_[5][6] += triangle_motif_counts_[5][6];

                /////////////////////////////////////////////////////////////////////
                //C2 Dir D8: M3
                triangle_motif_counts_[5][7] = countCaseC(t_graph,
                            start_pos_v_w, end_pos_v_w,
                            start_pos_w_u, end_pos_w_u,
                            start_pos_v_u, end_pos_v_u,
                            index_array[9], index_array[8], index_array[24], index_array[38]);
                // cout << "case C2 D8 : M3 addition: " << triangle_motif_counts_[5][7] << endl;
                motif_counts_[5][7] += triangle_motif_counts_[5][7];


                /////////////////////
                // adjustemnts for D4
                /////////////////////
                for (VertexEdgeId t_pos = start_pos_u_v; t_pos < end_pos_u_v; t_pos++)
                {
                    VertexEdgeId equals_w_u = t_graph.edgeTimeIntervalCount(start_pos_w_u, end_pos_w_u, t_graph.times_[t_pos], t_graph.times_[t_pos]);
                    if(equals_w_u > 0)
                    {
                        VertexEdgeId equals_v_w = t_graph.edgeTimeIntervalCount(start_pos_v_w, end_pos_v_w, t_graph.times_[t_pos], t_graph.times_[t_pos]);
                        if(equals_v_w > 0)
                        {
                            motif_counts_[2][3] -= equals_w_u * equals_v_w;
                            motif_counts_[3][3] -= equals_w_u * equals_v_w;
                            motif_counts_[4][3] -= equals_w_u * equals_v_w;
                            motif_counts_[5][3] -= equals_w_u * equals_v_w;
                        }
                    }
                }


                /////////////////////
                // adjustemnts for D7
                /////////////////////
                for (VertexEdgeId t_pos = start_pos_v_u; t_pos < end_pos_v_u; t_pos++)
                {
                    VertexEdgeId equals_u_w = t_graph.edgeTimeIntervalCount(start_pos_u_w, end_pos_u_w, t_graph.times_[t_pos], t_graph.times_[t_pos]);
                    
                    if(equals_u_w > 0)
                    {
                        VertexEdgeId equals_w_v = t_graph.edgeTimeIntervalCount(start_pos_w_v, end_pos_w_v, t_graph.times_[t_pos], t_graph.times_[t_pos]);
                        if(equals_w_v > 0)
                        {
                            motif_counts_[2][6] -= equals_u_w * equals_w_v;
                            motif_counts_[3][6] -= equals_u_w * equals_w_v;
                            motif_counts_[4][6] -= equals_u_w * equals_w_v;
                            motif_counts_[5][6] -= equals_u_w * equals_w_v;
                        }
                    }
                }
            }
        }
    }


    // // case A1 (\pi_1)
    // cout << endl;
    // cout << "case A1:" << endl;
    // cout << "D1: " << motif_counts_[0][0] << "  D2: " << motif_counts_[0][1] << "  D3: " << motif_counts_[0][2] << "  D4: " << motif_counts_[0][3] << 
    // "  D5: " << motif_counts_[0][4] << "  D6: " << motif_counts_[0][5] << "  D7: " << motif_counts_[0][6] << "  D8: " << motif_counts_[0][7] << endl;

    // // case A2 (\pi_2)
    // cout << endl;
    // cout << "case A2:" << endl;
    // cout << "D1: " << motif_counts_[1][0] << "  D2: " << motif_counts_[1][1] << "  D3: " << motif_counts_[1][2] << "  D4: " << motif_counts_[1][3] << 
    // "  D5: " << motif_counts_[1][4] << "  D6: " << motif_counts_[1][5] << "  D7: " << motif_counts_[1][6] << "  D8: " << motif_counts_[1][7] << endl;

    // // case B1 (\pi_3)
    // cout << endl;
    // cout << "case B1:" << endl;
    // cout << "D1: " << motif_counts_[2][0] << "  D2: " << motif_counts_[2][1] << "  D3: " << motif_counts_[2][2] << "  D4: " << motif_counts_[2][3] << 
    // "  D5: " << motif_counts_[2][4] << "  D6: " << motif_counts_[2][5] << "  D7: " << motif_counts_[2][6] << "  D8: " << motif_counts_[2][7] << endl;


    // // case B2 (\pi_4)
    // cout << endl;
    // cout << "case B2:" << endl;
    // cout << "D1: " << motif_counts_[3][0] << "  D2: " << motif_counts_[3][1] << "  D3: " << motif_counts_[3][2] << "  D4: " << motif_counts_[3][3] << 
    // "  D5: " << motif_counts_[3][4] << "  D6: " << motif_counts_[3][5] << "  D7: " << motif_counts_[3][6] << "  D8: " << motif_counts_[3][7] << endl;


    // // case C1 (\pi_5)
    // cout << endl;
    // cout << "case C1:" << endl;
    // cout << "D1: " << motif_counts_[4][0] << "  D2: " << motif_counts_[4][1] << "  D3: " << motif_counts_[4][2] << "  D4: " << motif_counts_[4][3] << 
    // "  D5: " << motif_counts_[4][4] << "  D6: " << motif_counts_[4][5] << "  D7: " << motif_counts_[4][6] << "  D8: " << motif_counts_[4][7] << endl;


    // // case C2 (\pi_6)
    // cout << endl;
    // cout << "case C2:" << endl;
    // cout << "D1: " << motif_counts_[5][0] << "  D2: " << motif_counts_[5][1] << "  D3: " << motif_counts_[5][2] << "  D4: " << motif_counts_[5][3] << 
    // "  D5: " << motif_counts_[5][4] << "  D6: " << motif_counts_[5][5] << "  D7: " << motif_counts_[5][6] << "  D8: " << motif_counts_[5][7] << endl;


    motif_type_counts_[0] = motif_counts_[0][7] + motif_counts_[1][5] + motif_counts_[2][2] + motif_counts_[3][0] + motif_counts_[4][1] + motif_counts_[5][4];
    motif_type_counts_[1] = motif_counts_[0][5] + motif_counts_[1][7] + motif_counts_[2][1] + motif_counts_[3][4] + motif_counts_[4][2] + motif_counts_[5][0];
    motif_type_counts_[2] = motif_counts_[0][4] + motif_counts_[1][1] + motif_counts_[2][0] + motif_counts_[3][2] + motif_counts_[4][5] + motif_counts_[5][7];
    motif_type_counts_[3] = motif_counts_[0][6] + motif_counts_[1][3] + motif_counts_[2][3] + motif_counts_[3][6] + motif_counts_[4][6] + motif_counts_[5][3];
    motif_type_counts_[4] = motif_counts_[0][2] + motif_counts_[1][0] + motif_counts_[2][7] + motif_counts_[3][5] + motif_counts_[4][4] + motif_counts_[5][1];
    motif_type_counts_[5] = motif_counts_[0][1] + motif_counts_[1][4] + motif_counts_[2][5] + motif_counts_[3][7] + motif_counts_[4][0] + motif_counts_[5][2];
    motif_type_counts_[6] = motif_counts_[0][0] + motif_counts_[1][2] + motif_counts_[2][4] + motif_counts_[3][1] + motif_counts_[4][7] + motif_counts_[5][5];
    motif_type_counts_[7] = motif_counts_[0][3] + motif_counts_[1][6] + motif_counts_[2][6] + motif_counts_[3][3] + motif_counts_[4][3] + motif_counts_[5][6];
}


void MotifCounter::printCounts()
{
    cout << endl;
    cout << "counts of temporal triangles for each motif type:" << endl;
    for(int i=0; i < 8; i++)
        cout << "type " << i+1 << ": " << motif_type_counts_[i] << endl;
}

void MotifCounter::printCountsFile(const char *path)
{   
    FILE* output_file = fopen(path, "w");
    for(int i=0; i < 8; i++)
    {
        fprintf(output_file, "%ld", motif_type_counts_[i]);
        fprintf(output_file, "\n");
    }
    fclose(output_file);
}

void MotifCounter::freeMemory()
{
    delete[] edge_count; 
    delete[] edge_count_cum;
}




