#include "ETTC/Graph.h"

#include <algorithm>

using namespace std;

// Graph class
Graph::Graph(VertexEdgeId num_vertices, VertexEdgeId num_edges): num_vertices_(num_vertices), num_edges_(num_edges)
{
    edges_ = new Edge[num_edges_];
}


void Graph::sortEdges()
{
  sort(edges_, edges_ + num_edges_); 
}

CSRGraph Graph::convertToCSR()
{
  sortEdges();
  CSRGraph csr_graph(num_vertices_, num_edges_);

  VertexEdgeId cur_vertex = 0;
  for(VertexEdgeId i = 0; i < num_edges_; i++)
  {
    csr_graph.nbrs_[i] = edges_[i].dst_;
    while(cur_vertex <= edges_[i].src_)
      csr_graph.offsets_[cur_vertex++] = i;
  }
  while(cur_vertex <= num_vertices_)
    csr_graph.offsets_[cur_vertex++] = num_edges_;

  return csr_graph;
}

void Graph::printGraph()
{
  for(VertexEdgeId i = 0; i < num_edges_; i++)
  {
    cout << edges_[i].src_ << " " << edges_[i].dst_ << endl;
  }
}

void Graph::deleteGraph()
{
  delete[] edges_;
}


// TemporalGraph class

TemporalGraph::TemporalGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges): num_vertices_(num_vertices), num_edges_(num_edges)
{
  temporal_edges_ = new TemporalEdge[num_edges_];
}

void TemporalGraph::sortEdges()
{
  sort(temporal_edges_, temporal_edges_ + num_edges_);
}


VertexEdgeId TemporalGraph::numberOfStaticDirectedEdges()
{
  VertexEdgeId last_src = -1, last_dsts = -1;
  VertexEdgeId num_static_edges = 0;
  for (int i=0; i < num_edges_; i++)
  {
    if(temporal_edges_[i].src_ != last_src || temporal_edges_[i].dst_ != last_dsts)
    {
      num_static_edges++; 
    }
    last_src = temporal_edges_[i].src_;
    last_dsts = temporal_edges_[i].dst_;
  }
  return num_static_edges;
}

Graph TemporalGraph::ExtractStaticGraph()
{
  sortEdges();
  VertexEdgeId num_middler_edges = numberOfStaticDirectedEdges();
  Graph middler(num_vertices_, num_middler_edges);
  VertexEdgeId last_src = -1, last_dsts = -1;
  VertexEdgeId edge_index = 0;
  for (VertexEdgeId i=0; i < num_edges_; i++)
  {
    if(temporal_edges_[i].src_ != last_src || temporal_edges_[i].dst_ != last_dsts)
    {
      middler.edges_[edge_index].src_ = min(temporal_edges_[i].src_, temporal_edges_[i].dst_);
      middler.edges_[edge_index].dst_ = max(temporal_edges_[i].src_, temporal_edges_[i].dst_);
      edge_index++;
    }
    last_src = temporal_edges_[i].src_;
    last_dsts = temporal_edges_[i].dst_;
  }
  middler.sortEdges();


  // find number of static edges
  VertexEdgeId num_graph_edges = 0;
  for (VertexEdgeId i=0; i < middler.num_edges_; i++)
  {
    if(middler.edges_[i].src_ != last_src || middler.edges_[i].dst_ != last_dsts)
    {
      num_graph_edges++;
    }
    last_src = middler.edges_[i].src_;
    last_dsts = middler.edges_[i].dst_;
  }


  // middler to graph
  Graph graph;
  graph.num_vertices_ = num_vertices_;
  graph.num_edges_ = 2 * num_graph_edges;
  graph.edges_ = new Edge[graph.num_edges_];
  edge_index = 0;
  for (VertexEdgeId i=0; i < middler.num_edges_; i++)
  {
    if(middler.edges_[i].src_ != last_src || middler.edges_[i].dst_ != last_dsts)
    {
      graph.edges_[edge_index].src_ = middler.edges_[i].src_;
      graph.edges_[edge_index].dst_ = middler.edges_[i].dst_;
      edge_index++;

      graph.edges_[edge_index].src_ = middler.edges_[i].dst_;
      graph.edges_[edge_index].dst_ = middler.edges_[i].src_;
      edge_index++;
    }
    last_src = middler.edges_[i].src_;
    last_dsts = middler.edges_[i].dst_;
  }
 
  graph.sortEdges();

  return graph;
}

CSRTemporalGraph TemporalGraph::convertToCSR()
{
  sortEdges();
  CSRTemporalGraph csr_temporal_graph(num_vertices_, num_edges_);

  VertexEdgeId cur_vertex = 0;
  for(VertexEdgeId i = 0; i < num_edges_; i++)
  {
    csr_temporal_graph.temporal_nbrs_[i].dst_ = temporal_edges_[i].dst_;
    csr_temporal_graph.temporal_nbrs_[i].time_ = temporal_edges_[i].time_;
    while(cur_vertex <= temporal_edges_[i].src_)
      csr_temporal_graph.offsets_[cur_vertex++] = i;
  }
  while(cur_vertex <= num_vertices_)
    csr_temporal_graph.offsets_[cur_vertex++] = num_edges_;

  return csr_temporal_graph;
}

void TemporalGraph::printGraph()
{
  for(VertexEdgeId i = 0; i < num_edges_; i++)
  {
    cout << temporal_edges_[i].src_ << " " << temporal_edges_[i].dst_ << " " << temporal_edges_[i].time_ << endl;
  }
}

void TemporalGraph::deleteGraph()
{
  delete[] temporal_edges_;
}


// CSRGraph
CSRGraph::CSRGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges): num_vertices_(num_vertices), num_edges_(num_edges)
{
  offsets_ = new VertexEdgeId[num_vertices_+1];
  nbrs_ = new VertexEdgeId[num_edges_];
}

VertexEdgeId CSRGraph::maxDegree()
{
  VertexEdgeId max_degree = 0;
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    max_degree = max(max_degree, offsets_[i+1] - offsets_[i]);
  }
  return max_degree;
}

void CSRGraph::sortNbrs()
{
  for (VertexEdgeId i=0; i < num_vertices_; i++)
        sort(nbrs_+offsets_[i], nbrs_+offsets_[i+1]);  
} 

struct DegreeComparator
{
    VertexEdgeId *degrees_;

    DegreeComparator(VertexEdgeId *degrees): degrees_(degrees) {}
    bool operator()(VertexEdgeId i1, VertexEdgeId i2)
    {
        return degrees_[i1] < degrees_[i2];
    }
};

void CSRGraph::findDegenOrdering()
{
  degen_order_ = new VertexEdgeId[num_vertices_];
  sort_by_degen_ = new VertexEdgeId[num_vertices_];

  VertexEdgeId *sort_by_deg = new VertexEdgeId[num_vertices_];
  VertexEdgeId *degrees = new VertexEdgeId[num_vertices_];
  VertexEdgeId *deg_order = new VertexEdgeId[num_vertices_];
  VertexEdgeId *deg_start_point = new VertexEdgeId[num_vertices_+1];
  bool *is_removed = new bool[num_vertices_];

  // sort_by_deg
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    sort_by_deg[i] = i;
    is_removed[i] = 0;
    degrees[i] = offsets_[i+1] - offsets_[i];
  }

  sort(sort_by_deg, sort_by_deg+num_vertices_, DegreeComparator(degrees));

  for (VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    deg_order[sort_by_deg[i]] = i;
  }

  // deg_start_point
  VertexEdgeId cur_vertex = sort_by_deg[0];
  VertexEdgeId min_degree = offsets_[cur_vertex+1] - offsets_[cur_vertex];
  for(VertexEdgeId i =0; i < num_vertices_; i++)
    deg_start_point[i] = -1;
  deg_start_point[min_degree] = 0;

  VertexEdgeId current_degree, running_degree;
  running_degree = min_degree;
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    cur_vertex = sort_by_deg[i];
    current_degree = offsets_[cur_vertex+1] - offsets_[cur_vertex];
    if (current_degree == running_degree)
      continue;
    deg_start_point[current_degree] = i;
    running_degree = current_degree;
  }


  VertexEdgeId to_remove;
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    to_remove = sort_by_deg[deg_start_point[min_degree]];
    degen_order_[to_remove] = i;
    is_removed[to_remove] = 1;

    if(i == num_vertices_ -1)
      continue;

    VertexEdgeId next = sort_by_deg[deg_start_point[min_degree]+1]; 
    sort_by_deg[deg_start_point[min_degree]] = -1;
    if(degrees[next] != degrees[to_remove])
    {
      deg_start_point[min_degree] = -1;
      min_degree = degrees[next];
    }
    else
    {
      deg_start_point[min_degree]++;
    }
    for (VertexEdgeId j = offsets_[to_remove]; j < offsets_[to_remove+1]; j++)
    {
      VertexEdgeId nbr = nbrs_[j];
      if (is_removed[nbr])
        continue;
      VertexEdgeId swap = sort_by_deg[deg_start_point[degrees[nbr]]];
      sort_by_deg[deg_start_point[degrees[nbr]]] = nbr;
      sort_by_deg[deg_order[nbr]] = swap;

      deg_order[swap] = deg_order[nbr];
      deg_order[nbr] = deg_start_point[degrees[nbr]];

      if (deg_start_point[degrees[nbr]]+1 < num_vertices_ &&  degrees[sort_by_deg[deg_start_point[degrees[nbr]]+1]] == degrees[nbr])
        deg_start_point[degrees[nbr]]++;
      else
        deg_start_point[degrees[nbr]] = -1;

      degrees[nbr]--;
      if(deg_start_point[degrees[nbr]] == -1)
        deg_start_point[degrees[nbr]] = deg_order[nbr];
      if(degrees[nbr] < min_degree)
        min_degree = degrees[nbr];
    }
  }

  for(VertexEdgeId i=0; i < num_vertices_; i++)
  {
    sort_by_degen_[degen_order_[i]] = i;
  }

  delete[] sort_by_deg;
  delete[]  degrees;
  delete[] deg_order;
  delete[] deg_start_point;
  delete[] is_removed;
}

void CSRGraph::relabelByDegenOrder()
{

  VertexEdgeId *new_offsets = new VertexEdgeId[num_vertices_+1];
  VertexEdgeId *new_nbrs = new VertexEdgeId[num_edges_];

  new_offsets[0] = 0;
  VertexEdgeId new_nbrs_index = 0;
  for (VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    VertexEdgeId inverse = sort_by_degen_[i];
    for(VertexEdgeId j = offsets_[inverse]; j < offsets_[inverse+1]; j++)
    {
      new_nbrs[new_nbrs_index] = degen_order_[nbrs_[j]];
      new_nbrs_index++;
    }
    new_offsets[i+1] = new_nbrs_index;
  }

  delete[] offsets_;
  delete[] nbrs_;

  offsets_ = new_offsets;
  nbrs_ = new_nbrs;

  sortNbrs();
}

VertexEdgeId CSRGraph::getEdgeIndx(VertexEdgeId u, VertexEdgeId v)
{
  VertexEdgeId low = offsets_[u];
  VertexEdgeId high = offsets_[u+1]-1;
  VertexEdgeId mid;

  while(low <= high)
  {
    mid = (low+high)/2;
    if (nbrs_[mid] == v)
      return mid;
    if (nbrs_[mid] > v)
      high = mid - 1;
    else
      low = mid + 1;
  }

  return -1;
}


void CSRGraph::printGraph()
{
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    for(VertexEdgeId j = offsets_[i]; j < offsets_[i+1]; j++)
      cout << i << " " << nbrs_[j] << endl;
  }
}

void CSRGraph::deleteGraph()
{
  delete[] offsets_;
  delete[] nbrs_;
}


// CSRTemporalGraph
CSRTemporalGraph::CSRTemporalGraph(VertexEdgeId num_vertices, VertexEdgeId num_edges): num_vertices_(num_vertices), num_edges_(num_edges)
{
  offsets_ = new VertexEdgeId[num_vertices_+1];
  temporal_nbrs_ = new CSRTemporalNbr[num_edges_];
}

VertexEdgeId CSRTemporalGraph::numberOfStaticDirectedEdges()
{
  VertexEdgeId num_static_dir_edges = 0;

  for (VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    VertexEdgeId running_nbr = -1;
    for (VertexEdgeId j = offsets_[i]; j < offsets_[i+1]; j++)
    {
      VertexEdgeId cur_nbr = temporal_nbrs_[j].dst_;
      if(cur_nbr == running_nbr)
        continue;
      num_static_dir_edges++;
      running_nbr = cur_nbr;
    }
  }

  return num_static_dir_edges;
}


void CSRTemporalGraph::sortNbrs()
{
  for (VertexEdgeId i=0; i < num_vertices_; i++)
        sort(temporal_nbrs_+offsets_[i], temporal_nbrs_+offsets_[i+1]);  
}

void CSRTemporalGraph::relabelByDegenOrder(VertexEdgeId *degen_order, VertexEdgeId *sort_by_degen)
{

  VertexEdgeId *new_offsets = new VertexEdgeId[num_vertices_+1];
  CSRTemporalNbr *new_temporal_nbrs = new CSRTemporalNbr[num_edges_];

  new_offsets[0] = 0;
  VertexEdgeId new_nbrs_index = 0;
  for (VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    VertexEdgeId inverse = sort_by_degen[i];
    for(VertexEdgeId j = offsets_[inverse]; j < offsets_[inverse+1]; j++)
    {
      new_temporal_nbrs[new_nbrs_index].dst_ = degen_order[temporal_nbrs_[j].dst_];
      new_temporal_nbrs[new_nbrs_index].time_ = temporal_nbrs_[j].time_;
      new_nbrs_index++;
    }
    new_offsets[i+1] = new_nbrs_index;
  }


  delete[] offsets_;
  delete[] temporal_nbrs_;

  offsets_ = new_offsets;
  temporal_nbrs_ = new_temporal_nbrs;

  sortNbrs();
}

CSRGraph CSRTemporalGraph::extractMultGraph()
{
  CSRGraph mult_graph(num_vertices_, numberOfStaticDirectedEdges());
  mult_graph.temporal_start_pos_ = new VertexEdgeId[num_edges_+1];


  VertexEdgeId edge_index = 0;
  mult_graph.offsets_[0] = 0;
  
  for (VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    VertexEdgeId running_nbr = -1;
    for(VertexEdgeId j = offsets_[i]; j < offsets_[i+1]; j++)
    {
      VertexEdgeId cur_nbr = temporal_nbrs_[j].dst_;
      if(cur_nbr == running_nbr)
      {
        continue;
      }
      mult_graph.nbrs_[edge_index] = cur_nbr;
      mult_graph.temporal_start_pos_[edge_index] = j;
      running_nbr = cur_nbr;
      edge_index++;
    }
    mult_graph.offsets_[i+1] = edge_index;
  }
  mult_graph.temporal_start_pos_[edge_index] = offsets_[num_vertices_];


  return mult_graph;
}

void CSRTemporalGraph::reloadTimes()
{
  times_ = new VertexEdgeId[num_edges_];
  for (VertexEdgeId i = 0; i < num_edges_; i++)
    times_[i] = temporal_nbrs_[i].time_;
}

VertexEdgeId CSRTemporalGraph::edgeTimeMinLimitSearch(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value)
{
  if(start_pos == -1)
    return -1;
  VertexEdgeId first_ind = lower_bound(times_ + start_pos, times_ + end_pos, min_value) - times_;
  if(first_ind == end_pos)
    return -1;
  return first_ind;
}

VertexEdgeId CSRTemporalGraph::edgeTimeMaxLimitSearch(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime max_value)
{
  if(start_pos == -1)
    return -1;
  VertexEdgeId last_ind = upper_bound(times_ + start_pos, times_ + end_pos, max_value) - times_;
  if (last_ind == start_pos)
    return -1;
  return last_ind; // one ahead of the actual last one
}

VertexEdgeId CSRTemporalGraph::edgeTimeMinLimitCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value)
{
  if(start_pos == -1)
    return 0;
  VertexEdgeId first_ind = lower_bound(times_ + start_pos, times_ + end_pos, min_value) - times_;
  if(first_ind == end_pos)
    return 0;
  return end_pos - first_ind;
}

VertexEdgeId CSRTemporalGraph::edgeTimeMaxLimitCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime max_value)
{
  if(start_pos == -1)
    return 0;
  VertexEdgeId last_ind = upper_bound(times_ + start_pos, times_ + end_pos, max_value) - times_;
  if (last_ind == start_pos)
    return 0;
  return last_ind - start_pos;
}

VertexEdgeId CSRTemporalGraph::edgeTimeIntervalCount(VertexEdgeId start_pos, VertexEdgeId end_pos, TemporalTime min_value, TemporalTime max_value)
{
  VertexEdgeId first_ind = edgeTimeMinLimitSearch(start_pos, end_pos, min_value);
  if(first_ind == -1)
    return 0;
  VertexEdgeId last_ind = edgeTimeMaxLimitSearch(start_pos, end_pos, max_value);
  if(last_ind == -1)
    return 0;

  return last_ind - first_ind;
}

void CSRTemporalGraph::printTimeSpan()
{
  TemporalTime min_time = times_[0];
  TemporalTime max_time = times_[0];
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    min_time = min(min_time, times_[i]);
    max_time = max(max_time, times_[i]);
  }
  cout << "Time span (days): " << (max_time - min_time) / (3600 * 24) << endl;
}

void CSRTemporalGraph::printGraph()
{
  for(VertexEdgeId i = 0; i < num_vertices_; i++)
  {
    for(VertexEdgeId j = offsets_[i]; j < offsets_[i+1]; j++)
      cout << i << " " << temporal_nbrs_[j].dst_ << " " << temporal_nbrs_[j].time_ << endl;
  }
}

void CSRTemporalGraph::deleteGraph()
{
  delete[] offsets_;
  delete[] temporal_nbrs_;
  delete[] times_;
}

CSRDAG::CSRDAG(CSRGraph &csr_graph)
{
  out_edge_dag_ = CSRGraph(csr_graph.num_vertices_, csr_graph.num_edges_/2);
  in_edge_dag_ = CSRGraph(csr_graph.num_vertices_, csr_graph.num_edges_/2);

  out_edge_dag_.offsets_[0] = 0;
  in_edge_dag_.offsets_[0] = 0;

  VertexEdgeId out_nbr_index = 0, in_nbr_index = 0;

  // VertexEdgeId max_out_deg = 0;

  for (VertexEdgeId i=0; i < csr_graph.num_vertices_; i++)
  {
    for (VertexEdgeId j = csr_graph.offsets_[i]; j < csr_graph.offsets_[i+1]; j++)
    {
      VertexEdgeId cur_nbr = csr_graph.nbrs_[j];

      if (i < cur_nbr)
      {
        // cout << "out" << endl;
        out_edge_dag_.nbrs_[out_nbr_index] = cur_nbr;
        out_nbr_index++;
      }
      else
      {
        // cout << "in" << endl;
        in_edge_dag_.nbrs_[in_nbr_index] = cur_nbr;
        in_nbr_index++;
      }
    }

    out_edge_dag_.offsets_[i+1] = out_nbr_index;
    in_edge_dag_.offsets_[i+1] = in_nbr_index;
  }
  out_edge_dag_.sortNbrs();
  in_edge_dag_.sortNbrs();
}

TemporalGraph loadTemporalGraph(const char *path)
{
  FILE* f = fopen(path, "r");

  TemporalGraph temporal_graph;
  temporal_graph.num_vertices_ = 0;
  temporal_graph.num_edges_ = 0;

  VertexEdgeId input_1, input_2;
  TemporalTime input_3;
  
  char line[1024];
  VertexEdgeId input_index = 0;
  while (fgets(line, sizeof(line), f))
  {
    if(temporal_graph.num_vertices_ == 0)
    {
      sscanf(line, "%ld%ld", &input_1, &input_2);
      temporal_graph.num_vertices_ = input_1;
      temporal_graph.num_edges_ = input_2;
      temporal_graph.temporal_edges_ = new TemporalEdge[temporal_graph.num_edges_];
    }
    else
    {
      sscanf(line, "%ld%ld%ld", &input_1, &input_2, &input_3);
      temporal_graph.temporal_edges_[input_index].src_ = input_1;
      temporal_graph.temporal_edges_[input_index].dst_ = input_2;
      temporal_graph.temporal_edges_[input_index].time_ = input_3;
      input_index++;
    }
  }

  fclose(f);
  return temporal_graph;
}


