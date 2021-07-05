#include "ETTC/Graph.h"
#include "ETTC/t_triangle_counting.h"
// #include <time.h>

using namespace std;

int main(int argc, char *argv[])
{
    TemporalTime delta = stoi(argv[2]);  // delta_{1,3} in the paper
    TemporalTime delta1 = stoi(argv[3]); // delta_{1,2} in the paper
    TemporalTime delta2 = stoi(argv[4]); // delta_{2,3} in the paper
    
    if(delta > delta1 + delta2)  // difference between first and third edge cannot be more than the sum
        delta = delta1 + delta2; // of difference between first and second edges and the difference between second and third edges

    TemporalGraph temporal_graph = loadTemporalGraph(argv[1]);  // read the input temporal graph

    cout << "Temporal graph loaded." << endl;
    
    Graph static_graph = temporal_graph.ExtractStaticGraph();   // extract the static graph

    cout << endl;
    cout << "Static graph extracted." << endl;

    CSRGraph static_csr_graph = static_graph.convertToCSR();   // convert the static graph into CSR format
    cout << endl;
    cout << "Static graph converted into CSR format." << endl;

    CSRTemporalGraph csr_temporal_graph = temporal_graph.convertToCSR();  // convert the input temporal graph into CSR format
    cout << endl;
    cout << "Temporal graph converted into CSR format." << endl;

    cout << endl;
    cout << "degeneracy ordering of vertices." << endl;
    static_csr_graph.findDegenOrdering();  // find the degeneracy ordering of the static graph (CSR format)

    cout << endl;
    cout << "Relabel vertices of the CSR temporal graph by degeneracy order of the static graph" << endl;
    csr_temporal_graph.relabelByDegenOrder(static_csr_graph.degen_order_, static_csr_graph.sort_by_degen_);  // relabel vertices of the temporal graph (CSR format) by the degeneracy ordering of the static graph

    cout << endl;
    cout << "Relabel vertices of the static CSR graph by degeneracy order" << endl;
    static_csr_graph.relabelByDegenOrder();   // relabel vertices of the static graph (CSR format) by the degenracy ordering of the static graph

    cout << endl;
    cout << "CSR static graph converted to CSRDAG" << endl;
    CSRDAG csr_dag(static_csr_graph);  // convert static graph (CSR format) to a DAG. In other words, we orient the static graph w.r.t. its degeneracy ordering.


    cout << endl;
    cout << "Degeneracy: " << csr_dag.out_edge_dag_.maxDegree() << endl; // max out-degree of the DAG oriented w.r.t. the degeneracy ordering
                                                                        //  is equal to the degeneracy of the graph

    // Counting temporal triangles for delta, delta_1, and delta_2
    MotifCounter motif_counter;
    motif_counter.countTemporalTriangle(csr_dag.out_edge_dag_, csr_temporal_graph, delta, delta1, delta2);

    cout << "\n///////////////// \nInfo\n/////////////////" << endl << endl;
    cout << "Number of vertices: " << csr_temporal_graph.num_vertices_ << endl;
    cout << "Number of edges: " << csr_temporal_graph.num_edges_ << endl;
    cout << "Number of static edges: " << csr_dag.out_edge_dag_.num_edges_ << endl;
    cout << "Number of static triangles: " << motif_counter.static_triangles_count_ << endl;
    cout << "Max multiplicity in one direction on a pair: " << motif_counter.highest_mult_ << endl;
    cout << "Avg multiplicity on a pair: " << csr_temporal_graph.num_edges_ / csr_dag.out_edge_dag_.num_edges_ << endl;
    csr_temporal_graph.printTimeSpan();

    // motif_counter.printCounts();
    motif_counter.printCountsFile(argv[5]);

    // free memory
    static_graph.deleteGraph();
    // cout << "a1" << endl;
    temporal_graph.deleteGraph();
    // cout << "a2" << endl;
    static_csr_graph.deleteGraph();
    // cout << "a3" << endl;
    csr_temporal_graph.deleteGraph();
    // cout << "a4" << endl;

    motif_counter.freeMemory();
    // cout << "a5" << endl;
    return 0;
 }

 