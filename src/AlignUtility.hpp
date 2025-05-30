#ifndef ALIGNUTILITY_H
#define ALIGNUTILITY_H


// std
#include <queue>
#include <limits>
#include <functional>

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

// Our code
#include "AlignParameters.hpp"
#include "PairToRangeAdaptor.hpp"
#include "FilteredGraph.hpp"
#include "GlobalVariables.hpp"

namespace GraphAlign {

// This file contains the following functions:
// getAnchorPreMap
// printResults
// printMatchTable
// getIndex
// printGraph


////////////////////////////////// DECLARATIONS //////////////////////////////////


/**
 * @brief Get a vector of vertex pairs corresponding to a premap of vertices between the two
 * underlying graph objects of a1 and a2.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type
 * @param a1 The first alignment graph
 * @param a2 The second alignment graph
 * @param anchor The anchor of vertices amongst the leaves (input graphs). anchor[i] is an ordered list of vertices mapped to each other from all the input graphs. 
 * anchor[i][j] is the i'th anchored vertex from graph j.
 * @return A list of (u, v) where u in G(a1) and v in G(a2) such that u -> v is a map
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::vector<std::pair<Vertex, Vertex>> getAnchorPreMap(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, const std::vector<std::vector<Vertex>> &anchor); 


template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void printResults(const std::vector<Graph> &inputGraphs, const Graph &alignmentGraph, const std::unordered_map<int, std::unordered_map<Vertex, Vertex>> &projection);

/**
 * Prints the match table of an alignment object in human readable manner.
 * 
 * @param a The alignment object containing the match table to be printed.
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void printMatchTable(const AlignObj<Graph> &a);

/**
 * @brief Get the id for the given vertex. Necessary when indexing into vectors with Vertex.
 * NOTE: Currently not very useful as we use vertex_descriptor for adjacency_list with vecS for vertices, i.e. vertex is an usigned long int.
 * The idea is that if we change this later, then updating this function and using it in the rest of the code should be a good idea.
 * 
 * @tparam Vertex Vertex type
 * @param v The vertex to get an id for
 * @return int Index for vertex v
 */
template<typename Vertex>
int getIndex(Vertex v);

/**
 * @brief Prints the content of the graph. 
 * 
 */

template <typename Graph>
void printGraph(const Graph &g);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template<typename Graph, typename Vertex>
std::vector<std::pair<Vertex, Vertex>> getAnchorPreMap(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, const std::vector<std::vector<Vertex>> &anchor) {
    std::vector<std::pair<Vertex, Vertex>> anchorPreMap = {};
    int fromA1Index = a1.contained[0]; // The index of some graph in a1, G_i
    int fromA2Index = a2.contained[0]; // The index of some graph in a2, G_j

    for(const std::vector<Vertex> &preMap : anchor) {
        Vertex u = preMap[fromA1Index]; // The anchored vertex in G_i
        Vertex v = preMap[fromA2Index]; // The anchored vertex in G_j
        // We now have u -> v mapping between the two leaf graphs from a1 and a2's subtrees respectively. 
        // Transform these into x -> y in a1 and a2 using matchTables.
        Vertex x, y;
        if(a1.contained.size() == 1) { // a1 is a leaf graph, so the extracted vertex u is already a vertex of a1.
            x = u;
        } else {
            for(Vertex m : asRange(vertices(a1.alignGraph))) {
                if(a1.matchTable[0][m] == u) {
                    x = m;
                    break;
                }
            }
        }
        
        if(a2.contained.size() == 1){
            y = v;
        } else {
            for(Vertex n : asRange(vertices(a2.alignGraph))) {
                if(a2.matchTable[0][n] == v) {
                    y = n;
                    break;
                }
            }
        }
        
        std::pair<Vertex, Vertex> foundMap = {x, y};
        anchorPreMap.push_back(foundMap);
    }
    return anchorPreMap;

} // getAnchorPreMap


template<typename Graph, typename Vertex>
void printResults(const std::vector<Graph> &inputGraphs, const Graph &alignmentGraph, const std::unordered_map<int, std::unordered_map<Vertex, Vertex>> &projection) {

    std::cout << std::endl << std::endl << "PRINTING ALIGNMENT" << std::endl;
    for(size_t i = 0; i < inputGraphs.size(); i++){
        Graph Gi = inputGraphs[i];
        std::cout << "INPUT GRAPH: " << i << std::endl;
        for(Vertex v : asRange(vertices(Gi))){
            Vertex x = projection.find(i)->second.find(v)->second;
            std::cout << "Vertex " << v << "(" << inputGraphs[i][v].label << ")" << " -> " << x << "(" << alignmentGraph[x].label << ")" << std::endl;
            // std::cout << "With labels " << inputGraphs[i][v].label << " and " << alignmentGraph[x].label << "\n";
            // std::cout << "and adjacency count " << out_degree(v, alignmentGraph) << " and " << out_degree(x, alignmentGraph) << "\n";
        }   
    }

    std::cout << "ALIGNMENT GRAPH:" << std::endl; 
    for(auto e : asRange(edges(alignmentGraph))){
        std::cout << source(e, alignmentGraph) << " - " << target(e, alignmentGraph) << std::endl;
    }

} // printResults

template <typename Graph, typename Vertex>
void printMatchTable(const AlignObj<Graph> &a){
    std::cout << "uV:\t";
    for(Vertex u : asRange(vertices(a.alignGraph))){
        std::cout << u << " ";
    }
    std::cout << "\n";
    for(size_t i = 0; i < a.contained.size(); i++){
        int graphIndex = a.contained[i];
        std::cout << "G_" << graphIndex << "\t";
        for(Vertex u : asRange(vertices(a.alignGraph))){
            Vertex domain = a.matchTable[i][u];
            if( domain == boost::graph_traits<Graph>::null_vertex()){
                std::cout << "NIL ";
            } else {
                std::cout << domain << " ";
            }
        }
        std::cout << "\n";
    }
} // printMatchTable


template<typename Vertex>
int getIndex(Vertex v) {
    return v;
} // getIndex

template <typename Graph>
void printGraph(const Graph &g){
    std::cout << "VERTICES AND NEIGHBOURS:\n";
    for(auto u : asRange(vertices(g))){
        std::cout << "Vertex: " << u << "\n";
        for(auto v : asRange(adjacent_vertices(u, g))){
            std::cout << "  -> " << v << "\n";
        }
    }
}

} // namespace GraphAlign

#endif // ALIGNUTILITY_H