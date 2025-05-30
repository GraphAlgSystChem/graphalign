#ifndef GENERATECANDIDATES_H
#define GENERATECANDIDATES_H

namespace GraphAlign {

/**
 * @brief Generate all candidate pairs based on considering the next vertex in g2 from the total order and only consider the 
 *        the neighbours of the match (in g1) of parent (from the parent map) of that vertex.
 *        Called by pseudoVF, the subroutine of trim.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type 
 * @param g1 The first graph
 * @param g2 FilteredGraph of the second graph
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @param order List of the vertices in g2 indicating a total order 
 * @param parent Map indicating the parent of any vertex in g2 according to the total order
 * @return List of pairs of candidate matches
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::vector<std::pair<Vertex, Vertex>>
    generateCandidatesTrim(const Graph &g1, const FilteredGraph<Graph> &g2, 
                        GlobalVariables<Vertex> &globalVars,
                        const std::vector<Vertex> &order,
                        const std::vector<Vertex> &parent);

/**
 * @brief Generates the candidates for matching vertices between two graphs, but instead of saving them to a vector and then returning
 *        them the callback is used to evaluate the candidate pair right away. 
 *        Called by expand. 
 *        
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type
 * @tparam Callback Type of the callback
 * @param g1 The first graph, the larger graph
 * @param g2 The second graph, the smaller graph
 * @param g1Match Vector representing the matches for g1
 * @param g2Match Vector representing the matches for g2
 * @param parent Parent of next vertex
 * @param nextVertex The vertex to be matched in this step
 * @param callback The callback function
 * @return true The program should continue, as no reason has been found to terminate
 * @return false If either isMaxColoredMatch or isSubIso returned true, should return false to indicate the program should stop 
 */
template<typename Graph, typename Vertex, typename Callback>
bool generateCandidatesRec(const Graph &g1, const Graph &g2, 
                        const std::vector<Vertex> &g1Match,
                        const std::vector<Vertex> &g2Match, 
                        const std::vector<Vertex> &parent,
                        Vertex nextVertex,
                        Callback callback);


////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////

template<typename Graph, typename Vertex>
std::vector<std::pair<Vertex, Vertex>>
    generateCandidatesTrim(const Graph &g1, const FilteredGraph<Graph> &g2, 
                        GlobalVariables<Vertex> &globalVars,
                        const std::vector<Vertex> &order,
                        const std::vector<Vertex> &parent) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    std::vector<std::pair<Vertex, Vertex>> P; // The list of candidate pairs to be returned

    // mathcSize vertices have been mapped (potentially 0). With 0-indexing, the matchSize'th vertex is the next in the order.
    Vertex nextVertex = order[globalVars.matchSize];
    Vertex parentVertex = parent[getIndex(nextVertex)];

    if(parentVertex == nullVertex) {
        // Handle vertex in a new connected component. This vertex can only be mapped to vertices
        // outside the neighbourhood of already matched vertices in the large graph.
        std::vector<int> occupied(num_vertices(g1), 0);
        for(Vertex u : asRange(vertices(g1))) {
            if(globalVars.getLargerMatch()[getIndex(u)] != nullVertex) {
                occupied[getIndex(u)] = 1;
                for(Vertex neighbour : asRange(adjacent_vertices(u, g1))) {
                    occupied[getIndex(neighbour)] = 1;
                }
            }
        }

        // All vertices in g1 that are not occupied (either already matched or adjacent to matched vertices) are potential candidates
        for(Vertex u : asRange(vertices(g1))) {
            if(!occupied[getIndex(u)]) {
                P.push_back({u, nextVertex});
            }
        }
    } 
    else {
        // Otherwise, we are still inside the same connected component

        // Find the vertex in g1 that is matched with parentVertex
        Vertex parentTarget = globalVars.getSmallerMatch()[getIndex(parentVertex)];
        

        // If there are fewer neighbour candidates in g1 than in g2 it means it is not possible to match
        // all the neighbour vertices of the parent in g2 to the neighbours of the parentTarget in g1.
        // Thus, we can just move on to the next recursive call in pseudoVF.
        const auto &adjG1 = adjacent_vertices(parentTarget, g1);
        const auto &adjG2 = adjacent_vertices(parentVertex, g2);
        if(std::distance(adjG2.first, adjG2.second) > std::distance(adjG1.first, adjG1.second)) {
            return  P;
        }

        // nextVertex can only be mapped to the neighbourhood of parent's matched vertex.
        // These must be unmapped.    
        for(Vertex u : asRange(adjacent_vertices(parentTarget, g1))) {
            if(globalVars.getLargerMatch()[getIndex(u)] == nullVertex) {
                P.push_back({u, nextVertex});
            }
        }
    }

    return P;
} //generateCandidatesTrim

template<typename Graph, typename Vertex, typename Callback>
bool generateCandidatesRec(const Graph &g1, const Graph &g2, 
                        const std::vector<Vertex> &g1Match,
                        const std::vector<Vertex> &g2Match, 
                        const std::vector<Vertex> &parent,
                        Vertex nextVertex,
                        Callback callback) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    bool computeOccupied = false;
    bool foundFull; // the callback returns foundIsoOrMaxColor. The recursive expansion algorithm should not continue in that case (nor should the candidates)

    Vertex parentVertex = parent[getIndex(nextVertex)]; // g2 vertex

    if(parentVertex != nullVertex){
        // Otherwise, we are still inside the same connected component
        // Need to determine whether the parent of the nextVertex was matched or not
        
        // parentVertex is non-null
        Vertex parentTarget = g2Match[getIndex(parentVertex)]; // g1 vertex
        // if the parentVertex is a part of the current match, we limit the nextVertex
        // to be mapped to the nhbourhood of parentTarget
        if(parentTarget != nullVertex) { 

            // nextVertex can only be mapped to the neighbourhood of parent's matched vertex.
            // These must be unmapped.
            for(Vertex u : asRange(adjacent_vertices(parentTarget, g1))) {
                if(g1Match[getIndex(u)] == nullVertex) {
                    foundFull = callback({u, nextVertex});
                    if(foundFull) {
                        return false;
                    }
                }
            }
        }
        else { // parentVertex has not been mapped
            bool matchedNeighbours = false;
            for(Vertex neigh : asRange(adjacent_vertices(nextVertex, g2))) {
                if(g2Match[getIndex(neigh)] != nullVertex) {
                    matchedNeighbours = true;
                    break;
                }   
            }

            // The parentVertex is not matched, but nextVertex has other vertices that are already matched
            // Find "substitute parent"
            if(matchedNeighbours) {
                std::vector<unsigned char> used(num_vertices(g1), 0);
                for(Vertex neigh : asRange(adjacent_vertices(nextVertex, g2))) {
                    Vertex neighG1 = g2Match[getIndex(neigh)];
                    // if the current neighbour is a matched neighbour add all its unmatched neighbours as candidates to nextVertex
                    if(neighG1 != nullVertex) {
                        for(Vertex adj : asRange(adjacent_vertices(neighG1, g1))) {
                            if(g1Match[getIndex(adj)] == nullVertex && !used[getIndex(adj)]) {
                                foundFull = callback({adj, nextVertex});
                                if(foundFull) {
                                    return false;
                                }
                                used[getIndex(adj)] = 1;
                            }
                        }
                    }
                }
            }
            // nextVertex has no neighbours that are matched, so we are in the same case if our parentVertex did not exist (new connected component)
            else {
                computeOccupied = true;
            }
        }
    }

    if(computeOccupied || parentVertex == nullVertex) {
        // Handle vertex in a new connected component. This vertex can only be mapped to vertices
        // outside the neighbourhood of already matched vertices in the large graph.
        std::vector<unsigned char> occupied(num_vertices(g1), 0);

        // Find all the vertices in g1 thar are matched vertices or their neighbours
        for(Vertex x : asRange(vertices(g1))) {
            if(g1Match[getIndex(x)] != nullVertex) {
                occupied[getIndex(x)] = 1;
                for(Vertex neigh : asRange(adjacent_vertices(x, g1))) {
                    occupied[getIndex(neigh)] = 1;
                }
            }
        } 

        // All vertices in g1 that are not occupied (either already matched or adjacent to matched vertices)
        for(Vertex u : asRange(vertices(g1))) {
            if(!occupied[getIndex(u)]) {
                foundFull = callback({u, nextVertex});
                if(foundFull) {
                    return false;
                }
            }
        }
    } 
    
    return true; // When all candidates have been considered, give recursiveExpand the opportunity to map u to nothing.

} // generateCandidatesRec


} // namespace GraphAlign

#endif // GENERATECANDIDATES_H