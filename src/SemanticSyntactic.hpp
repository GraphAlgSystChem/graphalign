#ifndef SEMANTICSYNTACTIC_H
#define SEMANTICSYNTACTIC_H


namespace GraphAlign {

// This file contains the relevant code for feasibility checks, both semantic and syntactic
// semanticCheck
// syntacticCheckExp
// syntacticCheckTrim

/**
 * @brief Checks the semantic constraints, which means that for a given vertex mapping (u,v) we need to ensure their labels match 
 * is allowed and the same for matched edges if we consider edge labels. 
 * Whether a label match is allowed depends on whether label mismatch is allowed by the user.
 *
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param g1 The first graph
 * @param g2 The second graph
 * @param g1Match The match vector for g1, namely g1Match[i] = j if v_i in G1 is mapped to v_j in G2.
 * @param g2Match The match vector for g2, namely g2Match[i] = j if v_i in G2 is mapped to v_j in G1.
 * @param match The current match that is being considered
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return true If the check is satisfied, namely that vertex and edge labels have to adhere to each other
 * @return false If the check is not satisfied
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool semanticCheck(const Graph &g1, const Graph &g2,
                const std::vector<Vertex> &g1Match,
                const std::vector<Vertex> &g2Match,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList);

/**
 * @brief Checks the semantic constraints in the same manner as above. Called by the trim algorithm.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Subgraph 
 * @tparam Vertex The vertex type 
 * @param g1 The first graph
 * @param g2 The second graph
 * @param match The current match that is being considered
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @return true If the check is satisfied, namely that vertex and edge labels have to adhere to each other
 * @return false If the check is not satisfied
 */
template<typename Graph, class SV, typename Subgraph, typename Vertex>
bool semanticCheck(const Graph &g1, const Subgraph &g2,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList, 
                GlobalVariables<Vertex> &globalVars);

/**
 * @brief Checks the syntactic constraints, i.e. for any two matches (u, v) and (u', v') that they need to agree on adjacency.
 * This means that if (u, u') is adjacent in g1 then (v,v') needs to be adjacent in g2 and vice versa.
 * By allowing the usage of ambiguous edges to obtain adjacency this becomes slightly more complicated 
 * by allowing adjacency constraints to be satisfied by these edges.
 * This function is called in the recursive expansion.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @tparam Vertex The vertex type
 * @param a1 First alignment object
 * @param a2 Second alignment object
 * @param a1Match The match vector for a1's underlying graph, g1, namely a1Match[i] = j if v_i in G1 is mapped to v_j in G2.
 * @param a2Match The match vector for a2's underlying graph, g2, namely a2Match[i] = j if v_i in G2 is mapped to v_j in G1.
 * @param match The potential match we are trying to check.
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return true If the check is satisfied
 * @return false If the check is not satisfied
 */
template<typename Graph, class SV, typename MatchTable, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool syntacticCheckExp(const AlignObj<Graph, MatchTable> &a1, const AlignObj<Graph, MatchTable> &a2,
                const std::vector<Vertex> &a1Match,
                const std::vector<Vertex> &a2Match,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList);

/**
 * @brief Checks the syntactic constraints, i.e. for any two matches (u, v) and (u', v') that they need to agree on adjacency.
 * This means that if (u, u') is adjacent in g1 then (v,v') needs to be adjacent and vice versa.
 * By allowing the usage of ambiguous edges to obtain adjacency this becomes slightly more complicated.
 * This function is called in the iterative trimming. Matches are now represented by two vectors stored in globalVars. 
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param g1 First input graph
 * @param g2 Second input graph
 * @param match The potential match we are trying to check (u, v)
 * @param ambG1 A list of pairs of vertices in g1, which becomes adjacenct through ambiguous edges
 * @param ambG2 A list of pairs of vertices in g2, which becomes adjacenct through ambiguous edges
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @return true If the check is satisfied
 * @return false If the check is not satisfied
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool syntacticCheckTrim(const Graph &g1, const FilteredGraph<Graph> &g2,
                const std::pair<Vertex,Vertex> &match, 
                const std::vector<std::pair<Vertex, Vertex>> &ambG1, 
                const std::vector<std::pair<Vertex, Vertex>> &ambG2,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////

template<typename Graph, class SV, typename Vertex>
bool semanticCheck(const Graph &g1, const Graph &g2,
                const std::vector<Vertex> &g1Match,
                const std::vector<Vertex> &g2Match,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList) {
    
    const Vertex u = match.first; // Vertex from g1
    const Vertex v = match.second; // Vertex from g2

    if(paramList.enforceVLabels) {
        if(g1[u].label != g2[v].label) {
            return false;
        }
    }
    if(paramList.enforceELabels) {
        // Compare loop-label (if any)
        auto edgeU = edge(u, u, g1);
        auto edgeV = edge(v, v, g2);
        if(edgeU.second && edgeV.second) {
            if(g1[edgeU.first].label != g2[edgeV.first].label) {
                return false;
            }
        }
        // Only consider edge labels between u - neighbour and neighbourMatch - v. If edges do not match up, return false, even though that is technically syntactic check.
        for(Vertex x : asRange(adjacent_vertices(u, g1))) {
            if(g1Match[getIndex(x)] != boost::graph_traits<Graph>::null_vertex()) {
                auto edgeV = edge(g1Match[getIndex(x)], v, g2);
                if(edgeV.second) {
                    if(g1[edge(u, x, g1).first].label != g2[edgeV.first].label) {
                        return false;
                    }
                }
                // If the corresponding edge does not exist in G2, it means syntactic check fails and no reason to spend time checking in that function as well
                else {
                    return false;
                }
            } 
        }
    }
    return true;

} // semanticCheck

template<typename Graph, class SV, typename Subgraph, typename Vertex>
bool semanticCheck(const Graph &g1, const Subgraph &g2,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList, 
                GlobalVariables<Vertex> &globalVars) {
    
    std::vector<Vertex> &largerMatch = globalVars.getLargerMatch();
    const Vertex u = match.first; // Vertex from g1
    const Vertex v = match.second; // Vertex from g2

    if(paramList.enforceVLabels) {
        if(g1[u].label != g2[v].label) {
            return false;
        }
    }
    if(paramList.enforceELabels) {
        // Compare loop-label (if any)
        auto edgeU = edge(u, u, g1);
        auto edgeV = edge(v, v, g2);
        if(edgeU.second && edgeV.second) {
            if(g1[edgeU.first].label != g2[edgeV.first].label) {
                return false;
            }
        }
        // Only consider edge labels between u - neighbour and neighbourMatch - v. If edges do not match up, return false, even though that is technically syntactic check.
        for(Vertex x : asRange(adjacent_vertices(u, g1))) {
            if(largerMatch[getIndex(x)] != boost::graph_traits<Graph>::null_vertex()) {
                auto edgeVal = edge(v, largerMatch[getIndex(x)], g2);
                if(edgeVal.second) {
                    if(g1[edge(u, x, g1).first].label != g2[edgeVal.first].label) {
                        return false;
                    }
                }
                // If the corresponding edge does not exist in G2, it means syntactic check fails and no reason to spend time checking in that function as well
                else {
                    return false;
                }
            } 
        }
    }
    return true;

} // semanticCheck


template<typename Graph, class SV, typename MatchTable, typename Vertex>
bool syntacticCheckExp(const AlignObj<Graph, MatchTable> &a1, const AlignObj<Graph, MatchTable> &a2,
                const std::vector<Vertex> &a1Match,
                const std::vector<Vertex> &a2Match,
                const std::pair<Vertex,Vertex> &match, 
                const AlignParameters<Graph, SV> &paramList) {
    
    const Vertex u = match.first; // The vertex from g1 in the match
    const Vertex v = match.second; // The vertex from g2 in the match

    // Loop consitency test
    if(edge(u, u, a1.alignGraph).second != edge(v, v, a2.alignGraph).second) return false;

    std::vector<unsigned char> ambN1(num_vertices(a1.alignGraph), 0);
    std::vector<unsigned char> ambN2(num_vertices(a2.alignGraph), 0);

    if(paramList.useAmbiguousEdges) {
        for(const std::pair<Vertex, Vertex> &pair : a1.ambiguousEdges) {
            if(u == pair.first) { 
                Vertex ambNId = getIndex(pair.second);
                ambN1[ambNId] = 1;
            }
            else if(u == pair.second) {
                Vertex ambNId = getIndex(pair.first);
                ambN1[ambNId] = 1;
            }
        }
        for(const std::pair<Vertex, Vertex> &pair : a2.ambiguousEdges) {
            if(v == pair.first) { 
                Vertex ambNId = getIndex(pair.second);
                ambN2[ambNId] = 1;
            }
            else if(v == pair.second) {
                Vertex ambNId = getIndex(pair.first);
                ambN2[ambNId] = 1;
            }
        }
    }
    // For all the neighbours of u in g1, check if those neighbours are in the match, 
    // if so check that the corresponding match is also a neighbour of v OR an ambiguous neighbour of v
    for(Vertex x: asRange(adjacent_vertices(u, a1.alignGraph))) {
        if(a1Match[getIndex(x)] != boost::graph_traits<Graph>::null_vertex()) {
            if(!edge(a1Match[getIndex(x)], v, a2.alignGraph).second && !ambN2[getIndex(a1Match[getIndex(x)])]) {
                return false;
            }
        }
    }
    for(Vertex x: asRange(adjacent_vertices(v, a2.alignGraph))) {
        if(a2Match[getIndex(x)] != boost::graph_traits<Graph>::null_vertex()) {
            if(!edge(a2Match[getIndex(x)], u, a1.alignGraph).second && !ambN1[getIndex(a2Match[getIndex(x)])]) {
                return false;
            }
        }
    }

    return true;

} // syntacticCheckRec

template<typename Graph, class SV, typename Vertex>
bool syntacticCheckTrim(const Graph &g1, const FilteredGraph<Graph> &g2,
                const std::pair<Vertex,Vertex> &match, 
                const std::vector<std::pair<Vertex, Vertex>> &ambG1, 
                const std::vector<std::pair<Vertex, Vertex>> &ambG2,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars) {
    
    const Vertex u = match.first;
    const Vertex v = match.second;
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();

    // Loop cconsistency, i.e. either both have or do not have loops.
    auto edgeU = edge(u, u, g1);
    auto edgeV = edge(v, v, g2);
    if(edgeU.second && !edgeV.second) {
        return false;
    }
    if(!edgeU.second && edgeV.second) {
        return false;
    }


    int candidateNeighG1 = 0;
    int candidateNeighG2 = 0;

    int exteriorNeighG1 = 0;
    int exteriorNeighG2 = 0;

    // Keeping track of the vertices that are neighbours to u and v through ambiguous edges
    std::vector<unsigned char> ambNeighG1(num_vertices(g1), 0);
    std::vector<unsigned char> ambNeighG2(num_vertices(g2), 0);

    if(paramList.useAmbiguousEdges) {
        // Computing set of ambiguous neighbours amongst u and v
        for(const std::pair<Vertex, Vertex> &e : ambG1) {
            Vertex x = e.first;
            Vertex y = e.second;
            Vertex z = u; // initial assignment to be able to check against
            if(u == x) {
                ambNeighG1[y] = 1;
                z = y;
            } else if (u == y) {
                ambNeighG1[x] = 1;
                z = x;
            }
            // Here, we distinguish between ambiguous neighbours that are exterior or candidates. 
            // Ambiguous neigbours that _have_ been matched cannot count towards the useful neighbourhood of g2 as these
            // are occupied.
            if(z != u && globalVars.getLargerMatch()[getIndex(z)] == nullVertex) {
                bool candidateNeighbour = false;
                // zNeigh is the neighbour of u's neighbour, z
                for(Vertex zNeigh: asRange(adjacent_vertices(z, g1))) {
                    if(globalVars.getLargerMatch()[zNeigh] != nullVertex) { // zNeigh is a matched vertex
                        candidateNeighbour = true;
                        break;
                    }
                }
                // if z is not matched, but neighbour to a matched vertex then it is a candidate
                if(candidateNeighbour) {
                    candidateNeighG1++;
                }
                // z is not matched and it has no matched neighbours
                else {
                    exteriorNeighG1++;
                }                
            }

        }
        for(const std::pair<Vertex, Vertex> &e : ambG2) {
            Vertex x = e.first;
            Vertex y = e.second;
            if(v == x) {
                ambNeighG2[y] = 1;
            } else if (v == y) {
                ambNeighG2[x] = 1;
            }
        }
    }

    // Check for all _matched_ neighbours x of u; x's image must necessarily be connected to v by an (ambiguous) edge. If not, u and v are incompatible.
    for(Vertex neigh: asRange(adjacent_vertices(u, g1))) {
        Vertex neighMatch = globalVars.getLargerMatch()[getIndex(neigh)];
        if(neighMatch != nullVertex) { // Neighbour has been matched.
            // Check if the matched vertex of neigh is adjacent to v, either through actual edge or through an ambiguous edge
            if(!edge(neighMatch, v, g2).second && !ambNeighG2[getIndex(neighMatch)]) return false;
        }
        else { // Neighbour has not been matched.
            bool candidateNeighbour = false;
            // nextNeigh is the neighbour of u's neighbour
            for(Vertex nextNeigh: asRange(adjacent_vertices(neigh, g1))) {
                if(globalVars.getLargerMatch()[nextNeigh] != nullVertex) {
                    candidateNeighbour = true;
                    break;
                }
            }
            // If neigh is not matched, but neighbour to a matched vertex then it is a candidate
            if(candidateNeighbour) {
                candidateNeighG1++;
            }
            // The neighbour is not matched and it has no matched neighbours
            else {
                exteriorNeighG1++;
            }
        }
    }
    // Symmetric, the other way around.
    for(Vertex neigh: asRange(adjacent_vertices(v, g2))) {
        Vertex neighMatch = globalVars.getSmallerMatch()[getIndex(neigh)];
        if(neighMatch != nullVertex) {
            if(!edge(neighMatch, u, g1).second && !ambNeighG1[getIndex(neighMatch)]) return false;
        }
        else { // neighMatch is null vertex
            bool candidateNeighbour = false;
            // nextNeigh is the neighbour of u's neighbour
            for(Vertex nextNeigh: asRange(adjacent_vertices(neigh, g2))) {
                if(globalVars.getSmallerMatch()[nextNeigh] != nullVertex) {
                    candidateNeighbour = true;
                    break;
                }
            }
            if(candidateNeighbour) {
                candidateNeighG2++;
            }
            else {
                exteriorNeighG2++;
            }
        }
    }

    // Look-ahead on only the neighbourhoods of u and v, considering the unmatched neighbours and dividing them into two categories
    // Either a neighbour is a candidate neighbour, which is when the neighbour of u or v is also a neighbour to a matched vertex
    // An exterior neighbour is a neighbour of u and v for which no of its other neighbours are matched
    // For both of these categories it is necessary that the number is greater or equal in g1 (larger graph) than in g2 (smaller graph)
    // It is a limited look-ahead as it does not consider the entire graphs, however, it requires less work to only consider two sets of neighbours
    // instead of all vertices in each graph
    if(candidateNeighG2 > candidateNeighG1) {
        // This would imply that the neighbourhood of v in the smaller graph is too large to be mapped to the neighbourhood of u in the larger graph.
        // This would leave neighbours dangling, and this mapping can therefore not result in a subgraph isomorphism.
        return false;
    }
    if(exteriorNeighG2 > exteriorNeighG1) {
        // Same argument as above. There should be enough completely blank u-neighbouring vertices in G1 to cover the blank v-neighbouring vertices in G2.
        return false;
    }

    return true;
    
} // syntacticCheckTrim

} // namespace GraphAlign

#endif // SEMANTICSYNTACTIC_H