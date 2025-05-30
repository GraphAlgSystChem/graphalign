#ifndef EXPAND_H
#define EXPAND_H

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

// our code
#include "AlignObj.hpp"
#include "AlignParameters.hpp"
#include "GlobalVariables.hpp"
#include "GenerateCandidates.hpp"

namespace GraphAlign {

// This file contains the following functions: 
// expand
// isMaxColoredMatch
// isSubIso

////////////////////////////////// DECLARATIONS //////////////////////////////////

/**
 * @brief The algorithm for pairwise graph alignment known as Expand. This a backtracking algorithm for finding maximal common subgraphs. 
 *        The search tree is guided using dynamic vertex order computation. 
 *        The search tree is pruned using branch-and-bound heuristics.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param a1 The first alignment object, which contains the necessary information for the first graph
 * @param a2 The second alignment object, which contains the necessary information for the second graph  
 * @param a1Match Vector representing the matched vertices in g1, i.e. if a1Match[x] = y, means that the vertex x in g1 is matched to vertex y in g2, is null vertex if not matched
 * @param a2Match Vector representing the matched vertices in g2. 
 * @param u The latest vertex added to the match from g1, can be null vertex if a vertex in the order is skipped
 * @param v The latest vertex added to the match from g2, can be null vertex if a vertex in the order is skipped
 * @param currentScore The score of the current match contained in a1Match, a2Match
 * @param allMatches Vector containing pairs of a1Match, a2Match, i.e. the best matches found so far.
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @param parent Vector of parents for vertices in g2
 * @return Boolean used by by the algorithm itself
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool expand(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
            std::vector<Vertex> &a1Match, std::vector<Vertex> &a2Match,
            Vertex u, Vertex v, 
            float currentScore,
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &allMatches,
            const AlignParameters<Graph, SV> &paramList,
            GlobalVariables<Vertex> &globalVars,
            std::vector<Vertex> &parent);

/**
 * @brief Check if unmatched vertices in g1 and g2 have disjoint sets of labels.
 * When vertex and edge labels are enforced such that mismatched labels cannot be matched
 * together this function allows for early termination. Two vectors are stored in globalVars
 * that keep track of how many vertices there are with each label in the respective graphs, 
 * thereby, for each label at least one of the two entries needs to be zero for the label sets to be disjoint. 
 * 
 * 
 * @tparam Vertex The vertex type 
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms 
 * @return true If the the unmatched vertices in g1 and g2 have disjoint sets of labels
 * @return false Otherwise
 */
template<typename Vertex>
bool isMaxColoredMatch(GlobalVariables<Vertex> &globalVars);

/**
 * @brief Checks if the given match between g1 and g2 is actually unlabelled (subgraph) isomorphism
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type 
 * @param g1 The first graph
 * @param g2 The second graph
 * @param g1Match The match vector for g1, namely g1Match[i] = j if v_i in G1 is mapped to v_j in G2.
 * @param g2Match The match vector for g2, namely g2Match[i] = j if v_i in G2 is mapped to v_j in G1.
 * @return true If an unlabelled isomorphism is found.
 * @return false Otherwise
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
bool isSubIso(const Graph &g1, const Graph &g2, const std::vector<Vertex> &g1Match, const std::vector<Vertex> &g2Match);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template<typename Graph, class SV, typename Vertex>
bool expand(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
            std::vector<Vertex> &a1Match, std::vector<Vertex> &a2Match,
            Vertex u, Vertex v, 
            float currentScore,
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &allMatches,
            const AlignParameters<Graph, SV> &paramList,
            GlobalVariables<Vertex> &globalVars,
            std::vector<Vertex> &parent) {
    
    // precondition, |V(a2)| <= |V(a1)|
    // Boolean indicating whether isMaxColoredMatch or isSubIso return true
    bool foundIsoOrMaxColor = false; 

    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();

    // Determine the number of vertices in the smaller graph
    size_t expOrder = num_vertices(a2.alignGraph);
    
    // If the match provided (u, v) is not null vertices, i.e. an actual match
    // calculate the value of said match
    if(u != nullVertex) {
        // Calculate the score achieved with the new match (u, v)
        currentScore = currentScore + scoreNewMatch(a1, a2, u, v, paramList);

        // If the current match is better than what is already stored in all matches
        // then the current match replaces the content of all matches
        if(currentScore > globalVars.bestScore) {
            allMatches.clear(); 
            allMatches.push_back({a1Match, a2Match});
            globalVars.bestScore = currentScore;
        }
    }

    // Check if the current match includes all the vertices in the smaller graph
    if(globalVars.matchSize < expOrder) {
        if(paramList.enforceVLabels && paramList.enforceELabels) {
            // When enforcing labels it is possible to terminate early if the label sets of the 
            // unmatched vertices in g1 and g2 is disjoint, i.e. if it is not possible to match
            // any of the remaining vertices due to labels
            if(isMaxColoredMatch<Vertex>(globalVars)) {
                foundIsoOrMaxColor = true;
                return foundIsoOrMaxColor;
            }
        }
    }
    // The smaller graph is completely contained in the current match
    else {
        if(paramList.enforceVLabels && paramList.enforceELabels) {
            // When allowing the usage of ambiguous edges it is necessary to check if the current match
            // is a subgraph isomorphism between g1 and g2
            if(paramList.useAmbiguousEdges) {
                foundIsoOrMaxColor = isSubIso(a1.alignGraph, a2.alignGraph, a1Match, a2Match);
            }
            // if ambiguous edges are not used then it is true as soon as all vertices in the smaller graph
            // is matched in the larger graph
            else {
                foundIsoOrMaxColor = true;
            }
        }
    }

    // The current match can potentially be extended (there are unmatched vertices in the smaller graph)
    int k = globalVars.degBagEnds.size() - 2;
    std::size_t numVerticesConsidered = globalVars.degBagEnds[k + 1] - globalVars.degBagEnds[k]; // including vertices that have been skipped
    if(numVerticesConsidered < expOrder) {
        // If the branch and bound check returns false, it means that the current path cannot lead to a better score
        // than bestScore
        if(!branchBoundCheck(a1, a2, currentScore, paramList, globalVars)) return foundIsoOrMaxColor;

        auto pairCallback = [&a1, &a2, &a1Match, &a2Match, currentScore, &allMatches, &paramList, &globalVars, &parent, numVerticesConsidered, expOrder, &nullVertex](std::pair<Vertex, Vertex> candidatePair) {
            bool semFeasibility = semanticCheck(a1.alignGraph, a2.alignGraph, a1Match, a2Match, candidatePair, paramList);

            if(semFeasibility) {
                // Check if syntactic requirements are satisfied, i.e. whether the neighbourhoods allow for the candidate pair to be matched
                bool synFeasibility = syntacticCheckExp(a1, a2, a1Match, a2Match, candidatePair, paramList);

                if(synFeasibility) {
                    // Both syntactic and semantic requirements are met resulting in the candidate pair being added to the current match

                    // Adding the candidates to the match vectors
                    a1Match[getIndex(candidatePair.first)] = candidatePair.second;
                    a2Match[getIndex(candidatePair.second)] = candidatePair.first;
                    globalVars.matchSize++;
                    
                    // The candidate pair should no longer be used for calculating hypothetical score, as they have been matched
                    if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
                        globalVars.notAllowedVerticesG1[getIndex(candidatePair.first)] = 1;
                        globalVars.notAllowedVerticesG2[getIndex(candidatePair.second)] = 1;

                        // Update matrices by decrementing the entries for the candidate pairs based on the label
                        // and since we are in the situation where we enforce labels we know they share the same label
                        int labelIndex = a1.alignGraph[candidatePair.first].labelIndex;
                        globalVars.g1LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG1[getIndex(candidatePair.first)]]--;
                        globalVars.g1LabelMatrix[labelIndex][0]--;

                        globalVars.g2LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG2[getIndex(candidatePair.second)]]--;
                        globalVars.g2LabelMatrix[labelIndex][0]--;
                    }
                    
                    // Maintain the label counters when adding the candidate pair to the current match
                    if(paramList.enforceVLabels && paramList.enforceELabels) {

                        globalVars.labelFrequencyG1[a1.alignGraph[candidatePair.first].labelIndex]--;
                        globalVars.labelFrequencyG2[a2.alignGraph[candidatePair.second].labelIndex]--;
                    }
                    
                    // Call recursively 
                    bool foundIsoOrMaxColor = expand(a1, a2, a1Match, a2Match, candidatePair.first, candidatePair.second, currentScore, allMatches, paramList, globalVars, parent);
                    // Resetting the previous candidates to null to try a different candidate pair from this state
                    a1Match[getIndex(candidatePair.first)] = nullVertex;
                    a2Match[getIndex(candidatePair.second)] = nullVertex;
                    globalVars.matchSize--;
                    
                    // Resetting the indicators for the vertices not allowed to be used for hypothetical score when doing branch and bound checking
                    if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
                        globalVars.notAllowedVerticesG1[getIndex(candidatePair.first)] = 0;
                        globalVars.notAllowedVerticesG2[getIndex(candidatePair.second)] = 0;

                        // Restore the matrices to the state before adding the candidate pair to the match
                        int labelIndex = a1.alignGraph[candidatePair.first].labelIndex;
                        globalVars.g1LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG1[getIndex(candidatePair.first)]]++;
                        globalVars.g1LabelMatrix[labelIndex][0]++;

                        globalVars.g2LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG2[getIndex(candidatePair.second)]]++;
                        globalVars.g2LabelMatrix[labelIndex][0]++;
                    }

                    // maintain the label counters when removing the candidate pair from the current match
                    if(paramList.enforceVLabels && paramList.enforceELabels) {
                        globalVars.labelFrequencyG1[a1.alignGraph[candidatePair.first].labelIndex]++;
                        globalVars.labelFrequencyG2[a2.alignGraph[candidatePair.second].labelIndex]++;
                    }

                    
                    // If labels are enforced and either a subgraph isomoprhism is found
                    // or it is not possible to extend the current match due to disjoint label sets between
                    // unmatched vertices in g1 and g2 break it is not needed to check other candidate pairs 
                    // for this match
                    if(paramList.enforceVLabels && paramList.enforceELabels && foundIsoOrMaxColor) {
                        return foundIsoOrMaxColor;
                    }

                    if(paramList.scoring == Scoring::ORDER) {
                        // If a recursive branch from this stage has achieved a better score, i.e. more vertices matched,
                        // or equals the best possible score from this point, return to parent. 
                        if(globalVars.matchSize + (expOrder - numVerticesConsidered) <= globalVars.bestScore) {
                            return foundIsoOrMaxColor;
                        }
                    }
                }
            }
            return false;
        };
        bool cont; // false if a maxcolor/subiso has been found, true otherwise.
        std::pair<Vertex, int> nextVertexBag = getNextVertex(a2.alignGraph, parent, globalVars);
        Vertex nextVertex = nextVertexBag.first;
        int oldBag = nextVertexBag.second;
        cont = generateCandidatesRec(a1.alignGraph, a2.alignGraph, a1Match, a2Match, parent, nextVertex, pairCallback);
        restoreDegBags(a2.alignGraph, nextVertex, oldBag, parent, globalVars);
        if(!cont) { // either isomorphism (without ambiguous edges) or max colored match.
            return true;
        }

        // If there is still at least two vertices left to match, we can choose to skip the chosen vertex
        // by calling recursively with null vertices. There is no need to skip a single remaining vertex
        // as the change in score is non-existing.
        if(numVerticesConsidered < expOrder - 1) {
            // Check if the number of (unconsidered) vertices left combined with the already matched vertices + skipping this vertex, 
            // is enough to achieve a potentially better score than the currently best score. 
            // If not (in which the condition is true), we simply return from this search branch.
            if(paramList.scoring == Scoring::ORDER) {
                if(globalVars.matchSize + (expOrder - (numVerticesConsidered + 1)) <= globalVars.bestScore) {
                    return false;
                }
            }
            // Maintain the notAllowedVerticesG2, as the vertex at order[orderIndex] is not to be mapped (it is being skipped for this iteration) 
            if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
                globalVars.notAllowedVerticesG2[getIndex(nextVertex)] = 1;

                // When skipping the skipped vertex is no longer valid for later matches 
                // therefore its corresponding entry in g2LabelMatrix needs to be decremented
                int labelIndex = a2.alignGraph[nextVertex].labelIndex;
                globalVars.g2LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG2[getIndex(nextVertex)]]--;
                globalVars.g2LabelMatrix[labelIndex][0]--;
            }
            oldBag = skipVertex(nextVertex, globalVars);
            foundIsoOrMaxColor = expand(a1, a2, a1Match, a2Match, nullVertex, nullVertex, currentScore, allMatches, paramList, globalVars, parent);
            restoreDegBagsSkipped(nextVertex, oldBag, globalVars);
            
            // Maintain the notAllowedVerticesG2 in regard to the previously skipped vertex from the order
            if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
                globalVars.notAllowedVerticesG2[getIndex(nextVertex)] = 0;

                // Restore the g2LabelMatrix entry for the vertex that is no longer being skipped.
                int labelIndex = a2.alignGraph[nextVertex].labelIndex;
                globalVars.g2LabelMatrix[labelIndex][globalVars.nInputGraphsForVertexG2[getIndex(nextVertex)]]++;
                globalVars.g2LabelMatrix[labelIndex][0]++;
            }
        }
    }
    
    return foundIsoOrMaxColor;

} // expand

template<typename Vertex>
bool isMaxColoredMatch(GlobalVariables<Vertex> &globalVars) {

    // The number of vertices matched is the maximum number of vertices that _can be matched_ given the label frequency in each graph.
    // Namely, this function returns false if there exists a vertex in the set of unmatched vertices in G1
    // whose label also exist in the set of unmatched vertices of G2.
    // If no such "label" exists, this matching is a max coloring (maximum sized matched).
    for(int i = 0; i < globalVars.nLabels; i++) {
        if(globalVars.labelFrequencyG1[i] > 0 && globalVars.labelFrequencyG2[i] > 0) {
            return false;
        }
    }

    return true;
}


template<typename Graph, typename Vertex, typename Edge>
bool isSubIso(const Graph &g1, const Graph &g2, const std::vector<Vertex> &g1Match, const std::vector<Vertex> &g2Match) {

    // potential check on vertices and edges up here, but we are not sure why this would make sense?
    // precondition: this is called if and only if one of the two input graphs has been fully matched.
    if(num_edges(g1) != num_edges(g2)) {
        return false;
    }


    // For all edges in g1, if there is not an edge in g2 between matched vertices of edge in g1 return false
    for(Edge e : asRange(edges(g1))) {
        Vertex s = source(e, g1);
        Vertex t = target(e, g1);

        Vertex sMapped = g1Match[getIndex(s)]; // s and t's images in g2, potentially non-existing.
        Vertex tMapped = g1Match[getIndex(t)];
        // Ensure that the edge considered has both its endpoint matched to vertices in the other graph. 
        // Otherwise, this edge cannot impact the subgraph isomorphism check.
        if(sMapped == boost::graph_traits<Graph>::null_vertex() || tMapped == boost::graph_traits<Graph>::null_vertex()) continue;
        std::pair<Edge, bool> edgePair = edge(sMapped, tMapped, g2);
        if(!edgePair.second) {
            return false;
        }
    }
    // For all edges in g2, if there is not an edge in g1 between matched vertices of edge in g2 return false
    for(Edge e : asRange(edges(g2))) {
        Vertex s = source(e, g2);
        Vertex t = target(e, g2);

        Vertex sMapped = g2Match[getIndex(s)]; // s and t's images in g1, potentially non-existing.
        Vertex tMapped = g2Match[getIndex(t)];
        // Note that this check will only fail for one of g1, g2 as one of the graphs has its entire vertex set matched.
        if(sMapped == boost::graph_traits<Graph>::null_vertex() || tMapped == boost::graph_traits<Graph>::null_vertex()) continue;
        std::pair<Edge, bool> edgePair = edge(sMapped, tMapped, g1);
        if(!edgePair.second) {
            return false;
        }
    }

    return true;
} // isSubIso

}

#endif // ALIGNRECURSIVE_H