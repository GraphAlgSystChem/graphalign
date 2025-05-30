#ifndef ALIGNTRIMMING_H
#define ALIGNTRIMMING_H

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

// our code
#include "AlignObj.hpp"
#include "AlignUtility.hpp"

namespace GraphAlign {

// This file include the following functions:
// trim
// genSubsets
// handleSubsets
// pseudoVF

////////////////////////////////// DECLARATIONS //////////////////////////////////

/**
 * @brief The algorithm for pairwise graph alignment known as Trim.
 *        The outer search tree consists of the powerset of all vertices. 
 *        The inner search is a VF-like approach for finding subgraph isomorphisms
 *        between a trimmed graph and the larger graph.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @param preMatch A list of pairs of vertices such that, for each pair the first element is a vertex from a1's underlying graph
 * and the second is a vertex from a2's underlying graph. 
 * @return A list of pairs of match vectors with the same maximal score, where each pair consists of the two vectors, a1Match and a2Match
 *         mapping matched vertices from a1 to the matched vertices in a2, and vice versa. Null vertex if the vertex is not matched.
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> trim(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                        const AlignParameters<Graph, SV> &paramList,
                        GlobalVariables<Vertex> &globalVars,
                        const std::vector<std::pair<Vertex, Vertex>> &preMatch = {});

/**
 * @brief Function to generate the needed subsets for each k value in iterativeTrimming, 
 *        which specify the size of the subsets of vertices to remove from the smaller graph. 
 *        Calls function handleSubset to preprocess the smaller graph based on the provided subset,
 *        and then call handleSubset and insert found matches in newMatches, most of the parameters are needed by handleSubset,
 *        not by this function itself.
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param f The filtered graph of the a2.alignGraph
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @param removableSmall List of vertices that are possible to remove from the smaller graph
 * @param k The required size of the subsets passed along
 * @param newMatches List to store the found matches in.
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void genSubsets(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                FilteredGraph<Graph> &f,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars,
                const std::vector<Vertex> &removableSmall,
                std::size_t k,
                std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches);

/**
 * @brief Takes a subset of vertices to be removed from the smaller graph between a1 and a2, prepares 
 *        said smaller graph correctly and compute order after which it calls pseudoVF. 
 *        Returns true in all cases except when scoring by order and paramList.oneMatch is true, in which case it returns false
 *        when at least one match has been added to newMatches.
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param f The filtered graph for the a2 graph
 * @param subset List of the vertices to remove from the smaller graph
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param GlobalVars The object representing global function parameters. Used to keep label IDs consistent between the labelScore vector and the matrix.
 * @param newMatches Reference to vector to add the found matches to
 * @return true Any time when not using order as scoring method or a match has not been found
 * @return false When using order and oneMatch and a match has been found 
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool handleSubset(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const FilteredGraph<Graph> &f,
                    const std::vector<Vertex> &subset,
                    const AlignParameters<Graph, SV> &paramList, 
                    GlobalVariables<Vertex> &globalVars,
                    std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches);

/**
 * @brief Finds the maximal MCs between g1 and the subgraph, based on VF2+ and VF3 optimizations to limit 
 *        the number of candidates considered with parent map. Match representation stored in globalVars
 * 
 * @tparam Graph Graph type
 * @tparam SV he type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param subgraph A candidate graph for being a subgraph in g1 (a subgrap of a2.alignGraph)
 * @param ambiguousSubgraph A list of the ambiguous edges between the vertices in the subgraph
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars Object of struct GlobalVariables containing information used by the the various algoritms
 * @param order A list consisting of the total order over the vertices in subgraph, 
 *                   such that the vertices are considered in the order of the lsit
 * @param parent A map indicating the parent of any vertex in the subgraph according to the total order
 * @param newMatches A list of the maximal size matches that have been found so far - this list is updated recursively.
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void pseudoVF(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const FilteredGraph<Graph> &subgraph,
                const std::vector<std::pair<Vertex, Vertex>> &ambiguousSubgraph,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars,
                const std::vector<Vertex> &order,
                const std::vector<Vertex> &parent,
                std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////

template<typename Graph, class SV, typename Vertex>
std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> trim(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                        const AlignParameters<Graph, SV> &paramList,
                        GlobalVariables<Vertex> &globalVars,
                        const std::vector<std::pair<Vertex, Vertex>> &preMatch) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    std::vector<Vertex> a1Match(num_vertices(a1.alignGraph), nullVertex);
    std::vector<Vertex> a2Match(num_vertices(a2.alignGraph), nullVertex);
    
    globalVars.anchor = preMatch;

    // Using the prematch to fill out the default match vectors of a1 and a2.
    for(const auto &pair : preMatch) {
        a1Match[getIndex(pair.first)] = pair.second;
        a2Match[getIndex(pair.second)] = pair.first;
    }
    globalVars.a1Match = a1Match;
    globalVars.a2Match = a2Match;

    // The set of all matches found.
    std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> allMatches;
    // If there is an anchor match to start, add this to all matches 
    // and save the score of the anchor match as the best score achieved so far
    
    if(!preMatch.empty()) {
        globalVars.bestScore = matchScore(a1, a2, a1Match, a2Match, paramList);
        globalVars.matchSize = preMatch.size();
        allMatches.push_back({a1Match, a2Match});
    }
    
    int smallestGraphWas; // Index telling us which graph is smallest. 
    Graph smallerGraph; // Need to save the smaller graph to create the filtered graph after
    
    if(num_vertices(a1.alignGraph) >= num_vertices(a2.alignGraph)){
        smallestGraphWas = 2;
        smallerGraph = a2.alignGraph;
        std::vector<Vertex> smallerAnchor;
        for(const auto &pair: preMatch) {
            smallerAnchor.push_back(pair.second);
        }
        globalVars.smallerAnchor = smallerAnchor;
    } else { // G2 > G1
        smallestGraphWas = 1;
        smallerGraph = a1.alignGraph;
        std::vector<Vertex> smallerAnchor;
        for(const auto &pair: preMatch) {
            smallerAnchor.push_back(pair.first);
        }
        globalVars.smallerAnchor = smallerAnchor;
    }
    // Creating the filtered graph for the smaller of the two graphs
    FilteredGraph<Graph> f = FilteredGraph(smallerGraph);

    // For branch and bound to be used we need to prepare various information to be used later
    if(paramList.enforceVLabels) {
        if(smallestGraphWas == 2) {
            prepareGlobalVars(a1, a2, preMatch, paramList, globalVars);
        }
        else {
            std::vector<std::pair<Vertex, Vertex>> flippedPreMatch;
            flippedPreMatch.reserve(preMatch.size());
            for(const auto &match : preMatch) {
                flippedPreMatch.push_back({match.second, match.first});
            }
            prepareGlobalVars(a2, a1, flippedPreMatch, paramList, globalVars);
        }
    }

    // Find the vertices that can be removed from the smaller graph
    std::vector<Vertex> removableSmall;
    if(smallestGraphWas == 2) {
        // Only consider vertices from removal that are not already prematched.
        for( Vertex v : asRange(vertices(a2.alignGraph)) ) {
            if(a2Match[getIndex(v)] == nullVertex) {
                removableSmall.push_back(v);
            }
        }
    }
    else {
        // Only consider vertices from removal that are not already prematched.
        for( Vertex v : asRange(vertices(a1.alignGraph))) {
            if(a1Match[getIndex(v)] == nullVertex) {
                removableSmall.push_back(v);
            }
        }
    }

    for(size_t k = 0; k < removableSmall.size(); k++) {
        std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> newMatches;
        
        // Check if largest hypothetical score achievable at this level is smaller than global best score
        // if so we can stop as it will only get worse by removing more vertices.
        // See documentation for evaluateKLevel to understand how this is done.
        if(smallestGraphWas == 2) {
            if(!evaluateKLevel(a1, a2, k, paramList, globalVars)) {
                return allMatches;
            }
        }
        else {
            if(!evaluateKLevel(a2, a1, k, paramList, globalVars)) {
                return allMatches;
            }
        }

        // (Precondition: The graph in second position is the donor for f)
        // genSubsets produces mapping from the first alignObj to the second alignObj and inserts them into newMatches.
        if(smallestGraphWas == 2){
            genSubsets(a1, a2, f, paramList, globalVars, removableSmall, k, newMatches);
        } else {
            genSubsets(a2, a1, f, paramList, globalVars, removableSmall, k, newMatches);
            // All subgraph isomorphisms are from G2 -> G1 in pseudoVF2, so they must be inverted.
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> reversedMatches;
            for(const auto &morphism : newMatches ){
                std::pair<std::vector<Vertex>, std::vector<Vertex>> reverse = {morphism.second, morphism.first};
                reversedMatches.push_back(reverse);
            }
            newMatches = reversedMatches;
        }
        
        if(newMatches.size() > 0) {
            
            if(paramList.oneMatch && paramList.scoring == Scoring::ORDER) {
                return newMatches;
            }

            for(const auto &match : newMatches) {
                float newScore = matchScore(a1, a2, match.first, match.second, paramList);
                // If a match is better than the best match found so far
                // empty allMatches and insert it as the new best match
                if(newScore > globalVars.bestScore) {
                    globalVars.bestScore = newScore;
                    allMatches.clear();
                    allMatches.push_back(match);
                }
                // Save all the matches giving the same best score
                else if(newScore == globalVars.bestScore) {
                    allMatches.push_back(match);
                }
            }
        }
    }

    return allMatches;
} 

template<typename Graph, class SV, typename Vertex>
void genSubsets(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                FilteredGraph<Graph> &f,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars,
                const std::vector<Vertex> &removableSmall,
                std::size_t k,
                std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches) {
    
    bool cont = true;
    std::size_t end = removableSmall.size();

    // When k == 0, the subset needed is the empty one, so this can be handled by itself
    if (k == 0) {
        std::vector<Vertex> emptySet;
        handleSubset(a1, a2, f, emptySet, paramList, globalVars, newMatches);
        return;
    }

    std::vector<Vertex> combination;
    combination.reserve(k);

    // Lambda function to build the subsets, whenever a subset is built it calls handleSubset
    auto combinationLambda = [&](std::vector<Vertex> &data, int start, std::size_t index, auto&& combinationLambda) {
        // Return type from handleSubset is boolean, used only for early termination when using order 
        if(!cont) return;

        // Ready to build the subset based on the elements in data
        if (index == k) {
            combination.clear();
            for (std::size_t j = 0; j < k; j++) {
                combination.push_back(data[j]);
            }
            f.remove(combination);
            cont = handleSubset(a1, a2, f, combination, paramList, globalVars, newMatches);
            f.insert(combination);
            return;
        }
        
        // while there are more elements not considered in data, 
        // and there is enough remaining elements to fill out the remaining needed spaces.
        for (std::size_t i = start; i < end && end - i >= k - index; i++) {
            data.push_back(removableSmall[i]);
            // New data vector containing the previously added elements, and the next element at position i in vector of removable vertices
            // Adding the next element to the data vector 
            combinationLambda(data, i + 1, index + 1, combinationLambda);
            // remove the last element in data to try a different vertex at this iteration
            data.pop_back();
        }
    };

    std::vector<Vertex> data;
    data.reserve(k);
    combinationLambda(data, 0, 0, combinationLambda);

    return;
} // genSubsets

template<typename Graph, class SV, typename Vertex>
bool handleSubset(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const FilteredGraph<Graph> &f,
                    const std::vector<Vertex> &subset,
                    const AlignParameters<Graph, SV> &paramList, 
                    GlobalVariables<Vertex> &globalVars,
                    std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches) {
    
    
    // If the branch bound check returns false it means the subset should not be explored
    // The reasoning depends on the scoring method:
    // if scoring by scheme it calculates the hypothetical score and if this value is smaller than the best score found, 
    // no reason to consider this subset. Also checks that the larger graph covers all the vertices in the smaller graph (based on labels).
    // If scoring by order the label frequencies are checked to ensure that the larger graph can cover all vertices in the smaller graph (based on labels).
    if(!branchBoundCheck(a1, a2, subset, paramList, globalVars)) {
        return true;
    }

    // Adjust the set of ambiguous edges s.t. only edges among visible vertices are considered.
    // Then, try to find an embedding of the smallest graph in the larger one.
    std::vector<std::pair<Vertex, Vertex>> ambiguousSmallerTrimmed;
    // Over allocating space to prevent reallocation being required.
    ambiguousSmallerTrimmed.reserve(a2.ambiguousEdges.size());
    for(const auto &edge : a2.ambiguousEdges ) {
        if(isVisible(edge.first, f) && isVisible(edge.second, f)) {
            ambiguousSmallerTrimmed.push_back(edge);
        }
    }

    // Compute the order
    auto res = computeOrder(f, globalVars.smallerAnchor);
    std::vector<Vertex> order = res.first;
    std::vector<Vertex> parent = res.second;

    // Vector to store the embeddings in pseudoVF
    std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> subgraphIsos;

    pseudoVF(a1, a2, f, ambiguousSmallerTrimmed, paramList, globalVars, order, parent, subgraphIsos);
    // Adding the found matches to the reference to newMatches, and if scoring by order and only wanting one match return
    newMatches.insert(newMatches.end(), subgraphIsos.begin(), subgraphIsos.end());


    if(paramList.oneMatch && !newMatches.empty() && paramList.scoring == Scoring::ORDER) {
        // Only return false when using order and oneMatch and a match has been found, as a "better" match will never be found
        return false;
    }

    return true;
} // handleSubset

template<typename Graph, class SV, typename Vertex>
void pseudoVF(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const FilteredGraph<Graph> &subgraph,
                const std::vector<std::pair<Vertex, Vertex>> &ambiguousSubgraph,
                const AlignParameters<Graph, SV> &paramList,
                GlobalVariables<Vertex> &globalVars,
                const std::vector<Vertex> &order,
                const std::vector<Vertex> &parent,
                std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> &newMatches) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    std::vector<Vertex> &largerMatch = globalVars.getLargerMatch();
    std::vector<Vertex> &smallerMatch = globalVars.getSmallerMatch();
    // A complete embedding has been found.
    if(globalVars.matchSize == num_visible_vertices(subgraph)) {
        newMatches.push_back({largerMatch, smallerMatch});
    }
    else {
        
        std::vector<std::pair<Vertex, Vertex>> candidatePairs = generateCandidatesTrim(a1.alignGraph, subgraph, globalVars, order, parent);

        for(const std::pair<Vertex, Vertex> &pair : candidatePairs) {

            bool semFeasibility = semanticCheck(a1.alignGraph, subgraph, pair, paramList, globalVars);

            if(semFeasibility) {

                bool synFeasibility = syntacticCheckTrim(a1.alignGraph, subgraph, pair, a1.ambiguousEdges, ambiguousSubgraph, paramList, globalVars);

                if(synFeasibility) {
                    
                    // Add the new candidate pair to the match and increase the match size
                    largerMatch[getIndex(pair.first)] = pair.second;
                    smallerMatch[getIndex(pair.second)] = pair.first;
                    globalVars.matchSize++;

                    pseudoVF(a1, a2, subgraph, ambiguousSubgraph, paramList, globalVars, order, parent, newMatches);

                    // Remove the candidate pair and decrease the match size to be able to try a different match at this stage
                    largerMatch[getIndex(pair.first)] = nullVertex;
                    smallerMatch[getIndex(pair.second)] = nullVertex;
                    globalVars.matchSize--;

                    // Only want the first embedding found
                    if(paramList.oneMatch && !newMatches.empty()) {
                        return;
                    }
                }
            }
        }
    }
    return;
} // pseudoVF

} // namespace GraphAlign

#endif // ALIGNTRIMMING_H