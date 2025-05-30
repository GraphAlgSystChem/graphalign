#ifndef MGA_ANALYSIS_H
#define MGA_ANALYSIS_H

// std
#include <vector>
#include <utility>
#include <string> 
#include <queue>
#include <map> 
#include <unordered_map>
#include <cmath>
#include <limits>
#include <memory> // for std::shared_ptr
#include <set> 
#include <algorithm> // set_intersection
#include <chrono> // Time checking

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

// us
#include "FatalException.hpp"
#include "PairToRangeAdaptor.hpp"
#include "Combinatorics.hpp"
#include "AlignObj.hpp"
#include "AlignParameters.hpp"
#include "AlignUtility.hpp"
#include "Trim.hpp"
#include "Expand.hpp"
#include "FilteredGraph.hpp"
#include "GlobalVariables.hpp"
#include "ComputeOrder.hpp"
#include "BuildAlignment.hpp"
#include "Scoring.hpp"
#include "GuideTree.hpp"
#include "SemanticSyntactic.hpp"
#include "GenerateCandidates.hpp"
#include "BranchBound.hpp"
#include "../test/GraphIO.hpp"

namespace GraphAlign {
/* 
    1) Guide Tree w. alignment vertices (parameterized to allow for customized construction order)
        e.g. with simillarity matrix (graph kernels), randomized etc.
    2) Array containing order of alignment cherries.
    3) Progressive alignment (recursive or trimming).
    4) Result should be an alignment graph with the final A-matrix.
*/

/**
 * @brief Function to be called by users, handles all the necessary function calls.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @param inputGraphs List of the input graphs to be aligned
 * @param scoreVertices Function to score the vertices, called based on scoringMethod
 * @param kernelType type of kernel to compute. Either delta (true) or linear (false).
 * @param isSparse Boolean indicating whether the input graphs are sparse (true) or dense (false).
 * @param measureType Boolean to determine whether to minimize distance (true) or maximize similarity (false).
 * @param enforceVLabels Boolean indicating whether vertices can only be mapped if they have the same label.
 * @param enforceELabels Boolean indicating whether edges can only be mapped if they have the same label.
 * @param useAmbiguousEdges Boolean indicating whether ambiguous edges can be used for connectivity.
 * @param oneMatch Boolean indicating whether only one match is wanted when finding MCSs or all maximal MCSs.
 * @param algorithm String stating what algorithm to use, either "trim" or "expand"
 * @param scoringMethod String indicating the scoring method, can be "order", "scheme"
 * @param anchor List of list of vectors consisting of the known mappings, where each list
 *               in the outer list represents an anchor vertex, where the inner lists maps the vertices in the 
 *               input graphs.
 * @return Final alignment graph, the final score, and a projection map which maps from each graph each vertex to the corresponding vertex
 *         in the alignment graph.
 */
template<typename Graph, class SV, 
        typename Vertex = boost::graph_traits<Graph>::vertex_descriptor,
        typename Edge = boost::graph_traits<Graph>::edge_descriptor>
std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> 
    runAlgorithm(const std::vector<Graph> &inputGraphs,
    SV scoreVertices,
    bool kernelType,
    bool isSparse,
    bool measureType,
    bool enforceVLabels = true,
    bool enforceELabels = true,
    bool useAmbiguousEdges = true,
    bool oneMatch = true,
    const std::string &algorithm = "trimming",
    const std::string &scoringMethod = "order",
    const std::vector<std::vector<Vertex>> &anchor = {});

// To reach this point we need input graphs to be in of the AlignObj<Graph> type, something the user needs to handle 
// before calling PGA (we should provide a default option for this).
// It is also necessary to have created the guide tree plus the queue.
// The score values should maybe be passed differently than this.
/**
 * @brief The function that carry out the entire alignment algorithm, returns the final AlignObj 
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type
 * @param Q A queue for the order in which alignments needs to be created, this should be returned by the getGuideTree function
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param anchor List of list of vertices, representing the anchor across the input graphs
 * @return Shared pointer to the final alignment object
 */
template<typename Graph, class SV,
        typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, 
        typename Edge = boost::graph_traits<Graph>::edge_descriptor>
std::shared_ptr<AlignObj<Graph>> PGA(std::queue<std::shared_ptr<AlignObj<Graph>>> &Q, 
                                    const AlignParameters<Graph, SV> &paramList, 
                                    GlobalVariables<Vertex> &globalVars,
                                    const std::vector<std::vector<Vertex>> &anchor = {}); 

/**
 * @brief The function that carries out the entire alignment algorithm. In this version, the guide tree
 * is not pre-computed but is instead constructed dynamically by computing the similarty scores for 
 * new alignment objects and the yet-to-be aligned alignment graphs.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type
 * @param paramList Object of struct AlignParameters that contains the input parameters psased from the user
 * @param globalVars The struct containing the global variables needed in the alignment algorithms.
 * @param anchor The set containing sets of anchored vertices in each graph.
 */
template<typename Graph, class SV, 
    typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, 
    typename Edge = boost::graph_traits<Graph>::edge_descriptor>
std::shared_ptr<AlignObj<Graph>> PGA_Dynamic(const AlignParameters<Graph, SV> &paramList, 
                                        GlobalVariables<Vertex> &globalVars,
                                        const std::vector<std::vector<Vertex>> &anchor = {});                                     

/**
 * @brief Does exactly the same as runAlgorithm except the kernelType is omitted
 * as all possible guide tree combinations will be computed and executed.
 */
template<typename Graph, class SV,
        typename Vertex = boost::graph_traits<Graph>::vertex_descriptor,
        typename Edge = boost::graph_traits<Graph>::edge_descriptor>
std::vector<std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>, std::chrono::milliseconds>>
    runAlgorithmChooseGuideTree(const std::vector<Graph> &inputGraphs,
    SV scoreVertices,
    bool isSparse,
    bool enforceVLabels = true,
    bool enforceELabels = true,
    bool useAmbiguousEdges = true,
    bool oneMatch = true,
    const std::string &algorithm = "trimming",
    const std::string &scoringMethod = "order",
    const std::vector<std::vector<Vertex>> &anchor = {},
    const std::string &treeString = "");



//////////////////// IMPLEMENTATION BELOW ///////////////////////////

template<typename Graph, class SV, typename Vertex, typename Edge>
std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>>
    runAlgorithm(const std::vector<Graph> &inputGraphs,
    SV scoreVertices,
    bool kernelType,
    bool isSparse,
    bool measureType,
    bool enforceVLabels,
    bool enforceELabels,
    bool useAmbiguousEdges,
    bool oneMatch,
    const std::string &algorithm,
    const std::string &scoringMethod,
    const std::vector<std::vector<Vertex>> &anchor) {
    
    Algorithm alg;
    Scoring scr;
    if(algorithm == "trim") {
        alg = Algorithm::TRIM;
    }
    else if(algorithm == "expand") {
        alg = Algorithm::REC;
    }
    else {
        std::string msgString = "Not a correct input for the parameter algorithm, given " + algorithm + "\n";
        const char *msg = msgString.c_str();
        throw FatalException(msg);
    }
    if(scoringMethod == "order") {
        scr = Scoring::ORDER;
    }
    else if(scoringMethod == "scheme") {
        scr  = Scoring::SCHEME;
    }
    else {
        std::string msgString = "Not a correct input for the parameter scoringMethod, given " + scoringMethod + "\n";
        const char *msg = msgString.c_str();
        throw FatalException(msg);
    }

    AlignParameters paramList(alg, scr, enforceVLabels, enforceELabels, useAmbiguousEdges, oneMatch, inputGraphs, scoreVertices, measureType);
    
    std::map<std::string, int> labelMap;
    std::vector<float> labelScore;
    int nLabels = 0;
    // Save information in regard to number of labels and label scores when enforcing vertex labels, as this is used in branch and bound check
    if(enforceVLabels) {
        // find the largest label index stored in .labelIndex
        for(const Graph &inputGraph : inputGraphs) {
            for(Vertex x: asRange(vertices(inputGraph))) {
                if(inputGraph[x].labelIndex > nLabels) {
                    nLabels = inputGraph[x].labelIndex;
                }
            }
        }
        // add one to the largest label index found to find the number different labels used
        // as first label is 0
        nLabels++;
        
        // Add dummy values to labelScore to be able to index into with the labelIndex value
        for(int i = 0; i < nLabels; i++) {
            labelScore.push_back(0);
        }
        
        // Find the label score of matching two vertices with the same label
        int seenLabels = 0;
        for(const Graph &inputGraph : inputGraphs) {
            // if the number of seen labels is equal to the number of labels, all label scores have been calculated, so break
            if(seenLabels == nLabels) break;
            for(Vertex x: asRange(vertices(inputGraph))) {
                if(labelScore[inputGraph[x].labelIndex] == 0) {
                    float score = std::invoke(scoreVertices, x, x, inputGraph, inputGraph);
                    labelScore[inputGraph[x].labelIndex] = score;
                    labelMap[inputGraph[x].label] = inputGraph[x].labelIndex;
                    seenLabels++;

                    if(seenLabels == nLabels) break;
                }
            }
        }
    }

    GlobalVariables<Vertex> globalVars(labelScore, nLabels);

    
    boost::numeric::ublas::matrix<float> simScore = kernelCalculation(inputGraphs, kernelType, isSparse);
    
    std::shared_ptr<AlignObj<Graph>> res;
    bool useDynamic = false;
    if(useDynamic){
        res = PGA_Dynamic<Graph>(paramList, globalVars, anchor);
    } else {
        std::queue<std::shared_ptr<AlignObj<Graph>>> guideTreeQueuePair = getGuideTree<Graph>(inputGraphs, simScore, measureType);
        std::queue<std::shared_ptr<AlignObj<Graph>>> Q = guideTreeQueuePair;
        res = PGA(Q, paramList, globalVars, anchor);
    }
    std::cout << "Guide tree used: " << res->tag << "\n";
    

    Graph alignmentGraph;
    float finalScore;
    std::unordered_map<int, std::unordered_map<Vertex, Vertex>> projections;
    // Unpack results
    // Mapping from an input graph -> vertex in that input graph -> vertex in final alignment graph
    alignmentGraph = res->alignGraph;
    finalScore = res->alignScore;

    // the input graphs and then the id of the vertices from the alignment graph 
    for(std::size_t i = 0; i < res->contained.size(); i++) {
        // For each input graph go through each vertex in the alignment graph and check if it exists in the 
        // given input graph, if it does add it to the projections map
        for(Vertex x: asRange(vertices(alignmentGraph))) {
            Vertex inputVertex = res->matchTable[i][x];
            if(inputVertex != boost::graph_traits<Graph>::null_vertex()) {
                projections[res->contained[i]][inputVertex] = x;
            }
        }
    }
    
    std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> toReturn = std::make_tuple(alignmentGraph, finalScore, projections);
    return toReturn;
} // runAlgorithm


template<typename Graph, class SV, typename Vertex, typename Edge>
std::shared_ptr<AlignObj<Graph>> PGA(std::queue<std::shared_ptr<AlignObj<Graph>>> &Q, 
                                    const AlignParameters<Graph, SV> &paramList,
                                    GlobalVariables<Vertex> &globalVars,
                                    const std::vector<std::vector<Vertex>> &anchor) {
                                        
    std::cout << "Progressive Alignment using Static Guide Tree\n";

    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    std::shared_ptr<AlignObj<Graph>> root = Q.back();

    while(Q.size() > 0){
        std::shared_ptr<AlignObj<Graph>> toAlign = Q.front();
        
        std::shared_ptr<AlignObj<Graph>> leftAlign = toAlign->left; // toAlign->left;
        std::shared_ptr<AlignObj<Graph>> rightAlign = toAlign->right; // toAlign->right;

        std::vector<std::pair<Vertex, Vertex>> preMap = getAnchorPreMap(*leftAlign, *rightAlign, anchor);
        std::vector<std::vector<std::pair<Vertex, Vertex>>> subgraphs;
        std::vector<Vertex> finalLeftMatch;
        std::vector<Vertex> finalRightMatch;

        globalVars.bestScore = 0.0;
        globalVars.matchSize = preMap.size();

        std::cout << "Working on alignment of " << leftAlign->tag << " and " << rightAlign->tag << std::endl;
        
        std::chrono::system_clock::time_point timerStart;
        std::chrono::system_clock::time_point timerEnd;

        if(paramList.algorithm == Algorithm::REC){
            timerStart = std::chrono::high_resolution_clock::now();
            
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> allMatches;
            // Create vectors containing the current matches going from one graph to the other and vice versa
            // These are maintained during execution so as to prevent constant memory usage in recursive calls
            std::vector<Vertex> leftMatch(num_vertices(leftAlign->alignGraph), boost::graph_traits<Graph>::null_vertex());
            std::vector<Vertex> rightMatch(num_vertices(rightAlign->alignGraph), boost::graph_traits<Graph>::null_vertex());
            // Add the anchor the leftMatch and rightMatch
            for(const std::pair<Vertex, Vertex> &match : preMap ) {
                leftMatch[getIndex(match.first)] = match.second;
                rightMatch[getIndex(match.second)] = match.first;
            }


            if(num_vertices(leftAlign->alignGraph) >= num_vertices(rightAlign->alignGraph)) {
                std::vector<Vertex> anchorVertices;
                for(auto &pair : preMap) {
                    anchorVertices.push_back(pair.second);
                }

                // Call function to prepare vectors used in branch and bound checking 
                prepareGlobalVars(*leftAlign, *rightAlign, preMap, paramList, globalVars);

                // parentMap cannot be const as it is updated dynamically
                std::vector<Vertex> parentMap(num_vertices(rightAlign->alignGraph), nullVertex);
                // initialise bags array needed for dynamic ordering over the smaller graph
                prepareDynamicOrder(rightAlign->alignGraph, rightMatch, parentMap, globalVars);

                
                // Let the startScore be the score from the anchor vertices, and only the anchor vertices
                // No scoring of previous matchings, if this is an alignment of alignments
                // If scoring by order we only want the anchor vertices
                float startScore = 0;
                if(paramList.scoring == Scoring::SCHEME) {
                    startScore = leftAlign->alignScore + rightAlign->alignScore;
                }
                for(const auto &pair: preMap) {
                    startScore = startScore + scoreNewMatch(*leftAlign, *rightAlign, pair.first, pair.second, paramList);
                }
                expand(*leftAlign, *rightAlign, leftMatch, rightMatch, nullVertex, nullVertex, startScore, allMatches, paramList, globalVars, parentMap);

                finalLeftMatch = allMatches[0].first;
                finalRightMatch = allMatches[0].second;
            }
            else {
                std::vector<Vertex> anchorVertices;
                for(auto &pair : preMap) {
                    anchorVertices.push_back(pair.first);
                }
                
                if(paramList.enforceVLabels) {
                    // flip the anchor mapping, as the function to prepare the vectors in globalVars expect the anchor to be from 
                    // the first graph to the second graph.
                    std::vector<std::pair<Vertex, Vertex>> flippedPreMap;
                    for(const auto &pair : preMap) {
                        flippedPreMap.push_back({pair.second, pair.first});
                    }

                    prepareGlobalVars(*rightAlign, *leftAlign, flippedPreMap, paramList, globalVars);
                }

                std::vector<Vertex> parentMap(num_vertices(leftAlign->alignGraph), nullVertex);
                prepareDynamicOrder(leftAlign->alignGraph, leftMatch, parentMap, globalVars);

                
                float startScore = 0;
                if(paramList.scoring == Scoring::SCHEME) {
                    startScore = leftAlign->alignScore + rightAlign->alignScore;
                }
                for(const auto &pair: preMap) {
                    startScore = startScore + scoreNewMatch(*leftAlign, *rightAlign, pair.first, pair.second, paramList);
                }
                expand(*rightAlign, *leftAlign, rightMatch, leftMatch, nullVertex, nullVertex, startScore, allMatches, paramList, globalVars, parentMap);

                finalLeftMatch = allMatches[0].second;
                finalRightMatch = allMatches[0].first;
            }
            
            timerEnd = std::chrono::high_resolution_clock::now();

        }
        else if (paramList.algorithm == Algorithm::TRIM){        
            timerStart = std::chrono::high_resolution_clock::now();
            
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> allMatches = trim(*leftAlign, *rightAlign, paramList, globalVars, preMap);
            // converting the result match vectors to vector of pairs to use old version of buildAlignment
            finalLeftMatch = allMatches[0].first;
            finalRightMatch = allMatches[0].second;

            timerEnd = std::chrono::high_resolution_clock::now();
        } 


        // If reached the final alignment, ensure we score by scheme, to get score of the alignment table, even when using order
        AlignParameters<Graph, SV> adjustedParamList = paramList;
        if(Q.size() == 1) {
            adjustedParamList.scoring = Scoring::SCHEME;
        }
        
        // Build alignment        
        buildAlignment(*toAlign, *leftAlign, *rightAlign, finalLeftMatch, finalRightMatch, adjustedParamList);
        
        std::cout << " taking time (sec): " << (std::chrono::duration_cast<std::chrono::microseconds>
                    (timerEnd - timerStart).count())/1000000.0 << "\n";

        // printMatchTable(*toAlign);
        Q.pop();
    }

    return root;
} // PGA

template<typename Graph, typename MatchTable, class SV, typename Vertex, typename Edge>
std::shared_ptr<AlignObj<Graph>> PGA_Dynamic(const AlignParameters<Graph, SV> &paramList, 
                                    GlobalVariables<Vertex> &globalVars,
                                    const std::vector<std::vector<Vertex>> &anchor){
    
    std::cout << "Progressive Alignment using Dynamic Guide Tree\n";
    // Compute the distance matrices for each of the input graphs.
    std::vector<std::vector<std::vector<int>>> floydMatrices = getDistanceMatrices(paramList.inputGraphs);
    size_t n = paramList.inputGraphs.size();
    float measureSum = 0;

    // Create cluster graph. Each vertex in the cluster graph represents an alignment object. 
    // The associated floyd matrix is stored with the vertex. Ensures that recomputation is unnecessary.
    struct CedgeProps{
        float measure;
    };
    struct CvertexProps{
        std::shared_ptr<AlignObj<Graph>> alignObj;
        std::shared_ptr<std::vector<std::vector<int>>> floydMatrix;
        float ownDistance;
    };
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, CvertexProps, CedgeProps> C;
    C clusterGraph;
    initialiseClusterGraph<Graph>(clusterGraph, paramList.inputGraphs, floydMatrices, paramList.measureType);


    // Extract the minimal (distance) or maximal (similarity) pair
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    using cVertex = typename C::vertex_descriptor;
    while (num_vertices(clusterGraph) > 1){
        // Extract next pair of alignment graphs to align
        std::pair<cVertex, cVertex> nextPair = getNextTreePair(clusterGraph, paramList.measureType);
        auto node1 = source(nextPair, clusterGraph);
        auto node2 = target(nextPair, clusterGraph);
        measureSum = measureSum + clusterGraph[edge(node1, node2, clusterGraph).first].measure;
        std::cout << "Working on alignment of " << clusterGraph[node1].alignObj->tag << " and " << clusterGraph[node2].alignObj->tag << "\n";

        // Preamble
        std::vector<std::pair<Vertex, Vertex>> preMap = getAnchorPreMap((*clusterGraph[node1].alignObj), (*clusterGraph[node2].alignObj), anchor);
        std::vector<Vertex> finalLeftMatch;
        std::vector<Vertex> finalRightMatch;
        globalVars.bestScore = 0.0;
        globalVars.matchSize = preMap.size();
        std::chrono::system_clock::time_point timerStart;
        std::chrono::system_clock::time_point timerEnd;

        // PGA
        if(paramList.algorithm == Algorithm::REC){
            timerStart = std::chrono::high_resolution_clock::now();
            
            // Vector of all matches (match representation is a vector)
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> allMatches;
            std::vector<Vertex> leftMatch(num_vertices(clusterGraph[node1].alignObj->alignGraph), boost::graph_traits<Graph>::null_vertex());
            std::vector<Vertex> rightMatch(num_vertices(clusterGraph[node2].alignObj->alignGraph), boost::graph_traits<Graph>::null_vertex());
        
            // Add the anchor the leftMatch and rightMatch
            for(const std::pair<Vertex, Vertex> &match : preMap ) {
                leftMatch[getIndex(match.first)] = match.second;
                rightMatch[getIndex(match.second)] = match.first;
            }


            if(num_vertices(clusterGraph[node1].alignObj->alignGraph) >= num_vertices(clusterGraph[node2].alignObj->alignGraph)) {
                std::vector<Vertex> anchorVertices;
                for(auto &pair : preMap) {
                    anchorVertices.push_back(pair.second);
                }

                // Call function to prepare vectors used in branch and bound checking 
                prepareGlobalVars(*clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, preMap, paramList, globalVars);
                prepareUnmatchedNeighbours(clusterGraph[node1].alignObj->alignGraph, globalVars);

                // parentMap cannot be const as it is updated dynamically
                std::vector<Vertex> parentMap(num_vertices(clusterGraph[node2].alignObj->alignGraph), nullVertex);
                // initialise bags array needed for dynamic ordering over the smaller graph
                prepareDynamicOrder(clusterGraph[node2].alignObj->alignGraph, rightMatch, parentMap, globalVars);

                
                // Let the startScore be the score from the anchor vertices, and only the anchor vertices
                // No scoring of previous matchings, if this is an alignment of alignments
                float startScore = clusterGraph[node1].alignObj->alignScore + clusterGraph[node2].alignObj->alignScore;
                for(const auto &pair: preMap) {
                    startScore = startScore + scoreNewMatch(*clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, pair.first, pair.second, paramList);
                }
                expand(*clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, leftMatch, rightMatch, nullVertex, nullVertex, startScore, allMatches, paramList, globalVars, parentMap);

                finalLeftMatch = allMatches[0].first;
                finalRightMatch = allMatches[0].second;
            }
            else {
                std::vector<Vertex> anchorVertices;
                for(auto &pair : preMap) {
                    anchorVertices.push_back(pair.first);
                }
                
                if(paramList.enforceVLabels) {
                    // flip the anchor mapping, as the function to prepare the vectors in globalVars expect the anchor to be from 
                    // the first graph to the second graph.
                    std::vector<std::pair<Vertex, Vertex>> flippedPreMap;
                    for(const auto &pair : preMap) {
                        flippedPreMap.push_back({pair.second, pair.first});
                    }

                    prepareGlobalVars(*clusterGraph[node2].alignObj, *clusterGraph[node1].alignObj, flippedPreMap, paramList, globalVars);
                }
                prepareUnmatchedNeighbours(clusterGraph[node2].alignObj->alignGraph, globalVars);

                std::vector<Vertex> parentMap(num_vertices(clusterGraph[node1].alignObj->alignGraph), nullVertex);
                prepareDynamicOrder(clusterGraph[node1].alignObj->alignGraph, leftMatch, parentMap, globalVars);

                float startScore = clusterGraph[node1].alignObj->alignScore + clusterGraph[node2].alignObj->alignScore;
                for(const auto &pair: preMap) {
                    startScore = startScore + scoreNewMatch(*clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, pair.first, pair.second, paramList);
                }
                expand(*clusterGraph[node2].alignObj, *clusterGraph[node1].alignObj, rightMatch, leftMatch, nullVertex, nullVertex, startScore, allMatches, paramList, globalVars, parentMap);

                finalLeftMatch = allMatches[0].second;
                finalRightMatch = allMatches[0].first;
            }

            timerEnd = std::chrono::high_resolution_clock::now();
            

        } else if(paramList.algorithm == Algorithm::TRIM){
            std::vector<std::pair<std::vector<Vertex>, std::vector<Vertex>>> allMatches = trim(*clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, paramList, globalVars, preMap);
            finalLeftMatch = allMatches[0].first;
            finalRightMatch = allMatches[0].second;
        }

       
        std::vector<int> toAlignContained = clusterGraph[node1].alignObj->contained;
        toAlignContained.insert(toAlignContained.end(), clusterGraph[node2].alignObj->contained.begin(), clusterGraph[node2].alignObj->contained.end());
        std::shared_ptr<AlignObj<Graph>> toAlign = std::make_shared<AlignObj<Graph>>(clusterGraph[node1].alignObj, clusterGraph[node2].alignObj, toAlignContained);
        toAlign->tag = "(" + clusterGraph[node1].alignObj->tag + "," + clusterGraph[node2].alignObj->tag + ")";

        if(num_vertices(clusterGraph) == 2){
            AlignParameters<Graph, SV> adjustedParamList = paramList; // insert the actual scheme score in the final alignment for inspection
            adjustedParamList.scoring = Scoring::SCHEME;
            buildAlignment(*toAlign, *clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, finalLeftMatch, finalRightMatch, adjustedParamList);
        } else {
            buildAlignment(*toAlign, *clusterGraph[node1].alignObj, *clusterGraph[node2].alignObj, finalLeftMatch, finalRightMatch, paramList);
        }
        std::cout << " taking time (sec): " << (std::chrono::duration_cast<std::chrono::microseconds>
            (timerEnd - timerStart).count())/1000000.0 << "\n";

        // Remove the previous nodes, add one node and all edges from this node to all the other nodes
        clear_vertex(node1, clusterGraph);
        clear_vertex(node2, clusterGraph);
        if(node1 < node2){
            remove_vertex(node2, clusterGraph);
            remove_vertex(node1, clusterGraph);
        } else {
            remove_vertex(node1, clusterGraph);
            remove_vertex(node2, clusterGraph);
        }
        std::vector<Graph> toAlignGraph = { toAlign->alignGraph };
        std::vector<std::vector<int>> alignDistMatrix = getDistanceMatrices(toAlignGraph)[0];
        auto alignVertex = add_vertex(clusterGraph);
        clusterGraph[alignVertex].alignObj = toAlign;
        clusterGraph[alignVertex].floydMatrix = std::make_shared<std::vector<std::vector<int>>>(alignDistMatrix);
        clusterGraph[alignVertex].ownDistance = getKernelScore(toAlign->alignGraph, toAlign->alignGraph, alignDistMatrix, alignDistMatrix, "delta");
        float measure = 0;
        float pairValue = 0;
        for(const auto u : asRange(vertices(clusterGraph))){
            if(!edge(u, alignVertex, clusterGraph).second && u != alignVertex){
                pairValue = getKernelScore(clusterGraph[u].alignObj->alignGraph, 
                                            toAlign->alignGraph,
                                            *clusterGraph[u].floydMatrix,
                                            *clusterGraph[alignVertex].floydMatrix, 
                                            "delta");
                if(paramList.measureType){
                    measure = sqrt( clusterGraph[u].ownDistance + clusterGraph[alignVertex].ownDistance - 2 * pairValue);
                } else {
                    measure = pairValue;
                }
                auto e = add_edge(u, alignVertex, clusterGraph);
                clusterGraph[e.first].measure = measure;
            }
        }

        // repeat until only one node in the cluster graph remains
    } // endwhile 

    return clusterGraph[0].alignObj;
} // PGA_Dynamic

template<typename Graph, class SV,  typename Vertex, typename Edge>
std::vector<std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>, std::chrono::milliseconds>>
    runAlgorithmChooseGuideTree(const std::vector<Graph> &inputGraphs,
    SV scoreVertices,
    bool isSparse,
    bool enforceVLabels,
    bool enforceELabels,
    bool useAmbiguousEdges,
    bool oneMatch,
    const std::string &algorithm,
    const std::string &scoringMethod,
    const std::vector<std::vector<Vertex>> &anchor,
    const std::string &treeString) {
    
    Algorithm alg;
    Scoring scr;
    if(algorithm == "trim") {
        alg = Algorithm::TRIM;
    }
    else if(algorithm == "expand") {
        alg = Algorithm::REC;
     }
    else {
        std::string msgString = "Not a correct input for the parameter algorithm, given " + algorithm + "\n";
        const char *msg = msgString.c_str();
        throw FatalException(msg);
    }
    if(scoringMethod == "order") {
        scr = Scoring::ORDER;
    }
    else if(scoringMethod == "scheme") {
        scr  = Scoring::SCHEME;
    }
    else {
        std::string msgString = "Not a correct input for the parameter scoringMethod, given " + scoringMethod + "\n";
        const char *msg = msgString.c_str();
        throw FatalException(msg);
    }

    std::vector<std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>, std::chrono::milliseconds>> allAlignments;
    AlignParameters paramList(alg, scr, enforceVLabels, enforceELabels, useAmbiguousEdges, oneMatch, inputGraphs, scoreVertices, false); // "dummy false" as guide tree has been established already.

    std::map<std::string, int> labelMap;
    std::vector<float> labelScore;
    int nLabels = 0;
    // Save information in regard to number of labels and label scores when enforcing vertex labels, as this is used in branch and bound check
    if(enforceVLabels) {
        // find the largest label index stored in .labelIndex
        for(const Graph &inputGraph : inputGraphs) {
            for(Vertex x: asRange(vertices(inputGraph))) {
                if(inputGraph[x].labelIndex > nLabels) {
                    nLabels = inputGraph[x].labelIndex;
                }
            }
        }
        // add one to the largest label index found to find the number different labels used
        // as first label is 0
        nLabels++;
        
        // Add dummy values to labelScore to be able to index into with the labelIndex value
        for(int i = 0; i < nLabels; i++) {
            labelScore.push_back(0);
        }
        
        // Find the label score of matching two vertices with the same label
        int seenLabels = 0;
        for(const Graph &inputGraph : inputGraphs) {
            // if the number of seen labels is equal to the number of labels, all label scores have been calculated, so break
            if(seenLabels == nLabels) break;
            for(Vertex x: asRange(vertices(inputGraph))) {
                if(labelScore[inputGraph[x].labelIndex] == 0) {
                    float score = std::invoke(scoreVertices, x, x, inputGraph, inputGraph);
                    labelScore[inputGraph[x].labelIndex] = score;
                    labelMap[inputGraph[x].label] = inputGraph[x].labelIndex;
                    seenLabels++;

                    if(seenLabels == nLabels) break;
                }
            }
        }
    }

    GlobalVariables<Vertex> globalVars(labelScore, nLabels);
    


    std::vector<std::queue<std::shared_ptr<AlignObj<Graph>>>> guideTrees;
    if(treeString.empty()){
        guideTrees = GraphAlignUtility::getAllGuideTrees<Graph>(inputGraphs);
        std::cout << "N GUIDETREES: " << guideTrees.size() << std::endl;
    } else{
        // compute guide tree here
        std::queue<std::shared_ptr<AlignObj<Graph>>> tree = getGuideTreeFromString<Graph>(inputGraphs, treeString);
        guideTrees.push_back(tree);
    }

    std::chrono::system_clock::time_point PGAtimerStart;
    std::chrono::system_clock::time_point PGAtimerEnd;
    int guideTreeNum = 0;
    for(auto &Q: guideTrees){
        std::shared_ptr<AlignObj<Graph>> root = Q.back();
        std::cout << "\tStarting GuideTree (" << guideTreeNum << "): " << root->tag << std::endl;
        guideTreeNum++;
        
        PGAtimerStart = std::chrono::high_resolution_clock::now();
        std::shared_ptr<AlignObj<Graph>> res = PGA(Q, paramList, globalVars, anchor);
        PGAtimerEnd = std::chrono::high_resolution_clock::now();
        
        // Unpack results
        // Mapping from an input graph -> vertex in that input graph -> vertex in final alignment graph
        Graph alignmentGraph = res->alignGraph;
        float finalScore = res->alignScore;
        std::unordered_map<int, std::unordered_map<Vertex, Vertex>> projections;
        // the input graphs and then the id of the vertices from the alignment graph 
        for(std::size_t i = 0; i < res->contained.size(); i++) {
            // For each input graph go through each vertex in the alignment graph and check if it exists in the 
            // given input graph, if it does add it to the projections map
            for(Vertex x: asRange(vertices(alignmentGraph))) {
                Vertex inputVertex = res->matchTable[i][x];
                if(inputVertex != boost::graph_traits<Graph>::null_vertex()) {
                    projections[res->contained[i]][inputVertex] = x;
                }
            }
        }
        std::chrono::milliseconds duration = std::chrono::duration_cast<std::chrono::milliseconds>(PGAtimerEnd - PGAtimerStart);
        std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>, std::chrono::milliseconds> toReturn = std::make_tuple(alignmentGraph, finalScore, projections, duration);
        allAlignments.push_back(toReturn);
        std::cout << std::endl;
        std::cout << "\tGuideTree took (sec): " << duration.count()/1000.0 << std::endl;
        std::cout << "\tFinal score: " << finalScore << std::endl;
        std::cout << "\t===================" << std::endl;
    }

    return allAlignments;
    
} // runAlgorithmChooseGuideTree


} // End namespace GraphAlign

#endif // MGA_ANALYSIS_H
