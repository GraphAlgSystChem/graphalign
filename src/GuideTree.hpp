#ifndef GUIDETREE_H
#define GUIDETREE_H


// 3rd party
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>


namespace GraphAlign {

    // kernelCalculation
    // getGuideTree
    // getSimDistMeasure
    // getGuideTreeFromString
    // getKernelScore

/**
 * @brief Kernel matrix calculation to help build the guide tree based on similarity measure of the input graphs are.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type
 * @param inputGraphs List of input graphs
 * @param kernelType type of kernel to compute. Either linear (0) or delta (1).
 * @param isSparse Boolean indicating whether the input graphs are sparse or dense
 * @return Matrix containing the kernel calculations, 
 * only the upper triangle of the matrix is filled out (M[j][i] contains kernel information iff j >= i). 
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
boost::numeric::ublas::matrix<float> kernelCalculation(const std::vector<Graph> &inputGraphs, const bool kernelType, const bool isSparse = true); 

/**
 * @brief Get the queue to guide the alignment process based on graph kernels.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam Vertex The vertex type
 * @param inputGraphs Vector of input graphs
 * @param kernelScores Matrix of the kernel scoring between the input graphs
 * @param measureType Boolean to determine whether we want to minimize distance or maximize similarity (1 = distance, 0 = similarity)
 * @return Queue determining the order the alignments are to be computed in, i.e. the guide tree
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::queue<std::shared_ptr<AlignObj<Graph>>> 
    getGuideTree(const std::vector<Graph> &inputGraphs, 
                const boost::numeric::ublas::matrix<float> &kernelScores,
                const bool measureType); 

template<typename Graph>
float getSimDistMeasure(const AlignObj<Graph> &g1, const AlignObj<Graph> &g2, const bool measure, const bool kernelType = true, const bool isSparse = true);

template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::queue<std::shared_ptr<AlignObj<Graph>>> 
    getGuideTreeFromString(const std::vector<Graph> &inputGraphs, 
                const std::string &treeString);

template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::shared_ptr<AlignObj<Graph>> _getGuideTreeFromStringAux(const std::vector<Graph> &inputGraphs, const std::string& treeString, std::queue<std::shared_ptr<AlignObj<Graph>>> &Q); 

std::pair<std::string, std::string> _splitGuideTreeString(const std::string &s);

// Computes the kernel measure of g1 and g2 w.r.t. kernelMeasure based on their shortest distance matrices.
template <typename Graph>
int getKernelScore(const Graph &g1, const Graph &g2, const std::vector<std::vector<int>> &g1Dist, const std::vector<std::vector<int>> &g2Dist, const std::string &kernelMeasure);

/**
 * @brief Computes the n x n matrices (n being the number of vertices for a given graph) of pairwise shortest distances for all pairs of vertices
 * 
 * @tparam Graph The graph type used to represent the graphs
 * @tparam Edge The edge descriptor used to represent the edges of the graphs
 * @param gs The grpahs to compute distance matrices for
 */
template <typename Graph, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
std::vector<std::vector<std::vector<int>>> getDistanceMatrices(const std::vector<Graph> &gs);

/**
 * @brief Initialises the (complete) cluser graph, each node representing a (trivial) alignment object along with its distance matrix. Edge weights are the similarity score of each underlying graphs.
 * 
 * @tparam Graph The graph type used to represent the graphs
 * @tparam ClusterGraph The graph type used to represent the ClusterGraph (i.e. with its associated vertex and edge props)
 * @param C The cluster graph
 * @param inputGraphs The graphs to be contained in the cluster graph vertices
 * @param floydMatrices The pairwise distance matrices used to compute the edge weights in C
 * @param measureType The type of similarity measure, e.g. distance or similarity
 */
template <typename Graph, typename ClusterGraph>
void initialiseClusterGraph(ClusterGraph &C, const std::vector<Graph> &inputGraphs, const std::vector<std::vector<std::vector<int>>> &floydMatrices, bool measureType);

/**
 * @brief Returns the next maximum/minimum edge pair in C
 * 
 * @tparam ClusterGraph The graph type used to represent the ClusterGraph (i.e. with its associated vertex and edge props)
 * @tparam Edge The edge type used to represent edges in the ClusterGraph
 * @param C The cluster graph
 * @param measureType Indicates which type of measure can be found on the edges of C, also indicates whether these should be maximised or minimised.
 */
template <typename ClusterGraph, typename Vertex = boost::graph_traits<ClusterGraph>::vertex_descriptor>
std::pair<Vertex, Vertex> getNextTreePair(const ClusterGraph &C, bool measureType);


////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template<typename Graph, typename Vertex, typename Edge>
boost::numeric::ublas::matrix<float> kernelCalculation(const std::vector<Graph> &inputGraphs, const bool kernelType, const bool isSparse) {

    // Floyd-Warshall transformation inspiration can be found here (p. 4): https://www.dbs.ifi.lmu.de/~borgward/papers/BorKri05.pdf
    // Delta/Linear kernel calculation inspiration can be found here (sl. 28): https://ethz.ch/content/dam/ethz/special-interest/bsse/borgwardt-lab/documents/slides/CA10_GraphKernels_intro.pdf

    // All the floyd-transformation matrices in order of the input graphs.
    // (i, j) in matrix k is shortest path distance between vertex i and j in graph K.
    std::vector<std::vector<std::vector<int>>> floydGraphs;
    size_t n = inputGraphs.size();
    
    // Creating a weight map for each of the input graphs as we are not requiring that the edges have weights
    // We map each edge descriptor to 1, to illustrate that traversing any edge costs 1.
    std::vector<std::map<Edge, int>> weightMaps;
    for(auto graphIt = inputGraphs.begin(); graphIt != inputGraphs.end(); graphIt++) {
        std::map<Edge, int> weightMap;
        for(auto edge : asRange(edges(*graphIt))){
            weightMap.insert(std::make_pair(edge, 1));
        }
        weightMaps.push_back(weightMap);
    }


    // Calculating all pairs shortest paths for all input graphs
    // Computed distance matrices correspond to the Floyd-Warshall transformed graphs.
    bool success;
    for(auto graphIt = inputGraphs.begin(); graphIt != inputGraphs.end(); graphIt++) {
        std::vector<std::vector<int>> distanceMatrix(num_vertices(*graphIt), std::vector<int>(num_vertices(*graphIt), 0)); 

        // Retrieving the map for the current graph
        std::map<Edge, int> edgeMap = weightMaps[graphIt - inputGraphs.begin()];
        // Transforming the std::map to a property map, which is what the Johnson and Floyd Warshall algorithms require
        boost::associative_property_map<std::map<Edge, int>> propMap = boost::make_assoc_property_map(edgeMap);

        if(isSparse){
            success = boost::johnson_all_pairs_shortest_paths( *graphIt, distanceMatrix, boost::weight_map(propMap));// Johnson algorithm on sparse graphs, uses Dijkstras and Bellman Ford, O(V^2*logV)
        } else{
            success = boost::floyd_warshall_all_pairs_shortest_paths( *graphIt, distanceMatrix, boost::weight_map(propMap) ); // Floyd Warshall on dense graphs. O(V^3) regardless of edges.
        }
        
        if(!success){
            std::string msgString = "Negative cycle detected in input graph. We do not accommodate this.";
            const char *msg = msgString.c_str();
            throw FatalException(msg);
        } else {
            
            floydGraphs.push_back(distanceMatrix);
            
        }
    }

    // Compute the (n x n) kernel matrix by considering all pairs of graphs and considering, for those graphs, all pairs of pairs of vertices.
    boost::numeric::ublas::matrix<float> kernelMatrix (n, n);
    int kernelValue;
    int toAdd;
    for(size_t i = 0; i < inputGraphs.size(); i++){
        for(size_t j = i; j < inputGraphs.size(); j++){
            kernelValue = 0;
            Graph g1 = inputGraphs[i];
            Graph g2 = inputGraphs[j];
            for(Vertex u : asRange(vertices(g1))){ // double sum over (u, v) and (u', v')
                for(Vertex v : asRange(vertices(g1))){
                    for(Vertex uP : asRange(vertices(g2))){
                        for(Vertex vP : asRange(vertices(g2))){
                            if(kernelType == false){ // Linear
                                toAdd = floydGraphs[i][u][v] * floydGraphs[j][uP][vP];
                            } else { // Delta
                                toAdd = floydGraphs[i][u][v] == floydGraphs[j][uP][vP] ? 1 : 0;
                            }
                            kernelValue += toAdd;
                        }
                    }
                }
            }
            kernelMatrix (i, j) = kernelValue;
        }
    }

    return kernelMatrix;
} // kernelCalculation


template<typename Graph, typename Vertex>
std::queue<std::shared_ptr<AlignObj<Graph>>>
    getGuideTree(const std::vector<Graph> &inputGraphs, 
                const boost::numeric::ublas::matrix<float> &kernelScores,
                const bool measureType) {
    // List of alignment vertices in the order they are expected to be aligned.
    std::queue<std::shared_ptr<AlignObj<Graph>>> Q;
    
    struct CedgeProps {
        float measure;
    };
    struct CvertexProps {
        std::string tag;
    };

    // Creating the cluster graph C where each vertex represents a graph (by index, tag).
    typedef boost::adjacency_list<boost::vecS, boost::listS, boost::undirectedS, CvertexProps, CedgeProps> C;
    C clusterGraph;
    float measure;
    for(auto graphIt = inputGraphs.begin(); graphIt != inputGraphs.end(); graphIt++){
        auto v = add_vertex(clusterGraph);
        size_t index = graphIt - inputGraphs.begin();
        clusterGraph[v].tag = std::to_string(index); 
    }

    // Adding edges between all pairs (no loops) of vertices in the graph with a measure given by the measure function
    for(auto v1 : asRange(vertices(clusterGraph))){
        int i = stoi(clusterGraph[v1].tag);
        for(auto v2 : asRange(vertices(clusterGraph))){
            int j = stoi(clusterGraph[v2].tag);

            if(!edge(v1, v2, clusterGraph).second && i != j){ // edge returns true in second position if an edge exists between v1 and v2
                float pairValue = j >= i ? kernelScores(i, j) : kernelScores(j, i); // matrix is only half filled
                if(measureType) { // Distance
                    measure = sqrt( kernelScores(i, i) + kernelScores(j, j) - 2*pairValue);
                } else { // Similarity
                    measure = pairValue;
                }
                auto e = add_edge(v1, v2, clusterGraph);
                clusterGraph[e.first].measure = measure;
            }
        }
    }

    std::vector<std::shared_ptr<AlignObj<Graph>>> Tsearch; // All alignment nodes that have yet to be aligned (for shorter searches).
    // Initialise alignObj leaves, each leaf having one of the inputgraphs
    for(auto graphIt = inputGraphs.begin(); graphIt != inputGraphs.end(); graphIt++){
        const Graph *constGraph = &(*graphIt); // Create Const Graph* by extracting a const graph and returning a reference to it.
        Graph *nonConstGraph = const_cast<Graph*>(constGraph); // Cast Const Graph* to Graph*
        std::shared_ptr<AlignObj<Graph>> g = std::make_shared<AlignObj<Graph>>(*nonConstGraph); // Dereference the now non-const pointer and give it to the constructor. (All in all, const removal).
        int index = graphIt - inputGraphs.begin();
        std::vector<int> contained = { index };
        g->createDefaultMatchTable(1, num_vertices(*constGraph));
        for( Vertex v : asRange(vertices(g->alignGraph))) {
            g->matchTable[0][v] = v;
        }

        g->tag = std::to_string(index);
        g->contained = contained;

        Tsearch.push_back(g);
    }


    // Concatenating the most similar vertices in C, two by two
    float optimalScore;
    using edgeDesc = typename C::edge_descriptor; // edge to extract later
    edgeDesc candidate; // One of the edges (graph pairs) with the optimal score.
    while(num_vertices(clusterGraph) > 1){
        // TODO: using edgeDesc = typename boost::graph_traits<boost::adjacency_list<boost::vecS, boost::listS, boost::undirectedS, CvertexProps, CedgeProps>>::edge_descriptor; // edge to extract later
        if(measureType){
            // Minimize distance
            optimalScore = std::numeric_limits<float>::max();
            for( auto e : asRange(edges(clusterGraph))){
                if(clusterGraph[e].measure < optimalScore){
                    optimalScore = clusterGraph[e].measure;
                    candidate = e;
                }
            }
        } else {
            // Maximize similarity
            optimalScore = std::numeric_limits<float>::lowest();
            for( auto e : asRange(edges(clusterGraph))){
                if(clusterGraph[e].measure > optimalScore){
                    optimalScore = clusterGraph[e].measure;
                    candidate = e;
                }
            }
        }

        // Updating all scores between the newly concatenated node and all other vertices
        auto u = source(candidate, clusterGraph);
        auto v = target(candidate, clusterGraph);
        auto concatVertex = add_vertex(clusterGraph);
        clusterGraph[concatVertex].tag = std::string( "(" + clusterGraph[u].tag + "," + clusterGraph[v].tag + ")");
        for(auto w : asRange(vertices(clusterGraph))){
            if(w != u && w != v && w != concatVertex){ // compute average of weight on (w, u), (w, v) and append it to (w, (u, v))
                std::pair<edgeDesc, bool> e = add_edge(concatVertex, w, clusterGraph);
                std::pair<edgeDesc, bool> wu = edge(w, u, clusterGraph);
                std::pair<edgeDesc, bool> wv = edge(w, v, clusterGraph);
                float measureWU = clusterGraph[wu.first].measure; // accessing weight on edge from w to u
                float measureWV = clusterGraph[wv.first].measure;
                clusterGraph[e.first].measure = (measureWU + measureWV) / 2;
            }
        }
        

        std::shared_ptr<AlignObj<Graph>> left, right;

        // Finding the alignment vertices corresponding to u and v in the current guide tree.
        // These will be aligned to, later, fill the guide tree out in a bottom-up manner.
        for(const std::shared_ptr<AlignObj<Graph>> &ptr : Tsearch){
            if(ptr->tag == clusterGraph[u].tag){
                // *tIt means dereferencing the iterator over Tsearch, here we get a shared_ptr back
                left = ptr;
            } else if(ptr->tag == clusterGraph[v].tag){
                right = ptr;
            }
        }

        // remove edges connecting the vertices u and v before removing the vertices themselves.
        clear_vertex(u, clusterGraph);
        clear_vertex(v, clusterGraph);
        remove_vertex(u, clusterGraph);
        remove_vertex(v, clusterGraph);

        std::vector<int> mergeContained = left->contained;
        mergeContained.insert(mergeContained.end(), right->contained.begin(), right->contained.end());

        std::shared_ptr<AlignObj<Graph>> alignVertex = std::make_shared<AlignObj<Graph>>(left, right, mergeContained);
        alignVertex->tag = std::string( "(" + left->tag + "," + right->tag + ")");

        // Clean tSearch to avoid looking at vertices that have already been aligned. Here it is left and right.
        for(auto tSearchIt = Tsearch.begin(); tSearchIt != Tsearch.end(); tSearchIt++){
            if((*tSearchIt)->tag == left->tag){
                Tsearch.erase(tSearchIt);
                break;
            }
        }
        for(auto tSearchIt = Tsearch.begin(); tSearchIt != Tsearch.end(); tSearchIt++){
            if((*tSearchIt)->tag == right->tag){
                Tsearch.erase(tSearchIt);
                break;
            }
        }
        
        Tsearch.push_back(alignVertex);
        Q.push(alignVertex);
    }

    return Q;
} // getGuideTree

template<typename Graph>
float getSimDistMeasure(const AlignObj<Graph> &g1, const AlignObj<Graph> &g2, const bool measure, const bool kernelType, const bool isSparse){
    std::vector<Graph> gs = {g1.alignGraph, g2.alignGraph};
    boost::numeric::ublas::matrix<float> kernelScores = kernelCalculation(gs, kernelType);
    // std::cout << "Graph i score: "  << kernelScores(0,0) << "\tGraph j score: " << kernelScores(1, 1) << "\t Pair score: " << kernelScores(0, 1) << "\n";
    if(measure){ // Distance
        return sqrt( kernelScores(0, 0) + kernelScores(1, 1) - 2*kernelScores(0, 1));
    } else { // Similarity
        return kernelScores(0, 1); 
    }
}

template<typename Graph, typename Vertex>
std::queue<std::shared_ptr<AlignObj<Graph>>> 
    getGuideTreeFromString(const std::vector<Graph> &inputGraphs, 
                const std::string &treeString){
    
    std::queue<std::shared_ptr<AlignObj<Graph>>> Q;
    _getGuideTreeFromStringAux(inputGraphs, treeString, Q);
    return Q;

} // getGuideTreeFromString

template<typename Graph, typename Vertex>
std::shared_ptr<AlignObj<Graph>> _getGuideTreeFromStringAux(const std::vector<Graph> &inputGraphs, const std::string& treeString, std::queue<std::shared_ptr<AlignObj<Graph>>> &Q){

    std::shared_ptr<AlignObj<Graph>> node;
    std::string leftString, rightString;

    if(std::count(treeString.begin(), treeString.end(), ',') == 0){
        // Node is a leaf, treestring = "x"
        int graphIndex = std::stoi(treeString);
        Graph g = inputGraphs[graphIndex];
        node = std::make_shared<AlignObj<Graph>>(g);
        // Create trivial match table for each leaf node
        
        node->createDefaultMatchTable(1, num_vertices(g));
        
        for( Vertex v : GraphAlign::asRange(vertices(g))){
            node->matchTable[0][v] = v;
        }
        node->contained = {graphIndex};
        node->tag = treeString;

        return node;
    }
    // Node is an internal node
    else if (std::count(treeString.begin(), treeString.end(), ',') == 1){
        // treestring = "(x,y)"
        auto sepPos = std::find(treeString.begin(), treeString.end(), ',');
        for(auto it = treeString.begin() + 1; it != sepPos; it++){
            leftString = leftString + *it;
        }
        for(auto it = sepPos + 1; it != treeString.end() - 1; it++){
            rightString = rightString + *it;
        }
    } else {
        std::pair<std::string, std::string> twohalves = _splitGuideTreeString(treeString);
        leftString = twohalves.first;
        rightString = twohalves.second;
    }

    std::shared_ptr<AlignObj<Graph>> leftChild = _getGuideTreeFromStringAux(inputGraphs, leftString, Q);
    std::shared_ptr<AlignObj<Graph>> rightChild = _getGuideTreeFromStringAux(inputGraphs, rightString, Q);

    std::vector<int> mergeContained = leftChild->contained;
    mergeContained.insert(mergeContained.end(), rightChild->contained.begin(), rightChild->contained.end());
    node = std::make_shared<GraphAlign::AlignObj<Graph>>(leftChild, rightChild, mergeContained);
    node->tag = std::string("(" + leftChild->tag + "," + rightChild->tag + ")");
    Q.push(node);
    return node;

} // _getGuideTreeFromStringAux

std::pair<std::string, std::string> _splitGuideTreeString(const std::string &s){
    std::string leftString, rightString;
    int stack = 0;
    int splitPos = 0;

    for(auto it = s.begin(); it != s.end(); it++){
        if (*it == '('){
            stack++;
        } else if (*it == ')'){
            stack--;
        }
        else if (*it == ',' && stack == 1){
            break;
        }
        splitPos++;
    }
    // splitPos holds the splitting position
    leftString = s.substr(1, splitPos - 1);
    rightString = s.substr(splitPos + 1, s.size() - splitPos + 1 - 3); // copy the rest excluding the string "(,)"

    return {leftString, rightString};
}

template <typename Graph>
int getKernelScore(const Graph &g1, const Graph &g2, const std::vector<std::vector<int>> &g1Dist, const std::vector<std::vector<int>> &g2Dist, const std::string &kernelMeasure){
    int kernelValue = 0;
    int toAdd = 0;
    for(auto u : asRange(vertices(g1))){ // double sum over (u, v) and (u', v')
        for(auto v : asRange(vertices(g1))){
            for(auto uP : asRange(vertices(g2))){
                for(auto vP : asRange(vertices(g2))){
                    if(kernelMeasure == "linear"){ // Linear
                        toAdd = g1Dist[u][v] * g2Dist[uP][vP];
                    } else if (kernelMeasure == "delta") { // Delta
                        toAdd = g1Dist[u][v] == g2Dist[uP][vP] ? 1 : 0;
                    }
                    kernelValue += toAdd;
                }
            }
        }
    }
    return kernelValue;
}

template <typename Graph, typename Edge> 
std::vector<std::vector<std::vector<int>>> getDistanceMatrices(const std::vector<Graph> &gs){
    std::vector<std::vector<std::vector<int>>> floydMatrices;
    floydMatrices.reserve(gs.size());

    for(auto it = gs.begin(); it != gs.end(); it++){
        std::map<Edge, int> weightMap;
        for(auto edge : asRange(edges(*it))){
            weightMap.insert(std::make_pair(edge, 1));
        }
        boost::associative_property_map<std::map<Edge, int>> propMap = boost::make_assoc_property_map(weightMap);
        std::vector<std::vector<int>> distanceMatrix(num_vertices(*it), std::vector<int>(num_vertices(*it), 0));
        // Only attempt this for SPARSE graphs!
        bool success = boost::johnson_all_pairs_shortest_paths( *it, distanceMatrix, boost::weight_map(propMap)); // Johnson algorithm on sparse graphs, uses Dijkstras and Bellman Ford, O(V^2*logV)

        if(!success){
            std::string msgString = "Negative cycle detected in input graph. We do not accommodate this.";
            const char *msg = msgString.c_str();
            throw FatalException(msg);
        } else {
            floydMatrices.push_back(distanceMatrix);
        }
    }

    return floydMatrices;
} // getDistanceMatrices

template <typename Graph, typename ClusterGraph>
void initialiseClusterGraph(ClusterGraph &C, const std::vector<Graph> &inputGraphs, const std::vector<std::vector<std::vector<int>>> &floydMatrices, bool measureType){
    for(auto it = inputGraphs.begin(); it != inputGraphs.end(); it++){
        const Graph *constGraph = &(*it); // Create Const Graph* by extracting a const graph and returning a reference to it.
        Graph *nonConstGraph = const_cast<Graph*>(constGraph); // Cast Const Graph* to Graph*
        int index = it - inputGraphs.begin();
        auto emptyAlign = std::make_shared<AlignObj<Graph>>(*nonConstGraph, index);;
        
        // Create vertex for the cluster graph
        auto v = add_vertex(C);
        C[v].alignObj = emptyAlign;
        C[v].floydMatrix = std::make_shared<std::vector<std::vector<int>>>(floydMatrices[index]);
        // compute the kernel(i,i) value and associate it with the node
        C[v].ownDistance = getKernelScore(*it, *it, floydMatrices[index], floydMatrices[index], "delta");
    }

    // Add edges and the similarity measure to the cluster graph
    float measure = 0;
    float pairValue = 0;
    for(const auto u : asRange(vertices(C))){
        for(const auto v : asRange(vertices(C))){
            if(!edge(u, v, C).second && u != v){
                // Compute the kernelscore for these two graphs only
                pairValue = getKernelScore(C[u].alignObj->alignGraph, 
                                                   C[v].alignObj->alignGraph, 
                                                   *C[u].floydMatrix, 
                                                   *C[v].floydMatrix, 
                                                   "delta");
                if(measureType){ // distance
                    measure = sqrt( C[u].ownDistance + C[v].ownDistance - 2*pairValue);
                } else { // similarity
                    measure = pairValue;
                }
                auto e = add_edge(u, v, C);
                C[e.first].measure = measure;
            }
        }
    }
} // initialiseClusterGraph

template <typename ClusterGraph, typename Vertex>
std::pair<Vertex, Vertex> getNextTreePair(const ClusterGraph &C, bool measureType){
    using cEdge = typename ClusterGraph::edge_descriptor;
    cEdge candidate;
    float optimalScore = 0;
    if(measureType){ // Distance
        optimalScore = std::numeric_limits<float>::max();
    } else { // similarity
        optimalScore = std::numeric_limits<float>::min();
    }
    for(const auto e : asRange(edges(C))){
        if(measureType){
            if(C[e].measure < optimalScore){
                optimalScore = C[e].measure;
                candidate = e;
            }
        } else {
            if(C[e].measure > optimalScore){
                optimalScore = C[e].measure;
                candidate = e;
            }
        }
    }
    return std::make_pair(source(candidate, C), target(candidate, C));
} // getNextTreePair


} // namespace GraphAlign

#endif 