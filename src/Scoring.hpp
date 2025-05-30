#ifndef SCORING_H
#define SCORING_H


namespace GraphAlign {

// This file include functions:
// matchScore
// scoreColumn
// scoreNewMatch


/**
 * @brief Computes the score for the given match for the two alignment graphs given as input (Does not score edges currently)
 *        (Efficient implementation that uses scoreNewMatch as subroutine for scoring a a pair of vertices, while 
 *         this function only responsibility is finding the pairs of vertices to score)
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type

 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param a1Match The match vector for a1's underlying graph, g1, namely a1Match[i] = j if v_i in G1 is mapped to v_j in G2.
 * @param a2Match The match vector for a2's underlying graph, g2, namely a1Match[i] = j if v_i in G2 is mapped to v_j in G1.
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return The score computed for the alignment achieved by the matched vertices in the two alignment graphs.
 */
template<typename Graph, class SV, 
    typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
float matchScore(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const std::vector<Vertex> &a1Match,
                const std::vector<Vertex> &a2Match,
                const AlignParameters<Graph, SV> &paramList);

/**
 * @brief Computes the score from matching two vertices, i.e. computing the sum-of-pairs score in the entire matchcolum
 *        created by the combining the entries for u and v in the match table from a1 and a2 
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param u Vertex from graph in a1
 * @param v Vertex from graph in a2
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return The score of matching u and v
 */
template<typename Graph, class SV, 
    typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
float scoreColumn(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const Vertex u, const Vertex v,
                const AlignParameters<Graph, SV> &paramList);


/**
 * @brief Computes the additional score achieved by matching two vertices (specifically for recursive expand), i.e. computing the sum-of-pairs score across the
 *        column (and not among vertices belonging to the same alignment graph) created by the combining the entries for u and v in the match table from a1 and a2 
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param u Vertex from graph in a1
 * @param v Vertex from graph in a2
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return The score of matching u and v
 */               
template<typename Graph, class SV,
    typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
float scoreNewMatch(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const Vertex u, const Vertex v,
                const AlignParameters<Graph, SV> &paramList);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template<typename Graph, class SV, typename Vertex, typename Edge>
float matchScore(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const std::vector<Vertex> &a1Match,
                const std::vector<Vertex> &a2Match,
                const AlignParameters<Graph, SV> &paramList) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    float currentScore = 0;
    if(paramList.scoring == Scoring::ORDER) {
        // checking which of the two graphs are smaller, as it is therefore quicker to 
        // go through the vertices of that graph and check if they are matched, rather than the larger of the two graphs
        if(num_vertices(a1.alignGraph) >= num_vertices(a2.alignGraph)) {
            for(std::size_t i = 0; i < a2Match.size(); i++) {
                if(a2Match[i] != nullVertex) {
                    currentScore++;
                }   
            }
        } else {
            for(std::size_t i = 0; i < a1Match.size(); i++) {
                if(a1Match[i] != nullVertex) {
                    currentScore++;
                }   
            }
        }
    }
    else { //paramList.scoring == Scoring::SCHEME

        for(Vertex x: asRange(vertices(a1.alignGraph))) {
            if(a1Match[getIndex(x)] != nullVertex) {
                currentScore = currentScore + scoreColumn(a1, a2, x, a1Match[getIndex(x)], paramList);
            }
            // Score a1 columns inserted onto empty a2 columns.
            else {
                currentScore = currentScore + scoreColumn(a1, a1, x, x, paramList);
            }
        }
        for(Vertex x: asRange(vertices(a2.alignGraph))) {
            // Score a2 columns inserted onto empty a1 columns.
            if(a2Match[getIndex(x)] == nullVertex) {
                currentScore = currentScore + scoreColumn(a2, a2, x, x, paramList);
            }
        }
    }
    return currentScore;
} // matchScore


template<typename Graph, class SV, typename Vertex, typename Edge>
float scoreColumn(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const Vertex u,
                const Vertex v,
                const AlignParameters<Graph, SV> &paramList) {
    float score = 0;
    if(paramList.scoring == Scoring::ORDER) {
        score++;
    }
    else {
        // Calculate the score of matching u to v, from g1 and g2
        Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
        std::vector<int> combinedContained = a1.contained; // The set of all graphs in this alignment.
        if(&a1 != &a2) { // If a1 == a2 refer to the same object, we are scoring a single column in an alignment graph.
            combinedContained.insert(combinedContained.end(), a2.contained.begin(), a2.contained.end());    
        }
        std::vector<Vertex> vertexMap(combinedContained.size(), nullVertex);

        // Find the vertices for each of the input graphs that a1 and a2 are created from
        for(std::size_t i = 0; i < a1.contained.size(); i++) {
            vertexMap[i] = a1.matchTable[i][u];
        }
        for(std::size_t j = a1.contained.size(); j < combinedContained.size(); j++) {
            vertexMap[j] = a2.matchTable[j - a1.contained.size()][v];
        }
        
        for(std::size_t i = 0; i < combinedContained.size(); i++) {
            for(std::size_t j = i + 1; j < combinedContained.size(); j++) {
                score = score + std::invoke(paramList.scoreVertices, vertexMap[i], vertexMap[j], 
                                paramList.inputGraphs[combinedContained[i]], paramList.inputGraphs[combinedContained[j]]);
            }
        }
    }

    return score;

} // scoreColumn

template<typename Graph, class SV, typename Vertex, typename Edge>
float scoreNewMatch(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2, 
                const Vertex u, const Vertex v,
                const AlignParameters<Graph, SV> &paramList) {
    float score = 0;
    if(paramList.scoring == Scoring::ORDER) {
        score++;
    }
    else {        
        for(std::size_t i = 0; i < a1.contained.size(); i++) {
            for(std::size_t j = 0; j < a2.contained.size(); j++) {
                score = score + std::invoke(paramList.scoreVertices, a1.matchTable[i][u], a2.matchTable[j][v],
                                paramList.inputGraphs[a1.contained[i]], paramList.inputGraphs[a2.contained[j]]);
            }
        }
    }

    return score;

} // scoreNewMatch


} // namespace GraphAlign

#endif // SCORING_H