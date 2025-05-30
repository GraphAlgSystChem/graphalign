#ifndef BUILDALIGNMENT_H
#define BUILDALIGNMENT_H

namespace GraphAlign {

/**
 * @brief Building the new alignment based on the two previous alignments and the matched vertices. 
 *        Updates parent accordingly.
 * 
 * @tparam Graph The graph type of the input graphs
 * @tparam SV The type of the function for scoring vertices
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @tparam Vertex The vertex type
 * @tparam Edge The edge type
 * @param parent The parent alignment being built based on the match between a1 and a2 
 * @param a1 The first alignment graph 
 * @param a2 The second alignment graph
 * @param a1Match Match vector for a1
 * @param a2Match Match vector for a2
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 */
template<typename Graph, class SV, 
         typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, 
         typename Edge = boost::graph_traits<Graph>::edge_descriptor>
void buildAlignment(AlignObj<Graph> &parent, 
                    const AlignObj<Graph> &a1, 
                    const AlignObj<Graph> &a2, 
                    std::vector<Vertex> &a1Match, 
                    std::vector<Vertex> &a2Match, 
                    const AlignParameters<Graph, SV> &paramList);


////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////

template<typename Graph, class SV, typename Vertex, typename Edge>
void buildAlignment(AlignObj<Graph> &parent, 
                        const AlignObj<Graph> &a1, 
                        const AlignObj<Graph> &a2, 
                        std::vector<Vertex> &a1Match, 
                        std::vector<Vertex> &a2Match, 
                        const AlignParameters<Graph, SV> &paramList) {
    
    Vertex nullVertex =  boost::graph_traits<Graph>::null_vertex();
    int matchSize = 0;
    for(Vertex x: asRange(vertices(a1.alignGraph))) {
        if(a1Match[getIndex(x)] != nullVertex) matchSize++;
    }
    Graph parentGraph;
    std::vector<std::pair<Vertex, Vertex>> parentAmbiguous;
    // Mapping from V(a1) -> V(alignGraph), V(a2) -> V(alignGraph)
    std::vector<Vertex> m1(num_vertices(a1.alignGraph), nullVertex); 
    std::vector<Vertex> m2(num_vertices(a2.alignGraph), nullVertex);

    size_t totalRows = a1.contained.size() + a2.contained.size(); // = totalGraphs
    size_t totalColumns = num_vertices(a1.alignGraph) + num_vertices(a2.alignGraph) - matchSize;

    parent.createDefaultMatchTable(totalRows, totalColumns);
    

    // Stacking a1 column onto a2 column. The a1 rows are first, followed by a2 rows.
    // For the match vertices.
    for(Vertex u : asRange(vertices(a1.alignGraph))) {
        Vertex v = a1Match[getIndex(u)];
        if (v == nullVertex) {
            continue;
        }
        Vertex x = add_vertex(parentGraph); // labels must exist in this case, otherwise boom
        if(paramList.enforceVLabels) {
            parentGraph[x].label = a1.alignGraph[u].label;
            parentGraph[x].labelIndex = a1.alignGraph[u].labelIndex;
        }
        m1[getIndex(u)] = x;
        m2[getIndex(v)] = x;
        for(size_t i = 0; i < totalRows; i++) {
            if(i < a1.contained.size()) {
                parent.matchTable[i][x] = a1.matchTable[i][u];
            } 
            // moved on to a2's rows
            else {
                size_t remappedIndex = i - a1.contained.size(); // offset into a2's contained vector. a1.contained.size() vertices have already been considered.
                parent.matchTable[i][x] = a2.matchTable[remappedIndex][v];
            }
        }
    }

    // Stacking a1 column onto "empty" a2 column.
    // Vertices that are in a1 but not in a2.
    for(Vertex u : asRange(vertices(a1.alignGraph))) {
        if (a1Match[getIndex(u)] != nullVertex) {
            continue;
        }
        else {
            Vertex x = add_vertex(parentGraph);
            m1[getIndex(u)] = x;
            for(size_t i = 0; i < totalRows; i++){
                if(i < a1.contained.size()){
                    parent.matchTable[i][x] = a1.matchTable[i][u];
                } else {
                    parent.matchTable[i][x] = nullVertex;
                }
            }
            if(paramList.enforceVLabels) {
                parentGraph[x].label = a1.alignGraph[u].label;
                parentGraph[x].labelIndex = a1.alignGraph[u].labelIndex;
            }
            // rows below are already initialised to null vertices
        }
    }

    // Stacking "empty" a1 column on top of a2 column.
    // For vertices in a2 that are not in a1.
    for(Vertex v: asRange(vertices(a2.alignGraph))) {
        if(a2Match[getIndex(v)] != nullVertex) {
            continue;
        }
        else {
            Vertex x = add_vertex(parentGraph);
            m2[getIndex(v)] = x;
            for(size_t i = 0; i < totalRows; i++){
                if(i < a1.contained.size()){
                    parent.matchTable[i][x] = nullVertex;
                } else {
                    size_t remappedIndex = i - a1.contained.size();
                    parent.matchTable[i][x] = a2.matchTable[remappedIndex][v];
                }
            }
            if(paramList.enforceVLabels) {
                parentGraph[x].label = a2.alignGraph[v].label;
                parentGraph[x].labelIndex = a2.alignGraph[v].labelIndex;
            }
        }
    }


    // Adding edges from a1 and a2 to the alignment graph.
    for(auto e : asRange(edges(a1.alignGraph))) {
        Vertex u = source(e, a1.alignGraph);
        Vertex v = target(e, a1.alignGraph);
        Vertex uPrime = m1[getIndex(u)];
        Vertex vPrime = m1[getIndex(v)];
        Edge eP = add_edge(uPrime, vPrime, parentGraph).first;
        if(paramList.enforceELabels) {
            parentGraph[eP].label = a1.alignGraph[e].label;
        }
    }
    for(auto e : asRange(edges(a2.alignGraph))) {
        Vertex u = source(e, a2.alignGraph);
        Vertex v = target(e, a2.alignGraph);
        Vertex uPrime = m2[getIndex(u)];
        Vertex vPrime = m2[getIndex(v)];
        // Avoid duplicate edges between matched vertices.
        if(!edge(uPrime, vPrime, parentGraph).second) {
            Edge eP = add_edge(uPrime, vPrime, parentGraph).first;
            if(paramList.enforceELabels) {
                parentGraph[eP].label = a2.alignGraph[e].label;
            }
        }
    }

    // Compute set of ambiguous edges
    if(paramList.useAmbiguousEdges) {
        bool ambiguous;
        std::vector<Vertex> parentVertices;
        for(Vertex v : asRange(vertices(parentGraph))) {
            parentVertices.push_back(v);
        }
        std::vector<std::vector<Vertex>> vertexPairs = GraphAlignUtility::getCombinations(parentVertices, parentVertices.size(), 2);
        for(const std::vector<Vertex> &pair : vertexPairs) {
            Vertex u = pair[0];
            Vertex v = pair[1];
            ambiguous = true;
            for(size_t i = 0; i < totalRows; i++) {
                Vertex existenceU = parent.matchTable[i][u];
                Vertex existenceV = parent.matchTable[i][v];
                if(existenceU != nullVertex && existenceV != nullVertex) {
                    ambiguous = false;
                    break;
                }
            }
            if(ambiguous) {
                parentAmbiguous.push_back({u, v});
            }
        }
    }

    float score = matchScore(a1, a2, a1Match, a2Match, paramList);
    parent.alignGraph = parentGraph;
    parent.ambiguousEdges = parentAmbiguous;
    parent.alignScore = score;
    std::cout << "\tAlignment of " << a1.tag << " and " << a2.tag << " with a match score of " << score;

} // buildAlignment

} // namespace GraphAlign

#endif // BUILDALIGNMENT_H