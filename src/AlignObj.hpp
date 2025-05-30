#ifndef ALIGNOBJ_H
#define ALIGNOBJ_H

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace GraphAlign {

template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
struct AlignObj {

public:
    AlignObj(Graph g) : alignGraph(g) {}

    // Constructor for a trivial alignment object representing a single graph.
    AlignObj(Graph g, int index) : alignGraph(g) {
        createDefaultMatchTable(1, num_vertices(g));
        contained = { index };
        for(Vertex v : asRange(vertices(g))){
            this->matchTable[0][v] = v;
        }
        tag = std::to_string(index);
    }
    AlignObj(std::shared_ptr<AlignObj> left, std::shared_ptr<AlignObj> right, std::vector<int> contained) : 
        left(left), right(right), contained(contained) {}
    AlignObj(Graph g, std::shared_ptr<AlignObj> left, std::shared_ptr<AlignObj> right, 
             std::vector<int> contained, std::string tag, std::vector<std::vector<Vertex>> matchTable, 
             std::vector<std::pair<Vertex,Vertex>> ambiguousEdges, float alignScore) :
        alignGraph(g), left(left), right(right), contained(contained), tag(tag), matchTable(matchTable), 
        ambiguousEdges(ambiguousEdges), alignScore(alignScore) {}
    
    // Destructor
    ~AlignObj() = default;

    // Copy constructor
    AlignObj(const AlignObj &other) {
        this->alignGraph = other.alignGraph;
        this->left = other.left;
        this->right = other.right;
        this->contained = other.contained;
        this->tag = other.tag;
        this->matchTable = other.matchTable;
        this->ambiguousEdges = other.ambiguousEdges;
        this->alignScore = other.alignScore;
    }
    
    // Copy assignment
    AlignObj &operator=(const AlignObj &other) {
        this->clear();
        AlignObj(other.alignGraph, other.left, other.right, other.contained,
                 other.tag, other.matchTable, other.ambiguousEdges, other.alignScore);
        return *this;
    }

    // Clears the alignment object
    void clear() {
        this->alignGraph = Graph();
        this->tag = "";
        this->matchTable = {};
        this->left = nullptr;
        this->right = nullptr;
        this->contained = {};
        this->ambiguousEdges = {};
        this->alignScore = 0;
    }

    // Creating a default match table of the right size.
    void createDefaultMatchTable(int nRows, int nColumns) {
        Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
        this->matchTable.reserve(nRows);
        for(int i = 0; i < nRows; i++) {
            std::vector<Vertex> emptyRow(nColumns, nullVertex);
            this->matchTable.push_back(emptyRow);
        }

    }

/*
*   alignGraph:     The underlying graph structure for representing this alignment graph.
*   tag:            A tag (x, y) stating that graph "x" and graph "y" have been aligned to create this alignment object. 
                    Used in getGuideTree.
*   matchTable:     A map from G_i -> u -> v, where 'i' in contained, 'u' in vertices(alignGraph) and 'v' is POTENTIALLY in vertices(G_i).
                    If matchTable[i][u] returns null vertex, it means that vertex 'u' in the alignment graph does not correspond to a vertex in
                    G_i.
                    Row i in the matchTable corresponds to the graph G_contained[i].
*   left:           The 'left' alignment object that was (or is) used to create this alignment graph.
*   right:          -||- 'right' -||-
*   contained:      The list of input graphs (by index in inputGraphs) that are included in this alignment graph.
*   ambiguousEdges: The set of ambiguous edges in the alignment graph.
*   alignScore:     The alignment score of this alignObj's underlying alignment.
*/
public:
    Graph alignGraph;
    std::string tag; 
    std::vector<std::vector<Vertex>> matchTable = {};
    std::shared_ptr<AlignObj<Graph>> left = nullptr; 
    std::shared_ptr<AlignObj<Graph>> right = nullptr; 
    std::vector<int> contained = {}; 
    std::vector<std::pair<Vertex, Vertex>> ambiguousEdges = {};
    float alignScore = 0;
}; // AlignObj


}


#endif // ALIGNOBJ_H