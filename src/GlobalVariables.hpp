#ifndef GLOBALVARS_H
#define GLOBALVARS_H

#include <vector>
#include <map>
#include <string>

namespace GraphAlign {

template<typename Vertex>
struct GlobalVariables{ 

    GlobalVariables(std::vector<float> labelScore) : labelScore(labelScore) {}
    GlobalVariables(std::vector<float> labelScore, int nLabels) : labelScore(labelScore), nLabels(nLabels) {}

    GlobalVariables(const GlobalVariables &other) : labelScore(other.labelScore), nLabels(other.nLabels) {}

    ~GlobalVariables() = default;

public:
    std::vector<Vertex>& getLargerMatch() {
        if(a1Match.size() >= a2Match.size()){
            return a1Match;
        }
        else {
            return a2Match;
        }
    }

    std::vector<Vertex>& getSmallerMatch() {
        if(a1Match.size() >= a2Match.size()) {
            return a2Match;
        }
        else {
            return a1Match;
        }
    }

public:
    std::vector<float> labelScore; // labelScore[i] = x if matching two vertices with identical labels 'i' yields a score of 'x'
    int nLabels = 0; // total number of labels across the graph sets
    float bestScore = 0; // the global best score achieved
    std::vector<Vertex> a1Match; // the current match vector for a1
    std::vector<Vertex> a2Match; // the current match vector for a2
    std::size_t matchSize = 0; // the number of non nullvertex entries in a1Match/a2Match

    std::vector<std::pair<Vertex,Vertex>> anchor; // the anchor vertex pairs amongst a1 and a2
    std::vector<Vertex> smallerAnchor; // the anchor vertices of the smaller graph of a1 and a2

    float initialScore = 0; // Initial score, meaning the score achieved from the anchor mapping, used as base value for hypothetical score.
    
    // Following vectors are used in branch and bound checking for expand algorithms, as they can be maintained as vertices are matched
    // instead of recreating them in each step.
    std::vector<unsigned char> notAllowedVerticesG1; // Vector indicating vertices that cannot be mapped in g1
    std::vector<unsigned char> notAllowedVerticesG2; // Vector indicating vertices that cannot be mapped in g2

    // Label matrices for branch and bound check when scoring by scheme, they are always assigned such that the larger
    // graph corresponds to g1LabelMatrix and the smaller to g2LabelMatrix.
    // The first vector is indexed into with a label id and the second is accessed by a number between 0 and size of contained for 
    // corresponding graph. 
    // The entry then states how many vertices with the given label appears in how many of the input graphs.
    // g1LabelMatrix[i][j] tells how many vertices in g1 that have label i appears in j input graphs.
    // Position 0 indicates how many vertices in the graph have the given label. 
    std::vector<std::vector<int>> g1LabelMatrix;
    std::vector<std::vector<int>> g2LabelMatrix;

    // Vectors indexed into with vertices from the alignment graphs returning how many input graphs
    // said vertices appear in.
    std::vector<int> nInputGraphsForVertexG1;
    std::vector<int> nInputGraphsForVertexG2;

    // label frequency vector for the graphs, used in trimming algorithm when scoring by order
    // and in expand in isMaxColoredMatch to check if the label sets are disjoint or not.
    // In same manner as the graph label matricces, the labelFrequencyG1 is always for the larger graph
    std::vector<int> labelFrequencyG1;
    std::vector<int> labelFrequencyG2;

    // For dynamic order computation
    std::vector<Vertex> degBags; // The vector containing vertices distributed into bags, each bag representing a deg-to-mapped class.
    std::vector<int> degBagEnds; // The vector containing information about the .end() position of the i'th degBag.
    std::vector<int> whichDegBag; // The vector indicating which degree bag a given vertex belongs to
    std::vector<int> bagPos; // The vector indicating exactly where in degBags a specific vertex is positioned.

    // For least-degree-distance first 
    std::vector<int> nUnmatchedNeighs; // Number of unmatched neighbours in the neighbour of a given vertex.
};

}

#endif // GLOBALVARS_H