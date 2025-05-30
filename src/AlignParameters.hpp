#ifndef ALIGNPARAMS_H
#define ALIGNPARAMS_H

namespace GraphAlign {


// This file includes the AlignParameters struct and the corresponding enum classes

enum class Algorithm {TRIM, REC};
enum class Scoring {ORDER, SCHEME};

/**
 * @brief Struct to store the input parameters, a const reference to this struct 
 * can then be passed around instead of passing all the input parameters around.
 * 
 * @tparam Graph Need the graph type for storign the input graphs
 * @tparam SV The type of the function that scores vertices
 */
template <typename Graph, class SV>
struct AlignParameters{ 

    AlignParameters(Algorithm algorithm, Scoring scoring,
                    bool enforceV, bool enforceE, 
                    bool useAmbiguous, bool oneMatch, 
                    const std::vector<Graph> &graphs, 
                    SV scoreVertices, 
                    bool measureType) : algorithm(algorithm), scoring(scoring), 
                                                enforceVLabels(enforceV), enforceELabels(enforceE), 
                                                useAmbiguousEdges(useAmbiguous), oneMatch(oneMatch), 
                                                inputGraphs(graphs), scoreVertices(scoreVertices), 
                                                measureType(measureType) {}
    // Copy constructor
    AlignParameters(const AlignParameters &other) : algorithm(other.algorithm), scoring(other.scoring), 
                                                    enforceVLabels(other.enforceVLabels), enforceELabels(other.enforceELabels),
                                                    useAmbiguousEdges(other.useAmbiguousEdges), oneMatch(other.oneMatch), 
                                                    inputGraphs(other.inputGraphs), scoreVertices(other.scoreVertices) {}

    ~AlignParameters() = default;

public:
    /*
    * Storing the input parameters:
    * algorithm: Enum specifying whether the pairwise alignment algorithm used is "trimming" (TRIM) or "expand" (REC).
    * scoring: Enum specifying the scoring method, either by number of vertices (ORDER) or label scoring on vertices (SCHEME)
    * enforceVLabels: Boolean to determine whether vertices are only allowed to be mapped if they have the same label.
    * enforceELabels:Boolean to determine whether edges are only allowed to be mapped if they have the same label.
    * useAmbiguousEdges: Boolean specifying whether the usage of ambiguous edges is allowed.
    * oneMatch: Boolean stating if only one match from TRIM / REC is wanted or all with the same score.
    * inputGraphs: Vector storing the input graphs
    * scoreVertices: The function to score vertices
    * measureType: true (distance) and false (similarity) w.r.t. guide tree calculations
    */
    Algorithm algorithm;
    Scoring scoring;
    bool enforceVLabels;
    bool enforceELabels;
    bool useAmbiguousEdges;
    bool oneMatch;
    const std::vector<Graph> &inputGraphs;
    const SV scoreVertices;
    bool measureType;
};

}

#endif // ALIGNPARAMS_H