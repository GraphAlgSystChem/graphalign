#ifndef BRANCHBOUND_H
#define BRANCHBOUND_H



namespace GraphAlign {

// This file includes the following:
// branchBoundCheck (trim)
// branchBoundCheck (expand)
// calculateHypotheticalScore
// evaluateKLevel
// prepareGlobalVars

/**
 * @brief Function to handle branch and bound check, called in trimming, different actions depending on order or scheme.
 *        Scheme: Update the label matrix for the smaller graph by decrementing values corresponding to the vertices in the subset
 *        and check both the necessary condition of enough vertices with each label in larger graph and the hypothetical score. 
 *        Order: The label frequencies for the two graphs were computed in prepareGlobalVars, so simply decrease the label frequencies
 *        for the vertices in the subset (i.e. the removed vertices) and check if there are enough vertices of each label 
 *        in g1 to cover all the ones in g2.
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the scoring function for scoring vertices
 * @tparam Vertex The type of the vertices
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param subset The subset of vertices that have been removed  (trimmed)
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars The object representing global function parameters. 
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool branchBoundCheck(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const std::vector<Vertex> &subset,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars);

/**
 * @brief Branch bound check for expand, calls calculateHypotheticalScore for the hypothetical score and adds to current score.
 *        If this value is smaller than or equal to the current best score return false otherwise return true. 
 *        The label matrices used for calculating hypothetical score should be updated before calling this function.
 *        Only relevant when scoring by scheme. 
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the scoring function for scoring vertices
 * @tparam Vertex The type of the vertices
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param currentScore The score achieved by the match so far
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars The object representing global function parameters. 
 * @return true If the the expand algorithm should continue to explore this path in the search tree
 * @return false If the expand algorithm should backtrack due to score smaller than or equal to best score
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool branchBoundCheck(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    float currentScore,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars);


/**
 * @brief Calculate the hypothetical score for the available vertices based on the label matrices.
 *        Score is calculated by "stacking" match columns onto one another, prioritising stacking columns with
 *        the most entries first. This provides the best possible score. 
 * 
 * @tparam Graph Graph type
 * @tparam Vertex The type of the vertices
 * @param a1 The first alignment object, representing the larger alignment graph
 * @param a2 The second alignment object, representing the smaller alignment graph
 * @param globalVars The object representing global function parameters
 * @return float The hypothetical score for the available vertices
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
float calculateHypotheticalScore(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                                GlobalVariables<Vertex> &globalVars);

/**
 * @brief Calculate the largest hypothetical score achieveable by removing k vertices from the smaller graph
 *        by saving the k smallest values added to the hypothetical and then subtracting them afterwards.
 *        This gives the largest possible hypothetical for all subsets when removing k vertices and as such
 *        if this hypothetical score is smaller than the best score it is not possible for any of the subsets
 *        created by removing k vertices to give a better score. 
 *        Therefore, stop execution of this k level and since the hypothetical score calculated at the next level 
 *        only gets worse the program can stop completely and return the best found match.
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the scoring function for scoring vertices
 * @tparam Vertex The type of the vertices
 * @param a1 The first alignment object
 * @param a2 The second alignment object
 * @param k The k level indicator
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars The object representing global function parameters. 
 * @return true If the program should continue, i.e. hypothetical score is greater than best score
 * @return false If the program should stop, i.e. hypothetical score is less than best score
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
bool evaluateKLevel(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    int k, const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars);

/**
 * @brief 
 * 
 * @tparam Graph Graph type
 * @tparam SV The type of the scoring function for scoring vertices
 * @param a1 The first alignment object, contains the larger graph
 * @param a2 The second alignment object, contains the smaller graph
 * @param g1LabelMatrix The label matrix for a1
 * @param g2LabelMatrix The label matrix for a2
 * @param k The k level indicator, telling how many values should be subtracted from hypothetical score
 * @param nLabels Number of labels possible
 * @param labelScore Vector providing the score achieved for matching two vertices with a given label for each label
 * @param initialScore Score provided by only the anchor mapping
 * @param bestScore The best score found so far
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @return true If the program should continue, i.e. hypothetical score is greater than best score
 * @return false If the program should stop, i.e. hypothetical score is less than best score
 */
template<typename Graph, class SV>
bool evaluateKLevel(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const std::vector<std::vector<int>> &g1LabelMatrix,
                    const std::vector<std::vector<int>> &g2LabelMatrix,
                    int k, 
                    int nLabels,
                    const std::vector<float> &labelScore,
                    float initialScore,
                    float bestScore,
                    const AlignParameters<Graph, SV> &paramList);

/**
 * @brief Prepare vectors in globalVars used primarily for branch bound checking, but also for isMaxColoredMatch in expand.
 *        The function also creates the initial state of the label matrices for the two graphs based on the anchor. 
 *     
 * @tparam Graph Graph type
 * @tparam SV The type of the scoring function for scoring vertices
 * @tparam Vertex The type of the vertices
 * @param a1 The first alignment object 
 * @param a2 The second alignment object, represents the smaller of the two graphs
 * @param preMap The anchor vertices in g1 and g2, always going from g1 to g2 
 * @param paramList Object of struct AlignParameters that contains the input parameters passed from the user
 * @param globalVars The object representing global function parameters
 */
template<typename Graph, class SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void prepareGlobalVars(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                        const std::vector<std::pair<Vertex, Vertex>> &preMap,
                        const AlignParameters<Graph, SV> &paramList,
                        GlobalVariables<Vertex> &globalVars);


////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template<typename Graph, class SV, typename Vertex>
bool branchBoundCheck(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const std::vector<Vertex> &subset,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars) {
    
    if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
        
        const std::vector<std::vector<int>> &largerLabelMatrix = globalVars.g1LabelMatrix;
        std::vector<std::vector<int>> &smallerLabelMatrix = globalVars.g2LabelMatrix;

        for(Vertex x : subset) {
            int labelId = a2.alignGraph[x].labelIndex;
            smallerLabelMatrix[labelId][globalVars.nInputGraphsForVertexG2[getIndex(x)]]--;
            smallerLabelMatrix[labelId][0]--;
        }

        bool valid = true;

        for(int i = 0; i < globalVars.nLabels; i++) {
            if(largerLabelMatrix[i][0] < smallerLabelMatrix[i][0]) {
                valid = false;
            }
        }

        // The inital score is the value provided by the anchor mapping, if such exists
        float hypotheticalScore = globalVars.initialScore;
        if(valid) {
            hypotheticalScore = hypotheticalScore + calculateHypotheticalScore(a1, a2, globalVars);
        }
        
        for(Vertex x : subset) {
            int labelId = a2.alignGraph[x].labelIndex;
            smallerLabelMatrix[labelId][globalVars.nInputGraphsForVertexG2[getIndex(x)]]++;
            smallerLabelMatrix[labelId][0]++;
        }

        // This subset is invalid as no embedding is possible. We thus need to return and continue with the next subset.
        if(!valid) {
            return false; 
        }
        // If the highest hypothetical score achievable from this subset is smaller than the current max score,
        // we do not even attempt the embedding of the smaller graph as the resulting score can never be better
        // than the current max.
        if(hypotheticalScore <= globalVars.bestScore) {
            return false;
        }
    }

    // Branch and bound checking when scoring by order in trimming algorithm is simply a matter of checking
    // whether the larger graph has enough vertices with each label to hypothetically 
    // cover the vertices of the given label in the smaller graph
    else if(paramList.enforceVLabels && paramList.scoring == Scoring::ORDER) {
    
        if(!globalVars.labelFrequencyG1.empty()) {
            std::vector<int> &labelFrequencyG2 = globalVars.labelFrequencyG2;
            for(Vertex x: subset) {
                labelFrequencyG2[a2.alignGraph[x].labelIndex]--;
            }

            // match the frequencies of the smaller graph with those of the larger graphs
            bool valid = true;
            for(int i = 0; i < globalVars.nLabels; i++) {
                // If there are more vertices with a given label in smallerTrimmed than biggerGraph then it is not possible 
                // to find a subgraph isomorphism between the two graphs and therefore, no reason to consider this subset.
                if(globalVars.labelFrequencyG1[i] < labelFrequencyG2[i]) {
                    valid = false;
                    break;
                }
            }
            
            for(Vertex x: subset) {
                labelFrequencyG2[a2.alignGraph[x].labelIndex]++;
            }
            if(!valid) {
                return false;
            }
        }
    }

    return true;
}

template<typename Graph, class SV, typename Vertex>
bool branchBoundCheck(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    float currentScore,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars) {
    
    if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
        
        float hypotheticalScore = calculateHypotheticalScore(a1, a2, globalVars);

        // If the hypothetical score is lower than the currently best achieved score, return
        if(hypotheticalScore + currentScore <= globalVars.bestScore) {
            return false;
        }
    }

    return true;
}


template<typename Graph, typename Vertex>
float calculateHypotheticalScore(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                                GlobalVariables<Vertex> &globalVars) {
    
    float hypotheticalScore = 0;

    const std::vector<std::vector<int>> &g1LabelMatrix = globalVars.g1LabelMatrix;
    const std::vector<std::vector<int>> &g2LabelMatrix = globalVars.g2LabelMatrix;

    for(int i = 0; i < globalVars.nLabels; i++) {

        int guardVal = std::min(g1LabelMatrix[i][0], g2LabelMatrix[i][0]);

        // If there are one of the two graphs that do not have any available vertices with label i it does not make sense
        // to do any work here, so we skip to the next label
        if(guardVal == 0) continue;

        // Find the starting point in the two label rows, starting from the right end
        // Save both the index of the entry found (representing the number of input graphs the vertices appear in)
        // and the actual number saved at that position (the number of vertices appearing in that number of input graphs with the given label)
        int containedValG1 = 0;
        int containedValG2 = 0;
        int remainingG1 = 0;
        int remainingG2 = 0;
        for(int j = a1.contained.size(); j > 0; j--) {
            if(g1LabelMatrix[i][j] > 0) {
                containedValG1 = j;
                remainingG1 = g1LabelMatrix[i][j];
                break;
            }
        }
        for(int j = a2.contained.size(); j > 0; j--) {
            if(g2LabelMatrix[i][j] > 0) {
                containedValG2 = j;
                remainingG2 = g2LabelMatrix[i][j];
                break;
            }
        }

        while(guardVal > 0) {
            // The scoring function is based on sum-of-pairs so the score is dependent on how many entries are added together
            // the pairings value represent the number of entries scored for a match.
            // Multiplying this value with the score for matching vertices of the given label gives the score for each match
            // multiply this by the number of matches done, i.e. the number of vertices for which the previously calculated score is valid
            
            float score = 0;
            int pairings = containedValG1 * containedValG2;
            
            // Boolean to guard locating the next values when calculated the hypothetical score for current values of
            // containedValG1 / G2 and remainingG1 / G2.
            bool foundNewSize = false;
            // When remainingG1 == remainingG2 both containedValG1 and containedValG2 needs to be updated after
            // calculating the addition to the hypothetical score
            if(remainingG1 == remainingG2) {
                score = globalVars.labelScore[i] * pairings * remainingG1;
                hypotheticalScore = hypotheticalScore + score;

                // Decrease the guard value (the number of vertices that can be matched with given label) 
                // by the number of vertices "matched" in this step
                guardVal = guardVal - remainingG1;

                foundNewSize = false;
                while(!foundNewSize) {
                    containedValG1--;
                    remainingG1 = g1LabelMatrix[i][containedValG1];
                    if(remainingG1 > 0) foundNewSize = true;
                }

                foundNewSize = false;
                while(!foundNewSize) {
                    containedValG2--;
                    remainingG2 = g2LabelMatrix[i][containedValG2];
                    if(remainingG2 > 0) foundNewSize = true;
                }
            }
            // Only containedValG1 needs to be updated after addition to the hypothetical score
            // remainingG2 needs to be decreased by the number of "matches", i.e. remainingG1
            else if(remainingG1 < remainingG2) {
                // When there are fewer vertices for the given number of input graphs in g1
                // calculate the score by using all of those and then decrease remainingG2 by remaininG1
                score = globalVars.labelScore[i] * pairings * remainingG1;
                hypotheticalScore = hypotheticalScore + score;

                // Need to decrease remainingG2 and guardVal before updating remainingG1
                remainingG2 = remainingG2 - remainingG1;
                guardVal = guardVal - remainingG1;

                // Find the next number of input graphs for which a vertex appears in g1
                foundNewSize = false;
                while(!foundNewSize) {
                    containedValG1--;
                    remainingG1 = g1LabelMatrix[i][containedValG1];
                    if(remainingG1 > 0) foundNewSize = true;
                }

            }
            // remainingG1 > remainingG2
            // Symmetric of above
            else {
                // When there are fewer vertices for the given number of input graphs in g2
                // calculate the score by using all of those and then decrease remainingG1 by remaininG2
                score = globalVars.labelScore[i] * pairings * remainingG2;
                hypotheticalScore = hypotheticalScore + score;

                // Need to decrease remainingG1 and guardVal before updating remainingG2
                remainingG1 = remainingG1 - remainingG2;
                guardVal = guardVal - remainingG2;

                // Find the next number of input graphs for which a vertex appears in g2
                foundNewSize = false;
                while(!foundNewSize) {
                    containedValG2--;
                    remainingG2 = g2LabelMatrix[i][containedValG2];
                    if(remainingG2 > 0) foundNewSize = true;
                }
            }
        }
    }
    return hypotheticalScore;
}

template<typename Graph, class SV, typename Vertex>
bool evaluateKLevel(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    int k,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars) {
    
    // Caller when access to globalVars
    return evaluateKLevel(a1, a2, globalVars.g1LabelMatrix, globalVars.g2LabelMatrix, k, globalVars.nLabels,
                    globalVars.labelScore, globalVars.initialScore, globalVars.bestScore, paramList);
}

template<typename Graph, class SV>
bool evaluateKLevel(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const std::vector<std::vector<int>> &g1LabelMatrix,
                    const std::vector<std::vector<int>> &g2LabelMatrix,
                    int k, 
                    int nLabels,
                    const std::vector<float> &labelScore,
                    float initialScore,
                    float bestScore,
                    const AlignParameters<Graph, SV> &paramList) {

    if(paramList.enforceVLabels && paramList.scoring == Scoring::SCHEME) {
        std::vector<float> smallestValues;
        smallestValues.reserve(k);
        float largest = 0;
        int saved = 0;
        int nRemove = k;
        // The inital score is the value provided by the anchor mapping, if such exists
        float hypotheticalScore = initialScore;
        // For each possible label check if there is enough vertices with given label in larger graph compared to the smaller graph
        // if so calculate the hypothetical score
        for(int i = 0; i < nLabels; i++) {
            int guardVal = std::min(g1LabelMatrix[i][0], g2LabelMatrix[i][0]);
            // If there are more vertices with label i in G2 they will always need to be removed 
            // and they do not contribute to the hypothetical score
            // We should only subtract the number of values that correspond to the vertices
            // that are removed beyond the ones that have to be removed for the subset to result in valid match
            if(g1LabelMatrix[i][0] < g2LabelMatrix[i][0]) {
                int diff = g2LabelMatrix[i][0] - g1LabelMatrix[i][0];
                nRemove = nRemove - diff;
            }
            // If there are one of the two graphs that do not have any available vertices with label i it does not make sense
            // to do any work here, so we skip to the next label
            if(guardVal == 0) continue;

            // Find the starting point in the two label rows, starting from the right end
            // Save both the index of the entry found (representing the number of input graphs the vertices appear in)
            // and the actual number saved at that position (the number of vertices appearing in that number of input graphs with the given label)
            int containedValG1 = 0;
            int containedValG2 = 0;
            int remainingG1 = 0;
            int remainingG2 = 0;
            for(int j = a1.contained.size(); j > 0; j--) {
                if(g1LabelMatrix[i][j] > 0) {
                    containedValG1 = j;
                    remainingG1 = g1LabelMatrix[i][j];
                    break;
                }
            }
            for(int j = a2.contained.size(); j > 0; j--) {
                if(g2LabelMatrix[i][j] > 0) {
                    containedValG2 = j;
                    remainingG2 = g2LabelMatrix[i][j];
                    break;
                }
            }

            while(guardVal > 0) {
                // The scoring function is based on sum-of-pairs so the score is dependent on how many entries are added together
                // the pairings value represent the number of entries scored for a match.
                // Multiplying this value with the score for matching vertices of the given label gives the score for each match
                // multiply this by the number of matches done, i.e. the number of vertices for which the previously calculated score is valid
                
                float score = 0;
                int pairings = containedValG1 * containedValG2;

                score = labelScore[i] * pairings;
                hypotheticalScore = hypotheticalScore + score;

                // Computed the score for one "match", and therefore an entry has been consumed from both
                // remaininG1 and remainingG2, and the guard should decrease by one
                guardVal--;
                remainingG1--;
                remainingG2--;
                
                // Boolean to guard locating the next values when calculated the hypothetical score for current values of
                // containedValG1 / G2 and remainingG1 / G2.
                bool foundNewSize = false;
                
                // After consuming an entry from remainingG1 and remainingG2 they may need to be updated if either reaches 0.
                // There are three cases determining which updates needs to be done. 
                if(remainingG1 == remainingG2 && remainingG1 == 0) { 
                    foundNewSize = false;
                    while(!foundNewSize) {
                        containedValG1--;
                        remainingG1 = g1LabelMatrix[i][containedValG1];
                        if(remainingG1 > 0) foundNewSize = true;
                    }

                    foundNewSize = false;
                    while(!foundNewSize) {
                        containedValG2--;
                        remainingG2 = g2LabelMatrix[i][containedValG2];
                        if(remainingG2 > 0) foundNewSize = true;
                    }
                }
                // Only containedValG1 needs to be updated after addition to the hypothetical score
                else if(remainingG1 < remainingG2 && remainingG1 == 0) { 

                    // This was the last vertex in remainingG1, but still some left in remainingG2
                    foundNewSize = false;
                    while(!foundNewSize) {
                        containedValG1--;
                        remainingG1 = g1LabelMatrix[i][containedValG1];
                        if(remainingG1 > 0) foundNewSize = true;
                    }
                }                
                // remainingG1 > remainingG2
                // Symmetric of above
                else if(remainingG1 > remainingG2 && remainingG2 == 0) {

                    // This was last entry in remainingG2, but still some left in remainingG1                    
                    foundNewSize = false;
                    while(!foundNewSize) {
                        containedValG2--;
                        remainingG2 = g2LabelMatrix[i][containedValG2];
                        if(remainingG2 > 0) foundNewSize = true;
                    }
                }

                // Updating the k smallest values is the same no matter which of the three cases from above was used
                // Need to always save k values, so if there have been saved less than k values so far
                // simply save the current score and increment the saved value and potentially update largest value. 
                if(saved < k) {
                    smallestValues.push_back(score);
                    saved++;
                    if(score > largest) {
                        largest = score;
                    }
                }
                // When k values have already been seen, switch out the largest one with new value, if new value is smaller than
                // the largest stored
                else if(score < largest) {
                    float max = 0;
                    int maxIndex = 0;
                    for(std::size_t i = 0; i < smallestValues.size(); i++) {
                        if(smallestValues[i] > max) {
                            max = smallestValues[i];
                            maxIndex = i;
                        } 
                    }
                    smallestValues[maxIndex] = score;
                    float maxRemaining = 0;
                    for(float value : smallestValues) {
                        if(value > maxRemaining) {
                            maxRemaining = value;
                        }
                    }
                    largest = maxRemaining;
                }
            }
        }
        
        // Subtracting the nRemove smallest values from the hypothetical score as this represents the
        // largest hypothetical score over all valid subsets generated when removing k vertices
        std::sort(smallestValues.begin(), smallestValues.end());
        for(int i = 0; i < nRemove; i++) {
            hypotheticalScore = hypotheticalScore - smallestValues[i];
        }

        // If the hypothetical score is not better than the best achieved score return false
        // to indicate the program should stop as it is impossible to find a better match at any later point
        if(hypotheticalScore <= bestScore) {
            return false;
        }
    }
    return true;
}


template<typename Graph, class SV, typename Vertex>
void prepareGlobalVars(const AlignObj<Graph> &a1, const AlignObj<Graph> &a2,
                    const std::vector<std::pair<Vertex, Vertex>> &preMap,
                    const AlignParameters<Graph, SV> &paramList,
                    GlobalVariables<Vertex> &globalVars) {
    
    // There is no need to do any preparations if we do not enforce vertex labels
    if(paramList.enforceVLabels) {
        Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();

        // First assign the anchor as already matched vertices
        std::vector<unsigned char> currentlyMatchedG1(num_vertices(a1.alignGraph), 0);
        std::vector<unsigned char> currentlyMatchedG2(num_vertices(a2.alignGraph), 0);
        for(const auto &pair : preMap) {
            currentlyMatchedG1[getIndex(pair.first)] = 1;
            currentlyMatchedG2[getIndex(pair.second)] = 1;
        }

        // Creating global currently matched vectors, so as to not create them in each recursive call
        globalVars.notAllowedVerticesG1 = currentlyMatchedG1;
        globalVars.notAllowedVerticesG2 = currentlyMatchedG2;


        // Prepare vectors used in isMaxColoredMatch and when branch bound in trimming and scoring by order.
        // The vectors maintain the frequency of each label amongst unmatched vertices.
        if(paramList.enforceELabels || (paramList.algorithm == Algorithm::TRIM && paramList.scoring == Scoring::ORDER)) {

            // increment the position for given label when encountering a vertex that is not a part of the anchor
            std::vector<int> labelCounterG1(globalVars.nLabels, 0);
            for(Vertex x : asRange(vertices(a1.alignGraph))) {
                if(!currentlyMatchedG1[getIndex(x)]) {
                    labelCounterG1[a1.alignGraph[x].labelIndex]++;
                }
            }

            // increment the position for given label when encountering a vertex that is not a part of the anchor
            std::vector<int> labelCounterG2(globalVars.nLabels, 0);
            for(Vertex x : asRange(vertices(a2.alignGraph))) {
                if(!currentlyMatchedG2[getIndex(x)]) {
                    labelCounterG2[a2.alignGraph[x].labelIndex]++;
                }
            }

            // Instead of calculating the frequencies of the labels of unmatched vertices constantly 
            // in isMaxColoredMatch maintain these counters when adding / removing from the current match
            globalVars.labelFrequencyG1 = labelCounterG1;
            globalVars.labelFrequencyG2 = labelCounterG2;

            // When the algorithm is trim and scoring by order the label frequencies are used to determine 
            // if it is even possible for the given trimmed graph to be covered by the vertices in the larger graph
            // however if it is possible to cover the smaller graph with the larger one without trimming any vertices
            // there are no reason to carry out the branch bound check for each subset as it will always be true
            // Let globalVars.labelFrequencyG1 be empty, when this is the case which is checked against in branchBoundCheck
            if(paramList.algorithm == Algorithm::TRIM && paramList.scoring == Scoring::ORDER) {
                bool allCovered = true;
                for(int i = 0; i < globalVars.nLabels; i++) {
                    if(labelCounterG1[i] < labelCounterG2[i]) {
                        allCovered = false;
                        break;
                    }
                }
                if(allCovered) globalVars.labelFrequencyG1 = {};
            }
        }

        // When scoring by scheme we can perform branch bound checking in terms of the hypothetical score.
        // To do so it is necessary to maintain information about how many unmatched vertices of each label there are
        // but also in how many input graphs each of those unmatched vertices appear in so they can be stacked in optimal fashion 
        // for the largest hypothetical score. 
        if(paramList.scoring == Scoring::SCHEME) {

            std::vector<std::vector<int>> g1LabelMatrix(globalVars.nLabels, std::vector<int>(a1.contained.size() + 1, 0));
            std::vector<std::vector<int>> g2LabelMatrix(globalVars.nLabels, std::vector<int>(a2.contained.size() + 1, 0));

            std::vector<int> nInputGraphsForVertexG1(num_vertices(a1.alignGraph), 0);
            std::vector<int> nInputGraphsForVertexG2(num_vertices(a2.alignGraph), 0);

            for(Vertex x : asRange(vertices(a1.alignGraph))) {
                if(!currentlyMatchedG1[getIndex(x)]) {
                    int appearanceCounter = 0;
                    for(std::size_t i = 0; i < a1.contained.size(); i++) {
                        if(a1.matchTable[i][x] != nullVertex) {
                            appearanceCounter++;
                        }
                    }
                    int labelId = a1.alignGraph[x].labelIndex;
                    nInputGraphsForVertexG1[getIndex(x)] = appearanceCounter;
                    g1LabelMatrix[labelId][appearanceCounter]++;
                    g1LabelMatrix[labelId][0]++;
                }
            }
            for(Vertex x : asRange(vertices(a2.alignGraph))) {
                if(!currentlyMatchedG2[getIndex(x)]) {
                    int appearanceCounter = 0;
                    for(std::size_t i = 0; i < a2.contained.size(); i++) {
                        if(a2.matchTable[i][x] != nullVertex) {
                            appearanceCounter++;
                        }
                    }
                    int labelId = a2.alignGraph[x].labelIndex;
                    nInputGraphsForVertexG2[getIndex(x)] = appearanceCounter;
                    g2LabelMatrix[labelId][appearanceCounter]++;
                    g2LabelMatrix[labelId][0]++;
                }
            }

            globalVars.g1LabelMatrix = g1LabelMatrix;
            globalVars.g2LabelMatrix = g2LabelMatrix;

            globalVars.nInputGraphsForVertexG1 = nInputGraphsForVertexG1;
            globalVars.nInputGraphsForVertexG2 = nInputGraphsForVertexG2;

            if(paramList.algorithm == Algorithm::TRIM) {
                globalVars.initialScore = matchScore(a1, a2, globalVars.getLargerMatch(), globalVars.getSmallerMatch(), paramList);
            }
        }
    }
}

} // namespace GraphAlign

#endif // BRANCHBOUND_H