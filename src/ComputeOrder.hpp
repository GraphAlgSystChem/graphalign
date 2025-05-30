#ifndef COMPUTEORDER_H
#define COMPUTEORDER_H

namespace GraphAlign {

/**
 * @brief Computes the order and parent list for vertices in g, prioritise adding vertices with a 
 *        larger degree to already added vertices first. 
 * 
 * @tparam Graph Graph type
 * @tparam Vertex Vertex type
 * @param g Graph to compute the order for
 * @param anchor A list containing the matched vertices, as the vertices thar are matched in g needs to be considered (the second vertex in each pair)
 * @return List of vertices in an ordered manner and list of parents for all vertices, 
 *         parent is either a vertex previously added to the order, or null_vertex
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::pair<std::vector<Vertex>, std::vector<Vertex>> 
    computeOrder(const FilteredGraph<Graph> &g, const std::vector<std::pair<Vertex, Vertex>> &anchor);

/**
 * @brief Prepares the degree bag and auxiliary vectors to maintain this bag given an initial (anchor mapping) in globalVars.
 * For a vertex v, we have globalVars.whichDegBag[v] = k if v is a neighbour to k already matched vertices.
 * 
 * @tparam Graph Graph type
 * @tparam Vertex Vertex type
 * @param g The domain graph
 * @param matchVector The vector indicating which vertices of g have been mapped already (i.e. if anchor vertex)
 * @param parentMap The parent map to be filled up. Will only be updated for vertices adjacent to the anchored vertices.
 * @param globalVars The structure containing the global variables, namely the one that will hold the four vectors needed to maintain the degree bags.
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void prepareDynamicOrder(const Graph &g, 
                         const std::vector<Vertex> &matchVector, 
                         std::vector<Vertex> &parentMap,
                         GlobalVariables<Vertex> &globalVars);

/**
 * @brief Extracts the next vertex whose degree to previously mapped vertices is maximum, namely in the largest non-empty bag. 
 * Returns the vertex and the bag it belonged to.
 * 
 * @tparam Graph Graph type
 * @tparam Vertex Vertex type
 * @param g The graph to extract a vertex from (needed to update the neighbours of the extracted vertex)
 * @param parentMap The parent map to be updated wrt. the neighbours of the extracted vertex.
 * @param globalVars The structure containing the global variables, namely the one that holds the four vectors needed to maintain the degree bags.
 */
template<typename Graph, typename Vertex>
std::pair<Vertex, int> getNextVertex(const Graph &g, std::vector<Vertex> &parentMap, GlobalVariables<Vertex> &globalVars);

/**
 * @brief Restores the degree bags to their original state before v was extracted (allowing v to be skipped).
 * 
 * @tparam Graph Graph type
 * @tparam Vertex Vertex type
 * @param g The graph from which v was extracted.
 * @param v The vertex that was previously extracted and must now be re-placed inside its original bag.
 * @param oldDegBag The integer returned from getNextVertex
 * @param parentMap The parent map to be updated wrt. the neighbours of the extracted vertex.
 * @param globalVars The structure containing the global variables, namely the one that holds the four vectors needed to maintain the degree bags.
 */
template<typename Graph, typename Vertex>
void restoreDegBags(const Graph &g, Vertex v, int oldDegBag, std::vector<Vertex> &parentMap, GlobalVariables<Vertex> &globalVars);

/**
 * @brief Update the degree bags, moving v into the "already considered" bag.
 * 
 * @tparam Vertex Vertex type
 * @param v The vertex to skip
 * @param globalVars The structure containing the global variables, namely the one that holds the four vectors needed to maintain the degree bags.
 */
template<typename Vertex>
int skipVertex(Vertex v, GlobalVariables<Vertex> &globalVars);

/**
 * @brief Inserts v into its original bag s.t. it can be considered again when backtracking from a skipping branch.
 * 
 * @tparam Graph Graph type
 * @tparam Vertex Vertex type
 * @param v The vertex that was previously skipped and must now be re-placed inside its original bag.
 * @param globalVars The structure containing the global variables, namely the one that holds the four vectors needed to maintain the degree bags.
 */
template<typename Vertex>
void restoreDegBagsSkipped(Vertex v, int oldDegBag, GlobalVariables<Vertex> &globalVars);

template<typename Graph, typename Vertex>
void prepareUnmatchedNeighbours(const Graph &g, GlobalVariables<Vertex> &globalVars);

template<typename Graph, typename Vertex>
void updateUnmatchedNeighbours(Vertex v, const Graph &g, GlobalVariables<Vertex> &globalVars);

template<typename Graph, typename Vertex>
void restoreUnmatchedNeighbours(Vertex v, const Graph &g, GlobalVariables<Vertex> &globalVars);

/**
 * @brief Swap function that swaps the position of the elements v1 and v2 in the vec array and the
 *        corresponding indexes in the pos array. 
 * 
 * @tparam T Type of the elements v1 and v2, important that this are elements valid as indexing into a vector
 * @param vec Vector containing elements of type T, must contain the elements v1 and v2
 * @param pos Positions vector such that pos[v1] gives the index of the position that v1 occupies in vec, i.e. vec[pos[v1]] returns v1
 * @param v1 First element
 * @param v2 Second element
 */
template<typename T>
void swap(std::vector<T> &vec, std::vector<int> &pos, T v1, T v2);

////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////

template<typename Graph, typename Vertex>
std::pair<std::vector<Vertex>, std::vector<Vertex>> 
    computeOrder(const FilteredGraph<Graph> &g, const std::vector<Vertex> &anchorVertices) {
    
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    size_t nVisibleVertices = num_visible_vertices(g); // The number of vertices remaining in the trimmed graph
    size_t totalVertices = num_vertices(g); // number of vertices in the graph when including the trimmed vertices
    
    // Find the maximum degree among the vertices in g
    int maxDegree = 0;
    for(Vertex x : asRange(vertices(g))) {
        const auto &pair = adjacent_vertices(x, g);
        int dist = std::distance(pair.first, pair.second);
        if(dist > maxDegree) {
            maxDegree = dist;
        }
    }

    // Bags are: 0, 1, 2, ..., maxDegree, consideredVertices (= maxDegree + 2 bags)
    
    std::vector<Vertex> bags(nVisibleVertices, nullVertex); // the placement of visible vertices in the bags
    std::vector<int> bagEndIndex(maxDegree + 2, nVisibleVertices); // bagEndIndex[i] = j if the end of bag i is at index j in the bags vector. 
                                                                   // The bagEndIndex is similar to .end() from std::vector as it points outside of the bag.
                                                                   // Initially, all bags except the 0-bag are empty, so all .end()-pointers point outside the bag array.
    std::vector<int> vertexBagIndex(totalVertices, 0); // vertexBagIndex[v] = i if v belongs to bag i
    std::vector<int> vertexBagsPosition(totalVertices, 0); // vertexBagsPosition[v] = i if bags[i] = v

    // The data structures used to maintain the order
    std::vector<Vertex> order; // the order of the vertices computed
    order.reserve(nVisibleVertices); // need to store all the visible vertices in the trimmed graph
    std::vector<Vertex> parent(totalVertices, nullVertex); // parent map, indexed with the vertices, therefore, needs room for all vertices
    std::size_t nUsed = 0;
    
    // Initially, all anchored vertices are in the last bag and all non-anchored vertices are in the first bag.
    for(auto v : anchorVertices){
        vertexBagIndex[v] = maxDegree + 1;
        order.push_back(v);
        nUsed++;
    }
    // Move the end of all bags (except the considered one) to the beginning of the last bag.
    for(int i = 0; i <= maxDegree; i++){
        bagEndIndex[i] = nVisibleVertices - anchorVertices.size();
    }

    // Fill up the bags. Initially all vertices are in bag 0, and all anchored vertices are in bags k + 1
    // Fill up the array left-><-right
    int leftPos = 0;
    int rightPos = nVisibleVertices - 1;
    for(auto v : asRange(vertices(g))){ // only traverses the visible vertices
        if(vertexBagIndex[v] == maxDegree + 1) { // anchor vertex
            bags[rightPos] = v;
            vertexBagsPosition[v] = rightPos;
            rightPos--;
        } else {
            bags[leftPos] = v;
            vertexBagsPosition[v] = leftPos;
            leftPos++;
        }
    }

    // All four vectors have been initialised and can now be used to move around vertices that are neighbours to anchored vertices.
    // The next step is to consider all anchored vertices and their neighbours, moving them up into higher bags and adjusting their parent maps.
    for(auto u : anchorVertices){
        for(auto v : asRange(adjacent_vertices(u, g))){
            int bag = vertexBagIndex[v]; // the bag v belongs to
            if(bag < maxDegree + 1){ // Only move neighbours if they are not anchored themselves.
                auto lastInBag = bags[bagEndIndex[bag] - 1]; // The last vertex in this bag.
                swap(bags, vertexBagsPosition, v, lastInBag); // swap v with the last verte in this bag.
                vertexBagIndex[v]++;
                bagEndIndex[bag]--;
                if(parent[v] == nullVertex){ // Set parent if not set previously
                    parent[v] = u;
                }
            }
        }
    }
    // The bags array has been prepared to compute the total order. Anchored vertices are already in the order, and the remaining vertices
    // have been put into bags corresponding to the number of edges they share with already mapped (anchored) vertices.
    

    // We now add all the vertices in the order.
    while(nUsed < nVisibleVertices) {
        Vertex nextVertex; // find the vertex stored at the first index of the greatest non-empty bag in the bags vector
        int oldDegBag;
        
        // If the zero-bag ends in the same position as the max-degree bag, it means that the 0-bag is the only non-empty bag (meaning that we are in a new connected component).
        // We thus must choose the vertex from this bag with the largest degree to start the new connected component.
        if(bagEndIndex[0] == bagEndIndex[maxDegree]) {
            nextVertex = bags[0];
            int bestDegree = std::distance(adjacent_vertices(nextVertex, g).first, adjacent_vertices(nextVertex, g).second);
            for(int i = 1; i < bagEndIndex[0]; i++) {
                auto potentialNext = bags[i];
                int degree = std::distance(adjacent_vertices(potentialNext, g).first, adjacent_vertices(potentialNext, g).second);
                if(degree > bestDegree) {
                    nextVertex = potentialNext;
                    bestDegree = degree;
                }
                else if(degree == bestDegree && getIndex(potentialNext) < getIndex(nextVertex)) {
                    nextVertex = potentialNext;
                }
            }
            // Once the next vertex is found, swap it such that it is positioned last in the bag
            swap(bags, vertexBagsPosition, nextVertex, bags[bagEndIndex[maxDegree] - 1]);
            oldDegBag = 0;
        }
        else {
            // Otherwise, just choose the next vertex as the "lexicographically largest" one in the largest non-empty bag
            nextVertex = bags[bagEndIndex[maxDegree] - 1];
            const auto &adjPair = adjacent_vertices(nextVertex, g);
            auto nextVertexDegree = std::distance(adjPair.first, adjPair.second);
            auto bag = vertexBagIndex[nextVertex];
            // choose lexicographically LARGEST according to overall degree
            for(int i = bagEndIndex[bag - 1]; i < bagEndIndex[bag]; i++){

                const auto &bagPair = adjacent_vertices(bags[i], g);
                auto bagVertexDegree = std::distance(bagPair.first, bagPair.second);
                if(bagVertexDegree > nextVertexDegree){
                    nextVertex = bags[i];
                    nextVertexDegree = bagVertexDegree;
                }
                // tie breaker is lexicographically smallest
                else if(bagVertexDegree == nextVertexDegree && bags[i] < nextVertex){
                    nextVertex = bags[i];
                    nextVertexDegree = bagVertexDegree;
                }
            }
            oldDegBag = bag;
            swap(bags, vertexBagsPosition, nextVertex, bags[bagEndIndex[maxDegree] - 1]);
        }

        // Shift all ends of vertex bags that were stacked on top of each other (to the right of nextVertex)
        // As nextVertex is in the right-most non-empty bag (excluding k + 1), the .ends() must be for all bags
        // starting in nextVertex's bag and up to k.
        for(int i = oldDegBag; i <= maxDegree; i++){
            bagEndIndex[i]--;
        }
        vertexBagIndex[nextVertex] = maxDegree + 1;
        order.push_back(nextVertex); // adding the next vertex to the order

        // Now, we just need to move all nextVertex's neighbours up a bag and add nextVertex as their parent if relevant.
        // Update neighbourhood bags
        for(Vertex neigh : asRange(adjacent_vertices(nextVertex, g))) {
            if(vertexBagIndex[neigh] < maxDegree + 1) { // Not considered yet
                int bag = vertexBagIndex[neigh];
                swap(bags, vertexBagsPosition, neigh, bags[bagEndIndex[bag] - 1]);
                vertexBagIndex[neigh]++;
                bagEndIndex[bag]--;
                if(parent[neigh] == nullVertex) {
                    parent[neigh] = nextVertex;
                }
            }
        }
        nUsed++;
    }
    return {order, parent};
} // computeOrder

template<typename Graph, typename Vertex>
void prepareDynamicOrder(const Graph &g, 
    const std::vector<Vertex> &matchVector, 
    std::vector<Vertex> &parentMap,
    GlobalVariables<Vertex> &globalVars){

    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();

    int numAnchored = 0;
    int k = 0; // maxDegree
    // Find max degree in g
    for(auto v : asRange(vertices(g))){
        auto rangePair = adjacent_vertices(v, g);
        auto deg = std::distance(rangePair.first, rangePair.second);
        if(deg > k){
            k = deg;
        }
    }
    // Initialise degree bags and positions references, one for each vertex
    std::vector<Vertex> degBags(num_vertices(g), nullVertex); 
    std::vector<int> degBagEnds(k + 2, 0); 
    std::vector<int> whichDegBag(num_vertices(g), 0);
    std::vector<int> bagPos(num_vertices(g), 0);

    // Inform all anchored vertices that they are positioned in the last bag
    for(size_t i = 0; i < matchVector.size(); i++){
        if(matchVector[i] != nullVertex){
            // anchored vertex, should later be located in the last bag
            whichDegBag[i] = k + 1;
            numAnchored++;
        }
    }
    // Move the end of all bags to the beginning of the last bag.
    for(int i = 0; i <= k; i++){
        degBagEnds[i] = num_vertices(g) - numAnchored;
    }
    degBagEnds[k + 1] = degBags.size(); // will never be used as this is trivially true

    // Fill up the bags. Initially all vertices are in bag 0, and all anchored vertices are in bags k + 1
    // Fill up the array left-><-right
    int leftPos = 0;
    int rightPos = num_vertices(g) - 1;
    for(auto v : asRange(vertices(g))){
        if(matchVector[v] != nullVertex){
            degBags[rightPos] = v;
            bagPos[v] = rightPos;
            rightPos--;
        } else {
            degBags[leftPos] = v;
            bagPos[v] = leftPos;
            leftPos++;
        }
    }

    // All four vectors have been initialised and can now be used to move around vertices that are neighbours to anchored vertices.
    // The next step is to consider all anchored vertices and their neighbours, moving them up into higher bags and adjusting their parent maps.
    for(auto u : asRange(vertices(g))){
        if(matchVector[u] != nullVertex){
            for(auto v : asRange(adjacent_vertices(u, g))){
                int bag = whichDegBag[v];
                if(bag < k + 1){ // Only move neighbours if they are not anchored themselves.
                    Vertex lastInBag = degBags[degBagEnds[bag] - 1]; // The last vertex in this bag.
                    swap(degBags, bagPos, v, lastInBag); // swap v with the last verte in this bag.
                    whichDegBag[v]++;
                    degBagEnds[bag]--;
                    if(parentMap[v] == nullVertex){ // Set parent if not set previously.
                        parentMap[v] = u;
                    }
                }
            }
        }
    }
    globalVars.degBags = degBags;
    globalVars.degBagEnds = degBagEnds;
    globalVars.whichDegBag = whichDegBag;
    globalVars.bagPos = bagPos;
} // prepareDynamicOrder

template<typename Graph, typename Vertex>
std::pair<Vertex, int> getNextVertex(const Graph &g, std::vector<Vertex> &parentMap, GlobalVariables<Vertex> &globalVars){
    
    int k = globalVars.degBagEnds.size() - 2; // k (max degree) is the position of k-bag. degBags.size() - 1 points to the (k + 1) bag.
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    // The next vertex (u) is always positioned in the largest k bag. All consecutive empty bags have their
    // .end() stacked on top of each other, so degBagEnds[k] - 1 always points to the last vertex in the largest
    // non-empty degree bag even if bag k is empty.
    // Note also that u's position in degBags is unchanged. Only the .end()'s are moved for the different bags.

    // The next vertex should furthermore have the largest overall degree and, for comparison reasons, be the lexicographically smallest vertex.
    Vertex u = globalVars.degBags[globalVars.degBagEnds[k] - 1];
    int oldDegBag = globalVars.whichDegBag[u];
    int previousBagEnds = oldDegBag == 0 ? 0 : globalVars.degBagEnds[oldDegBag - 1]; // otherwise oldDegBag - 1 = -1 which would give undefined behavior in degBagEnds.
    int bestDegree = std::distance(adjacent_vertices(u, g).first, adjacent_vertices(u, g).second);
    for(int i = previousBagEnds; i < globalVars.degBagEnds[k]; i++){
        Vertex v = globalVars.degBags[i];
        auto beginEnd = adjacent_vertices(v, g);
        auto vDegree = std::distance(beginEnd.first, beginEnd.second);
        if(vDegree > bestDegree){
            u = v;
            bestDegree = vDegree;
        } else if (vDegree == bestDegree && v < u){
            u = v;
        }
    }
    swap(globalVars.degBags, globalVars.bagPos, u, globalVars.degBags[globalVars.degBagEnds[k] - 1]);
    
    // Shift all ends of vertex bags that were stacked on top of each other (to the right of u)
    // As u is in the right-most non-empty bag (excluding k + 1), the .ends() must be for all bags
    // starting in u's bag and up to k.
    for(int i = oldDegBag; i <= k; i++){
        globalVars.degBagEnds[i]--;
    }
    globalVars.whichDegBag[u] = k + 1; 

    // Move all u's neighbour up a bag
    for(Vertex v : asRange(adjacent_vertices(u, g))){
        int bag = globalVars.whichDegBag[v];
        if(bag < k + 1){ // Do not move vertices that have already been considered
            Vertex lastInBag = globalVars.degBags[globalVars.degBagEnds[bag] - 1]; // The last vertex in this bag.
            swap(globalVars.degBags, globalVars.bagPos, v, lastInBag);
            globalVars.whichDegBag[v]++;
            globalVars.degBagEnds[bag]--;
            if(parentMap[v] == nullVertex){
                parentMap[v] = u;
            }
        }
    }

    return std::make_pair(u, oldDegBag);
} // getNextVertex

template<typename Graph, typename Vertex>
void restoreDegBags(const Graph &g, Vertex v, int oldDegBag, std::vector<Vertex> &parentMap, GlobalVariables<Vertex> &globalVars){
    Vertex nullVertex = boost::graph_traits<Graph>::null_vertex();
    int k = globalVars.degBagEnds.size() - 2;
    // We move v neighbours down a bag first by swapping them with the first vertex in their current bag and moving the .end() of the lower bag up.
    for(Vertex u : asRange(adjacent_vertices(v, g))){
        int bag = globalVars.whichDegBag[u];
        if(bag < k + 1){ // Do not move vertices that have already been considered
            Vertex firstInBag = globalVars.degBags[globalVars.degBagEnds[bag - 1]]; // The first vertex in the current bag.
            swap(globalVars.degBags, globalVars.bagPos, u, firstInBag);
            globalVars.whichDegBag[u]--;
            globalVars.degBagEnds[bag - 1]++; // grow smaller bag.
            if(parentMap[u] == v){ // revert parent maps
                parentMap[u] = nullVertex;
            }
        }
    }

    // When v was extracted, it was in the right-most non-empty bag, so we know it must be positioned in this bag again. (And all bags stacked on top must be moved right.)
    // v will be inserted into its old bag but not necessarily in the exact position as it was before as the ordering is irrelevant.
    // Because v is moved by shifting bag ends, v's position does not change.
    for(int i = oldDegBag; i <= k; i++){
        globalVars.degBagEnds[i]++;
    }
    globalVars.whichDegBag[v] = oldDegBag;

} // restoreDegBags

template<typename Vertex>
int skipVertex(Vertex v, GlobalVariables<Vertex> &globalVars){
    // Precondition: v is the last element in its degree back (which whill be the case after restore has been called)
    int k = globalVars.degBagEnds.size() - 2;

    // Shift all of the .ends() for bags after v down as v has been removed.
    for(int i = globalVars.whichDegBag[v]; i <= k; i++){
        globalVars.degBagEnds[i]--;
    }
    int oldDegBag = globalVars.whichDegBag[v];
    globalVars.whichDegBag[v] = k + 1;

    return oldDegBag;

} // skipVertex

template<typename Vertex>
void restoreDegBagsSkipped(Vertex v, int oldDegBag, GlobalVariables<Vertex> &globalVars){
    // Move v into its old bag (the end of it) and shift ends upwards for the bags stacked on top of each other.
    // v does not change position and therefore neither degBags nor bagsPos are changed.

    int k = globalVars.degBagEnds.size() - 2;

    for(int i = oldDegBag; i <= k; i++){
        globalVars.degBagEnds[i]++;
    }
    globalVars.whichDegBag[v] = oldDegBag;

} // restoreDegBagsSkipped

template<typename Graph, typename Vertex>
void prepareUnmatchedNeighbours(const Graph &g, GlobalVariables<Vertex> &globalVars){
    globalVars.nUnmatchedNeighs = std::vector<int>(num_vertices(g), 0);
    for(auto v : asRange(vertices(g))){
        globalVars.nUnmatchedNeighs[v] = std::distance(adjacent_vertices(v, g).first, adjacent_vertices(v, g).second);
    }
}

template<typename Graph, typename Vertex>
void updateUnmatchedNeighbours(Vertex v, const Graph &g, GlobalVariables<Vertex> &globalVars){
    for(auto u : asRange(adjacent_vertices(v, g))){
        globalVars.nUnmatchedNeighs[u]--;
    }
} // updateUnmatchedNeighbours

template<typename Graph, typename Vertex>
void restoreUnmatchedNeighbours(Vertex v, const Graph &g, GlobalVariables<Vertex> &globalVars){
    for(auto u : asRange(adjacent_vertices(v, g))){
        globalVars.nUnmatchedNeighs[u]++;
    }

} // restoreUnmatchedNeighbours

template<typename T>
void swap(std::vector<T> &vec, std::vector<int> &pos, T v1, T v2) {
    // save old values for the value and position, respectively in vec and pos
    auto tmpPos = pos[v2];
    auto tmpVal = vec[pos[v2]];

    // swap the values in vec
    vec[pos[v2]] = vec[pos[v1]];
    vec[pos[v1]] = tmpVal;

    // swap the position pointers 
    pos[v2] = pos[v1];
    pos[v1] = tmpPos;
}

} // namespace GraphAlign

#endif // COMPUTEORDER_H