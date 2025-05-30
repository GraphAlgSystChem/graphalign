// Library for solving combinatoric problems within alignment.

// std
#include <iostream>
#include <vector>
#include <queue>
#include <memory>

// ours
#include "AlignObj.hpp"
#include "PairToRangeAdaptor.hpp"

namespace GraphAlignUtility {

    // getCombinations
    // combinationUtil
    // getAllGuideTrees
    // _allRoots
    // allBipartitions
    // _allBipartitionsAux
    // makeQueue
    // _makeQueueAux

/**
 * Computes a vector of (n choose r) combinations of the values of arr, each combination being represented as a
 * a vector of length r.
 * 
 * Code adapted from: https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
 * C++ program to print all combination of size r in an array of size n. 
 * Code adapted to not print but return combinations instead.
 * Original code is contributed by rathbhupendra
 */
template <typename T>
std::vector<std::vector<T>> getCombinations(std::vector<T> &arr, int n, int r);

/**
 * arr ---> Input vector 
 * data[] ---> Temporary array to store current combination 
 * start & end ---> Starting and Ending indexes in arr[] 
 * index ---> Current index in data[] 
 * r ---> Size of a combination to be printed 
 */
template <typename T>
void combinationUtil(std::vector<T> arr, T data[], 
                            int start, int end, 
                            int index, int r, std::vector<std::vector<T>> &basket);

/**
 * @brief Produces a list of all possible guide trees, all guide trees together covering all combinations of aligning the graphs of Gs.
 *
 * @tparam Graph The type of input graphs.
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @param Gs The input graphs of which all guide trees should be produced.
 */
template <typename Graph>
std::vector<std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>>> getAllGuideTrees(const std::vector<Graph> &Gs);


/**
 * @brief Produces a list of all guide trees, namely the list of all root nodes in these.
 * 
 * @tparam Graph The type of input graphs.
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @tparam Vertex The vertex type
 * @param Gs The input graphs of which all roots should be produced.
 */
template <typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> _allRoots(const std::vector<Graph> &Gs, const std::vector<int> &indices);


/**
 * @brief Produces a list of all unique bipartitions of the list of elements. A bipartition is a pair of disjoint vectors whose union contains exactly the elements of elems.
 * 
 * @tparam T the type of elements in elems.
 * @param elems The original list of elements to partition.
 */
template <typename T>
std::vector<std::pair<std::vector<T>, std::vector<T>>> allBiPartitions(const std::vector<T> &elems);

/**
 * @brief Auxiliary function for creating all bipartitions.
 * 
 * @tparam T The type of elements wanting to create the bipartitions over.
 * @param index The next element to insert into partitions
 * @param elems The elements to create a bipartition
 * @param S1 The first part of the partition being built
 * @param S2 The second part of the partition being built
 * @param bag Bag containing the partitions
 */
template <typename T>
void _allBiPartitionsAux(size_t index, const std::vector<T> &elems, std::vector<T> &S1, std::vector<T> &S2, std::vector<std::pair<std::vector<T>, std::vector<T>>> &bag);

/**
 * @brief Produces a valid queue of alignments w.r.t. to the root of the guide tree.
 * 
 * @tparam Graph The type of input graphs.
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @param root The root of the guide tree of which a queue of alignment vertices must be created.
 */
template <typename Graph>
std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> makeQueue(const std::shared_ptr<GraphAlign::AlignObj<Graph>> &root);

/**
 * @brief Auxiliary function for filling up a queue recursively.
 * 
 * @tparam Graph The type of input graphs 
 * @tparam MatchTable the type of the match table used, either map of map or vector of vector
 * @param root The root of the subtree being built recursively
 * @param Q Queue containing the trees
 */
template <typename Graph>
void _makeQueueAux(const std::shared_ptr<GraphAlign::AlignObj<Graph>> &root, std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> &Q);


////////////////////////////////// IMPLEMENTATIONS //////////////////////////////////


template <typename T>
std::vector<std::vector<T>> getCombinations(std::vector<T> &arr, int n, int r){

    std::vector<std::vector<T>> combinations;
    
    if(r == 0){
        // choosing 0 yields one unique combination, namely the empty one.
        std::vector<T> oneCombination;
        combinations.push_back(oneCombination);
        return combinations;
    } else if(r == n){
        combinations.push_back(arr);
        return combinations;
    }
    size_t expSize = 1;
    for(int i = 1; i <= r; i++){
        expSize = expSize * (n + 1 - i)/(i);
    }
    combinations.reserve(expSize);
    // A temporary array to store
    // all combination one by one 
    T data[r];

    // Gathers all combinations in combinations
    combinationUtil(arr, data, 0, n-1, 0, r, combinations);

    return combinations;
} // getCombinations

template <typename T>
void combinationUtil(std::vector<T> arr, T data[], 
                    int start, int end, 
                    int index, int r, std::vector<std::vector<T>> &basket) 
{ 
    // Current combination is ready
    // to be printed, print it 
    if (index == r) 
    { 
        std::vector<T> combination;
        combination.reserve(r);
        for (int j = 0; j < r; j++) 
            combination.push_back(data[j]);
        basket.push_back(combination);
        return; 
    } 

    // replace index with all possible 
    // elements. The condition "end-i+1 >= r-index"
    // ensures that, upon inserting the index'th (index <= r) element in a combination,
    // there are at least (r - index) elements remaining in the input array 
    // to be inserted in later recursive calls. Skipping the tailing elements here
    // is acceptable as they will be considered again in the next layer of recursion.
    for (int i = start; i <= end && 
        end - i + 1 >= r - index; i++) 
    { 
        data[index] = arr[i]; 
        combinationUtil(arr, data, i+1, 
                        end, index+1, r, basket); 
    } 
} // combinationUtil

template <typename Graph>
std::vector<std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>>> getAllGuideTrees(const std::vector<Graph> &Gs){
    // create indexTable of graphs before passing them on to _allRoots
    std::vector<int> indices;
    for(std::size_t i = 0; i < Gs.size(); i++) {
        indices.push_back(i);
    }
    
    std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> bag;
    // Collecting all roots of all guide trees
    bag = _allRoots<Graph>(Gs, indices);
    // For each root, creating the queue representing the guide tree traversal
    std::vector<std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>>> Qs;
    for(const auto &root : bag){
        std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> Q;
        Q = makeQueue(root);
        Qs.push_back(Q);
    }

    return Qs;
} // getAllGuideTrees


template <typename Graph, typename Vertex>
std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> _allRoots(const std::vector<Graph> &Gs, const std::vector<int> &indices){
    std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> R;
    // Create a leaf node
    if(indices.size() == 1){
        int graphIndex = indices[0];
        Graph g = Gs[graphIndex];
        std::shared_ptr<GraphAlign::AlignObj<Graph>> r = std::make_shared<GraphAlign::AlignObj<Graph>>(g, graphIndex);
        
        // Create trivial match table for each leaf node
        r->contained = {graphIndex};
        R = {r};
        return R;
    } else{
        std::vector<std::pair<std::vector<int>, std::vector<int>>> P = allBiPartitions(indices);
        for(const auto &p : P){
            std::vector<int> p1 = p.first;
            std::vector<int> p2 = p.second;
            std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> S = _allRoots<Graph>(Gs, p1);
            std::vector<std::shared_ptr<GraphAlign::AlignObj<Graph>>> T = _allRoots<Graph>(Gs, p2);
            
            // S and T hold all guide trees corresponding to alignments of p1 and p2 respectively.
            // We now create all possible combinations of these guide trees (no symmetric pairs, however)
            for(const auto &s : S){
                for(const auto &t : T){
                    std::vector<int> mergeContained = s->contained;
                    mergeContained.insert(mergeContained.end(), t->contained.begin(), t->contained.end());
                    std::shared_ptr<GraphAlign::AlignObj<Graph>> r = std::make_shared<GraphAlign::AlignObj<Graph>>(s, t, mergeContained);
                    r->tag = std::string("(" + s->tag + "," + t->tag + ")");
                    R.push_back(r);
                }
            }
        }
        return R;
    }

}// _allRoots

template <typename T>
std::vector<std::pair<std::vector<T>, std::vector<T>>> allBiPartitions(const std::vector<T> &elems){
    assert(!elems.empty());
    std::vector<T> S1, S2;
    std::vector<std::pair<std::vector<T>, std::vector<T>>> B;
    S1.push_back(elems[0]);
    // Create all bipartitions where the first element is in the left set.
    // (As a means to avoid symmetric duplicates)
    _allBiPartitionsAux(1, elems, S1, S2, B);
    return B;

} // allBiPartitions

template <typename T>
void _allBiPartitionsAux(size_t index, const std::vector<T> &elems, std::vector<T> &S1, std::vector<T> &S2, std::vector<std::pair<std::vector<T>, std::vector<T>>> &bag){
    if (index == elems.size()){
        // A part in the partition might be empty as a result of the recursion
        // but is not allowed by definition, so we scrap the whole partition if that's the case.
        if(!S1.empty() && !S2.empty()){
            bag.push_back({S1, S2});
        } 
    } else {
        T e = elems[index];
        std::vector<T> S1_e = S1;
        std::vector<T> S2_e = S2;
        S1_e.push_back(e);
        S2_e.push_back(e);
        // The next element in line can either be in the part
        // containing the first element or in the other part.
        _allBiPartitionsAux(index + 1, elems, S1_e, S2, bag);
        _allBiPartitionsAux(index + 1, elems, S1, S2_e, bag);
    }
} // _allBiPartitionsAux

template <typename Graph>
std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> makeQueue(const std::shared_ptr<GraphAlign::AlignObj<Graph>> &root){
    std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> Q;
    _makeQueueAux(root, Q);
    return Q;
} // makeQueue


template <typename Graph>
void _makeQueueAux(const std::shared_ptr<GraphAlign::AlignObj<Graph>> &root, std::queue<std::shared_ptr<GraphAlign::AlignObj<Graph>>> &Q){
    if( root->left->contained.size() == 1 && root->right->contained.size() == 1){
        Q.push(root);
    } else if (root->left->contained.size() > 1 && root->right->contained.size() == 1){
        _makeQueueAux(root->left, Q);
        Q.push(root);
    } else if (root->left->contained.size() == 1 && root->right->contained.size() > 1){
        _makeQueueAux(root->right, Q);
        Q.push(root);
    } else{
        _makeQueueAux(root->left, Q);
        _makeQueueAux(root->right, Q);
        Q.push(root);
    }
} // _makeQueueAux


} // namespace GraphAlignUtility