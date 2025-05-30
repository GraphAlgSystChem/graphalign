#ifndef TESTSUITE_H
#define TESTSUITE_H

// std
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <string>
#include <sstream> // for replacing suffixes
#include <cassert> // for assertion
#include <thread> // for timeoutable thread
#include <mutex> // for the thread to lock on when waiting
#include <condition_variable> // variable indicating "wake-up" call
#include <chrono> // for time taking and variables

// 3rd party
#include <boost/graph/adjacency_list.hpp> // graph traits
#include <boost/algorithm/string.hpp> 
#include "../src/external/json.hpp" // for json parsing, header only by Niels Lohmann

// us
#include "../src/MGA_Analysis.hpp"
#include "../src/PairToRangeAdaptor.hpp"
#include "../src/FatalException.hpp"
#include "../src/AlignObj.hpp"
#include "GraphIO.hpp"

namespace GraphAlignTest {

    // readGraphs
    // writeAlignmentToFile
    // readJSON
    // _wrapper_runAlgorithm
    // testSingleWTimeout
    // testSingle
    // testSingleChooseGuideTree

    using json = nlohmann::json;

    

/** Taken from https://stackoverflow.com/questions/40550730/how-to-implement-timeout-for-function-in-c
 * @brief Wrapper function for runAlgorithm that waits for results or timeouts in the process.
 * @param timeoutWait is the only parameter that differens from runAlgorithm. This is the expected duration for waiting.
 */
template<typename Graph, typename SV,
        typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> 
    _wrapper_runAlgorithm(const std::chrono::seconds &timeoutWait,
    const std::vector<Graph> &inputGraphs,
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

/**
 * @brief Runs the PGA algorithm on the single .json graph instance with the parameters specified
 * in a config file. Uses the timeout parameter specified in the config file.
 * @param graphPath The filepath to the .json file containing the input graphs.
 * @param configPath The filepath to the .json file containing the PGA parameters.
 * @param reportPath The filepath to the .json file which will contain the final alignment report.
 */
template <typename Graph, typename SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void testSingleWTimeout(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::filesystem::path &reportPath="report.json");

/**
 * @brief Runs the PGA algorithm on the single .json graph instance with the parameters specified
 * in a config file. Prints |V(AG)|, |E(AG)| and the alignment score to std out, separated by tab characters.
 * @param graphPath The filepath to the .json file containing the input graphs.
 * @param configPath The filepath to the .json file containing the PGA parameters.
 * @param reportPath The filepath to the .json file which will contain the final alignment report.
 */
template <typename Graph, typename SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void testSingle(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::filesystem::path &reportPath="report.json");

/**
 * @brief Runs the PGA algorithm on the single .json graph instance with the parameters specified
 * in a config file - tries all guide tree configurations. 
 * Prints |V(AG)|, |E(AG)| and the alignment score to std out, separated by tab characters.
 * @param graphPath The filepath to the .json file containing the input graphs.
 * @param configPath The filepath to the .json file containing the PGA parameters.
 */
template <typename Graph, typename SV, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void testSingleChooseGuideTree(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::string &treeString = "", const std::filesystem::path &reportPath = "");


template<typename Graph, typename SV, typename Vertex>
std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> 
    _wrapper_runAlgorithm(const std::chrono::seconds &timeoutWait,
    const std::vector<Graph> &inputGraphs,
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
    const std::vector<std::vector<Vertex>> &anchor){
    
    std::mutex m;
    std::condition_variable cv;
    std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> retValue;


    std::thread t([&cv, &retValue, inputGraphs, &scoreVertices, &kernelType, 
                &isSparse, &measureType, &enforceVLabels, &enforceELabels, &useAmbiguousEdges, &oneMatch, algorithm, scoringMethod, anchor](){
        retValue = GraphAlign::runAlgorithm(inputGraphs, 
                                scoreVertices, 
                                kernelType,
                                isSparse,
                                measureType,
                                enforceVLabels,
                                enforceELabels,
                                useAmbiguousEdges,
                                oneMatch,
                                algorithm,
                                scoringMethod,
                                anchor);
        cv.notify_one();
    });

    t.detach(); // let thread run

    {
        std::unique_lock<std::mutex> l(m);
        if(cv.wait_for(l, timeoutWait) == std::cv_status::timeout) 
            throw TimeoutException("Timeout");
    }

    return retValue;

} // _wrapper_runAlgorithm

template <typename Graph, typename SV, typename Vertex>
void testSingleWTimeout(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::filesystem::path &reportPath){
    if(!std::filesystem::exists(graphPath)){
        std::cout << "File path " << graphPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::exists(configPath)){
        std::cout << "File path " << configPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::is_regular_file(configPath)){
        std::cout << "Path " << configPath.string() << " is not a regular file. Terminating.\n";
        return;
    } 

    std::ifstream readConfigFile(configPath);
    json val = GraphAlignIO::readJSON(readConfigFile);
    assert(val.contains("algorithm"));
    assert(val.contains("enforceV"));
    assert(val.contains("enforceE"));
    assert(val.contains("ambiguous"));
    assert(val.contains("score"));
    assert(val.contains("measure"));
    assert(val.contains("timeout"));
    assert(val.contains("useAnchor"));


    std::string algorithm = val["algorithm"];
    std::string score = val["score"];
    bool enforceV = val["enforceV"];
    bool enforceE = val["enforceE"];
    bool ambiguous = val["ambiguous"];
    bool measure = val["measure"];
    bool useAnchor = val["useAnchor"];
    int PGA_TIMEOUT_INT = val["timeout"];
    std::chrono::seconds PGA_TIMEOUT_SECONDS = std::chrono::seconds(PGA_TIMEOUT_INT);
    std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> graphsAndAnchor = GraphAlignIO::readGraphs<Graph>(graphPath);
    std::vector<Graph> inputGraphs = graphsAndAnchor.first;
    std::vector<std::vector<std::vector<Vertex>>> anchors = graphsAndAnchor.second;
    std::vector<std::vector<Vertex>> chosenAnchor;
    std::cout << "Read in " << inputGraphs.size() << " input graphs\n";
    std::cout << "Read in " << anchors.size() << " different anchor mappings \n";
    if(anchors.size() > 0 && useAnchor){
        chosenAnchor = anchors[0];
        std::cout << "Read in " << chosenAnchor.size() << " anchored vertices in first mapping\n";
    }
    std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> retValue;
    
    auto startTime = std::chrono::high_resolution_clock::now();
    bool timeouted = false;
    try{
        retValue = _wrapper_runAlgorithm(PGA_TIMEOUT_SECONDS,
                                        inputGraphs,
                                        scoringVertices,
                                        true, // kernelType
                                        false, // isSparse
                                        measure,
                                        enforceV,
                                        enforceE,
                                        ambiguous,
                                        true, // oneMatch
                                        algorithm,
                                        score,
                                        chosenAnchor);
    } catch (TimeoutException &e){
        // timeout
        timeouted = true;
    }
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count()/1000.0;
    if(timeouted){
        // no |V|, |E| or score.
        std::cout << "Timed out after " << duration << " seconds.\n";
    } else {
        std::cout << "Success after " << duration << " seconds \n";
        std::cout << "Score: " << get<1>(retValue) << "\n";
        GraphAlignIO::writeAlignmentToFile(inputGraphs, get<0>(retValue), get<2>(retValue), reportPath);
    }
} // testSingleWTimeout

template <typename Graph, typename SV, typename Vertex>
void testSingle(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::filesystem::path &reportPath){

    if(!std::filesystem::exists(graphPath)){
        std::cout << "File path " << graphPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::exists(configPath)){
        std::cout << "File path " << configPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::is_regular_file(configPath)){
        std::cout << "Path " << configPath.string() << " is not a regular file. Terminating.\n";
        return;
    } 

    std::ifstream readConfigFile(configPath);
    json val = GraphAlignIO::readJSON(readConfigFile);
    assert(val.contains("algorithm"));
    assert(val.contains("enforceV"));
    assert(val.contains("enforceE"));
    assert(val.contains("ambiguous"));
    assert(val.contains("score"));
    assert(val.contains("measure"));
    assert(val.contains("useAnchor"));

    std::string algorithm = val["algorithm"];
    std::string score = val["score"];
    bool enforceV = val["enforceV"];
    bool enforceE = val["enforceE"];
    bool ambiguous = val["ambiguous"];
    bool measure = val["measure"];
    bool useAnchor = val["useAnchor"];
    std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> graphsAndAnchor = GraphAlignIO::readGraphs<Graph>(graphPath);
    std::vector<Graph> inputGraphs = graphsAndAnchor.first;
    std::vector<std::vector<std::vector<Vertex>>> anchors = graphsAndAnchor.second;
    std::vector<std::vector<Vertex>> chosenAnchor = {};
    if(anchors.size() > 0 && useAnchor){
        chosenAnchor = anchors[0];
    }
    std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>> retValue;
    
    std::streambuf* original_buffer = std::cout.rdbuf();
    std::cout.rdbuf(NULL); // Disable prints from runAlgorithm
    auto startTime = std::chrono::high_resolution_clock::now();
    retValue = GraphAlign::runAlgorithm(inputGraphs,
                                        scoringVertices,
                                        true, // kernelType
                                        false, // isSparse
                                        measure,
                                        enforceV,
                                        enforceE,
                                        ambiguous,
                                        true, // oneMatch
                                        algorithm,
                                        score,
                                        chosenAnchor);
    std::cout.rdbuf(original_buffer);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count()/1000.0;

    Graph alignmentGraph = get<0>(retValue);
    float finalscore = get<1>(retValue);
    std::unordered_map<int, std::unordered_map<Vertex, Vertex>> mapping = get<2>(retValue);
    // V(AG), E(AG), SCORE, TIME
    std::cout << num_vertices(alignmentGraph) << "\t" << num_edges(alignmentGraph) << "\t" << finalscore << "\t" << duration << "\n";
    GraphAlignIO::writeAlignmentToFile(inputGraphs, get<0>(retValue), get<2>(retValue), reportPath);
} // testSingle


template <typename Graph, typename SV, typename Vertex>
void testSingleChooseGuideTree(const std::filesystem::path &graphPath, const std::filesystem::path &configPath, SV scoringVertices, const std::string &treeString, const std::filesystem::path &reportPath){

    if(!std::filesystem::exists(graphPath)){
        std::cout << "File path " << graphPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::exists(configPath)){
        std::cout << "File path " << configPath.string() << " does not exist. Terminating.\n";
        return;
    } 
    if(!std::filesystem::is_regular_file(configPath)){
        std::cout << "Path " << configPath.string() << " is not a regular file. Terminating.\n";
        return;
    } 

    std::ifstream readConfigFile(configPath);
    json val = GraphAlignIO::readJSON(readConfigFile);
    assert(val.contains("algorithm"));
    assert(val.contains("enforceV"));
    assert(val.contains("enforceE"));
    assert(val.contains("ambiguous"));
    assert(val.contains("score"));
    assert(val.contains("useAnchor"));

    std::string algorithm = val["algorithm"];
    std::string score = val["score"];
    bool enforceV = val["enforceV"];
    bool enforceE = val["enforceE"];
    bool ambiguous = val["ambiguous"];
    bool useAnchor = val["useAnchor"];
    std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> graphsAndAnchor = GraphAlignIO::readGraphs<Graph>(graphPath);
    std::vector<Graph> inputGraphs = graphsAndAnchor.first;
    std::vector<std::vector<std::vector<Vertex>>> anchors = graphsAndAnchor.second;
    std::vector<std::vector<Vertex>> chosenAnchor = {};
    if(anchors.size() > 0 && useAnchor){
        chosenAnchor = anchors[0];
    }
    std::vector<std::tuple<Graph, float, std::unordered_map<int, std::unordered_map<Vertex, Vertex>>, std::chrono::milliseconds>> retValue;
    auto startTime = std::chrono::high_resolution_clock::now();
    // No measure or kernel type as they are only used to create an "optimal" guide tree.
    retValue = GraphAlign::runAlgorithmChooseGuideTree(inputGraphs,
                                        scoringVertices,
                                        false, // isSparse
                                        enforceV,
                                        enforceE,
                                        ambiguous,
                                        true, // oneMatch
                                        algorithm,
                                        score,
                                        chosenAnchor,
                                        treeString);
    auto stopTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count()/1000.0;

    int guideTreeNum = 0;
    for(const auto &treeRes : retValue){
        Graph alignmentGraph = get<0>(treeRes);
        float finalscore = get<1>(treeRes);
        std::unordered_map<int, std::unordered_map<Vertex, Vertex>> mapping = get<2>(treeRes);
        std::chrono::milliseconds duration = get<3>(treeRes);
        // V(AG), E(AG), SCORE, TIME
        std::cout << "GT\t(" << guideTreeNum << ")\t";
        guideTreeNum++;
        std::cout << num_vertices(alignmentGraph) << "\t" << num_edges(alignmentGraph) << "\t" << finalscore << "\t" << duration.count()/1000.0 << "\n";
        
    }
    std::cout << "Overall duration (s): " << duration  << std::endl;
    if(!treeString.empty()){
        GraphAlignIO::writeAlignmentToFile(inputGraphs, get<0>(retValue[0]), get<2>(retValue[0]), reportPath);
    }
    

} // testSingleChooseGuideTree

} // namespace GraphAlignTest


#endif