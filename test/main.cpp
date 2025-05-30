
// std
#include <string>
#include <iostream>
#include <filesystem>

// 3rd party

// us
#include "../src/MGA_Analysis.hpp" // runAlgorithm
#include "../src/AlignUtility.hpp" // printResults
#include "TestSuite.hpp" // recursiveDescent 

using namespace GraphAlign;
using namespace GraphAlignTest;

struct VertexLabel {
    std::string label;
    int labelIndex;
};
struct EdgeLabel {
    std::string label;
};
using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexLabel, EdgeLabel>;
using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

float scoringVertices(const Vertex v1, const Vertex v2,
                    const Graph &g1, const Graph &g2 ) {
    
    if(v1 != boost::graph_traits<Graph>::null_vertex() && v2 != boost::graph_traits<Graph>::null_vertex()) {
        // Same label
        if(g1[v1].label == g2[v2].label) {
            return 1.0;
        }
        // mismatch labels
        else {
            return 0.1;
        }
    }
    // gap 
    if(v1 == boost::graph_traits<Graph>::null_vertex() || v2 == boost::graph_traits<Graph>::null_vertex()) {
        return 0.0;
    }
    // Double gap
    if(v1 == boost::graph_traits<Graph>::null_vertex() && v2 == boost::graph_traits<Graph>::null_vertex()) {
        return 0.0;
    }
    
    return 0.0;
}

enum Mode {SINGLE, SINGLETIMEOUT, UNSET, GUIDETREES};

/**
 * Call with a file path and a config path specifying parameters for a single execution.
 * OR
 * Call with a directory path and a config path specifying all wished combinations for execution.
 */
int main(int argc, char** argv){
    std::string path = "";
    std::string config = "";
    std::string report = "";
    std::string treeString = "";
    enum Mode m = UNSET;
    if(argc >= 4){
        if(strcmp(argv[1], "-s") == 0){
            m = SINGLE; 
        } else if (strcmp(argv[1], "-st") == 0){
            m = SINGLETIMEOUT;
        } else if (strcmp(argv[1], "-g") == 0){
            m = GUIDETREES;
        }
        path = std::string(argv[2]);
        config = std::string(argv[3]);
    if(argc == 5){
        if(m == SINGLE || m == SINGLETIMEOUT){
            report = std::string(argv[4]);
        } else if(m == GUIDETREES){
            treeString = std::string(argv[4]);
        }
    }
    } else {
        std::cout << "Wrong number of arguments given.\n";
        std::cout << "./test -MODE {FILE, DIR}_PATH CONFIG_PATH.json\n";  
        std::cout << "MODE:{-s, -st, -g}\n";
        std::cout << "./test -s FILE_PATH CONFIG_PATH.json OPTIONAL_REPORT_PATH.json\n";
        std::cout << "./test -st FILE_PATH CONFIG_PATH.json OPTIONAL_REPORT_FOLDER_PATH \n";
        std::cout << "./test -g FILE_PATH CONFIG_PATH.json OPTIONAL_GUIDETREE_STRING \n";
        std::cout << "-s:\tSingle instance w/o timeout, prints alignment statistics to std out (vertices, edges, score). Meant for python/bash scripting.\n";
        std::cout << "-st:\tSingle instance with timeout specified in CONFIG_PATH.json. Writes alignment to FILE_PATH_report.json and prints score to std out.\n";
        return 0;
    }

    std::string reportFileName;
    std::filesystem::path filePath = std::filesystem::path(path);
    std::filesystem::path configPath = std::filesystem::path(config);
    std::filesystem::path reportPath = filePath; // For single report w. timeout, produce filename_report.json
    std::filesystem::path logReportPath = report; // Output single report to this folder

    std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> res;
    std::vector<Graph> graphs;
    
    switch(m){
        case SINGLETIMEOUT:
            std::cout << "RUNNING SINGLE INSTANCE W. TIMEOUT ON FILE " << filePath.string() <<"\n";
            reportPath.replace_extension(""); // remove .json
            reportFileName = reportPath.filename().string() + "_report"; 
            reportPath.replace_filename(reportFileName); // append _report
            reportPath.replace_extension(".json"); // append .json
            logReportPath /= reportPath.filename();
            testSingleWTimeout<Graph>(filePath, config, scoringVertices, logReportPath);
            break;
        
        case SINGLE:
            if(report.compare("") == 0){
                testSingle<Graph>(filePath, config, scoringVertices);
            } else {
                testSingle<Graph>(filePath, config, scoringVertices, report);
            }
            break;
        
        case GUIDETREES:
            reportPath.replace_extension(""); // remove .json
            reportFileName = reportPath.filename().string() + "_report"; 
            reportPath.replace_filename(reportFileName); // append _report
            reportPath.replace_extension(".json"); // append .json
            logReportPath /= reportPath.filename();
            testSingleChooseGuideTree<Graph>(filePath, config, scoringVertices, treeString, logReportPath);
            break;

        default:
            std::cout << "No mode of operation was given.\n";
            std::cout << "./test -MODE ...\n";  
            return 0;    
    }

}