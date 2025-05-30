#ifndef GRAPHIO_H
#define GRAPHIO_H

// std
#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>

// 3rd party
#include <boost/graph/adjacency_list.hpp> // graph traits
#include "../src/external/json.hpp" // for json parsing, header only by Niels Lohmann


namespace GraphAlignIO {

	// for ease of use
    using json = nlohmann::json;

/**
 * @brief Reads a .json file with the format:
 *  { "inputGraphs" : 
 *      "graph_i" : {
 *          "vertices" : [
 *                      { 
 *                        "vid" : 
 *                        "label" : 
 *                      },
 *                      ]
 *          "edges" : [
 *                      { 
 *                        "source" :
 *                        "target" : 
 *                        "label" : 
 *                      },
 *                      ]
 *      },
 *  }
 * Labels may be empty.
 * 
 * @tparam Graph Graph type, Vertex descriptor type
 * @param filePath The path to the .json file containing the graphs to parse.
 * @return Vector of input graphs and the set of anchor mappings. E.g. A[i][j][k] is the j'th anchor vertex in graph k from the i'th anchor mapping.
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> readGraphs(const std::filesystem::path &filePath);

/**
 * Writes result from runAlgorithm to file for easy parsing into Python
 * @param inputGraphs Graphs given as input to the algorithm
 * @param alignmentGraph Final alignment graph
 * @param projection Map that, for each graph i, maps each vertex u to a vertex in alignmentGraph.
 * @param filePath The file to be written to.
 */

/**
 * @brief Writes result from runAlgorithm to file for easy parsing into Python
 * 
 * @tparam Graph The graph type
 * @tparam Vertex The vertex type, defaults to the vertex_descriptor for the graph type
 * @param inputGraphs Graphs given as input to the algorithm
 * @param alignmentGraph Final alignment graph
 * @param projection Map that, for each graph i, maps each vertex u to a vertex in alignmentGraph.
 * @param filePath The file to be written to.
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor>
void writeAlignmentToFile(const std::vector<Graph> &inputGraphs, const Graph &alignmentGraph, const std::unordered_map<int, std::unordered_map<Vertex, Vertex>> &projection, const std::string &filePath);

/** 
 * Taken from https://www.boost.org/doc/libs/1_84_0/libs/json/doc/html/json/input_output.html#json.input_output.parsing 
 * @brief Parses an input stream as an expected .json file and returns the corresponding json value.
 * 
 * @param is The inputstream from which a .json file is expected
 * @param ec An error code to be completed with information should the parsing fail
 * 
 * @return The json object corresponding to the content of the inputstream.
 */
json readJSON(std::istream& is);



template<typename Graph, typename Vertex>
std::pair<std::vector<Graph>, std::vector<std::vector<std::vector<Vertex>>>> readGraphs(const std::filesystem::path &filePath) {
    
    std::vector<std::vector<std::vector<Vertex>>> anchorMappings;
    std::vector<std::vector<Vertex>> legacyAnchoredVertices;
    std::vector<Graph> graphs;
    std::ifstream readFile(filePath);
    json obj = readJSON(readFile);
    assert(obj.size() >= 1);
    assert(obj.contains("inputGraphs"));

    // Adding each graph separately
    for(const json &graphObj : obj["inputGraphs"]){
        Graph currentGraph;
        int previousId = -1; // For checking sorted order of the vertices. Fails if they are not sorted
        // Adding vertices
        assert(graphObj.contains("vertices"));
        json verticesObj = graphObj["vertices"];
        assert(verticesObj.is_array());
        for(const auto &vertex : verticesObj){
            std::string label = vertex["label"].template get<std::string>();
            
            // Converting the string representation to integer representation.
            std::string vidString = vertex["vid"].template get<std::string>();
            std::stringstream vidStream;
            vidStream << vidString;
            int vidInt;
            vidStream >> vidInt;
            assert(vidInt == previousId + 1);
            auto v = add_vertex(currentGraph);
            currentGraph[v].label = label; // Under assumption that vertices come in sorted order
            previousId = vidInt;
        }

        // Adding edges
        assert(graphObj.contains("edges"));
        json edgesObj = graphObj["edges"];
        assert(edgesObj.is_array());
        for(const auto &edge : edgesObj){
            std::string u = edge["source"].template get<std::string>();
            std::string v = edge["target"].template get<std::string>();

            // Converting the string representation to integer representation
            std::stringstream uStream;
            std::stringstream vStream;
            uStream << u;
            vStream << v;
            int source, target;
            uStream >> source;
            vStream >> target;
            std::string label = edge["label"].template get<std::string>();
            auto e = add_edge(source, target, currentGraph);
            currentGraph[e.first].label = label;
        }
        graphs.push_back(currentGraph);
    }

    // Parse the list of anchor mappings embedded in the outer JSON object.
    if(obj.contains("anchor_mappings")){
        json anchorArray = obj["anchor_mappings"]; // the set of all anchor mappings
        assert(anchorArray.is_array());
        for(const auto &mapping : anchorArray){
            std::vector<std::vector<Vertex>> convertedMap;
            for(const auto &bag : mapping){ // bag of vertices mapped to each other
                std::vector<Vertex> vertexBag;
                for(const auto &uString : bag){
                    std::string u = uString.template get<std::string>();
                    std::stringstream uStream;
                    uStream << u;
                    int source;
                    uStream >> source;
                    vertexBag.push_back(source);
                }
                convertedMap.push_back(vertexBag);
            }
            anchorMappings.push_back(convertedMap);
        }
    }
    readFile.close();

    // Map each label to a labelIndex and assign this as a property to the vertices in each graph
    std::map<std::string, int> labelMap;
    int labelCounter = 0;
    for(Graph &graph : graphs) {
        for(Vertex x : GraphAlign::asRange(vertices(graph))) {
            if(!labelMap.contains(graph[x].label)) {
                labelMap[graph[x].label] = labelCounter;
                labelCounter++;
            }
        }
    }
    for(Graph &graph : graphs) {
        for(Vertex x : GraphAlign::asRange(vertices(graph))) {
            graph[x].labelIndex = labelMap[graph[x].label];
        }
    }

    return {graphs, anchorMappings};
} // readGraphs

template<typename Graph, typename Vertex>
void writeAlignmentToFile(const std::vector<Graph> &inputGraphs, const Graph &alignmentGraph, const std::unordered_map<int, std::unordered_map<Vertex, Vertex>> &projection, const std::string &filePath){

    std::cout << "Writing alignment to file " << filePath << "\n";
    std::ofstream writeFile(filePath);
    json obj;
    obj["inputGraphs"] = json::array();

    for(size_t i = 0; i < inputGraphs.size(); i++){
        Graph g = inputGraphs[i];
        json graphObj;
        json vertexArray = json::array();
        for(Vertex v : GraphAlign::asRange(vertices(g))){
            json vertexObj;
            vertexObj["vid"] = v;
            vertexObj["label"] = g[v].label;
            vertexArray.push_back(vertexObj);
        }
        graphObj["vertices"] = vertexArray;
        json edgeArray = json::array();
        for(auto e : GraphAlign::asRange(edges(g))){
            json edgeObj;
            edgeObj["source"] = source(e, g);
            edgeObj["target"] = target(e, g);
            edgeObj["label"] = g[e].label;
            edgeArray.push_back(edgeObj);
        }
        graphObj["edges"] = edgeArray;
        json projectionArray = json::array();
        for(Vertex v : GraphAlign::asRange(vertices(g))){
            json projectionObj;
            projectionObj["vid"] = v;
            projectionObj["alignid"] = projection.at(i).at(v);
            projectionArray.push_back(projectionObj);
            // writeFile << v << " " << projection.at(i).at(v) << "\n";
        }
        graphObj["projection"] = projectionArray;
        std::string graphIndex = "graph_" + std::to_string(i);
        
        // obj["inputGraphs"][graphIndex] = graphObj;
        graphObj["graph_index"] = graphIndex;
        obj["inputGraphs"].push_back(graphObj);
    }

    json alignObj;
    json alignVertexArray = json::array();
    for(Vertex v : GraphAlign::asRange(vertices(alignmentGraph))){
        json vertexObj;
        vertexObj["vid"] = v;
        vertexObj["label"] = alignmentGraph[v].label;
        alignVertexArray.push_back(vertexObj);
    }
    alignObj["vertices"] = alignVertexArray;
    json alignEdgeArray = json::array();
    for(auto e : GraphAlign::asRange(edges(alignmentGraph))){
        json edgeObj;
        edgeObj["source"] = source(e, alignmentGraph);
        edgeObj["target"] = target(e, alignmentGraph);
        edgeObj["label"] = alignmentGraph[e].label;
        alignEdgeArray.push_back(edgeObj);
    }
    alignObj["edges"] = alignEdgeArray;

    obj["alignmentGraph"] = alignObj;
    writeFile << std::setw(4) << obj;
    writeFile.close();
}

json readJSON(std::istream& is){
    json obj;
    obj = json::parse(is);
    return obj;
} // read_json

} // namespace GraphAlignIO



#endif // GRAPHIO_H