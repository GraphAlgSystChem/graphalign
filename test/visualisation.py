## Module for parsing an alignment report given by the WriteToFile function in
## TestSuite.hpp. The end result is a 2D visualization of the corresponding alignment.

## Code is taken from Marcos Laffitte's progressive alignment github and adapted to our format.
## Adapted by Kasper Halkjær Beider & Tobias Klink Lehn

## >>> python3 ./Visualization.py report_path.txt
import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import json # for report passing
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt, ceil

## GLOBAL DRAWING PARAMETERS
# drawing parameters
minX = -1.15
maxX = 1.15
minY = -1.15
maxY = 1.15
apart = 0.1
# figure-size offset and parameters
wOffset = 7
hOffset = 5
nodesPerInch = 5
graphsPerInch = 5

def visualize1D(graphs, alignmentGraph, projection, name1DPlotsPDF):
    numOfGraphs = len(graphs)
    nodesAlignment = list(alignmentGraph.nodes())
    numOfAlignmentNodes = alignmentGraph.order()
    print(nodesAlignment)
    # create file
    sanityCheckPDF = PdfPages(name1DPlotsPDF)
    # clear current figure
    plt.clf()
    # create figure
    vExtra = max(0, ceil((numOfAlignmentNodes-(wOffset*nodesPerInch))/nodesPerInch))
    gExtra = max(0, ceil((numOfGraphs-(hOffset*graphsPerInch))/graphsPerInch))
    fig = plt.figure(figsize = [wOffset + vExtra, hOffset + vExtra])
    fan = numOfAlignmentNodes
    valK = (1/sqrt(fan)) + apart
    posAlignment = nx.spring_layout(alignmentGraph, center = [0, 0], k = valK)
    # draw final alignment
    nx.draw_networkx(alignmentGraph, with_labels = True,
                     pos = posAlignment, node_size = 150,
                     font_size = 6, node_color = "silver",
                     edge_color = "k", width = 2)
    plt.ylabel("Alignment", fontsize = 14, weight = "light")
    plt.xlim([minX, maxX])
    plt.ylim([minY, maxY])
    plt.tight_layout()
    # save page
    sanityCheckPDF.savefig(fig)
    # clear current figure
    plt.clf()
    # create figure [width, height] in inches
    fig = plt.figure(figsize = [wOffset + vExtra, hOffset + gExtra])
    # build matrix
    nodePresence = dict()
    backMap = dict()
    for g in inputGraphs:
        backMap[g] = {projection[g][u]:u for u in list(g.nodes())}

    # Define label colors
    k = 0
    labelColor = dict()
    for g in graphs:
        for (v, nodeInfo) in g.nodes(data = True):
            nodeInfoTuple = ("label", nodeInfo["label"])
            if(not nodeInfoTuple in list(labelColor.keys())):
                k = k + 1
                labelColor[nodeInfoTuple] = k
    
    allValues = list(labelColor.values())
    if(len(allValues) == 1):
        myCmap = plt.cm.Blues
    else:
        myCmap = plt.cm.turbo
    myCmap.set_bad("white", 1.)
    # add space on border of colormap and map values
    minValue = min(allValues)-0.8
    maxValue = max(allValues)+0.8
    allValues = allValues + [minValue, maxValue]
    allValues.sort()
    myNorm = plt.Normalize(vmin = minValue, vmax = maxValue)
    colorCode = myCmap(myNorm(allValues))

    # Create existence dictionary, existence[v][g] == "-" if it doesn't exist, otherwise it's the vertex of that graph
    existence = {}
    for v in nodesAlignment:
        existence[v] = {}
        for g in inputGraphs:
            existence[v][g] = "-" # initialize to blank
            for u in g.nodes():
                if projection[g][u] == v:
                    existence[v][g] = u

    # Build columns
    for v in nodesAlignment:
        nodePresence[v] = []
        # Build rows
        for g in inputGraphs:
            # Check if the alignment vertex v is mapped to a vertex in g
            if(not existence[v][g] == "-"):
                tempDict = g.nodes[backMap[g][v]]["label"]
                nodeInfoTuple = ("label", g.nodes[backMap[g][v]]["label"])
                nodePresence[v].append(labelColor[nodeInfoTuple])
            else:
                nodePresence[v].append(np.nan)
    # draw heatmap of alignment
    finalArray = [nodePresence[v] for v in nodesAlignment]
    maskedArray = np.ma.array(finalArray, mask = np.isnan(finalArray))
    maskedArray = np.transpose(maskedArray)
    im = plt.imshow(maskedArray, cmap = myCmap, norm = myNorm, aspect = "equal")
    # set minor ticks
    ax = plt.gca()
    # set major tick positions
    ax.set_xticks(np.arange(len(nodesAlignment)))
    ax.set_yticks(np.arange(len(inputGraphs)))
    # set major tick labels
    inputGraph_ticks = [i for i in range(len(inputGraphs))]
    ax.set_xticklabels(nodesAlignment, fontsize = 5.5)
    ax.set_yticklabels(inputGraph_ticks, fontsize = 5.5)
    # set minor ticks
    ax.set_xticks(np.arange(-0.5, len(nodesAlignment), 1), minor = True)
    ax.set_yticks(np.arange(-0.5, len(inputGraphs), 1), minor = True)
    # set grid
    ax.grid(which = "minor", color = "k", linestyle = "-", linewidth = 1)
    # remove minor ticks
    ax.tick_params(which = "minor", bottom = False, left = False)
    # finish image
    plt.title("Matrix Representation of the Alignment\n", fontsize = 12, weight = "light")
    plt.xlabel("\nVertices of the alignment", fontsize = 10, weight = "light")
    plt.ylabel("Input Graphs\n", fontsize = 10, weight = "light")
    # save page
    sanityCheckPDF.savefig(fig)
    # clear current figure
    plt.clf()
    # create figure
    fig = plt.figure(figsize = [wOffset, hOffset])
    # plot color code information as legend
    legendInfo = dict()
    for labelTuple in list(labelColor.keys()):
        strLabels = []
        labelName, labelValue = labelTuple
        strLabels.append("(" + str(labelName) + ", " + str(labelValue) + ")")
        strLabelFinal = ", ".join(strLabels)
        legendInfo[strLabelFinal] = labelColor[labelTuple]
    for eachStr in list(legendInfo.keys()):
        if(eachStr == ""):
            plt.plot([0], marker = "s", linestyle = ":", color = colorCode[legendInfo[eachStr]], label = str(legendInfo[eachStr]) + ":  $\it{unlabeled}$")
        else:
            plt.plot([0], marker = "s", linestyle = ":", color = colorCode[legendInfo[eachStr]], label = str(legendInfo[eachStr]) + ": " + eachStr)
    plt.plot([0], marker = "s", linestyle = ":", color = "w", markersize = 10)
    # make legend box
    plt.legend(loc = "upper left", fontsize = 6, framealpha = 1)
    # set major tick positions
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    # drawing properties
    plt.title("Color code of vertex labels", fontsize = 8, weight = "light")
    plt.tight_layout()
    # save page
    sanityCheckPDF.savefig(fig)
    # save pdf
    sanityCheckPDF.close()
    plt.close()

def visualize2D(graphs, alignmentGraph, projection, name2DPlotsPDF):
    ## INITIALIZE VALUES FOR 2D PLOT
    indices = [i for i in range(len(graphs))]
    numOfGraphs = len(graphs)
    nodesAlignment = list(alignmentGraph.nodes())
    nodesAlignment.sort()
    numOfAlignmentNodes = alignmentGraph.order()
    vExtra = max(0, ceil((numOfAlignmentNodes-(wOffset*nodesPerInch))/nodesPerInch))
    gExtra = max(0, ceil((numOfGraphs-(hOffset*graphsPerInch))/graphsPerInch))
    # get positions for vertices in the alignment
    fan = numOfAlignmentNodes
    valK = (1/sqrt(fan)) + apart
    # *** DEFINE here custom positions for the alignment nodes
    # *** as pairs (a, b) with a and b in the interval [-1, 1],
    # *** e.g. posAlignment[0] = (0.75, -1)
    # *** and comment the following line
    posAlignment = nx.spring_layout(alignmentGraph, center = [0, 0], k = valK)
    k = 0
    nodeInfoTuple = []
    labelColor = dict()
    for g in graphs:
        for (v, nodeInfo) in g.nodes(data = True):
            nodeInfoTuple = ("label", nodeInfo["label"]) # nodeInfo["label"] is the label of the vertex, e.g. 'C', 'O', ...
            if(not nodeInfoTuple in list(labelColor.keys())):
                k = k + 1
                labelColor[nodeInfoTuple] = k
    
    allValues = list(labelColor.values())
    if(len(allValues) == 1):
        myCmap = plt.cm.Blues
    else:
        myCmap = plt.cm.turbo
    myCmap.set_bad("white", 1.)
    # add space on border of colormap and map values
    minValue = min(allValues)-0.8
    maxValue = max(allValues)+0.8
    allValues = allValues + [minValue, maxValue]
    allValues.sort()
    myNorm = plt.Normalize(vmin = minValue, vmax = maxValue)
    colorCode = myCmap(myNorm(allValues))

    sanityCheckPDF = PdfPages(name2DPlotsPDF)
    # draw input graphs
    graph_index = 0
    for g in graphs:
        # get position for vertices relative to alignment
        posLeaf = {v:posAlignment[projection[g][v]] for v in list(g.nodes())}
        # get color for vertices
        nodeList = []
        nodeColor = []
        for (v, nodeInfo) in list(g.nodes(data = True)):
            nodeInfoTuple = ("label", nodeInfo["label"])
            nodeColor.append(colorCode[labelColor[nodeInfoTuple]])
            nodeList.append(v)
        # create figure
        fig = plt.figure(figsize = [wOffset + vExtra, hOffset + vExtra])
        # plot alignment in background
        nx.draw_networkx(alignmentGraph, with_labels = False, pos = posAlignment, node_size = 15,
                         node_color = "lightgrey", edge_color = "lightgrey", width = 0.40)
        # plot input graph
        nx.draw_networkx(g, with_labels = True, labels = projection[g],
                         pos = posLeaf, node_size = 150, nodelist = nodeList, node_color = nodeColor,
                         font_size = 7, edge_color = "k", width = 2)
        plt.ylabel("Input Graph " + str(graph_index), fontsize = 14, weight = "light")
        plt.xlim([minX, maxX])
        plt.ylim([minY, maxY])
        plt.tight_layout()
        # save page
        sanityCheckPDF.savefig(fig)
        # clear current figure
        plt.clf()
        plt.close(fig)
        # save alignmet 3D information
        graph_index += 1
    # create figure
    fig = plt.figure(figsize = [wOffset + vExtra, hOffset + vExtra])
    # plot alignment
    nx.draw_networkx(alignmentGraph, with_labels = True,
                     pos = posAlignment, node_size = 150, node_color = "silver",
                     font_size = 7, edge_color = "k", width = 2)
    plt.ylabel("Alignment", fontsize = 14, weight = "light")
    plt.xlim([minX, maxX])
    plt.ylim([minY, maxY])
    plt.tight_layout()
    # save page
    sanityCheckPDF.savefig(fig)

    plt.clf() # Clear figure for match column 
    fig = plt.figure(figsize = [wOffset + vExtra, hOffset + vExtra])
    
    # plot the underlying alignment graph in the background
    nx.draw_networkx(alignmentGraph, with_labels = False, pos = posAlignment, node_size = 15,
                         node_color = "lightgrey", edge_color = "lightgrey", width = 0.40)
    matchColumns = findMatchColumnNodes(graphs, alignmentGraph, projection) # list of match color nodes
    matchPos = {v:posAlignment[v] for v in matchColumns}
    matchColors = []
    for v in matchColumns:
        matchLabel = alignmentGraph.nodes[v]["label"]
        matchInfoTuple = ("label", matchLabel)
        matchColors.append(colorCode[labelColor[matchInfoTuple]])
    matchGraph = nx.induced_subgraph(alignmentGraph, matchColumns)
    # plot the induced subgraph given by the match columns
    nx.draw_networkx(matchGraph, with_labels= True, nodelist=matchColumns,
                           pos = matchPos, node_size= 150, node_color=matchColors, edge_color = "k", width = 2, font_size = 7
                        )
    plt.ylabel(f"Match columns (of order {len(matchColumns)})", fontsize = 14, weight = "light")
    plt.xlim([minX, maxX])
    plt.ylim([minY, maxY])
    plt.tight_layout()
    sanityCheckPDF.savefig(fig)
    
    # save pdf
    sanityCheckPDF.close()
    plt.close()

def findMatchColumnNodes(graphs, alignmentGraph, projection):
    matchOccur = len(graphs) # number of graphs to be present in for a vertex to be 
    matchVertices = []
    occ_count = {v:0 for v in alignmentGraph.nodes}

    for g in graphs:
        for v in g.nodes:
            alignment_node = projection[g][v]
            occ_count[alignment_node] += 1
    
    for v in occ_count:
        if occ_count[v] == matchOccur:
            matchVertices.append(v)
    
    return matchVertices


if __name__ == "__main__":
    if len(sys.argv) == 2:
        inputFileName = sys.argv[1]
    else:
        print("No file argument given. Please provide a path to a .json file containing the alignment report as the first argument to the program.")
        exit()
    
    

    # settingVertices = False
    # settingEdges = False
    # settingProjection = False

    # Passing inputfile
    with open(inputFileName, "r") as f:

        parsedFile = json.load(f)
        
        inputGraphsArr = parsedFile["inputGraphs"]
        n = len(inputGraphsArr)
        inputGraphs = [0 for i in range(n)]
        projections = {} # G_i -> v -> x in alignment

        print("Number of graphs read in:", n)
        alignmentGraphObj = parsedFile["alignmentGraph"]

        for gObj in inputGraphsArr:
            g = nx.Graph()
            projections[g] = {}
            gVertices = gObj["vertices"]
            gEdges = gObj["edges"]
            gProjection = gObj["projection"]
            gIndex = int(gObj["graph_index"].split("_")[1])

            for v in gVertices:
                vid = v["vid"]
                vlabel = v["label"]
                g.add_node(vid, label=vlabel)

            for e in gEdges:
                source = e["source"]
                target = e["target"]
                elabel = e["label"]
                g.add_edge(source, target, label=elabel)
            
            for p in gProjection:
                vid = p["vid"]
                alignid = p["alignid"]
                projections[g][vid] = alignid
            inputGraphs[gIndex] = g
        
        alignmentGraph = nx.Graph()
        alignVertices = alignmentGraphObj["vertices"]
        alignEdges = alignmentGraphObj["edges"]
        for v in alignVertices:
            vid = v["vid"]
            vlabel = v["label"]
            alignmentGraph.add_node(vid, label=vlabel)
        for e in alignEdges:
            source = e["source"]
            target = e["target"]
            elabel = e["label"]
            alignmentGraph.add_edge(source, target, label=elabel)
    
    outFileName1D = inputFileName.replace(".json", "_1D_plot.pdf")
    outFileName2D = inputFileName.replace(".json", "_2D_plot.pdf")

    # for g in inputGraphs:
    #     print("GRAPH:")
    #     projection = projections[g]
    #     for v in projection:
    #         print(v, projection[v])

    visualize1D(inputGraphs, alignmentGraph, projections, outFileName1D)
    visualize2D(inputGraphs, alignmentGraph, projections, outFileName2D)