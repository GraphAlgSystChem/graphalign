# src

This folder contains the source code for the implementation to align multiple graphs in a progressive manner. A brief overview of the significant files is provided below. See each individual file for the documentation for each function. 

* **MGA_Analysis** is the entry point for the program. Either by providing a specific guide tree or  to execute a guide tree is constructed 

* **Expand** the first pairwise alignment algorithm. An algorithm that finds maximal common subgraphs that maximise the scoring function. It is a backtracking algorithm.

* **Trim** The second pairwise alignment algorithm. Functions by trying to remove subsets of vertices from the smaller graph and then use a subroutine inspired by VF algorithms to determine subgraph isomorphism.

* **BranchBound** contains the code for the branch and bound heuristic that is used to attempt to limit the search trees to improve performance. 

* **ComputeOrder** computes the vertex order.

* **BuildAlignment** builds the alignment graph based on a mapping between two (alignment) graphs that is returned from either **expand** or **trim**. 

* **GenerateCandidates** generates the candidate pairs for extending matches. 

* **GuideTree** computes the guide tree for the progressive alignment. 

* **Scoring** functions to score either complete alignments or a pair of vertices.

* **SemanticSyntacticCheck** carries out the feasibility checks, both the semantic check in terms of labels and the syntactic check in terms of satisfying the adjacency constriants. 

* **GlobalVariables** is a struct containing information used through large parts of the program that is passed around in this object in a reference format, to avoid having specifi
