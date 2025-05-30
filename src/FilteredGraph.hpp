#ifndef FILTEREDGRAPH_H
#define FILTEREDGRAPH_H

// 3rd party
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>

template<typename Vertex>
int getIndex(Vertex v) {
    return v;
}

// us
#include "PairToRangeAdaptor.hpp"

namespace GraphAlign {

/**
 * Class for representing a graph in which some vertices might have been "deleted".
 * That is, in the underlying graph structures, vertices really are not deleted, however,
 * 
 * Provides iteration over vertices, edges and adjacent vertices.
 */
template<typename Graph, typename Vertex = boost::graph_traits<Graph>::vertex_descriptor, typename Edge = boost::graph_traits<Graph>::edge_descriptor>
struct FilteredGraph{

	// JLA BEGIN (FilteredWrapper.hpp)
	using base_traits = boost::graph_traits<Graph>;
	// Graph
	using vertex_descriptor = typename base_traits::vertex_descriptor;
	using edge_descriptor = typename base_traits::edge_descriptor;
	using directed_category = typename base_traits::directed_category;
	using edge_parallel_category = typename base_traits::edge_parallel_category;
	using traversal_category = typename base_traits::traversal_category;
	// JLA END

public:
	
	using boostVertexIterator = boost::graph_traits<Graph>::vertex_iterator;
	using boostEdgeIterator = boost::graph_traits<Graph>::edge_iterator;
	using boostAdjacencyIterator = boost::graph_traits<Graph>::adjacency_iterator;
	
	/**
	 * Class for facading the vertex iterator from boost s.t.
	 * "deleted" vertices are ignored.
	 */
	struct vertexIterator : boost::iterator_facade<
			vertexIterator, // we adapt CRTP (CRTP WHAT IS THIS???)
			Vertex, // we want to return vertex descriptors
			std::random_access_iterator_tag // random access inherited from vertex iterator
			>
		{
			public:
				vertexIterator() = default;
				vertexIterator(const FilteredGraph &f) : begin(vertices(f.g).first), current(vertices(f.g).first), end(vertices(f.g).second), f(f) {}
				vertexIterator(const FilteredGraph &f, 
							   boostVertexIterator begin, 
							   boostVertexIterator current, 
							   boostVertexIterator end) : begin(begin), end(end), f(f) {
					while(current != end && !f.membership[(getIndex(*current))]){
						current++;
					}
					this->current = current;
				}

				Vertex operator*() const{
					return *current;
				}

				vertexIterator &operator++(){
					current++;
					while( current != end && !f.membership[getIndex(*current)]){
						current++;
					}
					return *this;
				}

				vertexIterator &operator--(){
					if(current == begin){
						current = end;
					} else{
						current--;
						while(!f.membership[getIndex(*current)]){
							if(current == begin){
								current = end;
								break;
							}
							current--;
						}
					}
					return *this;
				}

				// postfix here
				vertexIterator &operator++(int){
					vertexIterator tmp = this;
					++this;
					return tmp;
				}

				vertexIterator &operator--(int){
					vertexIterator tmp = this;
					--this;
					return tmp;
				}

				friend bool operator==(const vertexIterator &lhs, const vertexIterator &rhs){
					return (lhs.current == rhs.current);
				}

				friend bool operator!=(const vertexIterator &lhs, const vertexIterator &rhs){
					return (lhs.current != rhs.current);
				}
			private:
				friend class boost::iterator_core_access;
				boostVertexIterator begin;
				boostVertexIterator current;
				boostVertexIterator end;
				const FilteredGraph &f;

		};

	friend class vertexIterator; // to access membership
	
	// Constructor
	FilteredGraph(const Graph &g) : g(g) {
		membership.reserve(num_vertices(g));
		for(std::size_t i = 0; i < num_vertices(g); i++) {
			membership.push_back(true);
		}
		num_visible = num_vertices(g);
		cutted = g;
	}

	// Destructor
    ~FilteredGraph() = default;

	FilteredGraph(const FilteredGraph &other){
		FilteredGraph(other.g);
	}

	FilteredGraph &operator=(const FilteredGraph &other){
		FilteredGraph(other.g);
	}

	// Remove a single vertex
	void remove(const Vertex &v){
		if(membership[getIndex(v)]){
			clear_vertex(v, cutted);
			num_visible--;
		} 
		membership[getIndex(v)] = false;
		
	}

	// Remove a set of vertices
	void remove(const std::vector<Vertex> &vs){
		for(Vertex v : vs){
			remove(v);
		}
	}

	void insert(const Vertex &v){
		if(!membership[getIndex(v)]) {
			for(const auto u : asRange(adjacent_vertices(v, g))){
				bool edge_exists = edge(u, v, cutted).second;
				if (!edge_exists){ // do not add duplicate edges between neighbours
					auto e = add_edge(u, v, cutted);
					cutted[e.first].label = g[edge(u, v, g).first].label;
				}
			}
			num_visible++; 
		}
		membership[getIndex(v)] = true;
	}

	void insert(const std::vector<Vertex> &vs){
		for(Vertex v : vs){
			insert(v);
		}
	}

	// === VERTICES
	// The iterator pointing to the first vertex of G
	vertexIterator verticesBegin() const{
		return vertexIterator(*this, vertices(g).first, vertices(g).first, vertices(g).second);
	}

	// The iterator pointing to the last vertex of G
	vertexIterator verticesEnd() const{
		return vertexIterator(*this, vertices(g).first, vertices(g).second, vertices(g).second);
	}

	friend std::pair<vertexIterator, vertexIterator> vertices(const FilteredGraph &f){
		return {f.verticesBegin(), f.verticesEnd()};
	}

	// === EDGES
	friend std::pair<Edge, bool> edge(Vertex u, Vertex v, const FilteredGraph &f){
		return edge(u, v, f.cutted);
	} 

	friend std::pair<boostAdjacencyIterator, boostAdjacencyIterator> adjacent_vertices(Vertex v, const FilteredGraph &f){
		return adjacent_vertices(v, f.cutted);
	}

	// === Other
	friend size_t num_visible_vertices(const FilteredGraph &f){
		return f.num_visible;
	}

	friend size_t num_vertices(const FilteredGraph &f){
		return num_vertices(f.g);
	}

	friend bool isVisible(Vertex v, const FilteredGraph &f){
		return f.membership[getIndex(v)];
	}

	auto &operator[](Vertex v) const{
		assert(membership[v]);
		return g[v];
	}

	auto &operator[](Edge e) const{
		return g[e];
	}

private:
	const Graph &g;
	Graph cutted;
	size_t num_visible;
	std::vector<bool> membership;
};


}

#endif