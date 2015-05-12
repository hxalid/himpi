//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/graphviz.hpp>
#include <iostream>
#include <fstream>


int main()
{
  using namespace boost;
	typedef adjacency_list < vecS, vecS, undirectedS, property<vertex_name_t, std::string, property<vertex_color_t, std::string > >, property < edge_weight_t, double > > Graph;
	typedef graph_traits < Graph >::edge_descriptor Edge;
	typedef graph_traits < Graph >::vertex_descriptor Vertex;
	typedef std::pair<int, int> E;

	dynamic_properties dp;

	Graph graph(0);
	std::ifstream ifs;
	ifs.open("merged2.dot");
	if (!ifs) {fprintf(stderr, "Can't open DOT file\n"); exit(-1);}


	property_map<Graph, vertex_name_t>::type name = get(vertex_name, graph);
	dp.property("node_id",name);
	property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);
	dp.property("weight", weight);
	property_map<Graph, vertex_color_t>::type color = get(vertex_color, graph);
	dp.property("fillcolor", color);
	bool status = read_graphviz(ifs,graph,dp);
	ifs.close();

/*
  graph_traits<Graph>::edge_iterator eiter, eiter_end;
  for (boost::tie(eiter, eiter_end) = edges(graph); eiter != eiter_end; ++eiter) {
	std::cout << get(edge_weight,graph, *eiter) << std::endl;
  }
*/
  std::vector < Edge > spanning_tree;

  kruskal_minimum_spanning_tree(graph, std::back_inserter(spanning_tree));

  //std::ofstream fout("figs/kruskal-eg.dot");
  std::cout << "graph A {\n"
    << " rankdir=LR\n"
    << " size=\"3,3\"\n"
    << " ratio=\"filled\"\n"
    << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";
  graph_traits<Graph>::edge_iterator eiter, eiter_end;
  for (boost::tie(eiter, eiter_end) = edges(graph); eiter != eiter_end; ++eiter) {
    //std::cout << source(*eiter, graph) << " -- " << target(*eiter, graph);
    if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
        != spanning_tree.end()) {
		std::cout << "\"" <<  get(vertex_name, graph, source(*eiter, graph)) << "\" -- \"" << get(vertex_name, graph, target(*eiter, graph)) << "\"";
		std::cout << "[color=\"black\", label=\"" << get(edge_weight, graph, *eiter) << "\"];\n";
	}
    else {
      //std::cout << "[color=\"gray\", label=\"" << get(edge_weight, graph, *eiter)
       //    << "\"];\n";
	}
  }
  std::cout << "}\n";
  return EXIT_SUCCESS;
}
