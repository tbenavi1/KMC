#include <string>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <vector>

using namespace boost;

int main (int argc, char **argv)
{
  typedef property<vertex_name_t, std::string> VertexProperty;

  typedef adjacency_list <vecS, vecS, undirectedS, VertexProperty> Graph;

  typedef std::pair <std::string, std::string> Edge;
  //Edge edge_array[] = { Edge("aaa", "bbb"), Edge("bbb", "ccc"), Edge("aaa", "xxx"), Edge("bbb", "ccc")};
	std::vector<Edge> edge_array = { Edge("aaa", "bbb"), Edge("bbb", "ccc"), Edge("aaa", "xxx"), Edge("bbb", "ccc")};

  //const char* vertex_array[] = {"aaa", "bbb", "ccc", "xxx"};
	std::vector<std::string> vertex_array = {"aaa", "bbb", "ccc", "xxx"};
  std::map<std::string, Graph::vertex_descriptor> indexes;

  //const int nb_vertices = sizeof(vertex_array)/sizeof(vertex_array[0]);
	const int nb_vertices = vertex_array.size();

  // creates a graph with 4 vertices
  Graph g (nb_vertices); 

  // fills the property 'vertex_name_t' of the vertices
  for(int i = 0; i < nb_vertices; i++)
  {
    boost::put(vertex_name_t(), g, i, vertex_array[i]); // set the property of a vertex
    indexes[vertex_array[i]] = boost::vertex(i, g);     // retrives the associated vertex descriptor
  }

  // adds the edges
  // indexes[edges[0].first] maps "aaa" to the associated vertex index
  //for(int i = 0; i < sizeof(edge_array)/sizeof(edge_array[0]); i++)
	for(int i = 0; i < edge_array.size(); i++)
  {
    boost::add_edge(indexes[edge_array[i].first], indexes[edge_array[i].second], g);
  }
  //typedef property_map<Graph, vertex_index_t>::type IndexMap;
  //IndexMap index = get(vertex_index, g);
	typedef property_map<Graph, vertex_name_t>::type NameMap;
	NameMap name = get(vertex_name, g);

  std::cout << "vertices(g) = ";
  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  for (vp = vertices(g); vp.first != vp.second; ++vp.first)
    std::cout << name[*vp.first] <<  " ";
  std::cout << std::endl;
	std::cout << "edges(g) = ";
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	{
		std::cout << "(" << name[source(*ei, g)] << "," << name[target(*ei, g)] << ") ";
	}
	std::cout << std::endl;
  return 0;
}
