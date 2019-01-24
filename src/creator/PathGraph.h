//
// Created by makem on 17/01/19.
//

#ifndef TIGL_PATHGRAPH_H
#define TIGL_PATHGRAPH_H

#include <string>
#include <vector>

namespace tigl {

/**
 * @brief This class represent path graph object.
 * By path graph, we mean that the graph can be represented as simple directed line.
 *
 *
 */
class PathGraph {

public:
  /**
   * Create a empty path graph.
   */
  PathGraph();

  /**
   * Set the path for this graph.
   * The input represent the path as a vector of string.
   * The first element of the vector is the root,
   * the second the second node and so long.
   * We assume that each node is connected to the previous one by a edge.
   * @remark Each string of the vector should be unique.
   * @param newGraph the vector of sting
   */
  void set(std::vector<std::string> newGraph);

  /**
   * Reset the graph to is empty state.
   */
  virtual void clear();

  inline bool isEmpty() { return graph.size() == 0; };

  /**
   * @return the graph as a vector of string
   */
  inline std::vector<std::string> getGraphAsVector() { return graph; };

  /**
   * Returns the node between node 1 and node 2.
   * The order of the given nodes can be switched.
   * @remark: The two given nodes are also contained in the list.
   * @remark  If the graph is empty or it don't contain one of the node, it will
   * return a empty list
   * @param startNode
   * @param endNode
   * @return a vector of string
   */
  std::vector<std::string> getElementsInBetween(std::string uid1,
                                                std::string uid2);

  /**
   * Returns the nodes that are before the given node.
   * @remark The given node is not contain in the returned vector.
   * @remark If the graph is empty or it don't contain the given UID, it will
   * return a empty list
   * @param node
   * @return
   */
  std::vector<std::string> getElementsBefore(std::string uid);

  /**
   * Returns the nodes that are after the given node.
   * @remark The uid is not contains in the list.
   * @remark  If the graph is empty or it don't contain the given node, it will
   * return a empty list
   * @param uid
   * @return
   */
  std::vector<std::string> getElementsAfter(std::string uid);

  /**
   * Returns the first element of the graph.
   * @remark If the graph is empty, it will throw a exception
   * @return
   */
  std::string getRoot();

  /**
   * Return the last element of the graph path
   * @remark If the graph is empty, it will throw a exception
   * @return
   */
  std::string getLeaf();

  /**
   * @return  true if the node is in the graph, false otherwise.
   */
  bool contains(std::string uid);

  // todo implement [] operator

private:
  std::vector<std::string> graph;
};

} // namespace tigl
#endif // TIGL_PATHGRAPH_H
