//
// Created by makem on 17/01/19.
//

#include "PathGraph.h"
#include "CTiglError.h"
#include <algorithm>

tigl::PathGraph::PathGraph() {}

void tigl::PathGraph::set(std::vector<std::string> newGraph) {
  clear();
  // check if every nde is unique
  std::vector<std::string>::iterator it =
      std::unique(newGraph.begin(), newGraph.end());
  bool wasUnique = (it == newGraph.end());
  if (!wasUnique) {
    throw CTiglError("PathGraph::set(): The given graph has duplicated nodes.");
  }
  graph = newGraph;
}

void tigl::PathGraph::clear() { graph.clear(); }

std::vector<std::string>
tigl::PathGraph::getElementsInBetween(std::string uid1, std::string uid2) {

  std::vector<std::string> elementsBetween; // contain the start and the end

  if (!(contains(uid1) && contains(uid2))) {
    return elementsBetween;
  }

  bool afterStart = false;
  bool afterEnd = false;

  for (int i = 0; i < graph.size(); i++) {
    // cover the special case where the end and the start are the same
    if (graph[i] == uid1 && graph[i] == uid2) {
      afterStart = true;
      afterEnd = true;
      elementsBetween.push_back(graph[i]);
    }
    // we allow the start and end to be inverted
    else if (afterStart == false && (graph[i] == uid1 || graph[i] == uid2)) {
      afterStart = true;
      elementsBetween.push_back(graph[i]);
    } else if (afterStart == true && (graph[i] == uid1 || graph[i] == uid2)) {
      afterEnd = true;
      elementsBetween.push_back(graph[i]);
    } else if (afterStart == true && afterEnd == false) {
      elementsBetween.push_back(graph[i]);
    }
  }
  return elementsBetween;
}

std::vector<std::string> tigl::PathGraph::getElementsBefore(std::string uid) {

  std::vector<std::string> elementsBefore;

  if (!contains(uid)) {
    return elementsBefore;
  }

  bool after = false;

  for (int i = 0; i < graph.size(); i++) {
    if (graph[i] == uid) {
      after = true;
    }
    if (after == false) {
      elementsBefore.push_back(graph[i]);
    }
  }

  return elementsBefore;
}

std::vector<std::string> tigl::PathGraph::getElementsAfter(std::string uid) {
  std::vector<std::string> elementsAfter;

  bool after = false;

  std::vector<std::string>::iterator it;
  for (int i = 0; i < graph.size(); i++) {
    if (after) {
      elementsAfter.push_back(graph[i]);
    }
    if (graph[i] == uid) {
      after = true;
    }
  }
  return elementsAfter;
}

std::string tigl::PathGraph::getRoot() {

  if (!isEmpty()) {
    return graph[0];
  } else {
    throw CTiglError("PathGraph::getRoot called but graph is empty.");
  }
}

std::string tigl::PathGraph::getLeaf() {
  if (!isEmpty()) {
    return graph[graph.size() - 1];
  } else {
    throw CTiglError("PathGraph::getLeaf called but graph is empty.");
  }
}

bool tigl::PathGraph::contains(std::string uid) {

  if (std::find(graph.begin(), graph.end(), uid) != graph.end()) {
    return true;
  } else {
    return false;
  }
}
