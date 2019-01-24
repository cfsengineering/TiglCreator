//
// Created by makem on 17/01/19.
//

#include "PathGraph.h"
#include "CTiglError.h"


void tigl::PathGraph::set(std::vector<std::string> newGraph) {
    graph = newGraph;
}

void tigl::PathGraph::clear() {
    graph.clear();
}

std::vector<std::string> tigl::PathGraph::getElementsInBetween(std::string uid1, std::string uid2) {

    std::vector<std::string> elementsBetween; // contain the start and the end

    bool afterStart = false;
    bool afterEnd   = false;

    for (int i = 0; i < graph.size(); i++) {
        if(  graph[i] == uid1  && graph[i] == uid2 ){ // special case where the end and the start are the same
            afterStart = true;
            afterEnd = true;
            elementsBetween.push_back(graph[i]);
        }
        else if (afterStart == false && ( graph[i] == uid1  || graph[i] == uid2 ) ) { // we allow the start and end to be inverted
            afterStart = true;
            elementsBetween.push_back(graph[i]);
        }
        else if (afterStart == true && ( graph[i] == uid1  || graph[i] == uid2 ) ) {
            afterEnd = true;
            elementsBetween.push_back(graph[i]);
        }
        else if (afterStart == true && afterEnd == false) {
            elementsBetween.push_back(graph[i]);
        }
    }

    // the case where only one element where present in the list
    if (afterEnd != true ){
        elementsBetween.clear();
    }

    return elementsBetween;


}


std::vector<std::string> tigl::PathGraph::getElementsBefore(std::string uid) {

    std::vector<std::string> elementsBefore;

    bool after = false;


    for (int i = 0; i < graph.size(); i++) {
        if(graph[i] == uid){
            after = true;
        }
        if( after == false){
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
        if( after ){
            elementsAfter.push_back(graph[i]);
        }
        if(graph[i] == uid){
            after = true;
        }

    }
    return elementsAfter;
}

std::string tigl::PathGraph::getRoot() {

    if( ! isEmpty())
    {
        return graph[0];
    }
    else {
        throw CTiglError("PathGraph::getRoot called but graph is empty.");
    }
}

std::string tigl::PathGraph::getLeaf() {
    if( ! isEmpty())
    {
        return graph[graph.size() - 1];
    }
    else {
        throw CTiglError("PathGraph::getLeaf called but graph is empty.");
    }
}

