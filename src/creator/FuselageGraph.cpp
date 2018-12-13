//
// Created by cfse on 12/13/18.
//

#include "FuselageGraph.h"
#include "CCPACSFuselage.h"
#include "CCPACSFuselageSegment.h"

tigl::FuselageGraph::FuselageGraph() {
    fuselage = nullptr;
}

void tigl::FuselageGraph::build(tigl::CCPACSFuselage *inFuselage) {

    clear();

    if(inFuselage != nullptr) {

        fuselage = inFuselage;
        createGraph();
        centers = fuselage->GetElementsCenters();
        createSimpleGraph();
        circumferences = fuselage->GetCircumferenceOfElements();
    }
}

void tigl::FuselageGraph::createGraph() {

    graph.clear();

    CCPACSFuselageSegments& segments = fuselage->GetSegments();
    CCPACSFuselageSections& sections = fuselage->GetSections();

    std::string tempFromElmUID, tempToElmUID;

    for (int i = 1; i <= segments.GetSegmentCount(); i++) {
        CCPACSFuselageSegment& seg = segments.GetSegment(i);
        tempFromElmUID             = seg.GetStartSectionElementUID();
        tempToElmUID               = seg.GetEndSectionElementUID();

        if (graph.count(tempFromElmUID) == 0) {
            graph[tempFromElmUID] = std::vector<std::string>();
            graph[tempFromElmUID].push_back(tempToElmUID);
        }
        else {
            // check if the element exist
            if (std::find(graph[tempFromElmUID].begin(), graph[tempFromElmUID].end(), tempToElmUID) ==  graph[tempFromElmUID].end()) {
                graph[tempFromElmUID].push_back(tempToElmUID);
            }
        }
        if (graph.count(tempToElmUID) == 0) {
            graph[tempToElmUID] = std::vector<std::string>();
            graph[tempToElmUID].push_back(tempFromElmUID);
        }
        else {
            if (std::find(graph[tempToElmUID].begin(), graph[tempToElmUID].end(), tempFromElmUID) == graph[tempToElmUID].end()) {
                graph[tempToElmUID].push_back(tempFromElmUID);
            }
        }
    }
}

void tigl::FuselageGraph::createSimpleGraph() {

    simpleGraph.clear();

    root = "";
    leaf = "";

    bool stdIsPossible = true;

    std::map<std::string, std::vector<std::string>>::iterator it;
    for (it = graph.begin(); it != graph.end(); it++) {
        if (it->second.size() == 1) {
            if (root == "") {
                root = it->first;
            }
            else if (leaf == "") {
                leaf = it->first;
            }
            else {
                stdIsPossible = false;
            }
        }
        else if (it->second.size() < 1) {
            stdIsPossible = false;
        }
        else if (it->second.size() > 2) {
            stdIsPossible = false;
        }
    }

    if (leaf == "" || root == "") {
        stdIsPossible = false;
    }


    if( stdIsPossible){
        // the root is the extremity that is closest to the origin
        if (centers[leaf].norm2() < centers[root].norm2()) {
            std::string temp = root;
            root            = leaf;
            leaf            = temp;
        }

        // Now we have determine the root and we go though the graph starting from the root.
        std::string current = "";
        std::string next, nextA, nextB;
        std::set<std::string> usedElements; // to check previous and cycle

        next = root;
        while (next != "") {
            current = next;
            simpleGraph.push_back(current);
            usedElements.insert(current);

            if (graph[current].size() == 1) {
                next = graph[current][0];
                if (usedElements.count(next) != 0) {
                    next = ""; // end case
                }
            }
            else if (graph[current].size() == 2) {

                nextA = graph[current][0];
                nextB = graph[current][1];

                if (usedElements.count(nextA) == 0) {
                    next = nextA;
                }
                else if (usedElements.count(nextB) == 0) {
                    next = nextB;
                }
                else {
                    throw CTiglError("FuselageGraph::createSimpleGraph: unexpected graph");
                }
            }
            else {
                throw CTiglError("FuselageGraph::createSimpleGraph: unexpected graph");
            }
        }

    }else{
        //todo
        throw CTiglError("FuselageGraph::createSimpleGraph: Unstandard graph unsupported for the moment!");
    }



}

void tigl::FuselageGraph::clear() {

    fuselage = nullptr;
    graph.clear();
    simpleGraph.clear();
    circumferences.clear();
    centers.clear();
    root = "";
    leaf = "";


}

std::vector<std::string> tigl::FuselageGraph::getElementsInBetween(std::string startUID, std::string endUID) {

    std::vector<std::string> elementsBetween; // contain the start and the end

    bool afterStart = false;
    bool afterEnd   = false;

    std::vector<std::string>::iterator it;
    for (int i = 0; i < simpleGraph.size(); i++) {
        if(  simpleGraph[i] == startUID  && simpleGraph[i] == endUID ){ // special case where the end and the start are the same
            afterStart = true;
            afterEnd = true;
            elementsBetween.push_back(simpleGraph[i]);
        }
        else if (afterStart == false && ( simpleGraph[i] == startUID  || simpleGraph[i] == endUID ) ) { // we alow the start and end to be inverted
            afterStart = true;
            elementsBetween.push_back(simpleGraph[i]);
        }
        else if (afterStart == true && ( simpleGraph[i] == startUID  || simpleGraph[i] == endUID ) ) {
            afterEnd = true;
            elementsBetween.push_back(simpleGraph[i]);
        }
        else if (afterStart == true && afterEnd == false) {
            elementsBetween.push_back(simpleGraph[i]);
        }
    }

    return elementsBetween;
}

std::vector<std::string> tigl::FuselageGraph::getElementsBefore(std::string uid) {

    std::vector<std::string> elementsBefore;

    bool after = false;

    std::vector<std::string>::iterator it;
    for (int i = 0; i < simpleGraph.size(); i++) {
        if(simpleGraph[i] == uid){
            after = true;
        }
        if( after == false){
            elementsBefore.push_back(simpleGraph[i]);
        }
    }

    return elementsBefore;
}

std::vector<std::string> tigl::FuselageGraph::getElementsAfter(std::string uid) {

    std::vector<std::string> elementsAfter;

    bool after = false;

    std::vector<std::string>::iterator it;
    for (int i = 0; i < simpleGraph.size(); i++) {
        if( after ){
            elementsAfter.push_back(simpleGraph[i]);
        }
        if(simpleGraph[i] == uid){
            after = true;
        }

    }

    return elementsAfter;
}
