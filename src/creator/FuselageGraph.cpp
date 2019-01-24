//
// Created by cfse on 12/13/18.
//

#include "FuselageGraph.h"
#include "CCPACSFuselage.h"
#include "CCPACSFuselageSegment.h"

tigl::FuselageGraph::FuselageGraph() : PathGraph() {}

void tigl::FuselageGraph::set(tigl::CCPACSFuselage *fuselage) {

  clear();

  if (fuselage == nullptr) {
    throw CTiglError(
        "FuselageGraph::set(): Invalid fuselage pointer given as input.");
  }

  std::vector<std::string> newGraph;

  CCPACSFuselageSegments &segments = fuselage->GetSegments();
  CCPACSFuselageSections &sections = fuselage->GetSections();

  std::string tempFromElmUID, tempToElmUID;

  // We assume that the segments are order
  for (int i = 1; i <= segments.GetSegmentCount(); i++) {

    CCPACSFuselageSegment &seg = segments.GetSegment(i);
    tempFromElmUID = seg.GetStartSectionElementUID();
    newGraph.push_back(tempFromElmUID);

    // for the last one we need also to put the "to" element
    if (i == segments.GetSegmentCount()) {
      tempToElmUID = seg.GetEndSectionElementUID();
      newGraph.push_back(tempToElmUID);
    }
  }

  PathGraph::set(newGraph);
}
