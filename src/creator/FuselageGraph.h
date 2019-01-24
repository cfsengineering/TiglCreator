//
// Created by cfse on 12/13/18.
//

#ifndef TIGL_FUSELAGEGRAPH_H
#define TIGL_FUSELAGEGRAPH_H

#include "CTiglPoint.h"

#include "PathGraph.h"
#include <map>
#include <vector>

namespace tigl {

class CCPACSFuselage;
class FuselageGraph : public PathGraph {

public:
  FuselageGraph();

  void set(CCPACSFuselage *fuselageInput);

protected:
};
} // namespace tigl
#endif // TIGL_FUSELAGEGRAPH_H
