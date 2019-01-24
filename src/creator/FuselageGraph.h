//
// Created by cfse on 12/13/18.
//

#ifndef TIGL_FUSELAGEGRAPH_H
#define TIGL_FUSELAGEGRAPH_H

#include "CTiglPoint.h"


#include <map>
#include <vector>
#include "PathGraph.h"

namespace tigl
{

class CCPACSFuselage;
class FuselageGraph : public PathGraph
{

public:

    FuselageGraph();

    void set(CCPACSFuselage* fuselageInput);

protected:

};
}
#endif //TIGL_FUSELAGEGRAPH_H
