//
// Created by cfse on 12/13/18.
//

#ifndef TIGL_FUSELAGEGRAPH_H
#define TIGL_FUSELAGEGRAPH_H

#include "CTiglPoint.h"


#include <map>
#include <vector>

namespace tigl
{

class CCPACSFuselage;
class FuselageGraph
{

public:

    FuselageGraph();

    void build( CCPACSFuselage* inFuselage);

    inline bool isBuild() {return (fuselage != nullptr); };

    void clear();

    /**
     * Returns the element UIDs  between startUID and endUID.
     * The start and the end UID are also contains in the list if the they are present in the graph.
     *
     * @remark The graph should be build otherwise, it return a empty list
     * @param startUID
     * @param endUID
     * @return
     */
    std::vector<std::string> getElementsInBetween(std::string startUID, std::string endUID) ;

    /**
     * Returns the element UIDs before uid.
     * The uid is not contains in the list.
     *
     * @remark The graph should be build otherwise, it return a empty list
     * @param uid
     * @return
     */
    std::vector<std::string> getElementsBefore(std::string uid) ;



    /**
     * Returns the element UIDs after uid.
     * The uid is not contains in the list.
     *
     * @remark The graph should be build otherwise, it return a empty list
     * @param uid
     * @return
     */
    std::vector<std::string> getElementsAfter(std::string uid) ;


    inline std::string getNoseUID() {return root; };

    inline std::string getTailUID() {return leaf; };

    inline std::vector<std::string>& getSimpleGraph() { return simpleGraph; };

    inline std::map<std::string, std::vector<std::string> >& getGraph() { return graph; };


protected:

    void createGraph();
    void createSimpleGraph();


    CCPACSFuselage* fuselage;

    std::map<std::string, std::vector<std::string> > graph;
    std::vector<std::string> simpleGraph;
    std::string root;
    std::string leaf;



};
}
#endif //TIGL_FUSELAGEGRAPH_H
