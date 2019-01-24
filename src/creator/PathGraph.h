//
// Created by makem on 17/01/19.
//

#ifndef TIGL_PATHGRAPH_H
#define TIGL_PATHGRAPH_H

#include <vector>
#include <string>

namespace tigl {

    class PathGraph {

    public:
        PathGraph();

        void set( std::vector<std::string> newGraph );
        virtual void clear();

        inline bool isEmpty() { return graph.size() == 0; };

        /**
         * Return the graph as a vector of string
         *
         * @return
         */
        inline std::vector<std::string>& getGraph() { return graph; };


        /**
         * Returns the element UIDs  between uid1 and uid2.
         * The order of the given elements are not important, but the returned list will be ordered
         * The two given elements are also contained in the list.
         *
         * @remark  If the graph is empty or it don't contain one of the UID, it will return a empty list
         * @param startUID
         * @param endUID
         * @return
         */
        std::vector<std::string> getElementsInBetween(std::string uid1, std::string uid2) ;

        /**
         * Returns the element UIDs before uid.
         * The uid is not contains in the returned list.
         *
         * @remark If the graph is empty or it don't contain the given UID, it will return a empty list
         * @param uid
         * @return
         */
        std::vector<std::string> getElementsBefore(std::string uid) ;



        /**
         * Returns the element UIDs after uid.
         * The uid is not contains in the list.
         *
         * @remark  If the graph is empty or it don't contain the given UID, it will return a empty list
         * @param uid
         * @return
         */
        std::vector<std::string> getElementsAfter(std::string uid) ;


        /**
         * Returns the first element of the graph path
         *
         * @remark If the graph is empty, it will throw a exception
         * @return
         */
        std::string getRoot();

        /**
         * Return the last element of the graph path
         *
         *
         * @remark If the graph is empty, it will throw a exception
         * @return
         */
        std::string getLeaf();




    private:

        std::vector<std::string> graph;

    };

}
#endif //TIGL_PATHGRAPH_H
