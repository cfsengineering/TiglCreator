/*
* Copyright (C) 2018 CFS Engineering
*
* Created: 2018 Malo Drougard <malo.drougard@protonmail.com>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
//
// Created by makem on 18/04/18.
//

#include "CPACSTree.h"


#include <list>
#include <map>


#include "UniqueXPath.h"
#include "CTiglLogging.h"
#include "TIGLViewerException.h"

#include <tixicpp.h>



// create an over tree from tixi




namespace cpcr {

    CPACSTree::CPACSTree() {

        rootXPath = UniqueXPath();

        m_isBuild = false;

        m_root = nullptr;
    }




    void CPACSTree::build(TixiDocumentHandle handle, UniqueXPath rootXPath) {

        // TODO check if the given handle is valid

        tixiHandle = handle;
        // clean the old one if present
        if(isBuild()){
          clean();
        }

        this->rootXPath = rootXPath;

        if (!tixi::TixiCheckElement(tixiHandle, rootXPath.toString())) {
            throw TIGLViewerException("CPACSTree: The given unique xpath \"" + rootXPath.toString() + "\" seems not to be valid");
        }

        std::string rootType = rootXPath.getLastElementType();
        int rootIndex = rootXPath.getLastElementIndex();
        std::string rootUid = getUid(rootXPath, "");

        m_root = new CPACSTreeItem(this, rootXPath, rootType, rootIndex, rootUid) ;

        // recursive call to create element
        createChildrenRecursively(*m_root);
        m_isBuild = true;
    }


    void CPACSTree::createChildrenRecursively(CPACSTreeItem& parent) {

        UniqueXPath xpath = parent.getXPath();

        int childrenCount = getNumberOfChildren(xpath);

        int tixiRet = -1;

        // find out with type of children this elment have
        char *childType = NULL;
        std::list<std::string> types;
        for (int i = 1; i <= childrenCount; i++) {
            tixiRet = tixiGetChildNodeName(tixiHandle, xpath.toString().c_str(), i, &childType);
            // TODO: FIND HOW TIXI SHOULD MANAGE THE NONE NODE ELEMENT LIKE TEXT
            if (childType[0] == '#') {
                //std::cout << childType << std::endl;
            } else {
                types.push_back(childType);
            }

        }
        types.unique();

        // for each children create the associate CPACSItem and call this function
        UniqueXPath newXPath;
        int counterChildrenOfSameType = 0;
        for (std::string currentType : types) {

            tixiRet = tixiGetNamedChildrenCount(tixiHandle, xpath.toString().c_str(), currentType.c_str(),
                                                &counterChildrenOfSameType);
            if (counterChildrenOfSameType == 1) {
                newXPath = UniqueXPath(xpath.toString() ) ;
                newXPath.addParticleAtEnd( currentType ) ;
                CPACSTreeItem* newChildren  = parent.addChild(newXPath, currentType, 1, getUid(newXPath));
                createChildrenRecursively(*newChildren);    // Recursive call

            } else {
                for (int idx = 1; idx <= counterChildrenOfSameType; idx++) {
                    newXPath = UniqueXPath(xpath.toString() ) ;
                    newXPath.addParticleAtEnd( currentType + "[" + std::to_string(idx) + "]"  );
                    CPACSTreeItem* newChildren  = parent.addChild(newXPath, currentType, idx, getUid(newXPath) );
                    createChildrenRecursively(*newChildren);      // Recursive call
                }
            }
        }

        return;

    }


    CPACSTree::~CPACSTree() {
        clean();
    }

    void CPACSTree::clean() {
        m_isBuild = false;
        delete m_root;
        m_root = nullptr;
    }



    std::string CPACSTree::getUid(UniqueXPath target, std::string defaultRetrunedValue) {

            std::string r = defaultRetrunedValue;

            try {
                r = tixi::TixiGetTextAttribute(tixiHandle, target.toString().c_str(), "uID");
            } catch (tixi::TixiError error) {
                // we expect that some times this function is call on element that have no UID
                if( error.returnCode() != ATTRIBUTE_NOT_FOUND ){
                    throw error;
                }
            }

            return r;

    }

    int CPACSTree::getNumberOfChildren(UniqueXPath xpathObj) {
        std::string xpath  = xpathObj.toString();
        int childrenCount = -1;
        int tixiRet = -1;
        tixiRet = tixiGetNumberOfChilds(tixiHandle,xpath.c_str(),&childrenCount);
        return childrenCount;
    }

    void CPACSTree::reload() {

        UniqueXPath rootXPath = getRoot()->getXPath();
        this->build(tixiHandle, rootXPath);


    }


}   // end namespace cpcr



