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
// Created by makem on 29/05/18.
//

#include "ModificatorManager.h"
#include "CTiglUIDManager.h"


ModificatorManager::ModificatorManager( QTreeView* qView,
                                        //QWidget *applyInterface,
                                        ModificatorTransformationWidget* transformationModificator,
                                        ModificatorWingWidget* wingModificator,
                                        ModificatorPositioningsWidget* positioningsModificator,
                                        ModificatorFuselageWidget *fuselageModificator) {


    treeViewManager = new CPACSTreeView(qView);

    //this->widgetApply = widgetApply;
    //this->commitButton = widgetApply->findChild<QPushButton*>("commitButton");
    //this->cancelButton = widgetApply->findChild<QPushButton*>("cancelButton");

    this->transformationModificator = transformationModificator;
    this->wingModificator = wingModificator;
    this->positioningsModificator = positioningsModificator;
    this->fuselageModificator = fuselageModificator; // we use pointer for modificator widget because they are define in the GUI
    //this->fuselageModificator->init(this); needed?
    currentModificator = nullptr;
    this->config = nullptr;

    this->hideAll();

    connect(treeViewManager, SIGNAL( newSelectedTreeItem(cpcr::CPACSTreeItem *)), this, SLOT(dispatch(cpcr::CPACSTreeItem*)));

}


void ModificatorManager::setCPACSConfiguration(tigl::CCPACSConfiguration *newConfig) {
    this->config = newConfig;
    cpcr::UniqueXPath rootXPath("/cpacs/vehicles/aircraft/model[1]"); //TODO
    treeViewManager->displayNewTree(newConfig->GetTixiDocumentHandle(), rootXPath );

}



void ModificatorManager::applyCurrentModifications(){

    if(currentModificator != nullptr) {
        currentModificator->apply();
//        if (useCpacsStandard){
//            adapter->standardize();
//        }
    }else{
        LOG(WARNING) << "ModificatorManager::applyCurrentModifications: current modificator is null" <<std::endl;
    }
}



void ModificatorManager::applyCurrentCancellation() {
    if(currentModificator != nullptr){
        currentModificator->reset();
    }else{
        LOG(WARNING) << "ModificatorManager::applyCurrentCancellation: current modificator is null" <<std::endl;
    }
}


void ModificatorManager::dispatch(cpcr::CPACSTreeItem* item ) {


    if( ! item->isInitialized()) {
        LOG(ERROR) << "MODIFICATOR MANAGER GET A NULL ITEM";
    }
    else if(item->getType() == "transformation"){
        //currentModificator = transformationModificator;
        //this->setTransformationModificator(item);
    }
    else if(item->getType() == "wing"){
        //currentModificator = wingModificator;
        //this->setWingModificator(item);
    }
    else if (item->getType() == "positionings"){
        //currentModificator = positioningsModificator;
        //this->setPositioningsModificator(item);
    }
    else if (item->getType() == "fuselage"){
        currentModificator = fuselageModificator;
        this->setFuselageModificator(item);
    }
    else {
        currentModificator = nullptr;
        hideAll();
        LOG(INFO) <<  "MODIFICATOR MANAGER: item not suported";
    }
}


void ModificatorManager::hideAll() {
    bool visible = false;
    transformationModificator->setVisible(visible);
    wingModificator->setVisible(visible);
    positioningsModificator->setVisible(visible);
    fuselageModificator->setVisible(visible);
    //widgetApply->setVisible(visible);

}

void ModificatorManager::setFuselageModificator(cpcr::CPACSTreeItem* item) {
    hideAll();
    tigl::CTiglUIDManager& uidManager = config->GetUIDManager();
    tigl::CCPACSFuselage& fuselage = uidManager.ResolveObject<tigl::CCPACSFuselage>(item->getUid());
    fuselageModificator->setFuselage(fuselage);
    fuselageModificator->setVisible(true);
    //widgetApply->setVisible(true);
}

/*
void ModificatorManager::setTransformationModificator(cpcr::CPACSTreeItem * item ) {

    hideAll();
    transformationModificator->setTransformation(item);
    transformationModificator->setVisible(true);
    widgetApply->setVisible(true);

}


void ModificatorManager::setWingModificator(cpcr::CPACSTreeItem *item) {
    hideAll();
    wingModificator->setWing(item);
    wingModificator->setVisible(true);
    widgetApply->setVisible(true);
}


void ModificatorManager::setPositioningsModificator(cpcr::CPACSTreeItem *item) {
    hideAll();
    positioningsModificator->setPositionings(item);
    positioningsModificator->setVisible(true);
    widgetApply->setVisible(true);
}


void ModificatorManager::standardizeCurrentFile() {
    adapter->standardize();
    setUseCPACSStandard(true);
}

bool ModificatorManager::isStandardized() {
    adapter->isStandardized();
}

void ModificatorManager::setUseCPACSStandard(bool value) {
    useCpacsStandard = value;

}

void ModificatorManager::resetAdapter(CPACSCreatorAdapter* newAdapter) {
    adapter = newAdapter;
    currentModificator = nullptr;
    this->hideAll();
}

*/