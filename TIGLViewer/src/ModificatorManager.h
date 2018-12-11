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

#ifndef TIGL_MODIFICATORMANAGER_H
#define TIGL_MODIFICATORMANAGER_H

#include "ModificatorWidget.h"
#include "ModificatorTransformationWidget.h"
#include "ModificatorWingWidget.h"
#include "ModificatorPositioningsWidget.h"
#include "ModificatorFuselageWidget.h"
#include "CCPACSConfiguration.h"
#include "CPACSTreeItem.h"
#include "CPACSTreeView.h"

#include <QPushButton>

/**
 * @brief Main class to articulate the modification interface.
 *
 * This class holds four main components:
 * CCPACSConfiguration: that is the main tigl object to work with a cpacs configuration
 * CPACSTreeView: that is the main tree interface.
 * Different ModificatorWidgets: that are the interfaces for particular tigl object.
 * ApllyWidget: that hold the commit and cancel buttons
 *
 * To work with this class, first we need to set the CCPACSConfiguration. The CPACSTreeView will be updated and the
 * correct tree build and displayed.
 * Once the CPACSConfiguration is set, it will wait for a "newSelectedTreeItem" signal. If a "newSelectedTreeItem" signal
 * is emit, the "dispatch" function will be called. The dispatch function will look at the new selected element and if the
 * element has a associate ModificatorWidget, it will set this particular ModificatorWidget and display it.
 * Then if the user, click on the "apply" function, the modificator manager will apply the change onto the
 * CPACSConfiguration though the ModificatorWidget.
 *
 * @author Malo Drougard
 *
 */
class ModificatorManager: public QObject {
    Q_OBJECT



signals:
    void configurationEdited();

public slots:
    void dispatch(cpcr::CPACSTreeItem* item ) ;
    void applyCurrentModifications();
    // simply reset the values displayed by the current modificator from the file
    void applyCurrentCancellation();

public:
    ModificatorManager( QTreeView* qView,
                        QWidget* applyInterface,
                        ModificatorTransformationWidget* transformationModificator,
                        ModificatorWingWidget* wingModificator,
                        ModificatorPositioningsWidget* positioningsModificator,
                        ModificatorFuselageWidget * fuselageModificator
                       );

    void setCPACSConfiguration( tigl::CCPACSConfiguration* newConfig);

    //void setTransformationModificator(cpcr::CPACSTreeItem * item );
    //void setWingModificator(cpcr::CPACSTreeItem * item);
    //void setPositioningsModificator(cpcr::CPACSTreeItem * item);
    void setFuselageModificator(cpcr::CPACSTreeItem* item);
    void hideAll();

    //void standardizeCurrentFile();
    //bool isStandardized();
    //void setUseCPACSStandard(bool value);
    //inline bool useCPACSStandard() {return useCpacsStandard; };

    //inline ProfilesDBManager* getProfilesDB() {return profilesDB; }


protected:

    inline bool configurationIsSet() {return (config != nullptr); }


private:

    tigl::CCPACSConfiguration* config ;

    CPACSTreeView*  treeViewManager;


    ModificatorTransformationWidget* transformationModificator;
    ModificatorWingWidget* wingModificator;
    ModificatorPositioningsWidget* positioningsModificator;
    ModificatorFuselageWidget * fuselageModificator;

    ModificatorWidget* currentModificator;

    // cancel/apply interface
    QWidget* widgetApply;
    QPushButton* commitButton;
    QPushButton* cancelButton;

    //bool useCpacsStandard;


};


#endif //TIGL_MODIFICATORMANAGER_H
