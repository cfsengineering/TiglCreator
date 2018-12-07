//
// Created by cfse on 9/26/18.
//

#ifndef TIGL_TIGLVIEWERFUSELAGEWIDGET_H
#define TIGL_TIGLVIEWERFUSELAGEWIDGET_H

#include "ModificatorWidget.h"
#include "CCPACSFuselage.h"
//#include "CPACSCreatorLib/CPACSTreeItem.h"
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QComboBox>

class ModificatorFuselageWidget : public ModificatorWidget
{

    Q_OBJECT
    /*
public slots:

    void expendLengthDetails(bool checked);

    void setPartialLengthFromComboBoxes();
    // If a new partial length is set it will recompute the the global length
    void recomputeTotalLength(double newPartialLength);

    void setCircumferenceFromRadius(double newRadius);
    void setRadiusFromCircumference(double newCircumference);
*/

public:
    ModificatorFuselageWidget(QWidget* parent = 0);


    void setFuselage(tigl::CCPACSFuselage& fuselage);

    /*
       // Initialize the linked objects
       void init(ModificatorManager* associate) override ;

       void apply() override ;
       void reset() override ;



   private:

       cpcr::CPACSTreeItem * fuselageItem;

       // Length parameters
       QDoubleSpinBox* spinBoxLength;
       QPushButton*  btnExpendLengthDetails;
       QWidget* widgetLengthDetails;
       QComboBox* comboBoxLengthE1;
       QComboBox* comboBoxLengthE2;
       QDoubleSpinBox* spinBoxPartialLength;

       // Internal length parameters
       double internalLength;
       double internalPartialLength;


       // Circumference parameters:
       QDoubleSpinBox* spinBoxCircumference;
       QDoubleSpinBox* spinBoxRadius;

       // Internal circumference parameters:
       double internalCircumference;


   */

private:
    tigl::CCPACSFuselage* fuselage;
};

#endif //TIGL_TIGLVIEWERFUSELAGEWIDGET_H
