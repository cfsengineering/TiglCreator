//
// Created by cfse on 9/26/18.
//

#ifndef TIGL_MODIFICATORFUSELAGEWIDGET_H
#define TIGL_MODIFICATORFUSELAGEWIDGET_H

#include "ModificatorWidget.h"
#include "CCPACSFuselage.h"
#include "CPACSTreeItem.h"
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QComboBox>

class ModificatorFuselageWidget : public ModificatorWidget
{

    Q_OBJECT

public slots:

    void expendLengthDetails(bool checked);

    void setPartialLengthFromComboBoxes();
    // If a new partial length is set it will recompute the the global length
    void recomputeTotalLength(double newPartialLength);


    void expendRadiusDetails(bool checked);

    void setRadiusBetweenFromComboBoxes();
    //void setCircumferenceFromRadius(double newRadius);
    //void setRadiusFromCircumference(double newCircumference);



public:
    ModificatorFuselageWidget(QWidget* parent = 0);

    void setFuselage(tigl::CCPACSFuselage& fuselage);

    // Initialize the linked objects
    void init();

    void apply() override;

    void reset() override;

private:

    // Length parameters
    QDoubleSpinBox* spinBoxLength;
    QPushButton* btnExpendLengthDetails;
    QWidget* widgetLengthDetails;
    QComboBox* comboBoxLengthE1;
    QComboBox* comboBoxLengthE2;
    QDoubleSpinBox* spinBoxPartialLength;

    // Internal length parameters
    double internalLength;
    double internalPartialLength;

    // Circumference parameters:
    QDoubleSpinBox* spinBoxRadius;
    QPushButton* btnExpendRadiusDetails;
    QWidget* widgetRadiusDetails;
    QComboBox* comboBoxRadiusBE1;
    QComboBox* comboBoxRadiusBE2;
    QDoubleSpinBox* spinBoxRadiusBetween;


    // Internal circumference parameters:
    double internalRadius;
    double internalRadiusBetween;

private:
    tigl::CCPACSFuselage* fuselage;
};

#endif //TIGL_TIGLVIEWERFUSELAGEWIDGET_H
