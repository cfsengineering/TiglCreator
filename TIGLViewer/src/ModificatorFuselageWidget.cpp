
//
// Created by cfse on 9/26/18.
//

#include "ModificatorFuselageWidget.h"
#include "ModificatorManager.h"

ModificatorFuselageWidget::ModificatorFuselageWidget(QWidget* parent)
    : ModificatorWidget(parent)
{
}



void ModificatorFuselageWidget::setFuselage(tigl::CCPACSFuselage& newFuselage) {
    this->fuselage = &newFuselage;

    internalLength = fuselage->GetLength();
    spinBoxLength->setValue(internalLength);
    comboBoxLengthE1->clear();

    std::vector<std::string> fusleageGraph = fuselage->GetCreatorGraph();
    QStringList elementsUids;
    for (int i = 0 ; i < fusleageGraph.size(); i++){
        elementsUids.push_back( QString(fusleageGraph[i].c_str()) );
    }
    comboBoxLengthE1->addItems(elementsUids);
    comboBoxLengthE2->clear();
    comboBoxLengthE2->addItems(elementsUids);
    comboBoxLengthE2->setCurrentIndex(elementsUids.size() -1 ); // set the last element of the list
    setPartialLengthFromComboBoxes();
    // do total length after partial length, because changing partial can change total
    internalLength = fuselage->GetLength();
    spinBoxLength->setValue(internalLength);
    widgetLengthDetails->setVisible(false);

    // circumference
    //internalCircumference = associateManager->adapter->getFuselageMaximalCircumference(this->fuselageItem);
    //spinBoxCircumference->setValue(internalCircumference);

}


void ModificatorFuselageWidget::init() {

    ModificatorWidget::init(nullptr);

    spinBoxLength = this->findChild<QDoubleSpinBox*>("spinBoxLength");
    widgetLengthDetails = this->findChild<QWidget*>("widgetLengthDetails");
    btnExpendLengthDetails =  this->findChild<QPushButton*>("btnExpendLengthDetails");
    comboBoxLengthE1 = this->findChild<QComboBox*>("comboBoxE1");
    comboBoxLengthE2 = this->findChild<QComboBox*>("comboBoxE2");
    spinBoxPartialLength = this->findChild<QDoubleSpinBox*>("spinBoxPartialLength");

    //spinBoxCircumference = this->findChild<QDoubleSpinBox*>("spinBoxCircumference");
    //spinBoxRadius = this->findChild<QDoubleSpinBox*>("spinBoxRadius");

    widgetLengthDetails->hide();

    // connect the extend buttons with their slot
    connect(btnExpendLengthDetails, SIGNAL(clicked(bool)), this, SLOT(expendLengthDetails(bool)) );
    connect(comboBoxLengthE1,SIGNAL(currentIndexChanged(int )), this, SLOT(setPartialLengthFromComboBoxes()));
    connect(comboBoxLengthE2,SIGNAL(currentIndexChanged(int )), this, SLOT(setPartialLengthFromComboBoxes()));
    connect(spinBoxPartialLength, SIGNAL(valueChanged(double)), this, SLOT(recomputeTotalLength(double)));

    // connect spinBoxRadius with spinBoxCircumference
    //connect(spinBoxRadius, SIGNAL(valueChanged(double)), this, SLOT(setCircumferenceFromRadius(double)));
    //connect(spinBoxCircumference, SIGNAL(valueChanged(double)), this, SLOT(setRadiusFromCircumference(double)));



}


// inverse the visibility
void ModificatorFuselageWidget::expendLengthDetails(bool checked) {
    widgetLengthDetails->setVisible(! (widgetLengthDetails->isVisible() ));
    if(widgetLengthDetails->isVisible()){
        // Reset the values to the file values, avoid modifying from details and main at the same time
        internalLength = fuselage->GetLength();
        spinBoxLength->setValue(internalLength);
        setPartialLengthFromComboBoxes();
        spinBoxLength->setReadOnly(true);

    }else{
        // Reset the values to the file values, avoid modifying from details and main at the same time
        internalLength = fuselage->GetLength();
        spinBoxLength->setValue(internalLength);
        setPartialLengthFromComboBoxes();
        spinBoxLength->setReadOnly(false);

    }

}


void ModificatorFuselageWidget::setPartialLengthFromComboBoxes(){

    QString uid1 = comboBoxLengthE1->currentText();
    QString uid2 = comboBoxLengthE2->currentText();
    internalPartialLength = fuselage->GetLengthBetween(uid1.toStdString(), uid2.toStdString());
    // we reset the display value of total length, because the old displayed value can be have modified by recomputeTotalLength
    spinBoxLength->setValue(internalLength);
    spinBoxPartialLength->setValue(internalPartialLength);


}

// call when a new partial length is set
void ModificatorFuselageWidget::recomputeTotalLength(double newPartialLength){

    if( !( isApprox(newPartialLength, internalPartialLength) )){ // avoid diff between spin box implementation and double
        double diff = newPartialLength - internalPartialLength;
        spinBoxLength->setValue(internalLength + diff);
    }
}


/*
void ModificatorFuselageWidget::setCircumferenceFromRadius(double newRadius) {
    bool block = spinBoxCircumference->blockSignals(true); // to avoid infinite loop with setRadiusFromCircumference
    spinBoxCircumference->setValue(2.0* M_PI * newRadius);
    spinBoxCircumference->blockSignals(block);

}


void ModificatorFuselageWidget::setRadiusFromCircumference(double newCircumference) {
    bool block = spinBoxRadius->blockSignals(true); // to avoid infinite loop with setCircumferenceFromRadius
    spinBoxRadius->setValue( newCircumference/ ( M_PI * 2.0));
    spinBoxRadius->blockSignals(block);

}
*/

void ModificatorFuselageWidget::apply() {

   bool lengthHasChanged = ( (! isApprox(internalLength, spinBoxLength->value()) ) );
   bool partialLengthHasChanged = ( ! isApprox(internalPartialLength, spinBoxPartialLength->value() ) );
   bool isPartialCase = widgetLengthDetails->isVisible(); // if expend length details is shown, the details modifications prime on the main mofif
   //bool circumferenceHasChanged = ( (! isApprox(internalCircumference, spinBoxCircumference->value())) );
    bool circumferenceHasChanged = false;

    if(lengthHasChanged && (!isPartialCase)){
       internalLength = spinBoxLength->value();
       fuselage->SetLength(internalLength);
   }
   if(partialLengthHasChanged && isPartialCase){
        internalPartialLength = spinBoxPartialLength->value();
        QString uid1 = comboBoxLengthE1->currentText();
        QString uid2 = comboBoxLengthE2->currentText();
        fuselage->SetLengthBetween(uid1.toStdString(), uid2.toStdString(), internalPartialLength);
   }

   /*
   if( circumferenceHasChanged){
       internalCircumference = spinBoxCircumference->value();
       associateManager->adapter->setFuselageMaximalCircumference(fuselageItem, internalCircumference);
   }
*/

   //Todo: what we need to do XD
  /*
   if(lengthHasChanged || partialLengthHasChanged || circumferenceHasChanged ){
       associateManager->adapter->writeToFile();
   }
   */

}

void ModificatorFuselageWidget::reset() {
    if(fuselage != nullptr){
        this->setFuselage(*fuselage);
    }else{
        LOG(WARNING) << "ModificatorWingWidget: reset call but wing is not set!";
    }
}


