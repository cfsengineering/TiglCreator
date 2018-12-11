//
// Created by cfse on 11/5/18.
//

#ifndef TIGL_CPACSTREEVIEW_H
#define TIGL_CPACSTREEVIEW_H


#include "CPACSAbstractModel.h"
#include "CPACSTree.h"
#include <QTreeView>
#include <QObject>

/**
 * @brief Main class to manage the tree view.
 * It holds the cpcr::CPACSTree that contain the cpacs structure of the open file.
 * It holds the QTreeView that display the tree in the GUI.
 * It holds the CPACSAbstractModel that prepare the CPACSTree for the QTreeView.
 * It holds the QItemSelectionModel that emit a signal when a new element is selected in the QTreeView.
 *
 * The goal of this class is to articulate all the previous components together.
 * It has one function to display a new tree based on a TixiHandle.
 * It has one function to clear the tree.
 * It has one function to update the tree based on current TixiHandle.
 * It emit a signal when a new element in the tree is selected.
 *
 * @author Malo Drougard
 */
class CPACSTreeView: public QObject {

Q_OBJECT


signals:

    void newSelectedTreeItem(cpcr::CPACSTreeItem* );


private slots:

    void onSelectionChanged(const QItemSelection& newSelection, const QItemSelection& oldSelection);



public:

    CPACSTreeView(QTreeView* view);

    /**
     * Clear the displayed tree and delete the CPACSTree dat
     */
    void clear();


    /**
     * Build the new tree based on the TixiHandle and update the model and the display
     * @param handle : the TixiHandle used to retrieve the cpacs data
     * @param root :where the tree need to start
     */
    void displayNewTree(TixiDocumentHandle handle, cpcr::UniqueXPath root);

    /**
     * Rebuild the current tree based on the tixi data
     * @remark: the internal tixi handle remain the same
     */
    void refresh();



private:

    cpcr::CPACSTree tree;
    CPACSAbstractModel* model;
    QTreeView * treeView;
    QItemSelectionModel* selectionModel;

    //std::string backupSelectedUID;



};


#endif //TIGL_CPACSTREEVIEW_H
