//
// Created by cfse on 11/5/18.
//

#ifndef TIGL_CPACSTREEVIEW_H
#define TIGL_CPACSTREEVIEW_H


#include "CPACSAbstractModel.h"
#include "CPACSTree.h"
#include <QTreeView>
#include <QObject>


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


    void displayNewTree(TixiDocumentHandle handle, cpcr::UniqueXPath root);

    /**
     * Rebuild the actual tree build on the tixi data
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
