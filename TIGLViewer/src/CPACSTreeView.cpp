//
// Created by cfse on 11/5/18.
//

#include "CPACSTreeView.h"
#include "CPACSTreeItem.h"
#include "CTiglLogging.h"

CPACSTreeView::CPACSTreeView(QTreeView* view) {

    treeView = view;
    model =  new CPACSAbstractModel(nullptr);
    treeView->setModel(model);
    selectionModel = treeView->selectionModel();
    //backupSelectedUID = "";

    connect(selectionModel, SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this,
                   SLOT(onSelectionChanged( const QItemSelection&, const QItemSelection&)));


}

void CPACSTreeView::onSelectionChanged(const QItemSelection & newSelection, const QItemSelection & oldSelection) {

    if (model->isValid() ){
        cpcr::CPACSTreeItem* newSelectedItem = model->getItemFromSelection(newSelection);
        emit newSelectedTreeItem(newSelectedItem);
    }else{
        LOG(WARNING) <<  "CPACSTreeView: onSelectionChanged called but no valid model is set" << std::endl;
    }
}

void CPACSTreeView::clear() {
    model->disconnectInternalTree();
    tree.clean();
}

void CPACSTreeView::displayNewTree(TixiDocumentHandle handle,cpcr::UniqueXPath root) {
    tree.build(handle, root );
    model->resetInternalTree(&tree);
}

void CPACSTreeView::refresh() {
    model->disconnectInternalTree();
    tree.reload();
    model->resetInternalTree(&tree);
}

