/*
* Copyright (C) 2007-2014 German Aerospace Center (DLR/SC)
*
* Created: 2014-05-08 Martin Siggel <martin.siggel@dlr.de>
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

#ifndef TIGLVIEWERSCREENSHOTDIALOG_H
#define TIGLVIEWERSCREENSHOTDIALOG_H

#include "tigl_internal.h"
#include <QDialog>

namespace Ui
{
    class TiglViewerScreenshotDialog;
}

class TIGLViewerScreenshotDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit TIGLViewerScreenshotDialog(QString filename, QWidget *parent = 0);

    void setImageSize(int width, int height);
    void getImageSize(int&width, int& height) const;
    
    void setQualityValue(int quality);
    int  getQualityValue() const;
    
    bool getWhiteBGEnabled() const;

    ~TIGLViewerScreenshotDialog() OVERRIDE;
    
private:
    Ui::TiglViewerScreenshotDialog *ui;
};

#endif // TIGLVIEWERSCREENSHOTDIALOG_H
