#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_AnnotationTool_Oculus.h"

class AnnotationTool_Oculus : public QMainWindow
{
    Q_OBJECT

public:
    AnnotationTool_Oculus(QWidget *parent = Q_NULLPTR);

private:
    Ui::AnnotationTool_OculusClass ui;
};
