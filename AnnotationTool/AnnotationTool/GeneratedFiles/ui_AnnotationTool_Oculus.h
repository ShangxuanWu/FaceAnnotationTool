/********************************************************************************
** Form generated from reading UI file 'AnnotationTool_Oculus.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ANNOTATIONTOOL_OCULUS_H
#define UI_ANNOTATIONTOOL_OCULUS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_AnnotationTool_OculusClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *AnnotationTool_OculusClass)
    {
        if (AnnotationTool_OculusClass->objectName().isEmpty())
            AnnotationTool_OculusClass->setObjectName(QStringLiteral("AnnotationTool_OculusClass"));
        AnnotationTool_OculusClass->resize(600, 400);
        menuBar = new QMenuBar(AnnotationTool_OculusClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        AnnotationTool_OculusClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(AnnotationTool_OculusClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        AnnotationTool_OculusClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(AnnotationTool_OculusClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        AnnotationTool_OculusClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(AnnotationTool_OculusClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        AnnotationTool_OculusClass->setStatusBar(statusBar);

        retranslateUi(AnnotationTool_OculusClass);

        QMetaObject::connectSlotsByName(AnnotationTool_OculusClass);
    } // setupUi

    void retranslateUi(QMainWindow *AnnotationTool_OculusClass)
    {
        AnnotationTool_OculusClass->setWindowTitle(QApplication::translate("AnnotationTool_OculusClass", "AnnotationTool_Oculus", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class AnnotationTool_OculusClass: public Ui_AnnotationTool_OculusClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ANNOTATIONTOOL_OCULUS_H
