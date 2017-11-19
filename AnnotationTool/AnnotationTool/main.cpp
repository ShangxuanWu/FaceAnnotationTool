#define ENABLE_ASSERTS

#include "MainWindow.h"
#include <QStyleFactory>
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
	QStringList styleList = QStyleFactory::keys();
	qDebug() << styleList;
	//a.setStyle(QStyleFactory::create("Fusion"));
    MainWindow w;
    w.showMaximized();
    return a.exec();
}
