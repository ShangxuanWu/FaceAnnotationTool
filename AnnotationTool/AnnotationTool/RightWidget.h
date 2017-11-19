#ifndef RIGHTWIDGET_H
#define RIGHTWIDGET_H

#include <string>

#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>
#include <QVBoxLayout>
#include <QStatusBar>
#include <QLabel>
#include <QDebug>
#include <QPixmap>
#include <QString>
#include <QDebug>
//#include <QMainWindow>

#include "utils.h"
//#include "MainWindow.h"
#include "MeshVisibilityHelper.h"

//class MainWindow;

class RightWidget : public QWidget
{
	Q_OBJECT
public:
	RightWidget();
	void clearThumbnails(); // clear the thumbnail container
	void addThumbnailItem(QListWidgetItem* item); // add item to thumbnail container
	bool isEmpty();
	void clear();

	QListWidget* thumbnail_container;
	MeshVisibilityHelper* mesh_visibility_helper;

	QListWidgetItem* getCurrentItem();

	public slots:
	//void onThumbnailContainerItemClicked(QListWidgetItem*);
	void setHint(int point_id);
	void setHint3D(int point_id);
	void cancelHint();
	void cancelHint3D();
protected:

	void resizeEvent(QResizeEvent *event);
	QPixmap processImage(std::string fn);

	QLabel* warning_text;
	QLabel* hint_container;
	QLabel* hint_text;
	QVBoxLayout* vertical_layout;
};
#endif
