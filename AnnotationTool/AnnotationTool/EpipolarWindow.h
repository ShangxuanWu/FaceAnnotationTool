#ifndef EPIPOLARWINDOW_H
#define EPIPOLARWINDOW_H

#include <string>

#include <QDialog>
#include <QWidget>
#include <QListWidget>
#include <QListWidgetItem>
#include <QHBoxLayout>
#include <QLabel>
#include <QDebug>
#include <QPixmap>
#include <QString>
#include <QPushButton>
#include <QGraphicsView>
#include <QDialogButtonBox>
//#include <QMainWindow>

#include <opencv2/core/core.hpp>

#include <vector>

#include "EpipolarImageView.h"
#include "utils.h"
#include "../GeometryConstraint/EpipolarConstraint.h"
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/pts_3d_conf.h"
#include "../camera_propagation/camera_geometry.h"

//class EpipolarWindow : public QWidget
class EpipolarWindow : public QDialog
{
	Q_OBJECT
public:
	//EpipolarWindow();
	EpipolarWindow(std::vector<QString> image_paths, int point_id_, std::vector<mycamera> camera_cluster, QWidget *parent=0);
	void drawLine();
	QDialogButtonBox* button_box;

	public slots:
	void getLeftClickedPoints(double, double);
	void getRightClickedPoints(double, double);
	void accept();
signals:
	void epipolarLineOkSignal(int, cv::Point3d);
protected:
	//void ok();
	void getLeftAndRightImagePaths(std::vector<QString> image_paths);
	void getLeftAndRightCamera(std::vector<mycamera> camera_cluster);

	// 2 QLabels, no zoom-in / zoom-out function
	QString left_image_path;
	EpipolarImageView* left_panel;
	QString right_image_path;
	EpipolarImageView* right_panel;

	QPushButton* get_line_button;
	//QPushButton* ok_button;
	QPushButton* cancel_button;


	QHBoxLayout* horizontal_layout;

	mycamera camera_1;
	mycamera camera_2;

	cv::Point2d left_clicked_point;
	cv::Point2d right_clicked_point;
	cv::Point3d result_pt;

	// point_id is 0-based
	int point_id;

	// need mouse click function
	// EpipolarConstraint* epipolar_constraint;
};
#endif
