#ifndef EPIPOLARIMAGEVIEW_H
#define EPIPOLARIMAGEVIEW_H

#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>

#include <QGraphicsView>
#include <QTransform>
#include <QImage>
#include <QtGui>

#include "utils.h"

class EpipolarImageView : public QGraphicsView
{
	Q_OBJECT
public:
	EpipolarImageView(QWidget *parent, int type_);
	void setImage(const std::string& image_fpath);
	void enableClicking();
	void clearClickedPoint();
	void setEpipolarLine(cv::Mat line);

signals:
	void sendLeftPanelClickedPoint(double, double);
	void sendRightPanelClickedPoint(double, double);

protected:
	void paintEvent(QPaintEvent *event);

	void mouseMoveEvent(QMouseEvent* event);
	void mousePressEvent(QMouseEvent* event);
	void wheelEvent(QWheelEvent* event);

	cv::Mat m_img;

	cv::Mat m_img_plus_line;

	cv::Point2d clicked_point;

	double currScale_;
	double currRot_deg_;
	cv::Vec2d transVec_;

	QTransform m_matrix;
	QTransform m_matrix_inv;

	QPoint m_pntDownPos;

	int type;

	// used for the right panel, only be enabled after getting epipolar line
	bool enable;

	double line[3];
};
#endif
