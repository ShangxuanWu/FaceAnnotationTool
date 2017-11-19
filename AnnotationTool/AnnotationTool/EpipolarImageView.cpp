#include "EpipolarImageView.h"

EpipolarImageView::EpipolarImageView(QWidget *parent, int type_)
{
	type = type_;

	clicked_point = cv::Point2d(-1, -1);

	currScale_ = 0.15;
	currRot_deg_ = 0;
	transVec_[0] = 0;
	transVec_[1] = 0;

	enable = false;
}


void EpipolarImageView::setImage(const std::string & image_fpath)
{
	m_img = cv::imread(image_fpath);
	//imagePath_ = image_fpath;

	cv::cvtColor(m_img, m_img, CV_BGR2RGB);

	m_matrix.setMatrix(1, 0, 0, 0, 1, 0, 0, 0, 1);
	m_matrix_inv = m_matrix.inverted();
	viewport()->update();
}

void EpipolarImageView::paintEvent(QPaintEvent *event)
{
	// update m_matrix
	cv::Mat_<double> rot_mat(2, 2);
	//qDebug() << currRot_deg_;
	const double cos_angle = cos(PI * currRot_deg_ / 180.);
	const double sin_angle = sin(PI * currRot_deg_ / 180.);
	rot_mat(0, 0) = cos_angle;
	rot_mat(0, 1) = sin_angle;
	rot_mat(1, 0) = -sin_angle;
	rot_mat(1, 1) = cos_angle;

	m_matrix.setMatrix(
		currScale_ * rot_mat(0, 0), currScale_* rot_mat(0, 1), 0,
		currScale_ * rot_mat(1, 0), currScale_ * rot_mat(1, 1), 0,
		transVec_[0], transVec_[1], 1);

	m_matrix_inv = m_matrix.inverted();


	// draw image

	QPainter widget_painter(viewport());

	widget_painter.setWorldTransform(m_matrix);

	QImage qimg;

	if (type == EPIPOLAR_LEFT_WINDOW || (type == EPIPOLAR_RIGHT_WINDOW && enable == false))
	{
		qimg = QImage(m_img.ptr(), m_img.cols, m_img.rows, m_img.step, QImage::Format_RGB888);
	}
	else
	{
		// type == right window and enable == true
		qimg = QImage(m_img_plus_line.ptr(), m_img_plus_line.cols, m_img_plus_line.rows, m_img_plus_line.step, QImage::Format_RGB888);
	}

	widget_painter.drawImage(0, 0, qimg);

	// draw point

	const double curr_line_width = 2. / currScale_;

	if (clicked_point.x >= 0 && clicked_point.y >= 0 && clicked_point.x < m_img.cols && clicked_point.y < m_img.rows)
	{
		widget_painter.setPen(QPen(Qt::cyan, curr_line_width * 2));
		widget_painter.drawEllipse(clicked_point.x - 10, clicked_point.y - 10, 20, 20);
	}
}

void EpipolarImageView::mouseMoveEvent(QMouseEvent * event)
{
	QPoint pnt = event->pos();

	if (event->buttons() == Qt::LeftButton)
	{
		QPointF pntf = (pnt - m_pntDownPos) / currScale_;
		m_pntDownPos = event->pos();
		m_matrix.translate(pntf.x(), pntf.y());

		transVec_[0] += pntf.x() * currScale_;
		transVec_[1] += pntf.y() * currScale_;

	}
	viewport()->update();
}

void EpipolarImageView::mousePressEvent(QMouseEvent * event)
{
	m_pntDownPos = event->pos();
	if (event->buttons() == Qt::LeftButton && (type==EPIPOLAR_LEFT_WINDOW || (type == EPIPOLAR_RIGHT_WINDOW && enable == true)) && QGuiApplication::keyboardModifiers().testFlag(Qt::ShiftModifier) == true)
	{
		QPointF new_click_pos_in_image = m_matrix_inv.map(QPointF(event->pos().x(), event->pos().y()));
		clicked_point.x = new_click_pos_in_image.x();
		clicked_point.y = new_click_pos_in_image.y();
		if (type == EPIPOLAR_LEFT_WINDOW)
		{
			// send to parent window, calculate epipolar line, enable right panel clicking
			emit sendLeftPanelClickedPoint(clicked_point.x, clicked_point.y);
		}
		else
		{
			// right window && enable
			// calculate the y axis using x axis
			clicked_point.x = new_click_pos_in_image.x();
			clicked_point.y = -(line[2] + line[0] * new_click_pos_in_image.x()) / line[1];
			
			emit sendRightPanelClickedPoint(clicked_point.x, clicked_point.y);
		}
	}

	viewport()->update();
}

void EpipolarImageView::wheelEvent(QWheelEvent * event)
{
	QPoint angle_delta = event->angleDelta();

	QTransform curr_matrix = m_matrix;

	const double delta_scale = 0.05;
	if (angle_delta.y() > 0)
	{
		double curr_scale = m_matrix.m11();

		m_matrix.setMatrix(curr_matrix.m11() + delta_scale, 0, 0, 0, curr_matrix.m22() + delta_scale, 0, curr_matrix.m31(), curr_matrix.m32(), 1);


		currScale_ += delta_scale;
	}
	else if (angle_delta.y() < 0)
	{
		double curr_scale = m_matrix.m11();
		m_matrix.setMatrix(curr_matrix.m11() - delta_scale, 0, 0, 0, curr_matrix.m22() - delta_scale, 0, curr_matrix.m31(), curr_matrix.m32(), 1);

		if (currScale_ > 0.01)
			currScale_ -= delta_scale;
		else
			currScale_ = 0.01;
	}
	viewport()->update();
}

void EpipolarImageView::enableClicking()
{
	enable = true;
}

void EpipolarImageView::clearClickedPoint()
{
	clicked_point = cv::Point2d(-1, -1);
}

void EpipolarImageView::setEpipolarLine(cv::Mat line_)
{
	if (type == EPIPOLAR_RIGHT_WINDOW)
	{
		// copy an image
		m_img_plus_line = m_img.clone();

		// draw the line in the new image
		//double line[3];
		line[0] = line_.at<double>(0, 0);
		line[1] = line_.at<double>(1, 0);
		line[2] = line_.at<double>(2, 0);

		for (int x = 0; x < m_img_plus_line.cols; x++)
		{
			int y = -(line[2] + line[0] * x) / line[1];
			if (y >= 0 && y < m_img_plus_line.rows)
			{
				//std::cout << x << ' ' << y << std::endl;
				m_img_plus_line.at<cv::Vec3b>(y, x)[0] = 255;
				m_img_plus_line.at<cv::Vec3b>(y, x)[1] = 255;
				m_img_plus_line.at<cv::Vec3b>(y, x)[2] = 255;
			}
		}

	}
	viewport()->update();
}
