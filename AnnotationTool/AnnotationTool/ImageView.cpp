#include "ImageView.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <direct.h>

#include "utils.h"

using namespace std;

ImageView::ImageView(QWidget* parent)
	: QGraphicsView(parent), dataDirPrefix_("HairData_")
{
	this->setMouseTracking(true);

	this->initializeClickedPoints_();
	this_landmarks = std::vector<cv::Point>(TOTAL_POINTS_IN_SINGLE_VIEW, cv::Point(-1, -1));

	currScale_ = 1;
	currRot_deg_ = 0;
	transVec_[0] = 0;
	transVec_[1] = 0;

	pointIdxThatShowsAnnotation = -1;

	can_move_point_in_2D = false;

	visualizePastPoints_ = false;
	visualizePastCurves_ = false;
	visualizePastRegions_ = false;
	visualizeConnectedLines_ = true;
	visualze_moved_points = false;

	isSelectionMode_ = false;
	isSelectedOnePoint = false;
	selectedIndex = -1;

	isDragging = false;

	showAnnotations = true;
}

void ImageView::resetAll()
{
	m_img.release();
	this_landmarks = std::vector<cv::Point>(TOTAL_POINTS_IN_SINGLE_VIEW, cv::Point(-1, -1));
}

ImageView::~ImageView()
{

}

void ImageView::initializeClickedPoints_()
{
	currPointInImage_ = cv::Point2d(-1, -1);
	currCurveInImage_.clear();
	currRegionInImage_.clear();

	pastCurvesInImage_.clear();
	selectedPastAnnotation_iter_ = pastCurvesInImage_.end();

	viewport()->update();
}

// legacy function, useless.
void ImageView::deleteLatestAnnotation()
{
	if (!m_img.empty())
	{
		// if it is annotation mode...
		if (!isSelectionMode_)
		{
			if (!pastCurvesInImage_.empty())
			{
				//--currHairIndex_;
				//pastCurvesInImage_.erase(pastCurvesInImage_.end() - 1);

				string file_to_delete = pastCurvesInImage_.back().absFilePath;
				ifstream target_file_content(file_to_delete);
				if (target_file_content.fail())
				{
					std::cout << "WARNING: target annotation cannot be found..." << endl;
					target_file_content.close();
				}
				else
				{
					target_file_content.close();
					std::remove(file_to_delete.c_str());
					std::cout << "annotation data deleted: " << file_to_delete << endl;
				}
				pastCurvesInImage_.erase(--pastCurvesInImage_.end());
				viewport()->update();

			}
			else
				std::cout << "WARNING: No more curve data..." << endl;
		}
		else // if it is selection mode, delete selected data
		{
			// if any annotation data is selected
			if (selectedPastAnnotation_iter_ != pastCurvesInImage_.end())
			{
				string file_to_delete = selectedPastAnnotation_iter_->absFilePath;
				ifstream target_file_content(file_to_delete);
				if (target_file_content.fail())
				{
					cout << "WARNING: target annotation cannot be found..." << endl;
					target_file_content.close();
				}
				else
				{
					target_file_content.close();
					std::remove(file_to_delete.c_str());
					cout << "annotation data deleted: " << file_to_delete << endl;
				}
				pastCurvesInImage_.erase(selectedPastAnnotation_iter_);
				selectedPastAnnotation_iter_ = pastCurvesInImage_.end();

				viewport()->update();
			}

		}
	}
	else
		cout << "WARNING: Load image first..." << endl;
}


void ImageView::setImage(const std::string& image_fpath)
{
	//m_img = image.clone();
	m_img = cv::imread(image_fpath);
	imagePath_ = image_fpath;

	cv::cvtColor(m_img, m_img, CV_BGR2RGB);

	m_matrix.setMatrix(1, 0, 0, 0, 1, 0, 0, 0, 1);
	m_matrix_inv = m_matrix.inverted();

	// initialization for points
	this->initializeClickedPoints_();
	modified_points.clear();
	/*
	// get annotation directory name
	size_t slash_found = imagePath_.find_last_of("/");
	size_t dot_found = imagePath_.find_last_of(".");
	size_t image_name_size = dot_found - (slash_found + 1);
	imageNameNoExt_ = imagePath_.substr(slash_found + 1, image_name_size);
	//cout << "temp_dir_name: " << temp_dir_name << endl;
	dataDir_ = imagePath_.substr(0, slash_found) + "/" + dataDirPrefix_ + imageNameNoExt_;
	cout << "dataDir_: " << dataDir_ << endl;;

	_mkdir(dataDir_.c_str());
	*/

	std::string ID = getIDFromFilename(image_fpath);
	if (rotation_indicator[ID] == RIGHT)
	{
		resetRotation();
		resetToCenterAndScale(RIGHT);
		rotate90Clockwise();
	}
	else if (rotation_indicator[ID] == LEFT)
	{
		resetRotation();
		resetToCenterAndScale(LEFT);
		rotate90CounterClockwise();
	}
	else
	{
		resetRotation();
		resetToCenterAndScale(MIDDLE);
	}
	viewport()->update();
}

void ImageView::resetRotation()
{
	currRot_deg_ = 0;
}

void ImageView::resetToCenterAndScale(int orientation)
{
	if (orientation == LEFT)
	{
		transVec_[0] = 0;
		transVec_[1] = 1550;
	}
	else if (orientation == RIGHT)
	{
		transVec_[0] = 1150;
		transVec_[1] = -75;
	}
	else
	{
		transVec_[0] = 0;
		transVec_[1] = -100;
	}
	currScale_ = 0.3;
}

// legacy function, useless.
void ImageView::saveLatestAnnotationData()
{
	if (!m_img.empty())
	{
		// saving curves
		if (draw_mode_ == LINE_MODE && !currCurveInImage_.empty() && !isSelectionMode_)
		{
			int curve_index = this->getNextCurveIndex_();

			char out_file_path[1024];
			sprintf(out_file_path, "%s_curve%09d.txt", (dataDir_ + "/" + imageNameNoExt_).c_str(), curve_index);

			cout << "out_file_path: " << out_file_path << endl;

			ofstream out_curve_data(out_file_path);
			out_curve_data << currCurveInImage_.size() << endl;
			for (auto iter = currCurveInImage_.begin(); iter != currCurveInImage_.end(); ++iter)
			{
				out_curve_data << (*iter)[0] << " " << (*iter)[1] << endl;
			}

			out_curve_data.close();

			pastCurvesInImage_.push_back(
				OneCurveAnnotation(
					out_file_path,
					currCurveInImage_)
			);
			currCurveInImage_.clear();

			selectedPastAnnotation_iter_ = pastCurvesInImage_.end();

			viewport()->update();
		}

		// saving points
		if (draw_mode_ == POINT_MODE && currPointInImage_.x >= 0 && !isSelectionMode_)
		{

			int point_index = this->getNextPointIndex_();

			qDebug() << point_index << endl;

			char out_file_path[1024];
			sprintf(out_file_path, "%s_point%09d.txt", (dataDir_ + "/" + imageNameNoExt_).c_str(), point_index);

			qDebug() << "out_file_path: " << out_file_path << endl;

			ofstream out_curve_data(out_file_path);

			out_curve_data << currPointInImage_.x << " " << currPointInImage_.y << endl;

			out_curve_data.close();

			pastPointsInImage_.push_back(OnePointAnnotation(
				out_file_path,
				currPointInImage_));

			currPointInImage_ = cv::Point2d(-1, -1);

			//selectedPastAnnotation_iter_ = pastCurvesInImage_.end();

			viewport()->update();
		}
	}
}

// legacy function, useless.
void ImageView::loadPastCurveAnnotationData(const std::vector<std::string>& annotation_data_fpaths)
{
	//	cout << "Warning: Need to implement this function!" << endl;
	if (!m_img.empty())
	{

		pastCurvesInImage_.clear();
		currCurveInImage_.clear();


		const size_t curve_counts = annotation_data_fpaths.size();

		for (size_t c = 0; c < curve_counts; ++c)
		{
			//char curve_file_path[1024];
			//sprintf(curve_file_path, "%s_curve%09d.txt", (dataDir_ + "/" + imageNameNoExt_).c_str(), h);

			const string curve_file_path = annotation_data_fpaths[c];

			ifstream curve_file_content(curve_file_path);

			if (curve_file_content.fail())
			{
				cerr << "WARNING: annotation file does not exist: " << curve_file_path << endl;

			}
			else
			{
				cerr << "Now loading: " << curve_file_path;
				// load # of points
				int num_points;
				curve_file_content >> num_points;

				vector<cv::Vec2d> curr_curve_data;
				for (int p = 0; p < num_points; ++p)
				{
					cv::Vec2d new_point;
					curve_file_content >> new_point[0] >> new_point[1];

					curr_curve_data.push_back(new_point);
				}

				pastCurvesInImage_.push_back(OneCurveAnnotation(curve_file_path, curr_curve_data));

				std::cout << " DONE! - " << pastCurvesInImage_.back().curvePoints.size() << endl;

			}

		}

		viewport()->update();

	}
}

void ImageView::paintEvent(QPaintEvent *event)
{
	// if no content in this panel
	if (m_img.empty())
	{
		return;
	}

	// update m_matrix
	cv::Mat_<double> rot_mat(2, 2);
	//qDebug() << currRot_deg_;
	const double cos_angle = cos(M_PI * currRot_deg_ / 180.);
	const double sin_angle = sin(M_PI * currRot_deg_ / 180.);
	rot_mat(0, 0) = cos_angle;
	rot_mat(0, 1) = sin_angle;
	rot_mat(1, 0) = -sin_angle;
	rot_mat(1, 1) = cos_angle;

	m_matrix.setMatrix(
		currScale_ * rot_mat(0, 0), currScale_* rot_mat(0, 1), 0,
		currScale_ * rot_mat(1, 0), currScale_ * rot_mat(1, 1), 0,
		transVec_[0], transVec_[1], 1);

	m_matrix_inv = m_matrix.inverted();

	QPainter widgetpainter(viewport());

	widgetpainter.setWorldTransform(m_matrix);


	QImage qimg = QImage(m_img.ptr(), m_img.cols, m_img.rows, m_img.step, QImage::Format_RGB888);

	widgetpainter.drawImage(0, 0, qimg);

	// following are code for drawing annotations
	if (showAnnotations)
	{
		const double curr_line_width = 2. / currScale_;

		widgetpainter.setPen(QPen(Qt::red, curr_line_width));
		QFont font;
		font.setPixelSize(QPAINTER_FONT_SIZE);
		widgetpainter.setFont(font);

		// draw detected points
		if (!m_img.empty() && !this_landmarks.empty()) {
			for (int i = 0; i < this_landmarks.size(); i++) {
				if (visibility[i])
				{
					// draw point 
					widgetpainter.setPen(QPen(Qt::cyan, curr_line_width * 2));
					widgetpainter.drawEllipse(this_landmarks[i].x - 10, this_landmarks[i].y - 10, 20, 20);
				}
			}
		}

		// highlight moved points
		if (visualze_moved_points)
		{
			for (int i = 0; i < modified_points.size(); i++)
			{
				int idx = modified_points[i];
				widgetpainter.setPen(QPen(Qt::yellow, curr_line_width * 2));
				widgetpainter.drawEllipse(this_landmarks[idx].x - 10, this_landmarks[idx].y - 10, 20, 20);
			}
		}

		// highlight selected point
		if (selectedIndex >= 0) {
			if (visibility[selectedIndex])
			{
				widgetpainter.setPen(QPen(Qt::yellow, curr_line_width * 2));
				widgetpainter.drawEllipse(this_landmarks[selectedIndex].x - 10, this_landmarks[selectedIndex].y - 10, 20, 20);
			}
		}

		// draw text of nearest point
		if (pointIdxThatShowsAnnotation != -1)
		{
			if (visibility[pointIdxThatShowsAnnotation])
			{
				widgetpainter.setPen(QPen(Qt::yellow, curr_line_width * 2));
				widgetpainter.drawEllipse(this_landmarks[pointIdxThatShowsAnnotation].x - 10, this_landmarks[pointIdxThatShowsAnnotation].y - 10, 20, 20);

				widgetpainter.setPen(QPen(Qt::yellow, curr_line_width));
				widgetpainter.drawText(this_landmarks[pointIdxThatShowsAnnotation].x - 25, this_landmarks[pointIdxThatShowsAnnotation].y - 25, QString::fromStdString(to_string(pointIdxThatShowsAnnotation + 1)));
			}
		}

		// draw connecting lines
		if (visualizeConnectedLines_ && TOTAL_POINTS_IN_SINGLE_VIEW == 68)
		{
			widgetpainter.setPen(QPen(Qt::green, curr_line_width / 3));
			for (auto i = connected_points.begin(); i != connected_points.end(); i++)
			{
				int index_1 = i->first - 1;
				int index_2 = i->second - 1;
				if (visibility[index_1] && visibility[index_2])
				{
					widgetpainter.drawLine(this_landmarks[index_1].x, this_landmarks[index_1].y, this_landmarks[index_2].x, this_landmarks[index_2].y);
				}
			}
		}
	}
}

// legacy function, useless.
void ImageView::rotate90Clockwise()
{
	qDebug() << "Rotate 90 degree clockwise..." << endl;
	currRot_deg_ += 90;
	if (currRot_deg_ > 360)
		currRot_deg_ -= 360;
	viewport()->update();
}

// legacy function, useless.
void ImageView::rotate90CounterClockwise()
{
	qDebug() << "Rotate 90 degree counter-clockwise..." << endl;
	currRot_deg_ -= 90;

	if (currRot_deg_ < 0)
		currRot_deg_ += 360;
	viewport()->update();
}

// legacy function, useless.
void ImageView::draw_result(vector<vector<cv::Point2d>> points)
{
	for (int i = 0; i < 2; i++) {
		const string a = "C:\\Users\\shangxuanu\\Desktop\\MyAnnotationTool\\AnnotationTool_Oculus\\AnnotationTool_Oculus\\30.png";
		cv::Point2d point = points[15][i];
		std::vector<cv::Vec2d> tmpVec;
		tmpVec.push_back(cv::Vec2d(point));
		OneCurveAnnotation tmp(a, tmpVec);
		pastCurvesInImage_.push_back(tmp);
	}
	viewport()->update();
}

// used for draging points, and showing hints on thei right-hand window
void ImageView::mouseMoveEvent(QMouseEvent *event)
{
	QPoint pnt = event->pos();

	if (!m_img.empty())
	{

		//if (!showAnnotations)
		//{
		//	return;
		//}
		//qDebug() << "Moving!!!";
		// find nearest point label
		QPointF pnt_inversed = m_matrix_inv.map(QPointF(pnt.x(), pnt.y()));
		cv::Vec2d pnt_cv_vec2d(pnt_inversed.x(), pnt_inversed.y());
		pointIdxThatShowsAnnotation = searchClosestPointAnnotation(pnt_cv_vec2d, false);
		qDebug() << pointIdxThatShowsAnnotation;

		// for right button translation
		if (event->buttons() == Qt::LeftButton)
		{
		//{
			QPointF pntf = (pnt - m_pntDownPos) / currScale_;
			m_pntDownPos = event->pos();
			m_matrix.translate(pntf.x(), pntf.y());

			transVec_[0] += pntf.x() * currScale_;
			transVec_[1] += pntf.y() * currScale_;

		//}
		// for left button
		//if (event->buttons() == Qt::LeftButton)
		//{
			/*
			// for annotation mode
			if (!isSelectionMode_)
			{
				if (draw_mode_ == POINT_MODE) {
					;
				}
				else if (draw_mode_ == LINE_MODE)
				{
					QPointF new_click_pos_in_image = m_matrix_inv.map(QPointF(event->pos().x(), event->pos().y()));
					currCurveInImage_.push_back(cv::Vec2d(new_click_pos_in_image.x(), new_click_pos_in_image.y()));
				}
				else {
					;
				}
			}
			else {*/
			// for selection mode
			if (draw_mode_ == POINT_MODE) {
				QPointF new_click_pos_in_image = m_matrix_inv.map(QPointF(event->pos().x(), event->pos().y()));
				if (isDragging == false)
				{
					selectedIndex = searchClosestPointAnnotation(cv::Vec2d(new_click_pos_in_image.x(), new_click_pos_in_image.y()));
				}
				if (selectedIndex >= 0 && can_move_point_in_2D)
				{
					this_landmarks[selectedIndex].x = new_click_pos_in_image.x();
					this_landmarks[selectedIndex].y = new_click_pos_in_image.y();
				}
			}

		}
	}
	viewport()->update();

	QWidget::mouseMoveEvent(event);
}

void ImageView::mousePressEvent(QMouseEvent *event)
{
	m_pntDownPos = event->pos();

	// if there is no image in the current ImageView window
	if (m_img.empty())
	{
		return;
	}
	if (!showAnnotations)
	{
		return;
	}

	// if there is an image in the current ImageView window
	QPointF new_click_pos_in_image = m_matrix_inv.map(QPointF(event->pos().x(), event->pos().y()));

	if (!isSelectionMode_)
	{
		if (event->buttons() == Qt::LeftButton)
		{
			if (draw_mode_ == POINT_MODE) {
				QPointF new_click_pos_in_image = m_matrix_inv.map(QPointF(event->pos().x(), event->pos().y()));
				cv::Point2d tmpPoint(new_click_pos_in_image.x(), new_click_pos_in_image.y());
				currPointInImage_ = tmpPoint;
			}
			/*
			else if (draw_mode_ == LINE_MODE) {
				currCurveInImage_.clear();

				//QPointF new_click_pos_in_image = m_matrix_inv.map(event->pos());

				currCurveInImage_.push_back(cv::Vec2d(new_click_pos_in_image.x(), new_click_pos_in_image.y()));
			}
			else { // region mode
				;
			}
			*/
		}
	}
	else
	{
		if (event->buttons() == Qt::LeftButton) { // right button is for translation
			if (draw_mode_ == LINE_MODE) { // select line
				this->searchClosestLineAnnotation(cv::Vec2d(new_click_pos_in_image.x(), new_click_pos_in_image.y()));
			}
			if (draw_mode_ == POINT_MODE) { // select point
				/*if (!isSelectedOnePoint) { // didn't select a point
					selectedIndex = searchClosestPointAnnotation(cv::Vec2d(new_click_pos_in_image.x(), new_click_pos_in_image.y()));
					qDebug() << selectedIndex << endl;
					isSelectedOnePoint = true;
				}
				else {
					this_landmarks[selectedIndex].x = new_click_pos_in_image.x();
					this_landmarks[selectedIndex].y = new_click_pos_in_image.y();
					isSelectedOnePoint = false;
				}*/
			}
		}
	}
	QWidget::mousePressEvent(event);
}

void ImageView::wheelEvent(QWheelEvent* event)
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

void ImageView::mouseReleaseEvent(QMouseEvent* event)
{
	selectedIndex = -1;
	if (!showAnnotations)
	{
		return;
	}
	isDragging = false;
	viewport()->update();
}

int ImageView::searchClosestPointAnnotation(const cv::Vec2d& clicked_pos, bool isDraggingFlag)
{
	if (isDraggingFlag)
	{
		isDragging = true;
	}
	if (!m_img.empty() && isSelectionMode_)
	{
		if (draw_mode_ == POINT_MODE) {

			qDebug() << "now selecting points!" << endl;

			// iterate through all the points to find the closest one
			std::vector<float> distance(TOTAL_POINTS_IN_SINGLE_VIEW, -1);
			for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++) {
				if (this_landmarks[i].x != -1) {
					distance[i] = pow(clicked_pos[0] - this_landmarks[i].x, 2) + pow(clicked_pos[1] - this_landmarks[i].y, 2);
				}
			}
			// get the smallest distance
			float smallest = INT_MAX;
			int smallest_idx = 0;
			for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++) {
				if (distance[i] >= 0) {
					if (distance[i] < smallest) {
						smallest_idx = i;
						smallest = distance[i];
					}
				}
			}
			// judge whether drag distance is smaller than what we want
			if (smallest < MAX_SELECT_DISTANCE)
			{
				if (isDraggingFlag)
				{
					modified_points.push_back(smallest_idx);
				}
				return smallest_idx;
			}
			else
			{
				return -1;
			}
		}
	}

	viewport()->update();
}

void ImageView::searchClosestLineAnnotation(const cv::Vec2d& clicked_pos, const double thresh_dist2)
{
	if (isSelectionMode_)
	{
		cout << "now selecting!" << endl;

		selectedPastAnnotation_iter_ = pastCurvesInImage_.end();

		cv::Vec2d diff_vec = cv::Vec2d(9999999, 9999999);
		double min_diff2 = 99999999;
		for (auto iter = pastCurvesInImage_.begin(); iter != pastCurvesInImage_.end(); ++iter)
		{
			const vector<cv::Vec2d>& curr_curve_pts = iter->curvePoints;

			// point-to-line
			for (size_t p = 0; p < curr_curve_pts.size() - 1; ++p)
			{
				cv::Vec2d dir_vec = curr_curve_pts[p + 1] - curr_curve_pts[p];
				const double line_length2 = cv::norm(dir_vec, cv::NORM_L2);
				dir_vec = cv::normalize(dir_vec);

				diff_vec = clicked_pos - curr_curve_pts[p];

				cv::Vec2d parallel_component = diff_vec.ddot(dir_vec) * dir_vec;
				const double parallel_length2 = cv::norm(parallel_component, cv::NORM_L2);

				cv::Vec2d orthogonal_component = diff_vec - parallel_component;
				const double pt2line_dist2 = cv::norm(orthogonal_component, cv::NORM_L2);

				if (pt2line_dist2 < thresh_dist2
					&& pt2line_dist2 < min_diff2
					&& parallel_length2 < line_length2
					&& dir_vec.ddot(parallel_component) > 0)
				{
					min_diff2 = pt2line_dist2;
					selectedPastAnnotation_iter_ = iter;
				}
			}
		}

		// debug
		if (selectedPastAnnotation_iter_ == pastCurvesInImage_.end())
		{
			cout << "no close annotation data..." << endl;
		}
		else
		{
			cout << "FOUND: " << selectedPastAnnotation_iter_->absFilePath << endl;
		}
	}

	viewport()->update();
}

// legacy function
int ImageView::getNextCurveIndex_() const
{
	if (pastCurvesInImage_.empty())
	{
		return 0;
	}
	else
	{
		// get file name
		int curve_index;
		{
			const string last_data_fpath = pastCurvesInImage_.back().absFilePath;
			size_t slash_found = last_data_fpath.find_last_of('e'); // e of "curve"
			size_t dot_found = last_data_fpath.find_last_of(".");
			size_t image_name_size = dot_found - (slash_found + 1);
			const string last_curve_name = last_data_fpath.substr(slash_found + 1, image_name_size);

			cout << "file name: " << last_curve_name << endl;;

			cout << "curve index: " << stoi(last_curve_name) << endl;;
			curve_index = stoi(last_curve_name);
		}
		return curve_index + 1;
	}

}

// legacy code
int ImageView::getNextPointIndex_() const
{
	if (pastPointsInImage_.empty())
	{
		return 0;
	}
	else
	{
		// get file name
		int point_index;
		{

			const string last_data_fpath = pastPointsInImage_.back().absFilePath;
			size_t slash_found = last_data_fpath.find_last_of('n'); // e of "curve"
			size_t dot_found = last_data_fpath.find_last_of(".");
			size_t image_name_size = dot_found - (slash_found + 2);
			const string last_point_name = last_data_fpath.substr(slash_found + 2, image_name_size);

			qDebug() << "file name: " << last_point_name.c_str() << endl;;

			qDebug() << "curve index: " << stoi(last_point_name) << endl;;
			point_index = stoi(last_point_name);
		}
		return point_index + 1;
	}
}

void ImageView::set_draw_mode_(int draw_mode_number) {
	switch (draw_mode_number) {
	case POINT_MODE:
		draw_mode_ = POINT_MODE;
		qDebug() << "set to point mode" << endl;
		break;
	case LINE_MODE:
		draw_mode_ = LINE_MODE;
		qDebug() << "set to line mode" << endl;
		break;
	case REGION_MODE:
		draw_mode_ = REGION_MODE;
		qDebug() << "set to region mode" << endl;
		break;
	}
}

void ImageView::setDetectedKeypoint(std::vector<cv::Point> detected_landmarks) {
	for (int i = 0; i < detected_landmarks.size(); i++) {
		this_landmarks[i] = detected_landmarks[i];
	}
}

void ImageView::setDetectedKeypoint(std::map<int, pts_2d_tool> detected_landmarks)
{
	//this_landmarks.clear();

	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++) {
		this_landmarks[i] = cv::Point(detected_landmarks[i].x, detected_landmarks[i].y);
		if (detected_landmarks[i].conf == 1.0)
		{
			modified_points.push_back(i);
		}
	}
	viewport()->update();
}

void ImageView::toggleAnnotationsShown()
{
	if (showAnnotations)
	{
		showAnnotations = false;
	}
	else
	{
		showAnnotations = true;
	}
	viewport()->update();
}

void ImageView::setVisibility(vector<bool> this_view_visibility)
{
	visibility = this_view_visibility;
}

std::vector<bool> ImageView::getVisibility()
{
	return visibility;
}

std::string ImageView::getCameraName()
{
	std::string name = getFolderFromFilename(imagePath_);
	return name;
}
