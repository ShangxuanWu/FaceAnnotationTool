#ifndef _IMAGE_VIEW_H_
#define _IMAGE_VIEW_H_

#define POINT_MODE 0
#define LINE_MODE 1
#define REGION_MODE 2

#include <opencv2/opencv.hpp>

#include <QtGui>
#include <QGraphicsView>
#include <QDebug>
#include <QGraphicsTextItem>
#include <QFont>
#include <QString>

#include <string>
#include <list>
#include <vector>
#include <cmath>

#include "utils.h"
#include "../camera_propagation/pts_2d_tool.h"

struct OnePointAnnotation
{
	std::string absFilePath;
	//std::string visAbsFilePath;
	cv::Point2d Points;

	OnePointAnnotation(
		const std::string& abs_file_path,
		const cv::Point2d& points)
		:
		absFilePath(abs_file_path)
		, Points(points)
	{}
};

struct OneCurveAnnotation
{
	std::string absFilePath;
	//std::string visAbsFilePath;
	std::vector<cv::Vec2d> curvePoints;

	OneCurveAnnotation(
		const std::string& abs_file_path,
		const std::vector<cv::Vec2d>& curve_points)
		:
		absFilePath(abs_file_path)
		, curvePoints(curve_points)
	{}
};

class ImageView : public QGraphicsView
{
	Q_OBJECT

public:
	ImageView(QWidget *parent);
	~ImageView(void);
	void setImage(const std::string& image_fpath);

	void setDetectedKeypoint(std::vector<cv::Point> detected_landmarks);
	void ImageView::setDetectedKeypoint(std::map<int, pts_2d_tool> detected_landmarks);

	bool hasImage() const { return !(m_img.empty()); }

	const std::vector<cv::Vec2d>& currCurveInImage()  const { return currCurveInImage_; }

	void clearSelectedPoints()
	{
		this->initializeClickedPoints_();
	}

	void loadPastCurveAnnotationData(const std::vector<std::string>& annotation_data_fpaths);

	void deleteLatestAnnotation();

	void resetRotation();
	void rotate90Clockwise();
	void rotate90CounterClockwise();
	void resetToCenterAndScale(int orientation);

	void saveLatestAnnotationData();

	void setVisualizePastPoints(bool set)
	{
		visualizePastCurves_ = set;
		viewport()->update();
	}
	bool visualizePastPoints() const { return visualizePastCurves_; }
	
	bool isSelectionMode() const { return isSelectionMode_; }
	void setSelectionMode(bool set)
	{
		isSelectionMode_ = set;
		selectedIndex = -1;
		selectedPastAnnotation_iter_ = pastCurvesInImage_.end();
		if (isSelectionMode_)
			currCurveInImage_.clear();

		viewport()->update();
	}

	void searchClosestLineAnnotation(const cv::Vec2d& clicked_pos, const double thresh_dist2 = 2.5*2.5);

	int searchClosestPointAnnotation(const cv::Vec2d& clicked_pos, bool isDraggingFlag = true);

	void draw_result(std::vector<std::vector<cv::Point2d>> points);

	void set_draw_mode_(int draw_mode_number);

	std::vector<cv::Point> this_landmarks; // landmark equals -1 if not detected
	
	std::string folder_name;

	std::vector<int> modified_points;

	void resetAll();

	void toggleAnnotationsShown();

	void setVisibility(std::vector<bool> this_view_visibility);
	std::vector<bool> getVisibility();

	std::string getCameraName();

signals:
	void setHint(int point_id);
	void cancelHint();
	//void setHintOnMainStatusBar(int point_d)

protected:
	void paintEvent(QPaintEvent *event);
	
	cv::Mat m_img;
	std::string dataDirPrefix_;
	std::string imagePath_;
	std::string imageNameNoExt_;
	std::string dataDir_;

	double currScale_;
	double currRot_deg_;
	cv::Vec2d transVec_;

	QTransform m_matrix;
	QTransform m_matrix_inv;

	QPoint m_pntDownPos;

	void mouseMoveEvent(QMouseEvent* event);
	void mousePressEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);

	void wheelEvent(QWheelEvent* event);

	void initializeClickedPoints_();

	cv::Point2d currPointInImage_;
	std::vector<cv::Vec2d> currCurveInImage_;
	std::vector<cv::Vec2d> currRegionInImage_;

	bool can_move_point_in_2D;

	bool visualizePastPoints_;
	bool visualizePastCurves_;
	bool visualizePastRegions_;
	bool visualizeConnectedLines_;
	bool visualze_moved_points;

	std::list<OneCurveAnnotation> pastCurvesInImage_;
	std::list<OnePointAnnotation> pastPointsInImage_;

	int getNextCurveIndex_() const;
	int getNextPointIndex_() const;

	// for the selection
	bool isSelectionMode_;
	bool isSelectedOnePoint;
	int selectedIndex;
	// for the drawing modes
	int draw_mode_;
	
	bool showAnnotations;

	std::list<OneCurveAnnotation>::iterator selectedPastAnnotation_iter_;

	bool isDragging;

	int pointIdxThatShowsAnnotation;

	std::vector<bool> visibility;
	
};

#endif  // _IMAGE_VIEW_H_