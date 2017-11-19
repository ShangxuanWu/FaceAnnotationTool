#ifndef PCLVIEWER_H
#define PCLVIEWER_H

#include "vtkAutoInit.h" 
// following sentence is for suppressing the vtk stdout window
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingOpenGL);
VTK_MODULE_INIT(vtkRenderingFreeType);


#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include "utils.h"

// Qt
#include <QWidget>
#include <QPushButton>
#include <QSlider>
#include <QDebug>
#include <QSlider>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>

// Point Cloud Library
#include <pcl/io/ply_io.h>
#include <pcl/io/auto_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/conversions.h>
#include <pcl/filters/filter.h>

// Visualization Toolkit (VTK)
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <QVTKWidget.h>

#include <Eigen/Geometry>

#include "../camera_propagation/pts_3d_conf.h"

void pointPickingCallback(const pcl::visualization::PointPickingEvent& event, void* viewer_void);
void keyboardCallback(const pcl::visualization::KeyboardEvent& event, void* viewer_void);
//void mouseEventCallback(const pcl::visualization::MouseEvent& event, void* viewer_void);


class PCLViewer : public QWidget
{
	Q_OBJECT

public:
	explicit PCLViewer(QWidget *parent = 0);
	~PCLViewer();
	void hightlightPoint(int index_in_cloud, std::string name);
	void hightlightPoint(double x, double y, double z, std::string id);
	void setMesh(std::string path);

	// these functions are for the MainWindow.cpp
	//void show3DLandmarks(const std::vector<cv::Point3d>& triangulated_points);
	void set3DLandmarks(const std::vector<pts_3d_conf>& triangulated_points);
	void fitPropagatedLandmarks();
	int findNearestVertex(float x, float y, float z);
	
	// 
	void setPickedIndex(pcl::PointXYZRGBA current_point);
	int findNearestLandmark(float x, float y, float z);
	
	void onSwitchAnnotationShow();

	void onCancelPicking();

	void undo();

	void setHint(int index);

	public slots:
	void onChangeCameraView(int value);
	void onClear3DLandmarks();
	void onLoadDebugMesh();

protected:
	int render_mode;
	pcl::PLYReader Reader;

	// for MESH_MODE
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_for_mesh;
	pcl::PolygonMesh mesh;
	// for POINT_CLOUD_MODE
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_for_pointcloud;
	// for recording clicked points
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr clicked_points;

	// for the layout part
	QVBoxLayout* vertical_layout;
	QHBoxLayout* horizontal_layout;
	QVBoxLayout* button_layout;
	QVTKWidget* qvtkWidget;

	// for the resize of left panel
	void resizeEvent(QResizeEvent *event);

	// for the helper part
	//QSlider* horizontalSlider;
	//QLabel* instruction;
	QPushButton* undo_button_;
	QPushButton* cancel_button_;
	QPushButton* switch_display_button_;
	//QPushButton* find_nearest_button_;

	// a qlabel for what we want to show after mouse get close to some point
	QLabel* instruction_container;
	QLabel* image_container;
	const float R = 500.0;

	bool show_annotations;

	// for annotation
	bool has_picked_one_annotation;
	std::vector<pcl::PointXYZRGBA> original_landmarks;
	std::vector<pcl::PointXYZRGBA> fitted_landmarks;
	int now_chosen_landmark;
	//int picked_vertex_index;
	
	// for the undo function
	int last_changed_landmark_idx;
	pcl::PointXYZRGBA last_changed_landmark;


	void paintLandmarks(int width);
};

#endif // PCLVIEWER_H
