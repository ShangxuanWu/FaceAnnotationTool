#ifndef VTKVIEWER_H
#define VTKVIEWER_H

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
#include <QStatusBar>
//#include <QMainWindow>

// Visualization Toolkit (VTK)
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdTypeArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPlaneSource.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkCellLocator.h>
#include <vtkLine.h>
#include <vtkCubeSource.h>
#include <vtkCoordinate.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkHardwareSelector.h>
#include <QVTKWidget.h>

#include <Eigen/Geometry>

#include <pcl/io/ply_io.h>
#include <pcl/io/auto_io.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/conversions.h>
#include <pcl/filters/filter.h>
#include <pcl/surface/vtk_smoothing/vtk_utils.h>

#include "../camera_propagation/pts_3d_conf.h"

#include "MeshVisibilityHelper.h"
//#include "hintDialog.h"

// class forwarding
class MouseInteractorStyle;
//class MainWindow;

// real function
class VTKViewer : public QWidget
{
	Q_OBJECT

public:
	explicit VTKViewer(QWidget *parent = 0);
	~VTKViewer() {};
	//void closeDialog();
	void setMesh(std::string mesh_path);
	void pickPointCallback(vtkIdType idx);
	std::vector<pts_3d_conf> setAndShowLandmarks(std::vector<pts_3d_conf> points);
	vtkIdType findNearestLandmark(vtkIdType now_idx);
	void changePositionOfLastSelectedLandmark(vtkIdType new_position);
	std::vector<cv::Point3d> updateMainWindow();
	//void setHint(vtkIdType id);
	// visibility functions
	void refreshVisibility();
	void printVisibility();

	

	//void getAllVisiblePointsInOneView();
	//void countPoints();
	//void countPointsUsingHardwareSelector();

public slots:
	void debugSetMesh();
	//void popUpInstruction();
	void onCancel();
	void enableFineAnnotationButton();
	//void disableCancelButton();
	void toggleVisibility();

	void to2DButtonPressed_();

signals:
	void setHint(int point_id);
	void setHint3D(int point_id);
	void cancelHint();
	void cancelHint3D();
	void meshPointerSignal(vtkSmartPointer<vtkPolyData>);
	void to2DButtonPressed();
	void correctThisPointButtonPressed(int);

protected:
	//void initializePopUpInstruction();
	vtkIdType getNearestCell(float x, float y, float z);
	void highlightAllLandmarks();
	void highlightAllLandmarks_text();
	void highlightAllLandmarks_line();
	void removeRedundantActors();
	// helper functions
	void getCenterPointFromCellId(vtkIdType cell_id, double* p);

	std::vector<pts_3d_conf> setLandmarks(std::vector<pts_3d_conf> points);
	pcl::PLYReader Reader;

	// for the VTK 
	vtkSmartPointer<vtkPolyData> poly_data;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkActor> actor;
	pcl::PolygonMesh::Ptr helper_mesh;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
	vtkSmartPointer<MouseInteractorStyle> style;
	vtkSmartPointer<vtkRenderer> renderer;
	//vtkSmartPointer<vtkRenderWindow> renderWindow;
	// for rendering selected points
	vtkSmartPointer<vtkSelectionNode> selection_node;
	vtkSmartPointer<vtkSelection> selection;
	vtkSmartPointer<vtkExtractSelection> extract_selection;
	vtkSmartPointer<vtkUnstructuredGrid> selected;
	vtkSmartPointer<vtkDataSetMapper> selected_mapper;
	vtkSmartPointer<vtkActor> selected_actor;
	// for rendering lines
	vtkSmartPointer<vtkPolyData> lines_poly_data;
	vtkSmartPointer<vtkCellArray> lines;
	vtkSmartPointer<vtkPoints> pts;
	vtkSmartPointer<vtkUnsignedCharArray> colors;
	vtkSmartPointer<vtkPolyDataMapper> line_mapper;
	vtkSmartPointer<vtkActor> line_actor;
	// for rendering texts
	vtkSmartPointer<vtkFollower> follower;
	vtkSmartPointer<vtkVectorText> text_source;
	vtkSmartPointer<vtkPolyDataMapper> text_mapper;
	std::vector<vtkSmartPointer<vtkFollower>> followers;

	// for the layout part
	QVBoxLayout* vertical_layout;
	QHBoxLayout* horizontal_layout;
	QVBoxLayout* button_layout_left;
	QVBoxLayout* button_layout_right;
	QVTKWidget* qvtkWidget;

	std::string now_mesh_path;
	// for the resize of left panel

	// for the lower part	
	//QPushButton* undo_button;
	
	QPushButton* to_2d_button;
	QPushButton* correct_this_point_button;
	QPushButton* switch_num_button;
	QPushButton* switch_instruction_button;
	//QPushButton* propagate_back_to_2D_button;
	QStatusBar* status_bar;

	bool show_mesh;
	bool show_popup;

	// for annotation
	
	vtkSmartPointer<vtkIdTypeArray> all_landmarks_cell_id;
	//used for visibility
	vtkSmartPointer<vtkIdTypeArray> all_landmarks_point_id;
	
	bool has_chosen_one_annotation;
	vtkIdType now_chosen_cell_idx;
	int now_chosen_id_in_array;
	//int picked_vertex_index;

	// for the undo function
	int last_changed_landmark_idx;
	pcl::PointXYZRGBA last_changed_landmark;

	// for depth test
	
	// for buttons
	void correctThisPointButtonPressed_();

	
};
#endif