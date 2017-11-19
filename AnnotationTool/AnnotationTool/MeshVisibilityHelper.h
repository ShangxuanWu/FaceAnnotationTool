#ifndef MESHVISIBILITYHELPER_H
#define MESHVISIBILITYHELPER_H

//#include "vtkAutoInit.h" 
// following sentence is for suppressing the vtk stdout window
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingOpenGL);
//VTK_MODULE_INIT(vtkRenderingFreeType);

#include <vector>

#include <QWidget>
#include <QVBoxLayout>

#include <vtkSmartPointer.h>
#include <vtkCamera.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkObjectFactory.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCellLocator.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <QVTKWidget.h>

#include "utils.h"
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/pts_3d_conf.h"

// this is a QVTKWidget that is used to get the visibility of some specific point

// class forwarding
class NoInteractorStyle;

// real function
class MeshVisibilityHelper : public QWidget
{
	Q_OBJECT
public:
	MeshVisibilityHelper(QWidget *parent = 0);
	std::vector<bool> getThisViewVisibility(std::vector<pts_3d_conf> points);
	bool judgeIfPointIsNearPolyData(float x, float y, float z);

	public slots:
	void setCamera(mycamera camera_);
	void setMesh(vtkSmartPointer<vtkPolyData> poly_data_);

protected:
	vtkSmartPointer<vtkCamera> getNewVTKCameraByCameraMatrix();

	QHBoxLayout* layout;

	mycamera this_mycamera;

	QVTKWidget* qvtk_widget;
	vtkSmartPointer<vtkPolyData> poly_data;
	vtkSmartPointer<vtkPolyDataMapper> polydata_mapper;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor;
	//vtkSmartPointer<NoInteractorStyle> style;
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style;
	vtkSmartPointer<vtkSelectVisiblePoints> select_visible_points;
	vtkSmartPointer<vtkCellLocator> cellLocator;
};

#endif
