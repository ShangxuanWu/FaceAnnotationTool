#include "MeshVisibilityHelper.h"

#include <vtkInteractorStyleTrackballCamera.h>

// a vtk interactor that does nothing 
class NoInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:

	static NoInteractorStyle* New();
	vtkTypeMacro(NoInteractorStyle, vtkInteractorStyleTrackballCamera);

	NoInteractorStyle()
	{
		;
	}

	virtual void OnLeftButtonDown()
	{
		;
	}
	virtual void OnRightButtonDown()
	{
		;
	}
	virtual void 	OnMouseWheelForward()
	{
		;
	}

	virtual void 	OnMouseWheelBackward()
	{
		;
	}
};

vtkStandardNewMacro(NoInteractorStyle);
////////////////////////////////////////////////////

// setup the render process
MeshVisibilityHelper::MeshVisibilityHelper(QWidget *parent)
	: QWidget(parent)
{
	qvtk_widget = new QVTKWidget(this);
	poly_data = vtkSmartPointer<vtkPolyData>::New();

	layout = new QHBoxLayout(this);
	layout->addWidget(qvtk_widget);
	layout->setContentsMargins(0, 0, 0, 0);

	// Create an actor and mapper
	polydata_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	polydata_mapper->SetInputData(poly_data);

	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(polydata_mapper);

	// Create a renderer, render window, and interactor
	renderer = vtkSmartPointer<vtkRenderer>::New();

	render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	render_window_interactor->SetRenderWindow(qvtk_widget->GetRenderWindow());
	render_window_interactor->Initialize();

	//style = vtkSmartPointer<NoInteractorStyle>::New();
	style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	style->SetDefaultRenderer(renderer);
	render_window_interactor->SetInteractorStyle(style);
	

	// Add the actors to the scene
	renderer->AddActor(actor);

	renderer->SetBackground(0, 0, 0); // Black background
	renderer->ResetCamera();

	select_visible_points = vtkSmartPointer<vtkSelectVisiblePoints>::New();
	select_visible_points->SetRenderer(renderer);

	
	qvtk_widget->GetRenderWindow()->AddRenderer(renderer);
	qvtk_widget->update();
}


void MeshVisibilityHelper::setMesh(vtkSmartPointer<vtkPolyData> poly_data_)
{
	poly_data = poly_data_;
	polydata_mapper->SetInputData(poly_data);
	polydata_mapper->Update();

	renderer->ResetCamera();

	qvtk_widget->update();
	select_visible_points->SetInputData(poly_data);
	select_visible_points->Update();
}

void MeshVisibilityHelper::setCamera(mycamera camera_)
{
	this_mycamera = camera_;
	vtkSmartPointer<vtkCamera> now_camera = getNewVTKCameraByCameraMatrix();
	renderer->SetActiveCamera(now_camera);
	
	qvtk_widget->update();
	select_visible_points->Update();
}

/*
std::vector<bool> MeshVisibilityHelper::getVisibility(std::vector<vtkIdType> point_ids)
{
	select_visible_points->GetOutput()->GetNumberOfPoints();
	std::vector<bool> result;
	for (int i = 0; i < point_ids.size(); i++)
	{
		result.push_back(true);
	}
	return result;
}
*/

// get new vtkCamera by some maths
vtkSmartPointer<vtkCamera> MeshVisibilityHelper::getNewVTKCameraByCameraMatrix()
{
	vtkSmartPointer<vtkCamera> new_camera = vtkSmartPointer<vtkCamera>::New();

	/////////////////////////////////////////////////////// Start here ////////////////////
	double extrinsic_matrix[16] = {
		this_mycamera.extrinsic.at<double>(0,0),
		this_mycamera.extrinsic.at<double>(0,1),
		this_mycamera.extrinsic.at<double>(0,2),
		this_mycamera.extrinsic.at<double>(0,3),
		this_mycamera.extrinsic.at<double>(1,0),
		this_mycamera.extrinsic.at<double>(1,1),
		this_mycamera.extrinsic.at<double>(1,2),
		this_mycamera.extrinsic.at<double>(1,3),
		this_mycamera.extrinsic.at<double>(2,0),
		this_mycamera.extrinsic.at<double>(2,1),
		this_mycamera.extrinsic.at<double>(2,2),
		this_mycamera.extrinsic.at<double>(2,3),
		0,
		0,
		0,
		1
	};

	new_camera->SetModelTransformMatrix(extrinsic_matrix);
	// the camera can stay at the origin because we are transforming the scene objects
	new_camera->SetPosition(0, 0, 0);
	// look in the +Z direction of the camera coordinate system
	new_camera->SetFocalPoint(0, 0, 1);
	// the camera Y axis points down
	new_camera->SetViewUp(0, -1, 0);

	// ensure the relevant range of depths are rendered
	new_camera->SetClippingRange(1, 1000);

	// convert the principal point to window center (normalized coordinate system) and set it
	double wcx = -2 * (this_mycamera.intrinsic.at<double>(0, 2) - double(ORIGINAL_IMAGE_WIDTH ) / 2) / ORIGINAL_IMAGE_WIDTH;
	double wcy = 2 * (this_mycamera.intrinsic.at<double>(1, 2) - double(ORIGINAL_IMAGE_HEIGHT) / 2) / ORIGINAL_IMAGE_HEIGHT;
	new_camera->SetWindowCenter(wcx, wcy);

	// convert the focal length to view angle and set it
	double view_angle = 57.295 * (2.0 * std::atan2(ORIGINAL_IMAGE_HEIGHT / 2.0, (this_mycamera.intrinsic.at<double>(0,0) + this_mycamera.intrinsic.at<double>(1, 1))/2.0));
	//std::cout << "view_angle = " << view_angle << std::endl;
	new_camera->SetViewAngle(view_angle);

	return new_camera;
}

std::vector<bool> MeshVisibilityHelper::getThisViewVisibility(std::vector<pts_3d_conf> points)
{
	style->OnMouseWheelBackward();
	select_visible_points->Update();
	qDebug() << "Number of visible points: " << select_visible_points->GetOutput()->GetNumberOfPoints();
	std::vector<bool> result;
	for (int i = 0; i < points.size(); i++)
	{
		bool this_point_visibility = judgeIfPointIsNearPolyData(points[i].x, points[i].y, points[i].z);
		result.push_back(this_point_visibility);
		qDebug() << "Point: " << (i + 1);
	}
	return result;
}

bool MeshVisibilityHelper::judgeIfPointIsNearPolyData(float x, float y, float z)
{
	select_visible_points->Update();
	cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(select_visible_points->GetOutput());
	cellLocator->BuildLocator();

	double testPoint[3] = { x, y, z };

	//Find the closest points to TestPoint
	double closestPoint[3];//the coordinates of the closest point will be returned here
	double closestPointDist2; //the squared distance to the closest point will be returned here
	vtkIdType cell_id; //the cell id of the cell containing the closest point will be returned here
	int sub_id; //this is rarely used (in triangle strips only, I believe)
	cellLocator->FindClosestPoint(testPoint, closestPoint, cell_id, sub_id, closestPointDist2);
	qDebug() << closestPointDist2;
	//cellLocator->FindCell(testPoint, 0.1, );
	
	if ( closestPointDist2 < FIND_POINT_VISIBILE_TOLERANCE)
	{
		return true;
	}
	else
	{
		return false;
	}
}



