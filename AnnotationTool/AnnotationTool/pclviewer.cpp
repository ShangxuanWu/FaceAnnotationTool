#include "pclviewer.h"
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <QString>

PCLViewer::PCLViewer(QWidget *parent) :
	QWidget(parent)
{
	render_mode = NOW_RENDER_MODE;
	has_picked_one_annotation = false;
	show_annotations = true;
	now_chosen_landmark = -1;
	last_changed_landmark_idx = -1;

	vertical_layout = new QVBoxLayout(this);
	horizontal_layout = new QHBoxLayout(this);
	button_layout = new QVBoxLayout(this);

	qvtkWidget = new QVTKWidget(this);

	//qvtkWidget->setFixedHeight(1445);
	qvtkWidget->setFixedHeight(1255);

	/*horizontalSlider = new QSlider(Qt::Horizontal, this);
	horizontalSlider->setFixedWidth(200);
	horizontalSlider->setFixedHeight(30);
	horizontalSlider->setMinimum(-45);
	horizontalSlider->setValue(0);
	horizontalSlider->setMaximum(45);
	horizontalSlider->setVisible(false);*/

	/*instruction = new QLabel(this);
	instruction->setText("<a href=\"https://fb.quip.com/SFEqAGMhyyAP\">Labeling instructions</a>");
	//instruction ->setText("<a href=\"instruction.png\">Labeling instructions</a>");
	instruction->setTextFormat(Qt::RichText);
	instruction->setTextInteractionFlags(Qt::TextBrowserInteraction);
	instruction->setOpenExternalLinks(true);
	instruction->setFixedWidth(200);*/

	// Undo Button
	undo_button_ = new QPushButton("Undo", this);
	connect(undo_button_, &QPushButton::clicked, this, &PCLViewer::undo);
	
	// Cancel Button
	cancel_button_ = new QPushButton("Cancel Selection", this);
	connect(cancel_button_, &QPushButton::clicked, this, &PCLViewer::onCancelPicking);

	// Switch Display Button
	switch_display_button_ = new QPushButton("Switch Annotation Display", this);
	connect(switch_display_button_, &QPushButton::clicked, this, &PCLViewer::onSwitchAnnotationShow);

	// Find Nearest Annotation Button
	//find_nearest_button_ = new QPushButton("Find Nearest", this);
	//connect(find_nearest_button_, &QPushButton::clicked, this, &PCLViewer::onFindingNearestLandmark);


	// for the instructions of different point
	instruction_container = new QLabel(this);

	// set size for these layouts
	undo_button_->setFixedWidth(200);
	cancel_button_->setFixedWidth(200);
	switch_display_button_->setFixedWidth(200);
	instruction_container->setFixedWidth(200);


	image_container = new QLabel(this);
	// set the size of image_container to be 200 x 150
	//image_container->setFixedWidth(300);
	//image_container->setFixedHeight(150);
	// set the hint to initial state
	setHint(-1);

	// Set up the QVTK window
	viewer.reset(new pcl::visualization::PCLVisualizer("viewer", false));
	//viewer->registerMouseCallback(mouseEventOccurred, (void*)&viewer);
	//viewer->registerPointPickingCallback(pp_callback, (void*)&viewer);
	viewer->registerPointPickingCallback(pointPickingCallback, (void*)this);
	viewer->registerKeyboardCallback(keyboardCallback, (void*)this);
	//viewer->registerMouseCallback(mouseEventCallback, (void*)this);

	qvtkWidget->SetRenderWindow(viewer->getRenderWindow());
	viewer->setupInteractor(qvtkWidget->GetInteractor(), qvtkWidget->GetRenderWindow());
	viewer->setShowFPS(false);
	//viewer->addOrientationMarkerWidgetAxes(qvtkWidget->GetInteractor());
	qvtkWidget->update();

	vertical_layout->addWidget(qvtkWidget);
	vertical_layout->addLayout(horizontal_layout);
	//verticalLayout->addWidget(horizontalSlider);
	//verticalLayout->addWidget(instruction_container);
	horizontal_layout->addLayout(button_layout);
	horizontal_layout->addWidget(image_container);
	//verticalLayout->addWidget(instruction);
	//verticalLayout->addWidget(load_button_);
	//verticalLayout->addWidget(clear_button_);
	//verticalLayout->addWidget(find_nearest_button_);
	button_layout->addWidget(instruction_container);
	button_layout->addWidget(undo_button_);
	button_layout->addWidget(cancel_button_);
	button_layout->addWidget(switch_display_button_);
	vertical_layout->setAlignment(Qt::AlignTop);
	vertical_layout->setContentsMargins(0, 0, 0, 0);

	// qt slots
	//connect(horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(onChangeCameraView(int)));

	//onLoadDebugMesh();
}

//  this is an individual call-back function that is outside the PCLViewer class
void pointPickingCallback(const pcl::visualization::PointPickingEvent& event, void* viewer_void)
{
	if (event.getPointIndex() == -1)
	{
		return;
	}
	pcl::PointXYZRGBA current_point;
	event.getPoint(current_point.x, current_point.y, current_point.z);
	((PCLViewer*)viewer_void)->setPickedIndex(current_point);
}

void keyboardCallback(const pcl::visualization::KeyboardEvent& event, void* viewer_void)
{

	//std::string a =  event.getKeySym();
	// show or hide annotations
	if (event.getKeySym() == "c" && event.keyDown())
	{
		((PCLViewer*)viewer_void)->onSwitchAnnotationShow();
	}
	else if (event.getKeySym() == "Escape" && event.keyDown())
	{
		((PCLViewer*)viewer_void)->onCancelPicking();
	}
	else if (event.getKeySym() == "z" && event.keyDown())
	{
		((PCLViewer*)viewer_void)->undo();
	}
}

PCLViewer::~PCLViewer()
{
}

// called by pointPickingCallback
void PCLViewer::setPickedIndex(pcl::PointXYZRGBA current_point)
{
	// choose a clicked points
	if (!has_picked_one_annotation)
	{
		// find the nearest annotation
		int index = findNearestLandmark(current_point.x, current_point.y, current_point.z);
		qDebug() << index;
		if (index >= 0)
		{
			has_picked_one_annotation = true;
			// rendering, maybe adding one more cloud to display this point, and poping out the original point
			now_chosen_landmark = index;
			setHint(index);
		}
	}
	// or move the chosen point to a new position
	else
	{
		last_changed_landmark_idx = now_chosen_landmark;
		last_changed_landmark = fitted_landmarks[now_chosen_landmark];
		// TODO: set the picked point to the new position
		fitted_landmarks[now_chosen_landmark].x = current_point.x;
		fitted_landmarks[now_chosen_landmark].y = current_point.y;
		fitted_landmarks[now_chosen_landmark].z = current_point.z;
		// end
		has_picked_one_annotation = false;
		now_chosen_landmark = -1;
	}

	//clicked_points->points.push_back(current_point);

	paintLandmarks(POINT_SIZE_3D);
}

void PCLViewer::setMesh(std::string path) {

	if (render_mode == POINT_CLOUD_MODE)
	{
		// point cloud version
		qDebug() << "Start Loading Mesh" << endl;
		cloud_for_pointcloud.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
		pcl::io::load(path, mesh);

		pcl::fromPCLPointCloud2(mesh.cloud, *cloud_for_pointcloud);
		qDebug() << "Finish Loading Mesh" << endl;
		// show point cloud
		viewer->removePointCloud();
		viewer->addPointCloud(cloud_for_pointcloud, "cloud");
	}
	else if (render_mode == MESH_MODE)
	{
		// mesh version
		cloud_for_mesh.reset(new pcl::PointCloud<pcl::PointXYZ>);
		viewer->removePolygonMesh();
		pcl::io::load(path, mesh);
		viewer->addPolygonMesh(mesh, "meshes", 0);
		pcl::fromPCLPointCloud2(mesh.cloud, *cloud_for_mesh);
	}
	clicked_points.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
	viewer->resetCamera();
	qvtkWidget->update();
}

// the index should be regarding to the points on the mesh
void PCLViewer::hightlightPoint(int index_in_cloud, std::string name) {
	if (render_mode == MESH_MODE)
	{
		viewer->addSphere(cloud_for_mesh->points[index_in_cloud], 2, 1.0, 1.0, 0, name, 0);
	}
	else if (render_mode == POINT_CLOUD_MODE)
	{
		viewer->addSphere(cloud_for_pointcloud->points[index_in_cloud], 2, 1.0, 1.0, 0, name, 0);

	}
	qvtkWidget->update();
}

void PCLViewer::hightlightPoint(double x, double y, double z, std::string id) {
	pcl::ModelCoefficients sphere_coeff;
	sphere_coeff.values.resize(4);    // We need 4 values
	sphere_coeff.values[0] = x;
	sphere_coeff.values[1] = y;
	sphere_coeff.values[2] = z;
	// radius
	sphere_coeff.values[3] = 2;

	viewer->addSphere(sphere_coeff, id);
	pcl::PointXYZRGBA tmp;
	tmp.x = x - 10;
	tmp.y = y - 10;
	tmp.z = z - 10;
	tmp.r = 255;
	tmp.g = 255;
	tmp.b = 255;
	std::string text_id = "text" + id;
	qvtkWidget->update();
}

void PCLViewer::resizeEvent(QResizeEvent *event) {
	int wid = this->width();
	qDebug() << "MeshComtainer width adjusted! New width: " << wid;
	qvtkWidget->resize(wid, 1208);
	//viewer->resetCamera(); // here only resets zoom
}

void PCLViewer::onChangeCameraView(int value) {
	viewer->setCameraPosition(0, R * sin(value * PI / 180.0), R * cos(value * PI / 180.0), -1, 0, 0);
	viewer->setCameraClipDistances(10, 20000);
}

// please use this one
void PCLViewer::set3DLandmarks(const std::vector<pts_3d_conf>& triangulated_points)
{
	// remove previous shapes
	//viewer->removeAllShapes();
	for (int i = 0; i < triangulated_points.size(); i++)
	{
		//qDebug() << i;
		//hightlightPoint(triangulated_points[i].x, triangulated_points[i].y, triangulated_points[i].z, std::to_string(i));
		pcl::PointXYZRGBA now_point;
		now_point.x = triangulated_points[i].x;
		now_point.y = triangulated_points[i].y;
		now_point.z = triangulated_points[i].z;

		original_landmarks.push_back(now_point);
	}

	// fit all the points
	fitPropagatedLandmarks();
	// show all the landmarks
	paintLandmarks(POINT_SIZE_3D);
}

void PCLViewer::onLoadDebugMesh()
{
	std::string path = "C:\\Users\\shangxuanu\\Desktop\\image_original\\mesh.ply";
	setMesh(path);
}

void PCLViewer::onClear3DLandmarks()
{
	//viewer->removeAllShapes();
	//all_landmarks.clear();
	viewer->removePointCloud("clicked_points");
	qvtkWidget->update();
}

int PCLViewer::findNearestLandmark(float x, float y, float z)
{
	if (fitted_landmarks.size() == 0)
	{
		return -1;
	}

	std::vector<float> distance;
	for (int i = 0; i < fitted_landmarks.size(); i++)
	{
		float new_x = fitted_landmarks[i].x;
		float new_y = fitted_landmarks[i].y;
		float new_z = fitted_landmarks[i].z;
		distance.push_back(pow((x - new_x), 2) + pow((y - new_y), 2) + pow((z - new_z), 2));
	}
	float min_distance = INT_MAX;
	int min_index = -1;
	for (int i = 0; i < fitted_landmarks.size(); i++)
	{
		if (distance[i] < min_distance)
		{
			min_distance = distance[i];
			min_index = i;
		}
	}
	if (min_distance <= MAX_SELECT_DISTANCE_3D) {
		return min_index;
	}
	else
	{
		return -1;
	}
}

// fit the propagated landmarks (maybe outside the mesh) to the mesh
void PCLViewer::fitPropagatedLandmarks()
{
	if (original_landmarks.empty())
	{
		qDebug() << "No Landmarks Loaded !";
		return;
	}
	
	fitted_landmarks.clear();
	if (render_mode == MESH_MODE)
	{
		for (int i = 0; i < original_landmarks.size(); i++)
		{
			int idx = findNearestVertex(original_landmarks[i].x, original_landmarks[i].y, original_landmarks[i].z);
			pcl::PointXYZRGBA now_point;
			now_point.x = (*cloud_for_mesh)[idx].x;
			now_point.y = (*cloud_for_mesh)[idx].y;
			now_point.z = (*cloud_for_mesh)[idx].z;
			fitted_landmarks.push_back(now_point);
		}
	}
	// TODO: finish point cloud mode
}

// if no point is found, return -1
int PCLViewer::findNearestVertex(float x, float y, float z)
{
	if (render_mode == MESH_MODE)
	{
		// if no mesh loaded, return -1
		if (cloud_for_mesh->empty())
		{
			return -1;
		}

		srand(time(NULL));
		pcl::PointXYZ original_point(x, y, z);

		pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
		kdtree.setInputCloud(cloud_for_mesh);
		// K nearest neighbor search
		int K = 10;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);

		// need a search point
		if (kdtree.nearestKSearch(original_point, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
		{
			qDebug() << "    " << cloud_for_mesh->points[pointIdxNKNSearch[0]].x
			<< " " << cloud_for_mesh->points[pointIdxNKNSearch[0]].y
			<< " " << cloud_for_mesh->points[pointIdxNKNSearch[0]].z
			<< " (squared distance: " << pointNKNSquaredDistance[0] << ")";

			return pointIdxNKNSearch[0];
		}
		else // didn't get a nearest vertex
		{
			return -1;
		}
	}
	else if (render_mode == POINT_CLOUD_MODE)
	{
		// if no point cloud loaded, return -1
		if (cloud_for_pointcloud->empty())
		{
			return -1;
		}

		// also have to modify KdTreeFLANN

		return 0;
	}
	return -1;
}

void PCLViewer::onSwitchAnnotationShow()
{
	if (show_annotations)
	{
		show_annotations = false;
		paintLandmarks(POINT_MINIMAL_SIZE_3D);
	}
	else
	{
		show_annotations = true;
		paintLandmarks(POINT_SIZE_3D);
	}
	qvtkWidget->update();
}

// this is the paint function of this whole panel
void PCLViewer::paintLandmarks(int width)
{
	// put all the points now into clicked_points
	for (int i = 0; i < fitted_landmarks.size(); i++)
	{
		// green points
		fitted_landmarks[i].r = 0;
		fitted_landmarks[i].g = 255;
		fitted_landmarks[i].b = 0;
	}

	if (now_chosen_landmark >= 0)
	{
		// cyan points
		fitted_landmarks[now_chosen_landmark].r = 0;
		fitted_landmarks[now_chosen_landmark].g = 255;
		fitted_landmarks[now_chosen_landmark].b = 255;
	}

	clicked_points.reset(new pcl::PointCloud<pcl::PointXYZRGBA>);
	for (int i = 0; i < fitted_landmarks.size(); i++)
	{
		clicked_points->push_back(fitted_landmarks[i]);
	}
	// these are the same with both POINT_CLOUD and MESH mode
	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGBA> handler(clicked_points);
	viewer->removePointCloud("clicked_points");
	viewer->addPointCloud(clicked_points, handler, "clicked_points");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, width, "clicked_points");
	qvtkWidget->update();
}

void PCLViewer::onCancelPicking()
{
	has_picked_one_annotation = false;
	now_chosen_landmark = -1;
	// re-render
	paintLandmarks(POINT_SIZE_3D);
}

void PCLViewer::undo()
{
	if (last_changed_landmark_idx != -1)
	{
		fitted_landmarks[last_changed_landmark_idx] = last_changed_landmark;
		last_changed_landmark_idx = -1;
		paintLandmarks(POINT_SIZE_3D);
	}
	else
	{
		// TODO: a pop-up window
		;
	}
}

void PCLViewer::setHint(int index)
{
	if (index >= 0)
	{
		index = index + 1;
		// the upper text
		std::string now_string = "Index " + std::to_string(index);
		instruction_container->setText(QString::fromStdString(now_string));
		// the lower picture
		std::string file_name;
		//
		if (index >= 1 && index <= 3)
			file_name = "1-3";
		else if (index <= 9)
			file_name = "4-9";
		else if (index <= 14)
			file_name = "9-14";
		else if (index <= 17)
			file_name = "15-17";
		else if (index <= 22)
			file_name = "18-22";
		else if (index <= 27)
			file_name = "23-27";
		else if (index <= 29)
			file_name = "28-29";
		else if (index <= 36)
			file_name = "30-36";
		else if (index <= 42)
			file_name = "37-42";
		else if (index <= 48)
			file_name = "43-48";
		else if (index <= 68)
			file_name = "49-68";

		//
		std::string total_file_name = "C:\\Users\\shangxuanu\\Desktop\\MyAnnotationTool\\AnnotationTool_Oculus\\x64\\Release\\" + file_name + ".png";
		QPixmap image(QString::fromStdString(total_file_name));
		QPixmap scaled_image = image.scaled(400, 200);
		image_container->setPixmap(scaled_image);
	}
	else // index == -1
	{
		instruction_container->setText("No points selected now.");
		image_container->clear();
	}
}