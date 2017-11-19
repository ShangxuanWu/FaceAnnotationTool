#define ENABLE_ASSERTS

#include "MainWindow.h"

#include <QtCore>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QKeyEvent>
#include <QLabel>
#include <QFileDialog>
#include <QInputDialog>
#include <QtOpenGL>
#include <QIcon>
#include <QMatrix>
#include <QFileInfo>
#include <QStringRef>
#include <QMessageBox>
#include <iostream>
#include <fstream>
#include <tuple>
#include <omp.h>
#include <sstream>

#include "../camera_propagation/camera_geometry.h"
#include "../camera_propagation/IdDataParser.h"
#include "../camera_propagation/io_point.h"

#include "../camera_propagation/pts_2d_tool.h"
#include "../camera_propagation/pts_3d_conf.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent)
{
	check_mulitview_correction = false;

	omp_set_dynamic(0);
	omp_set_num_threads(12);

	rootFolder = "";
	isJustPropagated = false;
	isAnnotationShown = true;

	connect(this, SIGNAL(quitAppSelected()), qApp, SLOT(quit()));

	// data viewer
	middleContainer_ = new QSplitter(Qt::Horizontal);
	//this->setCentralWidget(dataViewWidget_);
	QHBoxLayout* hlayout = new QHBoxLayout(middleContainer_);
	hlayout->setContentsMargins(0, 0, 0, 0);

	////// image widget
	QWidget* image_widget = new QWidget(middleContainer_);
	QVBoxLayout* image_vlayout = new QVBoxLayout(image_widget);
	image_vlayout->setContentsMargins(0, 0, 0, 0);
	hlayout->addWidget(image_widget);

	// image_IO_widget is the toolbar
	QWidget* image_IO_widget = new QWidget(image_widget);
	QHBoxLayout* image_IO_hlayout = new QHBoxLayout(image_IO_widget);
	image_IO_hlayout->setContentsMargins(0, 0, 0, 0);
	image_vlayout->addWidget(image_IO_widget);
	image_IO_widget->setMaximumHeight(23);

	loadMeshButton_ = new QPushButton("Load 3D Mesh", image_IO_widget);
	image_IO_hlayout->addWidget(loadMeshButton_);
	connect(loadMeshButton_, &QPushButton::clicked, this, &MainWindow::onLoadMesh);
	loadMeshButton_->hide();

	PrecomputeButton_ = new QPushButton("Preprocess", image_IO_widget);
	image_IO_hlayout->addWidget(PrecomputeButton_);
	connect(PrecomputeButton_, &QPushButton::clicked, this, &MainWindow::onPrecompute);

	loadPastHairDataButton_ = new QPushButton(" Load Annotated Folder ", image_IO_widget);
	image_IO_hlayout->addWidget(loadPastHairDataButton_);
	connect(loadPastHairDataButton_, &QPushButton::clicked, this, &MainWindow::onLoadAnnotatedFolder);

	propagateButton_ = new QPushButton(" Propagate ", image_IO_widget);
	image_IO_hlayout->addWidget(propagateButton_);
	connect(propagateButton_, &QPushButton::clicked, this, &MainWindow::onPropagateImage);
	propagateButton_->hide();

	saveButton_ = new QPushButton(" Save ", image_IO_widget);
	image_IO_hlayout->addWidget(saveButton_);
	connect(saveButton_, &QPushButton::clicked, this, &MainWindow::onSaveAnnotationAndVisibility);
	saveButton_->setDisabled(true);
	
	imageFPathTextEdit_ = new QTextEdit(" (Current Image Folder Path Here) ", image_IO_widget);
	imageFPathTextEdit_->setReadOnly(true);
	image_IO_hlayout->addWidget(imageFPathTextEdit_);


	find_visibility_button = new QPushButton("Visibility", image_IO_widget);
	image_IO_hlayout->addWidget(find_visibility_button);
	connect(find_visibility_button, &QPushButton::clicked, this, &MainWindow::findVisibility);
	find_visibility_button->hide();
	
	// Point Annotation Mode Button -> useless
	pointAnnotateModeButton_ = new QPushButton("Point Mode", image_IO_widget);
	image_IO_hlayout->addWidget(pointAnnotateModeButton_);
	connect(pointAnnotateModeButton_, &QPushButton::clicked, this, &MainWindow::onSetPointMode);
	pointAnnotateModeButton_->hide();

	// Line Annotation Mode Button -> useless
	lineAnnotateModeButton_ = new QPushButton("Line Mode", image_IO_widget);
	image_IO_hlayout->addWidget(lineAnnotateModeButton_);
	connect(lineAnnotateModeButton_, &QPushButton::clicked, this, &MainWindow::onSetLineMode);
	lineAnnotateModeButton_->hide();

	// Annotation Mode Button -> useless
	annotationModeButton_ = new QPushButton("Annotation", image_IO_widget);
	image_IO_hlayout->addWidget(annotationModeButton_);
	connect(annotationModeButton_, &QPushButton::clicked, this, &MainWindow::onSetAnnotationMode);
	annotationModeButton_->hide();

	// Selection Mode Button -> useless
	selectionModeButton_ = new QPushButton("Selection", image_IO_widget);
	image_IO_hlayout->addWidget(selectionModeButton_);
	connect(selectionModeButton_, &QPushButton::clicked, this, &MainWindow::onSetSelectionMode);
	selectionModeButton_->hide();
	

	imageView_ = new ImageView(image_widget);
	image_vlayout->addWidget(imageView_);

	statusBar_ = new QStatusBar(image_widget);
	statusBar_->setFixedHeight(20);
	image_vlayout->addWidget(statusBar_);
	statusBar_->showMessage("Please precompute a folder or load an annotated folder.");

	right_container = new RightWidget();

	// for mesh view
	vtk_viewer = new VTKViewer();
	vtk_viewer->setFixedWidth(1000);

	// for the total layout
	totalLayout_ = new QSplitter(Qt::Horizontal, this);
	totalLayout_->addWidget(vtk_viewer);
	totalLayout_->addWidget(middleContainer_);
	//totalLayout_->addWidget(rightContainer_);
	totalLayout_->addWidget(right_container);
	totalLayout_->setContentsMargins(10, 10, 10, 10);
	totalLayout_->setCollapsible(0, false);
	totalLayout_->setCollapsible(1, false);
	totalLayout_->setCollapsible(2, false);
	//totalLayout_->setAttribute(Qt::WA_TransparentForMouseEvents);

	this->setCentralWidget(totalLayout_);

	this->setWindowTitle(tr("Multi-view Smart Annotation Tool @ Oculus Research Pittsburgh"));
	this->setFixedSize(2560, 1600);

	onSetSelectionMode();
	onSetPointMode();

	// this has to be modified
	connect(right_container->thumbnail_container, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(onRightPanelItemClicked(QListWidgetItem*)));

	setWindowIcon(QIcon("Resources/Oculus.ico"));

	// connecting set hint signal / slot
	connect(vtk_viewer, SIGNAL(setHint(int)), right_container, SLOT(setHint(int)));
	connect(vtk_viewer, SIGNAL(cancelHint()), right_container, SLOT(cancelHint()));
	connect(vtk_viewer, SIGNAL(setHint3D(int)), right_container, SLOT(setHint3D(int)));
	connect(vtk_viewer, SIGNAL(cancelHint3D()), right_container, SLOT(cancelHint3D()));
	connect(imageView_, SIGNAL(setHint(int)), right_container, SLOT(setHint(int)));
	connect(imageView_, SIGNAL(cancelHint()), right_container, SLOT(cancelHint()));

	// connecting status bar and their messages
	connect(this, SIGNAL(setStatusBarMessage(std::string)), this, SLOT(getStatusBarMessage(std::string)));

	// connecting the left vtkviewer and right vtkviewer
	connect(vtk_viewer, SIGNAL(meshPointerSignal(vtkSmartPointer<vtkPolyData>)), right_container->mesh_visibility_helper, SLOT(setMesh(vtkSmartPointer<vtkPolyData>)));

	// connect the left vtkviewer and the middle panel
	connect(vtk_viewer, SIGNAL(correctThisPointButtonPressed(int)), this, SLOT(onPopUpEpipolarWindow(int)));
	connect(vtk_viewer, SIGNAL(to2DButtonPressed()), this, SLOT(on3DAnnotationTo2D()));
}

void MainWindow::onSetHairCount()
{
	if (imageView_->hasImage())
	{
		cout << "Warning: no implementation here..." << endl;
	}
	else
	{
		cerr << "Error: set hair count AFTER you load image." << endl;
		cerr << "FILE: " << __FILE__ << endl;
		cerr << "LINE: " << __LINE__ << endl;
	}
}

void MainWindow::keyPressEvent(QKeyEvent* e)
{
	if (e->key() == Qt::Key_C)
	{
		onToggleAnnotationsShown();
	}
	else if (e->key() == Qt::Key_S)
	{
		this->onSaveImagePoints();

	}
	else if (e->key() == Qt::Key_Delete
		|| e->key() == Qt::Key_Backspace)
	{
		this->onDeleteLatestAnnotation();
	}
	else if (e->key() == Qt::Key_R)
	{
		imageView_->rotate90Clockwise();
	}
	else if (e->key() == Qt::Key_L)
	{
		imageView_->rotate90CounterClockwise();
	}
	else if (e->key() == Qt::Key_V)
	{
		imageView_->setVisualizePastPoints(!(imageView_->visualizePastPoints()));
	}

}

void MainWindow::onDeleteLatestAnnotation()
{
	imageView_->deleteLatestAnnotation();
}

void MainWindow::onSaveImagePoints()
{
	imageView_->saveLatestAnnotationData();
	statusBar_->showMessage("Recent annotation saved.");
}

void MainWindow::onLoadMesh()
{
	//qDebug() << "Now loading mesh" << endl;
	statusBar_->showMessage("Now loading mesh.");

	QString mesh_fpath = QFileDialog::getOpenFileName(this, "Open PLY mesh file", "../PerformanceMeshRenderingTool", tr("PLY Mesh (*.ply)"));

	if (mesh_fpath.isEmpty())
	{
		qDebug() << "Warning: mesh is not selected." << endl;
		return;
	}
	//meshContainer_->setMesh(mesh_fpath.toStdString());

	statusBar_->showMessage("Mesh loaded.");
}

void MainWindow::onPrecompute() {
	QString desktop_path = QDir::homePath() + "/Desktop";
	QString image_folder_path = QFileDialog::getExistingDirectory(this, tr("Choose Image Folder"), desktop_path, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

	if (image_folder_path.size() == 0) {
		statusBar_->showMessage("Please choose a folder.");
		return;
	}

	rootFolder = image_folder_path.toStdString();

	statusBar_->showMessage("Preprocessing image (might take 10 minutes). Please wait ...");

	QDirIterator it(image_folder_path, QDirIterator::NoIteratorFlags);

	std::vector<std::string> tmpImgPaths;

	while (it.hasNext()) {
		QString imageFolder = it.next();
		QFileInfo checkDir(imageFolder);
		if (imageFolder.endsWith(".") || !checkDir.isDir())
		{
			continue;
		}
		QString imagePath = imageFolder + "/image0000.png";
		tmpImgPaths.push_back(imagePath.toStdString());
	}

	for (int i = 0; i < tmpImgPaths.size(); i++) {

		// read the file number
		std::string imageFolderName = getFolderFromFilename(tmpImgPaths[i]);
		std::string imageID = getIDFromFilename(tmpImgPaths[i]);

		// read image
		cv::Mat img = cv::imread(tmpImgPaths[i]);

		int h = img.size().height;
		int w = img.size().width;

		// judge rotation
		if (rotation_indicator[imageID] == LEFT)
		{
			rot90(img, LEFT);
		}
		else if (rotation_indicator[imageID] == RIGHT)
		{
			rot90(img, RIGHT);
		}
		else
		{
			;
		}

		std::vector<cv::Point> detected_landmarks = detect_one_image_using_dlib_input_Mat(img, DLIB_MODEL_PATH);

		std::map<int, pts_2d_tool> hash_one_img_detected_landmarks;

		// adjust the detected points, and save all detected landmarks into our own format
		for (int j = 0; j < detected_landmarks.size(); j++)
		{
			int new_x, new_y;

			if (rotation_indicator[imageID] == LEFT)
			{
				// for all points
				new_x = w - detected_landmarks[j].y;
				new_y = detected_landmarks[j].x;
			}
			else if (rotation_indicator[imageID] == RIGHT)
			{
				// for all points
				new_x = detected_landmarks[j].y;
				new_y = h - detected_landmarks[j].x;
			}
			else
			{
				qDebug() << "MIDDLE";
				new_x = detected_landmarks[j].x;
				new_y = detected_landmarks[j].y;
			}
			hash_one_img_detected_landmarks[j] = pts_2d_tool(new_x, new_y, INITIAL_CONF, false);
		}

		// if dlib haven't detected any points, give 68 (0,0) points
		//if (detected_landmarks.size() == 0 || valid_dlib_view.find(imageFolderName) == valid_dlib_view.end()) {
		if (detected_landmarks.size() == 0)
		{
			for (int k = 0; k < TOTAL_POINTS_IN_SINGLE_VIEW; k++) {
				hash_one_img_detected_landmarks[k] = pts_2d_tool(0, 0, 0, false);
			}
		}

		hash_all_imgs_detected_landmarks[imageFolderName] = hash_one_img_detected_landmarks;
		std::string message = "Preprocessing: finished precomputing image " + imageFolderName + ". Got " + std::to_string(detected_landmarks.size()) + " keypoints. Please wait ...";
		statusBar_->showMessage(QString::fromStdString(message));
	}
	// propagate image once 
	convertAndPropagate(PREPROCESS_RANSAC_ITER);

	// set all confidence = 0 and anchor = 0
	for (auto it = hash_all_imgs_detected_landmarks.begin(); it != hash_all_imgs_detected_landmarks.end(); it++)
	{
		for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
		{
			(*it).second[i].anchor = false;
			(*it).second[i].conf = CONF_AFTER_FIRST_PROPAGATION_IN_PREPROCESS;
		}
	}

	// saving all the points
	onSaveAnnotation();

	// erase rootFolder variable and hash variable
	rootFolder = "";
	hash_all_imgs_detected_landmarks.clear();

	statusBar_->showMessage("Finished precomputing and saving images, you can now load the folder.");
}

// this function will perform a propagation after loading the annotated folder
void MainWindow::onLoadAnnotatedFolder()
{

	//QString desktop_path = QDir::homePath() + "/Desktop";
	QString now_path = QDir::currentPath() + "/../";
	QString image_folder_path = QFileDialog::getExistingDirectory(this, tr("Choose Image Folder"), now_path, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

	// if "cancel" is selected
	if (image_folder_path.size() == 0) {
		return;
	}
	try
	{
		hash_all_imgs_detected_landmarks.clear();
	
		right_container->clear();

		saveButton_->setDisabled(false);
		rootFolder = image_folder_path.toStdString();

		statusBar_->showMessage("Loading image and previous landmarks ... Please wait ...");

		std::string annotation_file_path = rootFolder + "/annotation.yml";
		std::string cali_file_path = rootFolder + "/calibration.txt";

		// if no annotated file, return
		QFileInfo check_anno_file(QString::fromStdString(annotation_file_path));
		QFileInfo check_cali_file(QString::fromStdString(cali_file_path));

		if (!check_anno_file.exists() || !check_cali_file.exists())
		{
			statusBar_->showMessage("Annotation / Calibration file didn't found! Please run \"Preprocess\" first .");
			return;
		}

		hash_all_imgs_detected_landmarks = load_annotated_points(annotation_file_path);

		imageFPathTextEdit_->setText(image_folder_path);

		imagePaths.clear();

		QDirIterator it(image_folder_path, QDirIterator::NoIteratorFlags);
		while (it.hasNext()) {
			QString imageFolder = it.next();
			QFileInfo checkDir(imageFolder);
			if (imageFolder.endsWith(".") || !checkDir.isDir())
			{
				continue;
			}
			qDebug() << imageFolder;
			QString imagePath = imageFolder + "/image0000.png";
			imagePaths.push_back(imagePath);
		}

#pragma omp parallel for schedule(static)
		for (int i = 0; i < imagePaths.size(); i++) {

			// read the file number
			std::string cameraNumber = imagePaths[i].toStdString();
			cameraNumber = getIDFromFilename(cameraNumber);

			// read image
			QImageReader reader(imagePaths[i]);
			QImage image = reader.read();

			// make the thumbnail
			QIcon* icon;
			if (rotation_indicator[cameraNumber] != MIDDLE) {
				QMatrix rm;
				int degree = 90 * rotation_indicator[cameraNumber];
				rm.rotate(degree);
				icon = new QIcon(QPixmap::fromImage(image.scaled(140, 105).transformed(rm)));
			}
			else {
				icon = new QIcon(QPixmap::fromImage(image.scaled(105, 140)));
			}

			// for ordering the images
			std::string order_str = find_camera_order_return_padded_string(cameraNumber);
			QListWidgetItem * now_item = new QListWidgetItem(*icon, QString::fromStdString(order_str));
#pragma omp critical
			{
				//rightContainer_->addItem(now_item);
				right_container->addThumbnailItem(now_item);
			}
			qDebug() << "Finished loading image " << i + 1 << endl;
		}
		// reset imageview
		imageView_->resetAll();

		// load mesh

		std::string mesh_fpath_std = rootFolder + "/mesh.ply";

		QString mesh_fpath = QString::fromStdString(mesh_fpath_std);

		if (mesh_fpath.isEmpty())
		{
			statusBar_->showMessage("Warning: no mesh found in this folder.");
		}
		else
		{
			//meshContainer_->setMesh(mesh_fpath_std);
			vtk_viewer->setMesh(mesh_fpath_std);
		}

		// finish
		statusBar_->showMessage("Finished loading images.");
	}
	catch (const std::exception&)
	{
		statusBar_->showMessage("Loading folder failed ! Please choose the correct folder.");
	}
	// propagate once
	onPropagateImage();
}

cv::Mat MainWindow::QImage2Mat(const QImage& src) {
	cv::Mat mat = cv::Mat(src.height(), src.width(), CV_8UC4, (uchar*)src.bits(), src.bytesPerLine());
	cv::Mat result = cv::Mat(mat.rows, mat.cols, CV_8UC3);
	int from_to[] = { 0,0,  1,1,  2,2 };
	cv::mixChannels(&mat, 1, &result, 1, from_to, 3);
	return result;
	//return mat;
}

void MainWindow::onAutoSegImage()
{
	cout << "now loading image" << endl;

	QString image_fpath = QFileDialog::getOpenFileName(this, "Open image file", ".", tr("Image (*.png)"));

	if (image_fpath.isEmpty())
	{
		cout << "Warning: image is not selected." << endl;
		return;
	}
	else
	{
		cout << "selected image name: " << image_fpath.toStdString() << endl;
		imageFilePath_ = image_fpath.toStdString();

		imageView_->setImage(image_fpath.toStdString());
		imageFPathTextEdit_->setText(imageFilePath_.c_str());

		std::string mask_path = "none";
		std::string output_yml_name = "./tmp.yml";
		std::string output_ori_png_name = "./tmp_ori.png";

		this->initializeSelectedPoints();

	}
	cout << "Image load done!" << endl;
}


// convert to the right format, do propagation, and get a set of 3D points
std::vector < pts_3d_conf> MainWindow::convertAndPropagate(int iteration)
{
	std::map<std::string, std::vector<pts_2d_conf>> pts_src;
	std::map<std::string, std::vector<pts_2d_conf>> pts_output;

	// read camera cluster matrix
	std::string calib_txt_path = rootFolder + "/calibration.txt";
	LoadIdCalibration(calib_txt_path, camera_cluster, false);
	statusBar_->showMessage("Reading calibration files and calculating landmarks ... Please wait ...");

	// convert hash... to pts_2d_conf
	for (auto it = hash_all_imgs_detected_landmarks.begin(); it != hash_all_imgs_detected_landmarks.end(); it++)
	{
		std::vector<pts_2d_conf> now_vec;
		for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
		{
			now_vec.push_back(it->second[i].convert_to_pts_2d_conf());
		}
		pts_src[it->first] = now_vec;
	}

	save_points_with_conf_multiview(rootFolder, 0, pts_src);

	// calculate
	multiview_optimization_multiple_points(pts_src, camera_cluster, pts_output, true, iteration);

	// convert back to hash...
	// loop for camera
	for (auto it2 = pts_output.begin(); it2 != pts_output.end(); it2++)
	{
		// loop for keypoint
		for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++) {
			// if is not user-moved points
			if (hash_all_imgs_detected_landmarks[it2->first][i].anchor == false)
			{
				hash_all_imgs_detected_landmarks[it2->first][i] = pts_2d_tool(it2->second[i]);
			}
		}
	}

	// show 3D landmarks
	std::vector < pts_3d_conf> keypoints_3d;

	// chose the first two cameras in camera_cluster
	triangulation_from_two_views(pts_output[camera_cluster[0].name], pts_output[camera_cluster[1].name], camera_cluster[0], camera_cluster[1], keypoints_3d);

	return keypoints_3d;
}

// this is the main propagation function
void MainWindow::onPropagateImage() {
	if (hash_all_imgs_detected_landmarks.size() != TOTAL_VIEW_COUNT)
	{
		statusBar_->showMessage("No annotated points loaded now.");
		return;
	}
	// get all the points from now imageView_
	// getPointsFromMiddlePanel();

	// check if user only corrected one view in a particular index
	if (check_mulitview_correction)
	{
		std::vector<int> lack_points;
		for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
		{
			int total_point_of_this_index = 0;
			for (auto it = hash_all_imgs_detected_landmarks.begin(); it != hash_all_imgs_detected_landmarks.end(); it++)
			{
				if (it->second[i].anchor == true)
				{
					total_point_of_this_index++;
				}
			}
			if (total_point_of_this_index == 1)
			{
				lack_points.push_back(i);
			}
		}
		if (lack_points.size() > 0)
		{
			QMessageBox msgBox;
			std::string point_string = "";
			for (int j = 0; j < lack_points.size(); j++)
			{
				point_string += std::to_string(lack_points[j] + 1);
				if (j != lack_points.size() - 1)
				{
					point_string += ", ";
				}
			}
			std::string message = "Please correct point " + point_string + " in at least two views !";
			msgBox.setText(QString::fromStdString(message));
			msgBox.exec();
			return;
		}
	}
	
	keypoints_3d = convertAndPropagate(PROPAGATE_RANSAC_ITER);

	// update shown points on current view
	isJustPropagated = true;
	//onRightPanelItemClicked(rightContainer_->currentItem());
	onRightPanelItemClicked(right_container->getCurrentItem());
	isJustPropagated = false;

	//meshContainer_->set3DLandmarks(keypoints_3d);
	keypoints_3d_mapped_to_mesh = vtk_viewer->setAndShowLandmarks(keypoints_3d);
	statusBar_->showMessage("Points propagated, 3D landmarks shown.");
}

void MainWindow::initializeSelectedPoints()
{
	imageView_->clearSelectedPoints();
}

void MainWindow::onSetSelectionMode()
{
	selectionModeButton_->setDisabled(true);
	annotationModeButton_->setDisabled(false);
	imageView_->setSelectionMode(true);
}

void MainWindow::onSetAnnotationMode()
{
	selectionModeButton_->setDisabled(false);
	annotationModeButton_->setDisabled(true);
	imageView_->setSelectionMode(false);
}

// this is also used as a refresh function for the main panel
void MainWindow::onRightPanelItemClicked(QListWidgetItem* item) {
	if (item == NULL) {
		return;
	}
	
	QString name_string = item->text();
	QStringRef sub_string(&name_string, 0, 2);
	int number = sub_string.toInt();
	std::string folder_name = getFolderFromID(camera_order[number-1]);
	std::string full_path = getFullPathFromFolder(rootFolder, folder_name);
	imageView_->setImage(full_path);
	imageView_->setDetectedKeypoint(hash_all_imgs_detected_landmarks[folder_name]);
	imageView_->folder_name = folder_name;

	// get camera from folder_name
	mycamera now_camera = searchCameraByFolderName(folder_name);
	
	right_container->mesh_visibility_helper->setCamera(now_camera);

	// calculate visibility and set it
	vector<bool> this_view_visibility = right_container->mesh_visibility_helper->getThisViewVisibility(keypoints_3d_mapped_to_mesh);
	imageView_->setVisibility(this_view_visibility);
}

mycamera MainWindow::searchCameraByFolderName(std::string folder_name)
{
	for (int i = 0; i < camera_cluster.size(); i++)
	{
		if (camera_cluster[i].name == folder_name)
		{
			return camera_cluster[i];
		}
	}

}


void MainWindow::onSetPointMode() {
	//statusBar_->showMessage("Set to Point Mode.");
	imageView_->set_draw_mode_(POINT_MODE);
	pointAnnotateModeButton_->setDisabled(true);
	lineAnnotateModeButton_->setDisabled(false);
}

void MainWindow::onSetLineMode() {
	//imageView_->set_draw_mode_(LINE_MODE);
	statusBar_->showMessage("Haven't implemented LINE MODE now.");
	pointAnnotateModeButton_->setDisabled(false);
	lineAnnotateModeButton_->setDisabled(true);
}

void MainWindow::onSaveAnnotation() 
{
	if (rootFolder.size() == 0) {
		statusBar_->showMessage("Please load a folder first.");
	}
	//getPointsFromMiddlePanel();
	std::string filename;
	filename = rootFolder + "/annotation.yml";
	save_annotated_points(filename, hash_all_imgs_detected_landmarks);
	statusBar_->showMessage("Annotation saved.");
}

void MainWindow::onSaveAnnotationAndVisibility()
{
	onSaveAnnotation();
	statusBar_->showMessage("Calculating final annotations. Please wait ...");

	// save visibility
	if (right_container->thumbnail_container->count() != TOTAL_VIEW_COUNT)
	{
		statusBar_->showMessage("There is some problem with saving. Please restart the program.");
		return;
	}
	std::map<std::string, std::vector<bool>> all_views_all_visibilities;

	// for all cameras
	for (int i = 0; i < right_container->thumbnail_container->count(); i++)
	{
		QListWidgetItem* item = right_container->thumbnail_container->takeItem(i);
		onRightPanelItemClicked(item);
		std::vector<bool> this_visibility = imageView_->getVisibility();
		std::string this_camera_name = imageView_->getCameraName();
		all_views_all_visibilities[this_camera_name] = this_visibility;
	}

	std::string filename = rootFolder + "/visibility.yml";
	save_points_visibility(filename, all_views_all_visibilities);
	statusBar_->showMessage("Annotation saved.");
}

// for getting the points from the image shown in the middle panel
void MainWindow::getPointsFromMiddlePanel()
{
	if (!imageView_->hasImage())
	{
		return;
	}
	// get points now
	std::string folder_name = imageView_->folder_name;
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		hash_all_imgs_detected_landmarks[folder_name][i].x = imageView_->this_landmarks[i].x;
		hash_all_imgs_detected_landmarks[folder_name][i].y = imageView_->this_landmarks[i].y;
	}
	for (int i = 0; i < imageView_->modified_points.size(); i++)
	{
		// confidence = 1 after user moving the point
		hash_all_imgs_detected_landmarks[folder_name][imageView_->modified_points[i]].conf = 1.0;
		// anchor = 1 after user moving the point
		hash_all_imgs_detected_landmarks[folder_name][imageView_->modified_points[i]].anchor = true;
		qDebug() << imageView_->modified_points[i];
	}
	imageView_->modified_points.clear();
	qDebug() << "POINTS UPDATED !!!";
}

int MainWindow::find_camera_order(std::string ID)
{
	for (int i = 0; i < TOTAL_VIEW_COUNT; i++)
	{
		if (camera_order[i] == ID)
		{
			return i;
		}
	}
	return TOTAL_VIEW_COUNT;
}

std::string MainWindow::find_camera_order_return_padded_string(std::string ID)
{
	// getting the right order and return something like "09 (Cam01)"
	for (int i = 0; i < TOTAL_VIEW_COUNT; i++)
	{
		if (camera_order[i] == ID)
		{
			std::ostringstream ss;
			ss << std::setw(2) << std::setfill('0') << (i+1);
			return ss.str() + " ( cam" + ID + " )";
		}
	}
}

void MainWindow::onToggleAnnotationsShown()
{
	if (isAnnotationShown)
	{
		isAnnotationShown = false;
		imageView_->toggleAnnotationsShown();
		statusBar_->showMessage("Annotation: hide.");
	}
	else
	{
		isAnnotationShown = true;
		imageView_->toggleAnnotationsShown();
		statusBar_->showMessage("Annotation: show.");
	}
}

void MainWindow::on3DAnnotationTo2D()
{
	if (right_container->isEmpty())
	{
		return;
	}

	// get 3D from the panel
	std::vector<cv::Point3d> new_3D_coordinates = vtk_viewer->updateMainWindow();

	// update keypoints_3d
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		keypoints_3d[i].x = new_3D_coordinates[i].x;
		keypoints_3d[i].y = new_3D_coordinates[i].y;
		keypoints_3d[i].z = new_3D_coordinates[i].z;
	}

	// propagate
	std::map<std::string, std::vector<cv::Point2d>> pts_dst;

	// void multi_view_projection(std::vector<cv::Point3d>& pts_src, std::vector<mycamera>& camera_cluster, std::map<std::string, std::vector<cv::Point2d>>& pts_dst, const bool consider_dist = default_consider_dist);
	multi_view_projection(new_3D_coordinates, camera_cluster, pts_dst);

	// update the 2D points to the middle window
	// pts_dst -> hash_all_imgs_detected_landmarks
	for (auto it = pts_dst.begin(); it != pts_dst.end(); it++)
	{
		for (int j = 0; j < TOTAL_POINTS_IN_SINGLE_VIEW; j++)
		{
			hash_all_imgs_detected_landmarks[it->first][j].x = pts_dst[it->first][j].x;
			hash_all_imgs_detected_landmarks[it->first][j].y = pts_dst[it->first][j].y;
			hash_all_imgs_detected_landmarks[it->first][j].conf = 1;
			hash_all_imgs_detected_landmarks[it->first][j].anchor = 0;
		}
	}

	// what if some image is selected in the middle panel? Update it.
	isJustPropagated = true;
	onRightPanelItemClicked(right_container->getCurrentItem());
	isJustPropagated = false;

	// update status bar
	statusBar_->showMessage("Updated 2D annotations from 3D panel.");
}

void MainWindow::setHint(int point_id)
{
	right_container->setHint(point_id);
}

void MainWindow::cancelHint()
{
	right_container->cancelHint();
}

void MainWindow::onPopUpEpipolarWindow(int point_id)
{
	if (imagePaths.empty())
	{
		setStatusBarMessage("No image folder loaded. Please load a folder first.");
		return;
	}

	epipolar_window = new EpipolarWindow(imagePaths, point_id, camera_cluster, this);
	epipolar_window->show();
	// connecting epipolar line constraint
	connect(this->epipolar_window, SIGNAL(epipolarLineOkSignal(int, cv::Point3d)), this, SLOT(getEpipolarLineOkSignal(int, cv::Point3d)));
	connect(this->epipolar_window->button_box, SIGNAL(rejected()), this->vtk_viewer, SLOT(onCancel()));
	connect(this->epipolar_window->button_box, SIGNAL(rejected()), this->epipolar_window, SLOT(reject()));
}

void MainWindow::getStatusBarMessage(std::string content)
{
	QString content_ = QString::fromStdString(content);
	statusBar_->showMessage(content_);
	return;
}

void MainWindow::getEpipolarLineOkSignal(int point_id, cv::Point3d point_3d)
{
	// get all 2D points from the 3D point
	std::map<std::string, cv::Point2d> pts_dst;
	multi_view_projection(point_3d, camera_cluster, pts_dst);
	
	// set the 2D points
	// for all cameras
	for (auto it = pts_dst.begin(); it != pts_dst.end(); it++)
	{
		// set single point id
		hash_all_imgs_detected_landmarks[it->first][point_id].x = it->second.x;
		hash_all_imgs_detected_landmarks[it->first][point_id].y = it->second.y;
		hash_all_imgs_detected_landmarks[it->first][point_id].conf = 1;
		hash_all_imgs_detected_landmarks[it->first][point_id].anchor = 1;
	}
	
	isJustPropagated = true;
	onRightPanelItemClicked(right_container->getCurrentItem());
	isJustPropagated = false;


	// cancel selection in left panel
	vtk_viewer->onCancel();

	// update 3D view
	// this one is not updated?
	keypoints_3d[point_id] = pts_3d_conf(point_3d.x, point_3d.y, point_3d.z, 1);
	keypoints_3d_mapped_to_mesh = vtk_viewer->setAndShowLandmarks(keypoints_3d);
	
	// update status bar
	statusBar_->showMessage("Updated 2D annotations from epipolar line constraint.");
}

// testing function for visibility button.
void MainWindow::findVisibility()
{
	vector<bool> this_view_visibility = right_container->mesh_visibility_helper->getThisViewVisibility(keypoints_3d_mapped_to_mesh);
	//imageView_->setVisibility(this_view_visibility);
	/*
	for (int i = 0; i < this_view_visibility.size(); i++)
	{
		qDebug() << i+1 <<  " : " <<  this_view_visibility[i];
	}
	*/
}