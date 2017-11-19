#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define ENABLE_ASSERTS

#include <opencv2/opencv.hpp>
#include "../facial_landmark_detection/main_function.h"

#include <map>
#include <string>

#include <QString>
#include <QMainWindow>
#include <QSplitter>
#include <QTextEdit>
#include <QPushButton>
#include <QSlider>
#include <QTextLine>
#include <QSpinBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QStatusBar>
#include <QListWidget>
#include <QListWidgetItem>
#include <QOpenGLWidget>
#include <pcl/visualization/pcl_visualizer.h>
#include <QSplitter>
#include <QList>
#include <QFileInfo>

#include "ImageView.h"
//#include "pclviewer.h"
#include "RightWidget.h"
#include "VTKViewer.h"
#include "utils.h"
#include "IO.h"
#include "EpipolarWindow.h"
#include "../camera_propagation/mycamera.h"
//#include "hintDialog.h"

// class forward declaration of "../camera_propagation/pts_2d_tool.h"
class pts_2s_tool;

class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		explicit MainWindow(QWidget *parent = 0);
		~MainWindow() {}

		void setHint(int point_id);
		void cancelHint();

	signals:
			void quitAppSelected();
			void setStatusBarMessage(std::string);
			void middlePanelCameraSetSignal(mycamera);

	public slots:
		void getStatusBarMessage(std::string);
	protected slots :
		void keyPressEvent(QKeyEvent* e);

		void onLoadMesh();
		void onPrecompute();
		void onLoadAnnotatedFolder();
		void onAutoSegImage();

		// Epipolar Button
		void onPopUpEpipolarWindow(int point_id);

		// key function for triangulating all the points
		void onPropagateImage();

		void onRightPanelItemClicked(QListWidgetItem* item);

		void onSetHairCount();

		void onSaveImagePoints();
		void initializeSelectedPoints();

		void onDeleteLatestAnnotation();

		void onSetAnnotationMode();
		void onSetSelectionMode();

		void onSetPointMode();
		void onSetLineMode();

		void onSaveAnnotation();
		void onSaveAnnotationAndVisibility();

		void getPointsFromMiddlePanel();

		void onToggleAnnotationsShown();

		// getting annotations from 3D panel
		void on3DAnnotationTo2D();

		void getEpipolarLineOkSignal(int, cv::Point3d);

protected:
	//void closeEvent(QCloseEvent * event);

	QSplitter* leftSplitter_;
	QSplitter* middleContainer_;
	QTextEdit* tmpTextEdit_;

	ImageView* imageView_;

	QPushButton* loadMeshButton_;
	QPushButton* PrecomputeButton_;
	QPushButton* propagateButton_;

	std::string imageFilePath_;

	QTextEdit* imageFPathTextEdit_;

	QPushButton* setHairCountButton_;
	QPushButton* loadPastHairDataButton_;

	QPushButton* annotationModeButton_;
	QPushButton* selectionModeButton_;

	QPushButton* pointAnnotateModeButton_;
	QPushButton* lineAnnotateModeButton_;
	QPushButton* saveButton_;

	// useless
	QPushButton* get_3D_annotation_button;
	QPushButton* epipolar_window_button;

	QPushButton* find_visibility_button;


	// camera
	// camera.name is just "330001"
	std::vector<mycamera> camera_cluster;

	// for Mesh 
	VTKViewer* vtk_viewer;

	// total layout splitter, horizontal
	QSplitter* totalLayout_;

	// right container and layout
	//QListWidget* rightContainer_;
	RightWidget* right_container;
	int find_camera_order(std::string ID);
	std::string MainWindow::find_camera_order_return_padded_string(std::string ID);

	// store the paths of all views from a person
	std::vector<QString> imagePaths;

	// status bar
	QStatusBar* statusBar_;

	//this is the key storage of all the points (2D and 3D)
	std::map<std::string, std::map<int, pts_2d_tool>> hash_all_imgs_detected_landmarks;
	std::vector < pts_3d_conf> keypoints_3d;
	// there will always be some slight difference between real 3D points and mapped 3D points
	std::vector<pts_3d_conf> keypoints_3d_mapped_to_mesh;

	cv::Mat MainWindow::QImage2Mat(const QImage& src);

	// root folder for all the images
	std::string rootFolder;

	// act as a lock
	bool isJustPropagated;

	// annotations shown control
	bool isAnnotationShown;

	std::vector < pts_3d_conf> convertAndPropagate(int iteration);

	EpipolarWindow* epipolar_window;

	bool check_mulitview_correction;

	// for emitting signal for the right panel
	mycamera searchCameraByFolderName(std::string folder_name);


	void findVisibility();
};

#endif // MAINWINDOW_H
