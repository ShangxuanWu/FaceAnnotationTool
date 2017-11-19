#include "EpipolarWindow.h"

EpipolarWindow::EpipolarWindow(std::vector<QString> image_paths, int point_id_, std::vector<mycamera> camera_cluster, QWidget *parent)
	: QDialog(parent)
{
	this->setWindowTitle(QString::fromStdString("Shift + Click point " + std::to_string(point_id_ + 1) + " in left image, then click in right image under white line constraint."));
	this->setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint & ~Qt::WindowCloseButtonHint);
	//this->setWindowFlags(windowFlags() | Qt::FramelessWindowHint | Qt::WindowTitleHint);

	point_id = point_id_;

	//camera_1 = camera1;
	//camera_2 = camera2;
	
	getLeftAndRightImagePaths(image_paths);
	getLeftAndRightCamera(camera_cluster);
	
	// make the layout
	left_panel = new EpipolarImageView(this, EPIPOLAR_LEFT_WINDOW);
	left_panel->setImage(left_image_path.toStdString());
	right_panel = new EpipolarImageView(this, EPIPOLAR_RIGHT_WINDOW);
	right_panel->setImage(right_image_path.toStdString());

	left_panel->setFixedHeight(EPIPOLAR_HEIGHT);
	right_panel->setFixedHeight(EPIPOLAR_HEIGHT);
	left_panel->setFixedWidth(EPIPOLAR_WIDTH);
	right_panel->setFixedWidth(EPIPOLAR_WIDTH);

	// get line button
	get_line_button = new QPushButton(">>", this);
	connect(get_line_button, &QPushButton::clicked, this, &EpipolarWindow::drawLine);
	get_line_button->hide();
	
	// okay button
	//ok_button = new QPushButton("OK",this);
	//connect(ok_button, &QPushButton::clicked, this, &EpipolarWindow::ok);
	//ok_button->hide();

	// try button box
	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Vertical);
	button_box->button(QDialogButtonBox::Ok)->setEnabled(false);
	
	// set layout
	horizontal_layout = new QHBoxLayout(this);
	horizontal_layout->addWidget(left_panel);
	//horizontal_layout->addWidget(get_line_button);
	horizontal_layout->addWidget(right_panel);
	//horizontal_layout->addWidget(ok_button);
	horizontal_layout->addWidget(button_box);
	horizontal_layout->setContentsMargins(10, 10, 10, 10);

	this->setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed));

	left_clicked_point = cv::Point2d(-1, -1);
	right_clicked_point = cv::Point2d(-1, -1);
	
	// should load image

	// should give epipolar_constraint: two cameras, one point
	//epipolar_constraint = new EpipolarConstraint(camera1, camera2);
	connect(left_panel, SIGNAL(sendLeftPanelClickedPoint(double, double)), this, SLOT(getLeftClickedPoints(double, double)));
	connect(right_panel, SIGNAL(sendRightPanelClickedPoint(double, double)), this, SLOT(getRightClickedPoints(double, double)));

	// these two are ok and cancel button
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	//connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));
}

// draw line in right panel
void EpipolarWindow::drawLine()
{
	if (left_clicked_point.x >= 0 && left_clicked_point.y >= 0)
	{
		cv::Mat epipolar_line = getEpipolarLine(camera_1, camera_2, left_clicked_point);
		
		right_panel->enableClicking();
		// draw a line on the right panel
		right_panel->setEpipolarLine(epipolar_line);

		right_panel->clearClickedPoint();
	}

}

void EpipolarWindow::accept()
{
	std::vector<cv::Point3d> pts_dst;

	std::vector<cv::Point2d> pts_src1;
	pts_src1.push_back(left_clicked_point);

	std::vector<cv::Point2d> pts_src2;
	pts_src2.push_back(right_clicked_point);

	// do something to send the 3D points
	triangulation_from_two_views(pts_src1, pts_src2, camera_1, camera_2, pts_dst);
	
	// emit the final result
	emit epipolarLineOkSignal(point_id, pts_dst[0]);

	// close the window
	QDialog::accept();
}

void EpipolarWindow::getLeftAndRightImagePaths(std::vector<QString> image_paths)
{
	// point_id is 0-based
	for (int i = 0; i < image_paths.size(); i++)
	{
		if (image_paths[i].contains(QString::fromStdString(point_to_epipolar_image_choice[point_id+1].first), Qt::CaseInsensitive))
		{
			left_image_path = image_paths[i];
		}
		else if (image_paths[i].contains(QString::fromStdString(point_to_epipolar_image_choice[point_id+1].second), Qt::CaseInsensitive))
		{
			right_image_path = image_paths[i];
		}
	}
	// Add Exception control when didn't find the folder

}

void EpipolarWindow::getLeftAndRightCamera(std::vector<mycamera> camera_cluster)
{
	for (int i = 0; i < camera_cluster.size(); i++)
	{
		if (camera_cluster[i].name == point_to_epipolar_image_choice[point_id+1].first)
		{
			camera_1 = camera_cluster[i];
		}
		else if (camera_cluster[i].name == point_to_epipolar_image_choice[point_id+1].second)
		{
			camera_2 = camera_cluster[i];
		}
	}

	// Add Exception control when didn't find the camera
}

/*
void EpipolarWindow::ok()
{
	// calculate 3d point

	// emit signal
}
*/

void EpipolarWindow::getLeftClickedPoints(double x, double y)
{
	left_clicked_point = cv::Point2d(x, y);
	drawLine();
}

void EpipolarWindow::getRightClickedPoints(double x, double y)
{
	right_clicked_point = cv::Point2d(x, y);
	button_box->button(QDialogButtonBox::Ok)->setEnabled(true);
}
