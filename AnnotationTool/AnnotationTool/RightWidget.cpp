#include "RightWidget.h"

RightWidget::RightWidget()
{
	//parent_window = (QMainWindow*)parentWidget();

	thumbnail_container = new QListWidget(this);
	thumbnail_container->setMinimumWidth(150);
	thumbnail_container->setViewMode(QListWidget::IconMode);
	thumbnail_container->setIconSize(QSize(150, 150));
	thumbnail_container->setResizeMode(QListWidget::Adjust);
	thumbnail_container->setMovement(QListView::Static);
	thumbnail_container->setSortingEnabled(true);
	// should we set the height of the QListWidget?

	mesh_visibility_helper = new MeshVisibilityHelper(this);
	mesh_visibility_helper->setFixedWidth(409);
	mesh_visibility_helper->setFixedHeight(409.0 / ORIGINAL_IMAGE_WIDTH * ORIGINAL_IMAGE_HEIGHT);

	hint_container = new QLabel(this);
	//hint_container->setFixedHeight(INSTRUCTION_IMG_HEIGHT);
	hint_container->setMinimumHeight(300);
	hint_container->setAlignment(Qt::AlignCenter);
	// should we set the height of the QHBoxLayout?
	
	//hint_text = new QStatusBar(this);
	hint_text = new QLabel(this);
	hint_text->setFixedHeight(20);

	warning_text = new QLabel(this);
	warning_text->setFixedHeight(20);
	warning_text->setText("Following is a helper window. Ignore it.");
	warning_text->setAlignment(Qt::AlignCenter);

	vertical_layout = new QVBoxLayout(this);
	vertical_layout->addWidget(thumbnail_container);
	vertical_layout->addWidget(warning_text);
	vertical_layout->addWidget(mesh_visibility_helper);
	vertical_layout->addWidget(hint_container);
	vertical_layout->addWidget(hint_text);
	vertical_layout->setContentsMargins(0, 0, 0, 0);

	// show nothing on the label page
	cancelHint();

	this->setFixedWidth(409);
}

QPixmap RightWidget::processImage(std::string fn)
{
	QPixmap image2(QString::fromStdString(fn));
	QRect rect((image2.width() - INSTRUCTION_IMG_WIDTH) / 2, (image2.height() - INSTRUCTION_IMG_HEIGHT) / 2, INSTRUCTION_IMG_WIDTH, INSTRUCTION_IMG_HEIGHT); // x, y, width, height
	QPixmap image2_crop = image2.copy(rect);
	return image2_crop;
}


void RightWidget::setHint(int point_id)
{
	if (point_id >= 0) {
		std::string file_name = "Resources\\" + std::to_string(point_id + 1) + ".png";
		QPixmap image2_crop = processImage(file_name);
		hint_container->setPixmap(image2_crop.scaled(hint_container->width(), hint_container->height(), Qt::KeepAspectRatio));
	}
	else
	{
		cancelHint();
		qDebug() << "point_id is smaller than zero! Not a valid point_id.";
	}
}

void RightWidget::setHint3D(int point_id)
{
	hint_text->setText(QString::fromStdString("Selected point " + std::to_string(point_id + 1)+ " in 3D panel."));
}

void RightWidget::cancelHint()
{
	//std::string file_name = "Resources\\empty.png";
	//QPixmap image2_crop = processImage(file_name);
	//hint_container->setPixmap(image2_crop);
	hint_container->setText("Instructions: \n\n Please correct all the annotations in 3D first. \n\n Then check annotations in 2D views. \n\n Select one specific point in 3D, \n\nand use \" Correct This Point \" button to make fine correction. \n\n 2D Views are just for checking, not annotating. ");
}

void RightWidget::cancelHint3D()
{
	hint_text->setText("No point selected in 3D panel.");
}

void RightWidget::clearThumbnails()
{
	thumbnail_container->clear();
}

void RightWidget::addThumbnailItem(QListWidgetItem * item)
{
	thumbnail_container->addItem(item);
}

/*
void RightWidget::onThumbnailContainerItemClicked(QListWidgetItem* item) {
	if (item == NULL) {
		return;
	}
	// synchronizing the previous points if not propogated just now
	if (isJustPropagated == false)
	{
		updatePointsFromImageView();
	}
	// updating the new folder view
	QString name_string = item->text();
	QStringRef sub_string(&name_string, 0, 2);
	int number = sub_string.toInt();
	std::string folder_name = getFolderFromID(camera_order[number]);
	std::string full_path = getFullPathFromFolder(rootFolder, folder_name);
	imageView_->setImage(full_path);
	imageView_->setDetectedKeypoint(hash_all_imgs_detected_landmarks[folder_name]);
	imageView_->folder_name = folder_name;
	//this->initializeSelectedPoints();
}
*/



void RightWidget::resizeEvent(QResizeEvent *event) {
	int wid = this->width();
	//qDebug() << "MeshComtainer width adjusted! New width: " << wid;
	//qvtkWidget->resize(wid, 1208);

}

bool RightWidget::isEmpty()
{
	if (thumbnail_container->count() == 0)
	{
		return true;
	}
	return false;
}

QListWidgetItem * RightWidget::getCurrentItem()
{
	return thumbnail_container->currentItem();
}

void RightWidget::clear()
{
	thumbnail_container->clear();
}
