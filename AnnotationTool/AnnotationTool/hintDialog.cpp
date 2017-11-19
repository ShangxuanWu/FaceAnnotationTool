#include <QtWidgets>

#include "hintDialog.h"

hintDialog::hintDialog(QWidget *parent)
	: QWidget(parent)
{
	// The window properties
	setWindowTitle("Instructions for 3D Panel");
	setWindowFlags(Qt::WindowStaysOnTopHint | Qt::WindowTitleHint | Qt::SubWindow);
	//setWindowFlags(Qt::WindowStaysOnTopHint);


	label_instruction_image = new QLabel(this);
	label_instruction_image->setGeometry(QRect(0, 0, INSTRUCTION_IMG_WIDTH, INSTRUCTION_IMG_HEIGHT));
	cancel();
	
	stylus_instruction = new QLabel(this);
	QPixmap image1("Resources\\stylus.png");
	stylus_instruction->setPixmap(image1);
	int h = image1.width();
	qDebug() << h;
	int w = image1.height();
	qDebug() << w;
	stylus_instruction->setGeometry(QRect(INSTRUCTION_IMG_WIDTH,  0, image1.width(), INSTRUCTION_IMG_HEIGHT));
	//label_instruction_image->setFixedSize(INSTRUCTION_IMG_WIDTH, INSTRUCTION_IMG_HEIGHT);

	//horizontal_layout = new QHBoxLayout(this);
	//horizontal_layout->addWidget(label_instruction_image);
	//horizontal_layout->addWidget(stylus_instruction);
	//horizontal_layout->setContentsMargins(10, 10, 10, 10);
}

void hintDialog::setHint(int id)
{
	if (id >= 0) {
		std::string file_name = "Resources\\" + std::to_string(id + 1) + ".png";
		QPixmap image2_crop = processImage(file_name);
		label_instruction_image->setPixmap(image2_crop);
	}
	else
	{
		cancel();
	}
}

QPixmap hintDialog::processImage(std::string fn)
{
	QPixmap image2(QString::fromStdString(fn));
	QRect rect((image2.width() - INSTRUCTION_IMG_WIDTH) / 2, (image2.height() - INSTRUCTION_IMG_HEIGHT) / 2, INSTRUCTION_IMG_WIDTH, INSTRUCTION_IMG_HEIGHT); // x, y, width, height
	QPixmap image2_crop = image2.copy(rect);
	return image2_crop;
}

void hintDialog::cancel()
{
	std::string file_name = "Resources\\empty.png";
	QPixmap image2_crop = processImage(file_name);
	label_instruction_image->setPixmap(image2_crop);
}