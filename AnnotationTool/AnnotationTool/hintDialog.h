#ifndef HINTDIALOG_H
#define HINTDIALOG_H

#include <string>

#include <QWidget>
#include <QHBoxLayout>
#include <QLabel>
#include <QRect>
#include <QDebug>

#include "utils.h"

class QLabel;

class hintDialog : public QWidget
{
	Q_OBJECT

public:
	hintDialog(QWidget *parent = 0);
	void setHint(int id);
	void cancel();

private:
	QHBoxLayout* horizontal_layout;
	QLabel* stylus_instruction;
	QLabel* label_instruction_image;

	QPixmap processImage(std::string fn);
};

#endif