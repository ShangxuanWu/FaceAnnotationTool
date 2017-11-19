// myglwidget.h

#ifndef MYGLWIDGET_H
#define MYGLWIDGET_H

#include <QGLWidget>
#include<pcl/io/ply_io.h>
#include <pcl/io/auto_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <gl/GLU.h>
#include <gl/GL.h>

class MyGLWidget : public QGLWidget
{
	Q_OBJECT
public:
	explicit MyGLWidget(QWidget *parent = 0);
	~MyGLWidget();
signals:

	public slots :

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);

	//QSize minimumSizeHint() const;
	//QSize sizeHint() const;
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);

	public slots:
	// slots for xyz-rotation slider
	void setXRotation(int angle);
	void setYRotation(int angle);
	void setZRotation(int angle);

signals:
	// signaling rotation from mouse movement
	void xRotationChanged(int angle);
	void yRotationChanged(int angle);
	void zRotationChanged(int angle);

private:
	void draw();

	int xRot;
	int yRot;
	int zRot;

	int tmp;

	QPoint lastPos;

	pcl::PointCloud<pcl::PointXYZRGB> cloud;
	pcl::PolygonMesh mesh;
	pcl::PLYReader Reader;

	float zoom_ratio = 400.0;
	float color_ratio = 255.0;
};

#endif // MYGLWIDGET_H