#pragma once

#include <QWidget>
#include <QDebug>
#include <QGLWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include <gl/GLU.h>
#include <gl/GL.h>
#include<pcl/io/ply_io.h>
#include <pcl/io/auto_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

class MeshRenderWidget : public QGLWidget
{
	Q_OBJECT
public:
	MeshRenderWidget(QWidget *parent = 0);
	~MeshRenderWidget();

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

protected:
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	pcl::PointCloud<pcl::PointXYZRGB> cloud;
	pcl::PolygonMesh mesh;
	pcl::PLYReader Reader;

	QOpenGLVertexArrayObject m_vao;
	QOpenGLBuffer m_logoVbo;
	QOpenGLShaderProgram *m_program;
	bool m_core;
	bool m_transparent;
	int m_projMatrixLoc;
	int m_mvMatrixLoc;
	int m_normalMatrixLoc;
	int m_lightPosLoc;
	float zoom_ratio = 400.0;
	float color_ratio = 255.0;
private:
	int xRot;
	int yRot;
	int zRot;
};