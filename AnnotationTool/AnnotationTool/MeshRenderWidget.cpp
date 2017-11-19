#include "MeshRenderWidget.h"
#include <QMouseEvent>
#include <QOpenGLShaderProgram>
#include <QCoreApplication>
#include <QGLFormat>

#include <QtWidgets>
#include <QtOpenGL>

MeshRenderWidget::MeshRenderWidget(QWidget *parent)
	: QGLWidget(parent)
{
	
	// load mesh
	qDebug() << "Start Loading Mesh" << endl;
	pcl::io::load("C:\\Users\\shangxuanu\\Desktop\\MyAnnotationTool\\AnnotationTool_Oculus\\PerformanceMeshRenderingTool\\meshcolor10.ply", mesh);
	//Reader.read("C:\\Users\\shangxuanu\\Desktop\\MyAnnotationTool\\AnnotationTool_Oculus\\PerformanceMeshRenderingTool\\simple.ply", cloud);
	pcl::fromPCLPointCloud2(mesh.cloud, cloud);
	qDebug() << "Finish Loading Mesh" << endl;
}

MeshRenderWidget::~MeshRenderWidget()
{

}

void MeshRenderWidget::initializeGL()
{
	glClearColor(0, 0, 0, 1);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
}

void MeshRenderWidget::paintGL()
{
	qDebug() << "Painting!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	glPushMatrix();
	//gluPerspective(45, (float)w / h, 1.0, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//glEnable(GL_DEPTH_TEST);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glEnable(GL_DEPTH_TEST);
}

void MeshRenderWidget::resizeGL(int w, int h)
{
	qDebug() << "Resizing!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glEnable(GL_DEPTH_TEST);
	
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < mesh.polygons.size(); i++) {
		//for (int i = mesh.polygons.size()-1; i >=0; i--) {
		for (int j = 0; j < 3; j++) {
			pcl::PointXYZRGB data = cloud.points[mesh.polygons[i].vertices[j]];
			glColor3f(data.r / color_ratio, data.g / color_ratio, data.b / color_ratio);
			glVertex3f(data.y / zoom_ratio, data.x / zoom_ratio, data.z / zoom_ratio);
		}
	}
	glEnd();

	//glPushMatrix();
	//glEnable(GL_CULL_FACE);
	glMatrixMode(GL_PROJECTION);
	gluPerspective(45, (float)w / h, 1.0, 100.0);
	glViewport(5, 5, w, h);
	gluLookAt(0, 0, 5, -0.1, 0, 0, 0, -1, 0); //  (eye_x, eye_y, eye_z, center_x, center_y, center_z, up_x, up_y, up_z)
}