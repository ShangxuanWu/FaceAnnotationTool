// myglwidget.cpp

#include <QtWidgets>
#include <QtOpenGL>

#include "myglwidget.h"

MyGLWidget::MyGLWidget(QWidget *parent)
	: QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
	xRot = 0;
	yRot = 0;
	zRot = 0;
	tmp = 0;

	// load mesh
	qDebug() << "Start Loading Mesh" << endl;
	pcl::io::load("C:\\Users\\shangxuanu\\Desktop\\MyAnnotationTool\\AnnotationTool_Oculus\\PerformanceMeshRenderingTool\\meshcolor10.ply", mesh);
	pcl::fromPCLPointCloud2(mesh.cloud, cloud);
	qDebug() << "Finish Loading Mesh" << endl;

	// normalize
	double total_x = 0;
	double total_y = 0;
	double total_z = 0;

	for (int i = 0; i < mesh.polygons.size(); i++) {
		//for (int i = mesh.polygons.size()-1; i >=0; i--) {
		for (int j = 0; j < 3; j++) {
			pcl::PointXYZRGB data = cloud.points[mesh.polygons[i].vertices[j]];
			total_x += data.x;
			total_y += data.y;
			total_z += data.z;
		}
	}

	total_x /= (mesh.polygons.size() * 3);
	total_y /= (mesh.polygons.size() * 3);
	total_z /= (mesh.polygons.size() * 3);

	qDebug() << total_x << endl;
	qDebug() << total_y << endl;
	qDebug() << total_z << endl;


	for (int i = 0; i < mesh.polygons.size(); i++) {
		//for (int i = mesh.polygons.size()-1; i >=0; i--) {
		for (int j = 0; j < 3; j++) {
			pcl::PointXYZRGB data = cloud.points[mesh.polygons[i].vertices[j]];
			data.x -= total_x;
			data.y -= total_y;
			data.z -= total_z;
		}
	}
}

MyGLWidget::~MyGLWidget()
{
}

/*QSize MyGLWidget::minimumSizeHint() const
{
	return QSize(50, 50);
}

QSize MyGLWidget::sizeHint() const
{
	return QSize(1500, 700);
}*/

static void qNormalizeAngle(int &angle)
{
	while (angle < 0)
		angle += 360 * 2;
	while (angle > 360)
		angle -= 360 * 2;
}

void MyGLWidget::setXRotation(int angle)
{
	tmp += angle;
	qNormalizeAngle(angle);
	if (angle != xRot) {
		xRot = angle;
		emit xRotationChanged(angle);
		updateGL();
	}
}

void MyGLWidget::setYRotation(int angle)
{
	qNormalizeAngle(angle);
	if (angle != yRot) {
		yRot = angle;
		emit yRotationChanged(angle);
		updateGL();
	}
}

void MyGLWidget::setZRotation(int angle)
{
	qNormalizeAngle(angle);
	if (angle != zRot) {
		zRot = angle;
		emit zRotationChanged(angle);
		updateGL();
	}
}

void MyGLWidget::initializeGL()
{
	//qglClearColor(Qt::white);
	glClearColor(0, 0, 0, 1);

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);
	//glShadeModel(GL_SMOOTH);
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	//
	//static GLfloat lightPosition[4] = { 0, 0, 10, 1.0 };
	//glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
}

void MyGLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	//glTranslatef(0.0, 0.0, -10.0);
	//glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
	//glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
	//glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
	draw();
}

void MyGLWidget::resizeGL(int width, int height)
{
	int side = qMin(width, height);
	//glViewport((width - side) / 2, (height - side) / 2, side, side);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(-2, +2, -2, +2, 1.0, 15.0);
	// views
	gluPerspective(45, (float)width / height, 1.0, 100.0);
	glViewport(5, 5, width, height);
	gluLookAt(tmp, 0, 5, -0.1, 0, 0, 0, -1, 0); //  (eye_x, eye_y, eye_z, center_x, center_y, center_z, up_x,  up_y , up_z)
	// views
	glMatrixMode(GL_MODELVIEW);
}

void MyGLWidget::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();
}

void MyGLWidget::mouseMoveEvent(QMouseEvent *event)
{
	int dx = event->x() - lastPos.x();
	int dy = event->y() - lastPos.y();

	if (event->buttons() & Qt::LeftButton) {
		setXRotation(xRot + 1 * dy);
		setYRotation(yRot + 1 * dx);
	}
	else if (event->buttons() & Qt::RightButton) {
		setXRotation(xRot + 1 * dy);
		setZRotation(zRot + 1 * dx);
	}

	lastPos = event->pos();
}

void MyGLWidget::draw()
{
	/*qglColor(Qt::red);
	glBegin(GL_QUADS);
	glNormal3f(0, 0, -1);
	glVertex3f(-1, -1, 0);
	glVertex3f(-1, 1, 0);
	glVertex3f(1, 1, 0);
	glVertex3f(1, -1, 0);

	glEnd();
	glBegin(GL_TRIANGLES);
	glNormal3f(0, -1, 0.707);
	glVertex3f(-1, -1, 0);
	glVertex3f(1, -1, 0);
	glVertex3f(0, 0, 1.2);
	glEnd();
	glBegin(GL_TRIANGLES);
	glNormal3f(1, 0, 0.707);
	glVertex3f(1, -1, 0);
	glVertex3f(1, 1, 0);
	glVertex3f(0, 0, 1.2);
	glEnd();
	glBegin(GL_TRIANGLES);
	glNormal3f(0, 1, 0.707);
	glVertex3f(1, 1, 0);
	glVertex3f(-1, 1, 0);
	glVertex3f(0, 0, 1.2);
	glEnd();
	glBegin(GL_TRIANGLES);
	glNormal3f(-1, 0, 0.707);
	glVertex3f(-1, 1, 0);
	glVertex3f(-1, -1, 0);
	glVertex3f(0, 0, 1.2);
	glEnd();*/
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
}