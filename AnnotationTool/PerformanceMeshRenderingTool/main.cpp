#include <stdlib.h>//should before glew.h
#include <GL/glew.h>

//#ifndef BYTE
//#define BYTE unsigned char
//#endif

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

#include "GLSL/GLSLProgram.hpp"

#include <math.h>
#include <stdio.h>

#include "opencv2/opencv.hpp"


#include "Mesh/TriMesh.h"
#include "Mesh/TriMesh_algo.h"
#include "Mesh/GLCamera.h"
#include "Mesh/XForm.h"

#include "CVec.cpp"

#include "CCameraArray.h"


#ifdef WIN32
#include <conio.h>
#endif

#include "GlobalAppState.h"

#define OFF_SCREEN_RENDER false

#define DYNAMIC_SCENE
#define USEOPENMP

#define READ_MATRIX_FROM_CALIBRATION

#define TEXTURE_MAPPING

#define VIEW_DIST -5.0f

#define WRITE_TO_IMAGE


struct flag_mesh {
	GLuint vertex_buffer, element_buffer,texture_buffer;
	GLsizei element_count;
};

struct flag_vertex {
	GLfloat position[4];
	GLfloat normal[4];
	GLfloat diffuse[4];
	GLfloat shininess;
	GLfloat specular[4];
	GLfloat uvcoord[2];
};


struct CProj_mat
{
	 GLfloat p_matrix[16], mv_matrix[16];
	 GLfloat light_position[7][4];
};

struct Cprogram
{
	GLuint vertex_shader, fragment_shader, program;

	GLint position, normal, diffuse, shininess, specular,uvcoord;


	GLint p_matrix,mv_matrix;
	GLint light_position[7];
	GLint texture;
};

//for iccv submission
float g_array_proj[3][4];
float g_array_mv[16];
//end

//reference point
point g_staticpt;
point g_refpt;
point g_meshctr;


Cprogram g_prog;

static int g_mode = 0;


static int g_width = 480;		// The width of the texture we'll be rendering to
static int g_height = 640;//for canon seq
static int g_winwidth = 480;
static int g_winheight = 640;



int viewportscl = 1;
int viewportsclx = 8;
int viewportscly = 8;

static int fbo_width = g_width*viewportsclx;//(g_width/16)*4;
static int fbo_height = g_height*viewportscly;///4;


TriMesh::BSphere global_bsph;
xform global_xf;//the extrsinic for global camera
GLCamera g_camera;

char g_prefixname[256];
int g_mouseindx = 0;
bool g_start_mouse_cap = false;
bool g_pointcloud_render = false;


//TriMesh *g_pdebugmesh;



xform g_xforms;//the extrinsic for local camera,say rotation or translation
xform g_extrinsic;
xform g_intrinsic;
TriMesh *g_pmesh;
std::string g_xffilenames;

static int g_index = 0;
static int g_rotindex = 0;

static int g_imgindex = 0;
static int g_nummesh = 50;
static int g_viewind = 0;
static int g_initial_index = 100;
static int g_end_index = 100;
static int g_rotate_init = 0;
static int g_frameinterval = 1;


//for brdf control setting 0
float g_brdf_sp_scl = 1.0f;
float g_brdf_shy_scl = 1.0f;

//for brdf control: setting 1
//float g_brdf_sp_scl = 0.91f;
//float g_brdf_shy_scl = 10.0f;


//for brdf control setting 2
//float g_brdf_sp_scl = 2.0f;
//float g_brdf_shy_scl = 3.0f;


float g_brdf_diffuse_col = 1.0f;
bool g_b_texture = true;



GLuint fbo;					// Our handle to the FBO
GLuint depthBuffer;			// Our handle to the depth render buffer
GLuint img, img2;			// Our handle to a textures



// Used for drawing the 3D cube with our rendered texture on it
GLfloat	xrot = 0;			// X Rotation
GLfloat	yrot = 0;			// Y Rotation
GLfloat xspeed = 0.2f;		// X Rotation Speed
GLfloat yspeed = 0.1f;		// Y Rotation Speed

GLSLProgram * mrtprogram;	// The shader we will be using for MRT rendering output
flag_vertex *g_pvertex_array;
flag_mesh g_task_mesh;
CProj_mat g_projmat;

float lit_pos[4] = {0.408248, -0.816497, 0.408248,0.0};


vector<TriMesh*> g_vec_mesh;
FILE *g_pcalfile, *g_logfile;
std::string resultdir,g_outputdir;

std::string camfilename;

vector<vector<CVec3f>> gvec_hairmap;

void update_mesh(TriMesh * pmesh,flag_mesh *pout_mesh)
{	
	int numvertices = pmesh->vertices.size();

	flag_vertex *vertex_data
		= (flag_vertex*) malloc(numvertices * sizeof(flag_vertex));


	vector<Color> vec_diffuse_alb(numvertices);
	if(g_b_texture)
	{
		for (int i=0;i<numvertices;i++)
		{
			vec_diffuse_alb[i] = pmesh->colors[i];
		}
	}
	else
	{
		//float m_col_H = g_brdf_diffuse_col;
		//float m_col_S = 0.75;
		//float m_col_V = 0.75;

		Color tmpdiffuse_alb(123,128,164) ;//= Color::hsv(m_col_H,m_col_S,m_col_V);
		for (int i=0;i<numvertices;i++)
		{
			vec_diffuse_alb[i] =tmpdiffuse_alb;
		}
	}


	//pmesh->need_adjacentfaces();
	for (int i=0;i<numvertices;i++)
	{
		vertex_data[i].position[0] = pmesh->vertices[i][0];
		vertex_data[i].position[1] = pmesh->vertices[i][1];
		vertex_data[i].position[2] = pmesh->vertices[i][2];
		vertex_data[i].position[3] = 1.0f;

		vertex_data[i].normal[0] = pmesh->normals[i][0];
		vertex_data[i].normal[1] = pmesh->normals[i][1];
		vertex_data[i].normal[2] = pmesh->normals[i][2];
		//int faceind = pmesh->adjacentfaces[i][0];
		//vertex_data[i].normal[0] = pmesh->faces[faceind].norfac[0];
		//vertex_data[i].normal[1] = pmesh->faces[faceind].norfac[1];
		//vertex_data[i].normal[2] = pmesh->faces[faceind].norfac[2];
		//vertex_data[i].normal[3] = 1.0f;

		vertex_data[i].diffuse[0] = vec_diffuse_alb[i][0];
		vertex_data[i].diffuse[1] = vec_diffuse_alb[i][1];
		vertex_data[i].diffuse[2] = vec_diffuse_alb[i][2];
		vertex_data[i].diffuse[3] = 1.0f;

		//Color m_tmpcol;
		//m_tmpcol[0] = (pmesh->normals[i][0]+1.0)/2.0;
		//m_tmpcol[1] = (pmesh->normals[i][1]+1.0)/2.0;
		//m_tmpcol[2] = (pmesh->normals[i][2]+1.0)/2.0;
		//vertex_data[i].diffuse[0] = m_tmpcol[0];
		//vertex_data[i].diffuse[1] = m_tmpcol[1];
		//vertex_data[i].diffuse[2] = m_tmpcol[2];
		//vertex_data[i].diffuse[3] = 1.0f;

		vertex_data[i].shininess   = g_brdf_shy_scl;
		vertex_data[i].specular[0] = g_brdf_sp_scl;
		vertex_data[i].specular[1] = g_brdf_sp_scl;
		vertex_data[i].specular[2] = g_brdf_sp_scl;
		vertex_data[i].specular[3] = 1.0;

#ifdef TEXTURE_MAPPING
		vertex_data[i].uvcoord[0] = pmesh->uvcoord[i][0];
		vertex_data[i].uvcoord[1] = pmesh->uvcoord[i][1];
#endif
	}


	GLsizei numfaces = g_pmesh->faces.size();

	GLsizei element_count = 3 * numfaces;

	GLuint *element_data
		= (GLuint*) malloc(element_count * sizeof(GLuint));
	GLuint index;

	for (int i=0;i<numfaces;i++)
	{
		element_data[i*3] = pmesh->faces[i][0];
		element_data[i*3+1] = pmesh->faces[i][1];
		element_data[i*3+2] = pmesh->faces[i][2];
	}

	pout_mesh->element_count = element_count;

	glBindBuffer(GL_ARRAY_BUFFER, pout_mesh->vertex_buffer);
	glBufferData(
		GL_ARRAY_BUFFER,
		numvertices * sizeof(flag_vertex),
		vertex_data,
		GL_STREAM_DRAW
		);


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pout_mesh->element_buffer);
	glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		element_count * sizeof(GLuint),
		element_data,
		GL_STATIC_DRAW
		);
	

	//texture
	int m_tmpwidth = 2000;
	int m_tmpheight = 2000;
	int m_tmpsizegrid = 40;
	BYTE *pixels = (BYTE*)malloc(3*m_tmpwidth*m_tmpheight * sizeof(BYTE));
	for(int i=0;i<m_tmpwidth;i++)
	{
		for(int j=0;j<m_tmpheight;j++)
		{
			int m_tmpind = j*m_tmpwidth+i; 
			BYTE tmpval = 0;

			int m_x = i/m_tmpsizegrid;
			int m_y = j/m_tmpsizegrid;

			bool b_x= true;
			bool b_y = true;
			if(m_x%2==0)
				b_x=false;
			if(m_y%2==0)
				b_y = false;

			if(b_x&&b_y)
			{
				pixels[m_tmpind*3] = 0;			
				pixels[m_tmpind*3+1] = 0;			
				pixels[m_tmpind*3+2] = 0;			
			}else
			{
				pixels[m_tmpind*3] = 255;
				pixels[m_tmpind*3+1] = 255;
				pixels[m_tmpind*3+2] = 255;
			}

			//pixels[m_tmpind*3] = 255;
			//pixels[m_tmpind*3+1] = 255;
			//pixels[m_tmpind*3+2] = 255;
		}
	}
	glBindTexture(GL_TEXTURE_2D, pout_mesh->texture_buffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,     GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,     GL_CLAMP_TO_EDGE);
    glTexImage2D(
        GL_TEXTURE_2D, 0,           /* target, level */
        GL_RGB8,                    /* internal format */
        m_tmpwidth, m_tmpheight, 0,           /* width, height, border */
        GL_BGR, GL_UNSIGNED_BYTE,   /* external format, type */
        pixels                      /* pixels */
    );

	free((void*)element_data);
	free((void*)vertex_data);
	free(pixels);

}


void update_mesh(TriMesh * pmesh, flag_mesh *pout_mesh, const std::string &textureimgname)
{
	int numvertices = pmesh->vertices.size();

	flag_vertex *vertex_data
		= (flag_vertex*)malloc(numvertices * sizeof(flag_vertex));


	vector<Color> vec_diffuse_alb(numvertices);
	if (g_b_texture)
	{
		for (int i = 0; i<numvertices; i++)
		{
			vec_diffuse_alb[i] = pmesh->colors[i];
		}
	}
	else
	{
		//float m_col_H = g_brdf_diffuse_col;
		//float m_col_S = 0.75;
		//float m_col_V = 0.75;

		Color tmpdiffuse_alb(123, 128, 164);//= Color::hsv(m_col_H,m_col_S,m_col_V);
		for (int i = 0; i<numvertices; i++)
		{
			vec_diffuse_alb[i] = tmpdiffuse_alb;
		}
	}


	for (int i = 0; i<numvertices; i++)
	{
		vertex_data[i].position[0] = pmesh->vertices[i][0];
		vertex_data[i].position[1] = pmesh->vertices[i][1];
		vertex_data[i].position[2] = pmesh->vertices[i][2];
		vertex_data[i].position[3] = 1.0f;

		vertex_data[i].normal[0] = pmesh->normals[i][0];
		vertex_data[i].normal[1] = pmesh->normals[i][1];
		vertex_data[i].normal[2] = pmesh->normals[i][2];
		vertex_data[i].normal[3] = 1.0f;

		vertex_data[i].diffuse[0] = vec_diffuse_alb[i][0];
		vertex_data[i].diffuse[1] = vec_diffuse_alb[i][1];
		vertex_data[i].diffuse[2] = vec_diffuse_alb[i][2];
		vertex_data[i].diffuse[3] = 1.0f;

		//Color m_tmpcol;
		//m_tmpcol[0] = (pmesh->normals[i][0]+1.0)/2.0;
		//m_tmpcol[1] = (pmesh->normals[i][1]+1.0)/2.0;
		//m_tmpcol[2] = (pmesh->normals[i][2]+1.0)/2.0;
		//vertex_data[i].diffuse[0] = m_tmpcol[0];
		//vertex_data[i].diffuse[1] = m_tmpcol[1];
		//vertex_data[i].diffuse[2] = m_tmpcol[2];
		//vertex_data[i].diffuse[3] = 1.0f;

		vertex_data[i].shininess = g_brdf_shy_scl;
		vertex_data[i].specular[0] = g_brdf_sp_scl;
		vertex_data[i].specular[1] = g_brdf_sp_scl;
		vertex_data[i].specular[2] = g_brdf_sp_scl;
		vertex_data[i].specular[3] = 1.0;

#ifdef TEXTURE_MAPPING
		vertex_data[i].uvcoord[0] = pmesh->uvcoord[i][0];
		vertex_data[i].uvcoord[1] = pmesh->uvcoord[i][1];
#endif
	}


	GLsizei numfaces = g_pmesh->faces.size();

	GLsizei element_count = 3 * numfaces;

	GLuint *element_data
		= (GLuint*)malloc(element_count * sizeof(GLuint));
	GLuint index;

	for (int i = 0; i<numfaces; i++)
	{
		element_data[i * 3] = pmesh->faces[i][0];
		element_data[i * 3 + 1] = pmesh->faces[i][1];
		element_data[i * 3 + 2] = pmesh->faces[i][2];
	}

	pout_mesh->element_count = element_count;

	glBindBuffer(GL_ARRAY_BUFFER, pout_mesh->vertex_buffer);
	glBufferData(
		GL_ARRAY_BUFFER,
		numvertices * sizeof(flag_vertex),
		vertex_data,
		GL_STREAM_DRAW
		);


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pout_mesh->element_buffer);
	glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		element_count * sizeof(GLuint),
		element_data,
		GL_STATIC_DRAW
		);



	cv::Mat textimg = cv::imread(textureimgname.c_str(),1);

	int tmpwidth = textimg.cols;
	int tmpheight = textimg.rows;

	std::cout << tmpwidth << std::endl;
	std::cout << tmpheight << std::endl;
	
	BYTE *pixels = textimg.data;
	
	glBindTexture(GL_TEXTURE_2D, pout_mesh->texture_buffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(
		GL_TEXTURE_2D, 0,           /* target, level */
		GL_RGB8,                    /* internal format */
		tmpwidth, tmpheight, 0,           /* width, height, border */
		GL_BGR, GL_UNSIGNED_BYTE,   /* external format, type */
		pixels                      /* pixels */
		);

	free((void*)element_data);
	free((void*)vertex_data);
	//free(pixels);

}



flag_vertex *update_mesh(const char * prefix,flag_mesh *pout_mesh)
{
	char filename[256];
	sprintf(filename,"%s.off",prefix);
	printf("%s\n",filename);
	g_pmesh = TriMesh::read(filename);
	GLsizei numvertices = g_pmesh->vertices.size();
	g_pmesh->need_normals();
	g_pmesh->need_bsphere();

	// 	bilateral_smooth_mesh(pmesh, 20);
	// 	pmesh->write("smooth.ply");

	// 	subdiv(pmesh, SUBDIV_LOOP);
	// 	pmesh->write("subdiv.ply");

	//sprintf(filename,"albedo_%s.txt",prefix);
	sprintf(filename,"%s_color.txt",prefix);
	FILE *palbedofile = fopen(filename,"r"); 
	if (g_pmesh->colors.size()!=numvertices)
	{
		g_pmesh->colors.resize(numvertices);
	}

	for (int i=0;i<numvertices;++i)
	{
		Color vec_color(1.0,1.0,1.0);

		if(palbedofile!=NULL)
			fscanf(palbedofile,"%f %f %f\n",&vec_color[0],&vec_color[1],&vec_color[2]);

		vec_color = vec_color/180.0f;
		g_pmesh->colors[i] = vec_color;
	}
	if(palbedofile!=NULL)
		fclose(palbedofile);


	flag_vertex *vertex_data
		= (flag_vertex*) malloc(numvertices * sizeof(flag_vertex));

	for (int i=0;i<numvertices;i++)
	{
		vertex_data[i].position[0] = g_pmesh->vertices[i][0];
		vertex_data[i].position[1] = g_pmesh->vertices[i][1];
		vertex_data[i].position[2] = g_pmesh->vertices[i][2];
		vertex_data[i].position[3] = 1.0f;

		vertex_data[i].normal[0] = g_pmesh->normals[i][0];
		vertex_data[i].normal[1] = g_pmesh->normals[i][1];
		vertex_data[i].normal[2] = g_pmesh->normals[i][2];
		vertex_data[i].normal[3] = 1.0f;

		vertex_data[i].diffuse[0] = g_pmesh->colors[i][0];
		vertex_data[i].diffuse[1] = g_pmesh->colors[i][1];
		vertex_data[i].diffuse[2] = g_pmesh->colors[i][2];
		vertex_data[i].diffuse[3] = 1.0f;

		vertex_data[i].shininess   = 0.0f;
		vertex_data[i].specular[0] = 0;
		vertex_data[i].specular[1] = 0;
		vertex_data[i].specular[2] = 0;
		vertex_data[i].specular[3] = 0;
	}



	GLsizei numfaces = g_pmesh->faces.size();

	GLsizei element_count = 3 * numfaces;

	GLuint *element_data
		= (GLuint*) malloc(element_count * sizeof(GLuint));
	GLuint index;

	for (int i=0;i<numfaces;i++)
	{
		element_data[i*3] = g_pmesh->faces[i][0];
		element_data[i*3+1] = g_pmesh->faces[i][1];
		element_data[i*3+2] = g_pmesh->faces[i][2];
	}

	pout_mesh->element_count = element_count;

	glBindBuffer(GL_ARRAY_BUFFER, pout_mesh->vertex_buffer);
	glBufferData(
		GL_ARRAY_BUFFER,
		numvertices * sizeof(flag_vertex),
		vertex_data,
		GL_STREAM_DRAW
		);



	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pout_mesh->element_buffer);
	glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		element_count * sizeof(GLuint),
		element_data,
		GL_STATIC_DRAW
		);

	free((void*)element_data);


	return vertex_data;
}

static void argConvGLcpara2_anothertype( float cparam[3][4], int width, int height, float
							gnear, float gfar, float m[16] )
{
	float   p[3][3], q[4][4];


	/* Camera parameters are converted openGL representation. */
	/* Camera parameter represents the transformation from camera coordinates
	to screen coordinates[pixel]. OpenGL projection matrix represents the
	transformation from the camera coordinates to normalized view volume. */
	for(int  i = 0; i < 3; i++ ) {
		for(int j = 0; j < 3; j++ ) {
			p[i][j] = cparam[i][j];// / cparam[2][2];
		}
	}
	q[0][0] = (2.0 * p[0][0] / width);
	q[0][1] = (2.0 * p[0][1] / width);
	q[0][2] = ((2.0 * p[0][2] / width)  - 1.0);
	q[0][3] = 0.0;

	q[1][0] = 0.0;
	q[1][1] = (2.0 * p[1][1] / height);
	q[1][2] = ((2.0 * p[1][2] / height) - 1.0);
	q[1][3] = 0.0;

	q[2][0] = 0.0;
	q[2][1] = 0.0;
	q[2][2] = -(gfar + gnear)/(gfar - gnear);
	q[2][3] = -2.0 * gfar * gnear / (gfar - gnear);

	q[3][0] = 0.0;
	q[3][1] = 0.0;
	q[3][2] = -1.0;
	q[3][3] = 0.0;

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			m[j*4+i] = q[i][j];
		}
	}


}


static void argConvGLcpara2( float cparam[3][4], int width, int height, float
							gnear, float gfar, float m[16] )
{
	float   p[3][3], q[4][4];


	/* Camera parameters are converted openGL representation. */
	/* Camera parameter represents the transformation from camera coordinates
	to screen coordinates[pixel]. OpenGL projection matrix represents the
	transformation from the camera coordinates to normalized view volume. */
	for(int  i = 0; i < 3; i++ ) {
		for(int j = 0; j < 3; j++ ) {
			p[i][j] = cparam[i][j];// / cparam[2][2];
		}
	}
	q[0][0] = (2.0 * p[0][0] / width);
	q[0][1] = (2.0 * p[0][1] / width);
	q[0][2] = ((2.0 * p[0][2] / width)  - 1.0);
	q[0][3] = 0.0;

	q[1][0] = 0.0;
	q[1][1] = -(2.0 * p[1][1] / height);
	q[1][2] = -((2.0 * p[1][2] / height) - 1.0);
	q[1][3] = 0.0;

	q[2][0] = 0.0;
	q[2][1] = 0.0;
	q[2][2] = (gfar + gnear)/(gfar - gnear);
	q[2][3] = -2.0 * gfar * gnear / (gfar - gnear);

	q[3][0] = 0.0;
	q[3][1] = 0.0;
	q[3][2] = 1.0;
	q[3][3] = 0.0;

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			m[j*4+i] = q[i][j];
		}
	}
}

Color GetImgInten(IplImage * pinimg, CvPoint2D32f imgpos,bool IsBilinear )
{
	Color res;
	res = 0.0f; 

	int m_width = pinimg->width;
	int m_height = pinimg->height;

	if(IsBilinear)
	{
		//do linear interpolation
		int xpos1 = cvFloor(imgpos.x);
		int ypos1 = cvFloor(imgpos.y);
		int xpos  = xpos1+1;
		int ypos  = ypos1+1;

		float cof1,cof2,cof3,cof4;
		cof1=imgpos.x-xpos1;
		cof2=xpos-imgpos.x;
		cof3=imgpos.y-ypos1;
		cof4=ypos-imgpos.y;

		int nInd1 = ((int)ypos1 * m_width + (int)(xpos1)) * pinimg->nChannels;
		int nInd2 = ((int)ypos1 * m_width + (int)(xpos)) * pinimg->nChannels;
		int nInd3 = ((int)ypos * m_width + (int)(xpos1)) * pinimg->nChannels;
		int nInd4  =((int)ypos * m_width + (int)(xpos)) * pinimg->nChannels;

		if(pinimg->nChannels==3)
		{
			int nB1 = (BYTE)pinimg->imageData[nInd1];
			int nG1 = (BYTE)pinimg->imageData[nInd1 + 1];
			int nR1 = (BYTE)pinimg->imageData[nInd1 + 2];
			//int nGr1 = 0.299 * nR1 + 0.587 * nG1 + 0.114 * nB1;


			int nB2 = (BYTE)pinimg->imageData[nInd2];
			int nG2 = (BYTE)pinimg->imageData[nInd2 + 1];
			int nR2 = (BYTE)pinimg->imageData[nInd2 + 2];
			//int nGr2 = 0.299 * nR2 + 0.587 * nG2 + 0.114 * nB2;

			int nB3 = (BYTE)pinimg->imageData[nInd3];
			int nG3 = (BYTE)pinimg->imageData[nInd3 + 1];
			int nR3 = (BYTE)pinimg->imageData[nInd3 + 2];
			//int nGr3 = 0.299 * nR3 + 0.587 * nG3 + 0.114 * nB3;

			int nB4 = (BYTE)pinimg->imageData[nInd4];
			int nG4 = (BYTE)pinimg->imageData[nInd4 + 1];
			int nR4 = (BYTE)pinimg->imageData[nInd4 + 2];
			//int nGr4 = 0.299 * nR4 + 0.587 * nG4 + 0.114 * nB4;

			float tmp1=cof1*nB2+cof2*nB1;
			float tmp2=cof1*nB4+cof2*nB3;
			float nB=cof3*tmp2+cof4*tmp1;

			tmp1=cof1*nG2+cof2*nG1;
			tmp2=cof1*nG4+cof2*nG3;
			float nG=cof3*tmp2+cof4*tmp1;

			tmp1=cof1*nR2+cof2*nR1;
			tmp2=cof1*nR4+cof2*nR3;
			float nR=cof3*tmp2+cof4*tmp1;

			res[0] = nR/255.0f;
			res[1] = nG/255.0f;
			res[2] = nB/255.0f;

			//res = 0.299 * nR + 0.587 * nG + 0.114 * nB;
		}
		else
		{
			int nGr1 = (BYTE)pinimg->imageData[nInd1];
			int nGr2 = (BYTE)pinimg->imageData[nInd2];
			int nGr3 = (BYTE)pinimg->imageData[nInd3];
			int nGr4 = (BYTE)pinimg->imageData[nInd4];

			float tmp1=cof1*nGr2+cof2*nGr1;
			float tmp2=cof1*nGr4+cof2*nGr3;
			float tmpval=cof3*tmp2+cof4*tmp1;

			res = tmpval/255.0f;

		}


	}
	else
	{
		int xpos = cvRound(imgpos.x);
		int ypos = cvRound(imgpos.y);

		int nInd = ((int)ypos * m_width + (int)(xpos)) * pinimg->nChannels;

		if(pinimg->nChannels==3)
		{

			int nB = (BYTE)pinimg->imageData[nInd];
			int nG = (BYTE)pinimg->imageData[nInd + 1];
			int nR = (BYTE)pinimg->imageData[nInd + 2];

			res[0] = nR / 255.0f;
			res[1] = nG / 255.0f;
			res[2] = nB / 255.0f;

			//res = 0.299 * nR + 0.587 * nG + 0.114 * nB;
		}
		else
		{
			int nGr = (BYTE)pinimg->imageData[nInd];

			res=(float)nGr;

		}
	}

	return res;
}

static void update_matrix_read_test (float array_proj[3][4],float array_mv[16], CProj_mat * resource)
{
	
// 	for (int j=0;j<3;j++)
// 	{
// 		array_proj[0][j] = array_proj[0][j]*fbo_width/g_width;
// 		array_proj[1][j] = array_proj[1][j]*fbo_height/g_height;
// 
// 	}

	//calculate the  near plane and far plane respectively
	xform mat_mv(array_mv);
	xform mat_tmp= transpose(rot_only(mat_mv))*trans_only(mat_mv);


	point cam_ctr(-mat_tmp[12], -mat_tmp[13], -mat_tmp[14]);
	float dist_cam2near =len(cam_ctr - g_pmesh->bsphere.center);

	float m_scale = 1.5;
	float nearplane = dist_cam2near - m_scale*g_pmesh->bsphere.r;
	float farplane = dist_cam2near + m_scale*g_pmesh->bsphere.r;

//	printf("near = %f, far = %f\n",nearplane,farplane);	
	

 	//array_proj[0][2] = g_width/2;
 	//array_proj[1][2] = g_height/2;

	//  	array_proj[0][2] *= -1.0;
	//  	array_proj[1][2] *= -1.0;
	//  	array_proj[2][2] *= -1.0;

// 	nearplane = 100;
// 	farplane = 10000;

	argConvGLcpara2_anothertype( array_proj, g_width, g_height,nearplane, farplane, resource->p_matrix );

//  	for (int i=0;i<4;i++)
//  	{
//  		printf("pre: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
 


// another way to generate the instrinsic parameter, however the y and z axles should re reverted
// 	double m_top = array_proj[1][2]*nearplane/array_proj[1][1];
// 	double m_bottom = (array_proj[1][2]-g_height)*nearplane/array_proj[1][1];
// 	double m_left = -array_proj[0][2]*nearplane/array_proj[0][0];
// 	double m_right = (g_width-array_proj[0][2])*nearplane/array_proj[0][0];
// 	glMatrixMode(GL_PROJECTION);
// 	glLoadIdentity();
// 	glFrustum(m_left,m_right,m_bottom,m_top,nearplane,farplane);
//  	glGetFloatv(GL_PROJECTION_MATRIX,resource->p_matrix);
//   	for (int i=0;i<4;i++)
//  	{
//  		printf("latter: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
// 	//-z and -y if generate intrinsic by glFrustum
 	for(int i = 0; i<4; ++i)
 	{
 		array_mv[i*4+1] *= -1;
 		array_mv[i*4+2] *= -1; 
 	}
//end of glFrustum 

	for (int i=0;i<16;i++)
	{
		resource->mv_matrix[i] = array_mv[i];
	}
}





static void update_matrix_read (float array_proj[3][4],float array_mv[16], CProj_mat * resource, int width, int height)
{
	
// 	for (int j=0;j<3;j++)
// 	{
// 		array_proj[0][j] = array_proj[0][j]*fbo_width/g_width;
// 		array_proj[1][j] = array_proj[1][j]*fbo_height/g_height;
// 
// 	}

	//calculate the  near plane and far plane respectively
	xform mat_mv(array_mv);
	xform mat_tmp= transpose(rot_only(mat_mv))*trans_only(mat_mv);


	point cam_ctr(-mat_tmp[12], -mat_tmp[13], -mat_tmp[14]);
	float dist_cam2near =len(cam_ctr - g_pmesh->bsphere.center);

	float m_scale = 1.5;
	float nearplane = dist_cam2near - m_scale*g_pmesh->bsphere.r;
	float farplane = dist_cam2near + m_scale*g_pmesh->bsphere.r;

	nearplane = 600;
	farplane = 1200.0f;

	//float nearplane = dist_cam2near - m_scale*g_pdebugmesh->bsphere.r;
	//float farplane = dist_cam2near + m_scale*g_pdebugmesh->bsphere.r;
	

//	printf("near = %f, far = %f\n",nearplane,farplane);	
	

 	//array_proj[0][2] = g_width/2;
 	//array_proj[1][2] = g_height/2;

	//  	array_proj[0][2] *= -1.0;
	//  	array_proj[1][2] *= -1.0;
	//  	array_proj[2][2] *= -1.0;

// 	nearplane = 100;
// 	farplane = 10000;

	argConvGLcpara2(array_proj, width, height, nearplane, farplane, resource->p_matrix);

//  	for (int i=0;i<4;i++)
//  	{
//  		printf("pre: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
 


// another way to generate the instrinsic parameter, however the y and z axles should re reverted
// 	double m_top = array_proj[1][2]*nearplane/array_proj[1][1];
// 	double m_bottom = (array_proj[1][2]-g_height)*nearplane/array_proj[1][1];
// 	double m_left = -array_proj[0][2]*nearplane/array_proj[0][0];
// 	double m_right = (g_width-array_proj[0][2])*nearplane/array_proj[0][0];
// 	glMatrixMode(GL_PROJECTION);
// 	glLoadIdentity();
// 	glFrustum(m_left,m_right,m_bottom,m_top,nearplane,farplane);
//  	glGetFloatv(GL_PROJECTION_MATRIX,resource->p_matrix);
//   	for (int i=0;i<4;i++)
//  	{
//  		printf("latter: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
// 	//-z and -y if generate intrinsic by glFrustum
// 	for(int i = 0; i<4; ++i)
// 	{
// 		array_mv[i*4+1] *= -1;
// 		array_mv[i*4+2] *= -1; 
// 	}
//end of glFrustum 

	for (int i=0;i<16;i++)
	{
		resource->mv_matrix[i] = array_mv[i];
	}
}




static void update_matrix_read (const char * filename,CProj_mat * resource)
{
	FILE *infile = fopen(filename, "r");
	if (infile == NULL)
	{
		fprintf(stderr, "Can't open '%s'! Exiting...\n", filename);
		exit(-1);
	}

	float dump_proj[3][4];

	fscanf(infile, "P:\n");
	fscanf(infile, "%f %f %f %f\n", &dump_proj[0][0], &dump_proj[0][1], &dump_proj[0][2], &dump_proj[0][3]);
	fscanf(infile, "%f %f %f %f\n", &dump_proj[1][0], &dump_proj[1][1], &dump_proj[1][2], &dump_proj[1][3]);
	fscanf(infile, "%f %f %f %f\n", &dump_proj[2][0], &dump_proj[2][1], &dump_proj[2][2], &dump_proj[2][3]);

// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[0][0], dump_proj[0][1], dump_proj[0][2], dump_proj[0][3]);
// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[1][0], dump_proj[1][1], dump_proj[1][2], dump_proj[1][3]);
// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[2][0], dump_proj[2][1], dump_proj[2][2], dump_proj[2][3]);
// 	fprintf(g_pcalfile,"\n");


	float array_proj[3][4];
	fscanf(infile, "\nK:\n");
	fscanf(infile, "%f %f %f\n", &array_proj[0][0], &array_proj[0][1], &array_proj[0][2]);
	fscanf(infile, "%f %f %f\n", &array_proj[1][0], &array_proj[1][1], &array_proj[1][2]);
	fscanf(infile, "%f %f %f\n", &array_proj[2][0], &array_proj[2][1], &array_proj[2][2]);

	
// 	for (int j=0;j<3;j++)
// 	{
// 		array_proj[0][j] = array_proj[0][j]*fbo_width/g_width;
// 		array_proj[1][j] = array_proj[1][j]*fbo_height/g_height;
// 
// 	}

	float array_mv[16];
	fscanf(infile, "\nM:\n");
	fscanf(infile, "%f %f %f %f\n", &array_mv[0], &array_mv[4], &array_mv[8], &array_mv[12]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[1], &array_mv[5], &array_mv[9], &array_mv[13]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[2], &array_mv[6], &array_mv[10], &array_mv[14]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[3], &array_mv[7], &array_mv[11], &array_mv[15]);

	fprintf(g_pcalfile,"0\n");
	fprintf(g_pcalfile,"%f %f %f\n", array_proj[0][0], array_proj[0][1], array_proj[0][2]);
	fprintf(g_pcalfile,"%f %f %f\n", array_proj[1][0], array_proj[1][1], array_proj[1][2]);
	fprintf(g_pcalfile,"%f %f %f\n", array_proj[2][0], array_proj[2][1], array_proj[2][2]);
	fprintf(g_pcalfile,"0.0 0.0\n");


	fprintf(g_pcalfile,"%f %f %f %f\n", array_mv[0], array_mv[4], array_mv[8], array_mv[12]);
	fprintf(g_pcalfile,"%f %f %f %f\n", array_mv[1], array_mv[5], array_mv[9], array_mv[13]);
	fprintf(g_pcalfile,"%f %f %f %f\n", array_mv[2], array_mv[6], array_mv[10], array_mv[14]);
	fprintf(g_pcalfile,"\n");

	fclose(infile);

	//calculate the  near plane and far plane respectively
	xform mat_mv(array_mv);
	xform mat_tmp= transpose(rot_only(mat_mv))*trans_only(mat_mv);


	point cam_ctr(-mat_tmp[12], -mat_tmp[13], -mat_tmp[14]);
	float dist_cam2near =len(cam_ctr - g_pmesh->bsphere.center);

	float m_scale = 1.0;
	float nearplane = dist_cam2near - m_scale*g_pmesh->bsphere.r;
	float farplane = dist_cam2near + m_scale*g_pmesh->bsphere.r;

//	printf("near = %f, far = %f\n",nearplane,farplane);	
	

 	//array_proj[0][2] = g_width/2;
 	//array_proj[1][2] = g_height/2;

	//  	array_proj[0][2] *= -1.0;
	//  	array_proj[1][2] *= -1.0;
	//  	array_proj[2][2] *= -1.0;

// 	nearplane = 100;
// 	farplane = 10000;

	argConvGLcpara2( array_proj, g_width, g_height,nearplane, farplane, resource->p_matrix );

//  	for (int i=0;i<4;i++)
//  	{
//  		printf("pre: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
 


// another way to generate the instrinsic parameter, however the y and z axles should re reverted
// 	double m_top = array_proj[1][2]*nearplane/array_proj[1][1];
// 	double m_bottom = (array_proj[1][2]-g_height)*nearplane/array_proj[1][1];
// 	double m_left = -array_proj[0][2]*nearplane/array_proj[0][0];
// 	double m_right = (g_width-array_proj[0][2])*nearplane/array_proj[0][0];
// 	glMatrixMode(GL_PROJECTION);
// 	glLoadIdentity();
// 	glFrustum(m_left,m_right,m_bottom,m_top,nearplane,farplane);
//  	glGetFloatv(GL_PROJECTION_MATRIX,resource->p_matrix);
//   	for (int i=0;i<4;i++)
//  	{
//  		printf("latter: %f %f %f %f\n",resource->p_matrix[i],resource->p_matrix[i+4*1],resource->p_matrix[i+4*2],resource->p_matrix[i+4*3]);
//  	}
// 	//-z and -y if generate intrinsic by glFrustum
// 	for(int i = 0; i<4; ++i)
// 	{
// 		array_mv[i*4+1] *= -1;
// 		array_mv[i*4+2] *= -1; 
// 	}
//end of glFrustum 

	for (int i=0;i<16;i++)
	{
		resource->mv_matrix[i] = array_mv[i];
	}
}

void *file_contents(const char *filename, GLint *length)
{
	FILE *f = fopen(filename, "r");
	void *buffer;

	if (!f) {
		fprintf(stderr, "Unable to open %s for reading\n", filename);
		return NULL;
	}

	fseek(f, 0, SEEK_END);
	*length = ftell(f);
	fseek(f, 0, SEEK_SET);

	buffer = malloc(*length+1);
	*length = fread(buffer, 1, *length, f);
	fclose(f);
	((char*)buffer)[*length] = '\0';

	return buffer;
}

void show_info_log(
				   GLuint object,
				   PFNGLGETSHADERIVPROC glGet__iv,
				   PFNGLGETSHADERINFOLOGPROC glGet__InfoLog
				   )
{
	GLint log_length;
	char *log;

	glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
	log = (char*)malloc(log_length);
	glGet__InfoLog(object, log_length, NULL, log);
	fprintf(stderr, "%s", log);
	free(log);
}


GLuint make_shader(GLenum type, const char *filename)
{
	GLint length;
	GLchar *source = (GLchar*)file_contents(filename, &length);
	GLuint shader;
	GLint shader_ok;

	if (!source)
		return 0;

	shader = glCreateShader(type);
	glShaderSource(shader, 1, (const GLchar**)&source, &length);
	free(source);
	glCompileShader(shader);

	glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_ok);
	if (!shader_ok) {
		fprintf(stderr, "Failed to compile %s:\n", filename);
		show_info_log(shader, glGetShaderiv, glGetShaderInfoLog);
		glDeleteShader(shader);
		return 0;
	}
	return shader;
}

GLuint make_program(GLuint vertex_shader, GLuint fragment_shader)
{
	GLint program_ok;

	GLuint program = glCreateProgram();

	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	glGetProgramiv(program, GL_LINK_STATUS, &program_ok);
	if (!program_ok) {
		fprintf(stderr, "Failed to link shader program:\n");
		show_info_log(program, glGetProgramiv, glGetProgramInfoLog);
		glDeleteProgram(program);
		return 0;
	}
	return program;
}



static int make_flag_program(
							 GLuint *vertex_shader,
							 GLuint *fragment_shader,
							 GLuint *program
							 ) 
{
	*vertex_shader = make_shader(GL_VERTEX_SHADER, "perf.v.glsl");
	if (*vertex_shader == 0)
		return 0;
	*fragment_shader = make_shader(GL_FRAGMENT_SHADER, "perf.f.glsl");
	if (*fragment_shader == 0)
		return 0;

	*program = make_program(*vertex_shader, *fragment_shader);
	if (*program == 0)
		return 0;

	return 1;
}

// static void enact_flag_program(
// 							   GLuint vertex_shader,
// 							   GLuint fragment_shader,
// 							   GLuint program
// 							   ) 
// {
// 	g_resources.flag_program.vertex_shader = vertex_shader;
// 	g_resources.flag_program.fragment_shader = fragment_shader;
// 
// 	g_resources.flag_program.program = program;
// 
// 	g_resources.flag_program.uniforms.p_matrix
// 		= glGetUniformLocation(program, "p_matrix");
// 	g_resources.flag_program.uniforms.mv_matrix
// 		= glGetUniformLocation(program, "mv_matrix");
// 
// 	g_resources.flag_program.attributes.position
// 		= glGetAttribLocation(program, "position");
// 	g_resources.flag_program.attributes.normal
// 		= glGetAttribLocation(program, "normal");
// 	g_resources.flag_program.attributes.diffuse
// 		= glGetAttribLocation(program, "diffuse");
// 	g_resources.flag_program.attributes.shininess
// 		= glGetAttribLocation(program, "shininess");
// 	g_resources.flag_program.attributes.specular
// 		= glGetAttribLocation(program, "specular");
// }


// Update global bounding sphere.
void update_bsph()
{
	point boxmin(1e38, 1e38, 1e38);
	point boxmax(-1e38, -1e38, -1e38);

	point c = g_xforms * g_pmesh->bsphere.center;
	float r = g_pmesh->bsphere.r;
	for (int j = 0; j < 3; j++) {
		boxmin[j] = min(boxmin[j], c[j]-r);
		boxmax[j] = max(boxmax[j], c[j]+r);
	}

	point &gc = global_bsph.center;
	float &gr = global_bsph.r;
	gc = 0.5f * (boxmin + boxmax);
	gr = 0.0f;

	c = g_xforms * g_pmesh->bsphere.center;
	r = g_pmesh->bsphere.r;
	gr = max(gr, dist(c, gc) + r);

	//debug
	gc = g_xforms * g_pmesh->bsphere.center;
	gr = g_pmesh->bsphere.r;


	
}



 void resetview()
 {
 	g_camera.stopspin();
 	g_xforms = xform();
 
 	update_bsph();
 	global_xf = xform::trans(0, 0, VIEW_DIST * global_bsph.r) *
 		xform::trans(-global_bsph.center);
 }




static int make_resources1(void)
{
	GLuint vertex_shader, fragment_shader, program;

	char filename[256];
	sprintf(filename,"org.ply");
	g_pvertex_array = update_mesh(filename,&g_task_mesh);

	if (!make_flag_program(&vertex_shader, &fragment_shader, &program))
		return 0;

	g_prog.vertex_shader = vertex_shader;
	g_prog.fragment_shader = fragment_shader;
	g_prog.program = program;

	g_prog.p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
	g_prog.mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");

	g_prog.position = glGetAttribLocation(g_prog.program, "position");
	g_prog.normal = glGetAttribLocation(g_prog.program, "normal");
	g_prog.diffuse = glGetAttribLocation(g_prog.program, "diffuse");
	g_prog.shininess = glGetAttribLocation(g_prog.program, "shininess");
	g_prog.specular = glGetAttribLocation(g_prog.program, "specular");

	for (int i=0;i<4;i++)
	{
		g_projmat.light_position [0][i] = lit_pos[i];
	}
	
	//update_matrix(&g_resources, 8);
	update_matrix_read ("proj06.dat",&g_projmat);

	return 1;
}


int make_resources(void)
{
	char filename[256];

	//update_matrix_read ("proj08.dat",&g_projmat);


	std::string vetshdname, frgshdname;
	std::string glshaderpath = "GLSL/shaders/";
	if (GlobalAppState::getInstance().s_glShaderPath != "")
		glshaderpath = GlobalAppState::getInstance().s_glShaderPath;
	vetshdname = glshaderpath + "/perf.v.glsl";
	frgshdname = glshaderpath + "/perf.f.glsl";

	mrtprogram = new GLSLProgram(vetshdname, frgshdname);

 	g_prog.program = mrtprogram->getHandle();
 
 	g_prog.p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
 	g_prog.mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");
	g_prog.light_position[0] = glGetUniformLocation(g_prog.program,"light0_position");
	g_prog.light_position[1] = glGetUniformLocation(g_prog.program, "light1_position");
	g_prog.light_position[2] = glGetUniformLocation(g_prog.program, "light2_position");
	g_prog.light_position[3] = glGetUniformLocation(g_prog.program, "light3_position");
	g_prog.light_position[4] = glGetUniformLocation(g_prog.program, "light4_position");
	g_prog.light_position[5] = glGetUniformLocation(g_prog.program, "light5_position");
	g_prog.texture = glGetUniformLocation(g_prog.program,"texture");
 
 	g_prog.position = glGetAttribLocation(g_prog.program, "position");
 	g_prog.normal = glGetAttribLocation(g_prog.program, "normal");
 	g_prog.diffuse = glGetAttribLocation(g_prog.program, "diffuse");
 	g_prog.shininess = glGetAttribLocation(g_prog.program, "shininess");
 	g_prog.specular = glGetAttribLocation(g_prog.program, "specular");
	g_prog.uvcoord = glGetAttribLocation(g_prog.program, "uvcoord");



	//g_projmat.light_position[0] = ;

	for (int i=0;i<4;i++)
	{
		g_projmat.light_position [0][i] = lit_pos[i];
	}

	return 1;
}


void init_gl_state()     
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	// Setup our FBO
	glGenFramebuffersEXT(1, &fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Create the render buffer for depth	
	glGenRenderbuffersEXT(1, &depthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, fbo_width, fbo_height);
	// Attach the depth render buffer to the FBO as it's depth attachment
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthBuffer);


	// Now setup the first texture to render to
  	glGenTextures(1, &img);
  	glBindTexture(GL_TEXTURE_2D, img);
  	glTexImage2D(GL_TEXTURE_2D, 0, 3,  fbo_width, fbo_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, img, 0);


	// Create the render buffer for depth	
// 	glGenRenderbuffersEXT(1, &img);
// 	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, img);
// 	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGB, fbo_width, fbo_height);
//  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, img);


	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(status != GL_FRAMEBUFFER_COMPLETE_EXT)
		exit(1);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	// Unbind the FBO for now


}


void init()     
{
	glShadeModel(GL_SMOOTH);
	glClearColor(0.0f, 0.0f, 0.2f, 0.5f);
	glClearDepth(1.0f);					
	glEnable(GL_DEPTH_TEST);			
	glDepthFunc(GL_LEQUAL);				
	glViewport(0,0,800,600);

	// Setup our FBO
	glGenFramebuffersEXT(1, &fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Create the render buffer for depth	
	glGenRenderbuffersEXT(1, &depthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, g_width, g_height);

	// Now setup the first texture to render to
	glGenTextures(1, &img);
	glBindTexture(GL_TEXTURE_2D, img);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,  g_width, g_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, img, 0);

	// after that setup the second texture to render to
	glGenTextures(1, &img2);
	glBindTexture(GL_TEXTURE_2D, img2);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,  g_width, g_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, img2, 0);

	// Attach the depth render buffer to the FBO as it's depth attachment
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depthBuffer);

	// Define an array which contains the targets we wish to render to...
	GLenum mrt[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
	// ... then inform OpenGL that we wish to render to these two targets
	glDrawBuffers(2,mrt);

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(status != GL_FRAMEBUFFER_COMPLETE_EXT)
		exit(1);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	// Unbind the FBO for now

	mrtprogram = new GLSLProgram("perf.v.glsl","perf.f.glsl");
}

void ShutDown()
{
	glDeleteFramebuffersEXT(1, &fbo);
	glDeleteRenderbuffersEXT(1, &depthBuffer);
	glDeleteTextures(1,&img);
	glDeleteTextures(1,&img2);

	delete mrtprogram;
}

void reshape(int w,int h)			
{
	glViewport( 0, 0, w, h );
	glMatrixMode(GL_PROJECTION);	
	glLoadIdentity();					
	if ( h==0 )							
		gluPerspective(80,(float)w,1.0,5000.0);
	else
		gluPerspective(80,(float)w/(float)h,1.0,5000.0);
	glMatrixMode(GL_MODELVIEW);	
	glLoadIdentity();					
}


void idle(void)
{
	xform tmp_xf = global_xf;
	tmp_xf = global_xf * g_xforms;

	if (g_camera.autospin(tmp_xf))
		glutPostRedisplay();
	else
		usleep(10000);

	g_xforms = inv(global_xf) * tmp_xf;
	update_bsph();

	

	
	//if()


	
}

void RenderVisibilityTriangle(void)
{
	vector<vector<Color>> vec_faceintenlist(g_pmesh->faces.size());

	char filename[256];
	char norname[256];

	glShadeModel(GL_FLAT);
	//glShadeModel(GL_SMOOTH);
	//  	glEnable(GL_POLYGON_OFFSET_FILL); 
	//	glPolygonOffset(1.1f, 4.0f);
	//glDisable(GL_BLEND);


	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);

	for (int view = 0; view <8; view++)
	{
		sprintf(norname,"proj%02d.dat",view);
		sprintf(filename,resultdir.c_str(),norname);
		//printf("%s\n",filename);
		update_matrix_read (norname/*filename*/,&g_projmat);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMultMatrixf(g_projmat.p_matrix);

		glMatrixMode ( GL_MODELVIEW ); 
		glLoadIdentity ();
		glMultMatrixf(g_projmat.mv_matrix);


		// Then render as normal
		// Today's scene is a wonderful multi-coloured spinning cube ;)
		// The destination of the data is controlled by the fragment shader we are using
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

		int numfaces = g_pmesh->faces.size();
		if (numfaces>16777215)
		{
			printf("Unable to generate visibility!\n");
		}

		for (GLuint i=0;i<g_pmesh->faces.size();i++)
		{
			int indtmp0 = g_pmesh->faces[i][0];
			int indtmp1 = g_pmesh->faces[i][1];
			int indtmp2 = g_pmesh->faces[i][2];

			point p0 = g_pmesh->vertices[indtmp0];
			point p1 = g_pmesh->vertices[indtmp1];
			point p2 = g_pmesh->vertices[indtmp2];

			GLubyte m_R = ((i+1)>>16)& 0xff;//i/65535;
			GLubyte m_G = ((i+1)>>8)& 0xff;//(i%65535)/255;
			GLubyte m_B = (i+1) & 0xff;//(i%65535)%255;



			//glColor3ub(m_R, m_G, m_B);
			glBegin (GL_TRIANGLES);

			glColor3ub(m_R, m_G, m_B);
			glVertex3f (p0[0],p0[1],p0[2]);
			glVertex3f (p1[0],p1[1],p1[2]);
			glVertex3f (p2[0],p2[1],p2[2]);

			glEnd();

		}

		GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

		glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

		IplImage *pimgout = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);
		sprintf(norname,"Movie2_Scene4_camera%d_0065.png",view);
		sprintf(filename,resultdir.c_str(),norname);
		printf("%s\n",filename);
		IplImage *pimg = cvLoadImage(norname/*filename*/,1);

		for (int i=0;i<fbo_height;i++)
		{
			for (int j=0;j<fbo_width;j++)
			{
				int nindgl = (fbo_height-i-1)*fbo_width+j;
				GLubyte m_B = poutbuffer[nindgl*3+2];
				GLubyte m_G = poutbuffer[nindgl*3+1];
				GLubyte m_R = poutbuffer[nindgl*3];

				int facind = ((GLuint)m_R)*65536 + ((GLuint)m_G)*256 + (GLuint)m_B - 1;

				int nind = i*fbo_width +j;

				pimgout->imageData[nind*3] = m_B;
				pimgout->imageData[nind*3+1] = m_G;
				pimgout->imageData[nind*3+2] = m_R;

				if (facind<0/*||facind>=numfaces*/)
				{
					continue;
				}

				CvPoint2D32f imgpos;
				imgpos.x = (float)j/(float)fbo_width*(float)g_width;
				imgpos.y = (float)i/(float)fbo_height*(float)g_height;

				Color tmpcol;
				tmpcol = GetImgInten(pimg,imgpos,true);

				tmpcol.grey();
				vec_faceintenlist[facind].push_back(tmpcol);

			}
		}

		sprintf(norname,"%d.png",view);
		sprintf(filename,resultdir.c_str(),norname);


		cvSaveImage(norname/*filename*/,pimgout);
		cvReleaseImage(&pimg);
		cvReleaseImage(&pimgout);

	}

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	//glDisable(GL_POLYGON_OFFSET_FILL); 
	glShadeModel(GL_SMOOTH);


	//filter the intenlist
	vector<Color> vec_filterinten(g_pmesh->faces.size());
	for (int f=0;f<g_pmesh->faces.size();f++)
	{

		if (vec_faceintenlist[f].size())
		{
			Color tmpfacecol(0.0,0.0,0.0);
			int m_tmplen = vec_faceintenlist[f].size();

#ifdef MEDIAN_TEXTURE
			for(int k =0;k<m_tmplen;k++)
			{
				tmpfacecol[0] += vec_faceintenlist[f][k][0];
				tmpfacecol[1] += vec_faceintenlist[f][k][1];
				tmpfacecol[2] += vec_faceintenlist[f][k][2];
			}

			vec_filterinten[f] = tmpfacecol / (float)m_tmplen;
#else
			sort(vec_faceintenlist[f].begin(),vec_faceintenlist[f].end());

			tmpfacecol = vec_faceintenlist[f][m_tmplen/2];
			vec_filterinten[f] = tmpfacecol;

#endif
		}


	}

	g_pmesh->need_adjacentfaces();
	for (int i=0;i<g_pmesh->vertices.size();i++)
	{
		int m_numadj = g_pmesh->adjacentfaces[i].size();
		int m_divfactor = m_numadj;
		Color m_tmpval(0.0f,0.0f,0.0f);

		for (int j=0;j<m_numadj;j++)
		{
			int facind = g_pmesh->adjacentfaces[i][j];

			if (vec_faceintenlist[facind].size())
			{
				m_tmpval= m_tmpval + vec_filterinten[facind];
			}
			else
			{
				m_divfactor -= 1;
			}
		}

		if (m_divfactor)
		{
			m_tmpval = m_tmpval / (float)m_divfactor;
		}
		else
		{
			m_tmpval = 0.0f;
		}

		g_pmesh->colors[i] = m_tmpval;

	}

	g_pmesh->write("D://color.ply");



}





void RenderVisibilityVertex(void)
{
	vector<vector<Color>> vec_vetintenlist(g_pmesh->vertices.size());
	int numvertices = g_pmesh->vertices.size();

	char filename[256];
	char norname[256];


	for (int i=0;i<numvertices;++i)
	{
 		GLubyte m_R = ((i+1)>>16)& 0xff;//i/65535;
 		GLubyte m_G = ((i+1)>>8)& 0xff;//(i%65535)/255;
 		GLubyte m_B = (i+1) & 0xff;//(i%65535)%255;

		g_pmesh->colors[i][0] = m_R/255.0f;
		g_pmesh->colors[i][1] = m_G/255.0f;
		g_pmesh->colors[i][2] = m_B/255.0f;

	}



	update_mesh(g_pmesh,&g_task_mesh);

	
	//glShadeModel(GL_SMOOTH);
//  	glEnable(GL_POLYGON_OFFSET_FILL); 
//	glPolygonOffset(1.1f, 4.0f);
	//glDisable(GL_BLEND);


	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);


	mrtprogram->use();

	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;

	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[0][0], 
		g_projmat.light_position[0][1],
		g_projmat.light_position[0][2],
		g_projmat.light_position[0][3]
	);



	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);



	for (int view=0;view<1;view++)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer


		sprintf(norname,"proj%02d.dat",view);
		sprintf(filename,resultdir.c_str(),norname);
		update_matrix_read (filename,&g_projmat);

		glUniformMatrix4fv(
			p_matrix,
			1, GL_FALSE,
			g_projmat.p_matrix
			);
		glUniformMatrix4fv(
			mv_matrix,
			1, GL_FALSE,
			g_projmat.mv_matrix
			);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
		glDrawElements(
			GL_TRIANGLES,
			g_task_mesh.element_count,
			GL_UNSIGNED_INT,
			(void*)0
			);

		GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

		glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

		IplImage *pimgout = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);
		sprintf(norname,"%d.bmp",view);
		sprintf(filename,resultdir.c_str(),norname);
		printf("%s\n",filename);
		IplImage *pimg = cvLoadImage(filename,1);

		for (int i=0;i<fbo_height;i++)
		{
			for (int j=0;j<fbo_width;j++)
			{
				int nindgl = (fbo_height-i-1)*fbo_width+j;
				GLubyte m_B = poutbuffer[nindgl*3+2];
				GLubyte m_G = poutbuffer[nindgl*3+1];
				GLubyte m_R = poutbuffer[nindgl*3];

				int vetind = ((GLuint)m_R)*65536 + ((GLuint)m_G)*256 + (GLuint)m_B - 1;

				int nind = i*fbo_width +j;

				pimgout->imageData[nind*3] = m_B;
				pimgout->imageData[nind*3+1] = m_G;
				pimgout->imageData[nind*3+2] = m_R;

				if (vetind<0||vetind>=numvertices)
				{
					continue;
				}

				CvPoint2D32f imgpos;
				imgpos.x = (float)j/(float)fbo_width*(float)g_width;
				imgpos.y = (float)i/(float)fbo_height*(float)g_height;

				Color tmpcol;
				tmpcol = GetImgInten(pimg,imgpos,true);

				tmpcol.grey();
				vec_vetintenlist[vetind].push_back(tmpcol);

			}
		}

		sprintf(norname,"%d.png",view);
		sprintf(filename,resultdir.c_str(),norname);


		cvSaveImage(filename,pimgout);
		cvReleaseImage(&pimg);
		cvReleaseImage(&pimgout);
	}



	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	mrtprogram->disable();

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	//glShadeModel(GL_SMOOTH);

	//filter the intenlist
	for (int v=0;v<g_pmesh->vertices.size();v++)
	{

		if (vec_vetintenlist[v].size())
		{
			Color tmpfacecol(0.0,0.0,0.0);
			int m_tmplen = vec_vetintenlist[v].size();
			
#ifdef MEDIAN_TEXTURE
			for(int k =0;k<m_tmplen;k++)
			{
				tmpfacecol[0] += vec_vetintenlist[v][k][0];
				tmpfacecol[1] += vec_vetintenlist[v][k][1];
				tmpfacecol[2] += vec_vetintenlist[v][k][2];
			}

			g_pmesh->colors[v] = tmpfacecol / (float)m_tmplen;

#else
			sort(vec_vetintenlist[v].begin(),vec_vetintenlist[v].end());

			tmpfacecol = vec_vetintenlist[v][m_tmplen/2];
			g_pmesh->colors[v] = tmpfacecol;			
#endif
		}
		else
		{
			g_pmesh->colors[v] = 0.0f;
		}


	}



	//for (int i=0;i<numvertices;++i)
	//{
	//	GLubyte m_R = ((i+1)>>16)& 0xff;//i/65535;
	//	GLubyte m_G = ((i+1)>>8)& 0xff;//(i%65535)/255;
	//	GLubyte m_B = (i+1) & 0xff;//(i%65535)%255;

	//	g_pmesh->colors[i][0] = ((float)m_R)/255.0f;
	//	g_pmesh->colors[i][1] = ((float)m_G)/255.0f;
	//	g_pmesh->colors[i][2] = ((float)m_B)/255.0f;

	//}

	g_pmesh->write("D://colorVet.ply");



}

void renderICCV(void)//from the same camera viewpoint
{
	glShadeModel(GL_SMOOTH);

	//glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, g_winwidth, g_winheight);//**control the rendering resolution


	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

#ifndef READ_MATRIX_FROM_CALIBRATION

	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	glPushMatrix();
// 	glMultMatrixd(global_xf);
// 	glMultMatrixd(g_xforms);
	g_extrinsic = global_xf*g_xforms;
	//glMultMatrixd(g_extrinsic);

	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);



#else
	glPushMatrix();

	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);

	//int viewnum = 7;
	//char norname[256];
	//char filename[256];
	//sprintf(norname,"proj%02d.dat",viewnum);
	//sprintf(filename,resultdir,norname);
	//update_matrix_read (filename,&g_projmat);
	update_matrix_read(g_array_proj, g_array_mv, &g_projmat, g_width,g_height);

	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	//glMultMatrixf(g_projmat.mv_matrix);

		//glPushMatrix();
// 	glMultMatrixd(global_xf);
// 	glMultMatrixd(g_xforms);
	//g_extrinsic = global_xf*g_xforms;
	//glMultMatrixf(g_xforms);


 	//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	

	// 	glMultMatrixd(global_xf);
	// 	glMultMatrixd(g_xforms);


#endif

	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_projmat.mv_matrix
		);

	//fprintf(g_logfile,"%f %f %f %f\n",g_projmat.light_position[1][0],
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]);

	for (int i = 0; i < 3;i++)
	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[i][0], 
		g_projmat.light_position[i][1],
		g_projmat.light_position[i][2],
		g_projmat.light_position[i][3]
	);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);


	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//glPolygonMode(GL_FRONT,GL_LINE);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);



	//for element array buffer
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);


	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	glPopMatrix();
	//glPopAttrib();
	glutSwapBuffers();

}

void render3(void)//automatic rotate
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	update_matrix_read(g_array_proj,g_array_mv,&g_projmat,g_width,g_height);


	//for(int i = 0; i<4; ++i)
	//{
	//	prarray_mv[i*4+1] *= -1;
	//	array_mv[i*4+2] *= -1; 
	//}
	
	global_xf = xform(g_projmat.mv_matrix);
	//end of adding

	//point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	//g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
	glLoadIdentity();
 	//glMultMatrixf(global_xf);
 	//glMultMatrixf(g_xforms);
	//glTranslatef(g_pmesh->bsphere.center[0],g_pmesh->bsphere.center[1],g_pmesh->bsphere.center[2]);
	g_extrinsic = global_xf * rot_only(g_xforms);

	//xform mat_trans = xform::trans(-g_pmesh->bsphere.center[0],-g_pmesh->bsphere.center[1],-g_pmesh->bsphere.center[2]);
	//g_extrinsic = mat_trans*g_extrinsic ;

	//for volker canon rotation axle: 0.395764 -0.083446 -0.914553
	//point rot_axle(0.395764, -0.083446, -0.914553);//old
	//point rot_axle(-0.914553, -0.083446,0.395764);

	//point rot_axle(1, 0,0);
	//xform mat_rot = xform::rot((-45+g_index)/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);

	//glRotatef(45.0f,rot_axle[0],rot_axle[1],rot_axle[2]);
	//glTranslatef(-g_pmesh->bsphere.center[0],-g_pmesh->bsphere.center[1],-g_pmesh->bsphere.center[2]);
	
	//g_extrinsic = mat_rot*g_extrinsic ;
	//xform mat_negtrans = xform::trans(g_pmesh->bsphere.center[0],g_pmesh->bsphere.center[1],g_pmesh->bsphere.center[2]);
	//g_extrinsic = mat_negtrans*g_extrinsic ;

	//glMultMatrixf(g_extrinsic);

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	//for debug
	////canon sequence setting 0
	//g_projmat.light_position[1][0] = 0.569416;
	//g_projmat.light_position[1][1] = -0.068804;
	//g_projmat.light_position[1][2] = 0.819165;
	//g_projmat.light_position[1][3] = 0.000000;

	point litdir(g_projmat.light_position[1][0],
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2]);

	//xform mat_litrot = xform::rot(-(-45+g_index)/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	//point newlitdir = mat_litrot * litdir;
	point newlitdir = litdir;
	

	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic
		//g_projmat.mv_matrix
		);

	glUniform4f(
		g_prog.light_position[0],
		newlitdir[0],//g_projmat.light_position[1][0], 
		newlitdir[1],//g_projmat.light_position[1][1],
		newlitdir[2],//g_projmat.light_position[1][2],
		0.0//g_projmat.light_position[1][3]
		);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	glPopMatrix();
	glutSwapBuffers();

}






void render2(void)//automatic rotate
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	update_matrix_read(g_array_proj, g_array_mv, &g_projmat,  g_width, g_height);



	//for(int i = 0; i<4; ++i)
	//{
	//	prarray_mv[i*4+1] *= -1;
	//	array_mv[i*4+2] *= -1; 
	//}
	
	global_xf = xform(g_projmat.mv_matrix);
	//end of adding

	//point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	//g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
	glLoadIdentity();
 	//glMultMatrixf(global_xf);
 	//glMultMatrixf(g_xforms);
	//glTranslatef(g_pmesh->bsphere.center[0],g_pmesh->bsphere.center[1],g_pmesh->bsphere.center[2]);
	g_extrinsic = global_xf;//*g_xforms;

	xform mat_trans = xform::trans(-g_pmesh->bsphere.center[0],-g_pmesh->bsphere.center[1],-g_pmesh->bsphere.center[2]);

	g_extrinsic = mat_trans*g_extrinsic ;

	//for volker canon rotation axle: 0.395764 -0.083446 -0.914553
	//point rot_axle(0.395764, -0.083446, -0.914553);//old
	//point rot_axle(-0.914553, -0.083446,0.395764);

	//point rot_axle(1, 0,0);//for canon seq
	point rot_axle(0, 1,0);//for gopro seq

	xform mat_rot = xform::rot(g_rotindex/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	//glRotatef(45.0f,rot_axle[0],rot_axle[1],rot_axle[2]);
	//glTranslatef(-g_pmesh->bsphere.center[0],-g_pmesh->bsphere.center[1],-g_pmesh->bsphere.center[2]);
	
	g_extrinsic = mat_rot*g_extrinsic ;

	xform mat_negtrans = xform::trans(g_pmesh->bsphere.center[0],g_pmesh->bsphere.center[1],g_pmesh->bsphere.center[2]);

	g_extrinsic = mat_negtrans*g_extrinsic ;

	//glMultMatrixf(g_extrinsic);

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	//for debug
	////canon sequence setting 0
	//g_projmat.light_position[1][0] = 0.569416;
	//g_projmat.light_position[1][1] = -0.068804;
	//g_projmat.light_position[1][2] = 0.819165;
	//g_projmat.light_position[1][3] = 0.000000;

	point litdir(g_projmat.light_position[1][0],
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2]);

	xform mat_litrot = xform::rot(-(-45+g_index)/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	//point newlitdir = mat_litrot * litdir;
	point newlitdir = litdir;
	

	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic
		//g_projmat.mv_matrix
		);

	glUniform4f(
		g_prog.light_position[0],
		newlitdir[0],//g_projmat.light_position[1][0], 
		newlitdir[1],//g_projmat.light_position[1][1],
		newlitdir[2],//g_projmat.light_position[1][2],
		0.0//g_projmat.light_position[1][3]
		);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);


	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	glPopMatrix();
	glutSwapBuffers();

}


void render1rotate_writeout(int mview,int imgindex)
{
	char filename[256];
	char norname[256];

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	//update_matrix_read(g_array_proj,g_array_mv,&g_projmat);
	//global_xf = xform(g_projmat.mv_matrix);
	//end of adding


	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);

	glPushMatrix();
	g_extrinsic = global_xf*g_xforms;
	

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	
	point scenectr = global_xf * global_bsph.center;
	

	xform mat_trans = xform::trans(-scenectr[0],-scenectr[1],-scenectr[2]);
	g_extrinsic = mat_trans*g_extrinsic ;
	//for volker canon rotation axle: 0.395764 -0.083446 -0.914553
	//point rot_axle(0.395764, -0.083446, -0.914553);//old
	//point rot_axle(-0.914553, -0.083446,0.395764);
	
	//point rot_axle(1, 0,0);//for canon seq
	point rot_axle(0, 1,0);//for gopro seq

	int mrotangle = -g_rotindex;	

	xform mat_rot = xform::rot(mrotangle/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	g_extrinsic = mat_rot*g_extrinsic ;	

	xform mat_negtrans = xform::trans(scenectr[0],scenectr[1],scenectr[2]);
	g_extrinsic = mat_negtrans*g_extrinsic ;
	


	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic		
		);

	
	//point initlit(-0.868953,0.467307,0.162925);
	point initlit(-0.000985, -0.000173, 1.000000);
	point newlit = mat_rot * initlit;

	g_projmat.light_position[1][0] = newlit[0];
	g_projmat.light_position[1][1] = newlit[1];
	g_projmat.light_position[1][2] = newlit[2];
	g_projmat.light_position[1][3] = 0.000000;
		

		
	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[1][0], 
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2],
		g_projmat.light_position[1][3]
	);

	


	//fprintf(g_logfile,"%f %f %f %f\n",
	//	g_projmat.light_position[1][0],
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]);

	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, fbo_width, fbo_height);
	

    glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//glPolygonMode(	GL_FRONT,GL_LINE);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);



	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);


	GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

	IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);

	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			int nind = i*fbo_width+j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}

	//sprintf(filename,"%s_%d.png",g_prefixname,i);
	//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
	//sprintf(filename,g_outputdir.c_str(),mview,imgindex);
	//sprintf(norname,resultdir.c_str(),filename);

	sprintf(filename, "render_%04d.png", imgindex);

	sprintf(norname, resultdir.c_str(), filename);	
	cvSaveImage(norname, pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);


	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	glPopMatrix();
	

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

}

void render1_writeout(int mview,int imgindex)
{
	char filename[256];
	char norname[256];

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	//update_matrix_read(g_array_proj,g_array_mv,&g_projmat);
	//global_xf = xform(g_projmat.mv_matrix);
	//end of adding


	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
// 	glMultMatrixd(global_xf);
// 	glMultMatrixd(g_xforms);
	g_extrinsic = global_xf*g_xforms;
	//glMultMatrixf(g_extrinsic);

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	//for (int i=0;i<6;i++)
	//{
	//	glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	//}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic
		);
		


	vector<CVec4f> vec_lits(6);
	vec_lits[0] = CVec4f( lit_pos[0], lit_pos[1], lit_pos[2], 0 );
	vec_lits[1] = CVec4f( -lit_pos[0], -lit_pos[1], -lit_pos[2], 0 );
	vec_lits[2] = CVec4f( lit_pos[2], 0, -lit_pos[0], 0 );
	vec_lits[3] = CVec4f( -lit_pos[2], 0, lit_pos[0], 0 );
	vec_lits[4] = CVec4f( 0, lit_pos[2], -lit_pos[1], 0 );
	vec_lits[5] = CVec4f( 0, -lit_pos[2], lit_pos[1], 0 );

	for (int i = 0; i < 3;i++)
		glUniform4f(
		g_prog.light_position[i],
		vec_lits[i][0],
		vec_lits[i][1],
		vec_lits[i][2],
		vec_lits[i][3]
	);

	//glUniform4f(
	//	g_prog.light_position[0],
	//	lit_pos[0], 
	//	lit_pos[1], 
	//	lit_pos[2], 
	//	lit_pos[3]
	//);

	
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, fbo_width, fbo_height);
	


	//fprintf(g_logfile,"%f %f %f %f\n",
	//	g_projmat.light_position[1][0],
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]);

	

    glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//glPolygonMode(	GL_FRONT,GL_LINE);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);



	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);


	GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

	IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);

	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			int nind = i*fbo_width+j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}

	//sprintf(filename,"%s_%d.png",g_prefixname,i);
	//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
	sprintf(filename,g_outputdir.c_str(),mview,imgindex);
		
	cvSaveImage(filename, pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);


	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	
	
	glPopAttrib();
	glPopMatrix();

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	

}


void render_PointCloud_writeout(int mview, int imgindex)
{
	char filename[256];
	char norname[256];

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
	// 	glMultMatrixd(global_xf);
	// 	glMultMatrixd(g_xforms);
	g_extrinsic = global_xf*g_xforms;
	//glMultMatrixf(g_extrinsic);

	glGetFloatv(GL_MODELVIEW_MATRIX, g_projmat.mv_matrix);
	glGetFloatv(GL_PROJECTION_MATRIX, g_projmat.p_matrix);

	//for (int i=0;i<6;i++)
	//{
	//	glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	//}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix, g_width, g_height, g_intrinsic);

	//glUniformMatrix4fv(
	//	g_prog.p_matrix,
	//	1, GL_FALSE,
	//	g_projmat.p_matrix
	//	);

	//glUniformMatrix4fv(
	//	g_prog.mv_matrix,
	//	1, GL_FALSE,
	//	g_extrinsic
	//	);


	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixf(g_projmat.p_matrix);
	//glMultMatrixf(vec_glrenderCam[view].p_matrix);


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixf(g_extrinsic);



	//vector<CVec4f> vec_lits(6);
	//vec_lits[0] = CVec4f(lit_pos[0], lit_pos[1], lit_pos[2], 0);
	//vec_lits[1] = CVec4f(-lit_pos[0], -lit_pos[1], -lit_pos[2], 0);
	//vec_lits[2] = CVec4f(lit_pos[2], 0, -lit_pos[0], 0);
	//vec_lits[3] = CVec4f(-lit_pos[2], 0, lit_pos[0], 0);
	//vec_lits[4] = CVec4f(0, lit_pos[2], -lit_pos[1], 0);
	//vec_lits[5] = CVec4f(0, -lit_pos[2], lit_pos[1], 0);

	float *lightdir = lit_pos;

	GLfloat light0_position[] = { lightdir[0], lightdir[1], lightdir[2], 0 };
	GLfloat light1_position[] = { -lightdir[0], -lightdir[1], -lightdir[2], 0 };
	GLfloat light2_position[] = { lightdir[2], 0, -lightdir[0], 0 };
	GLfloat light3_position[] = { -lightdir[2], 0, lightdir[0], 0 };
	GLfloat light4_position[] = { 0, lightdir[2], -lightdir[1], 0 };
	GLfloat light5_position[] = { 0, -lightdir[2], lightdir[1], 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	glLightfv(GL_LIGHT3, GL_POSITION, light3_position);
	glLightfv(GL_LIGHT4, GL_POSITION, light4_position);
	glLightfv(GL_LIGHT5, GL_POSITION, light5_position);


	//for (int i = 0; i < 3; i++)
	//	glUniform4f(
	//	g_prog.light_position[i],
	//	vec_lits[i][0],
	//	vec_lits[i][1],
	//	vec_lits[i][2],
	//	vec_lits[i][3]
	//	);

	//glUniform4f(
	//	g_prog.light_position[0],
	//	lit_pos[0], 
	//	lit_pos[1], 
	//	lit_pos[2], 
	//	lit_pos[3]
	//);


	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, fbo_width, fbo_height);



	//fprintf(g_logfile,"%f %f %f %f\n",
	//	g_projmat.light_position[1][0],
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]);

	//glPointSize(10);

	glBegin(GL_POINTS);
	for (int v = 0; v < g_pmesh->vertices.size(); v++)
	{
		point vetpos = g_pmesh->vertices[v];
		point vetnor = g_pmesh->normals[v];

		glColor3f((vetnor[0] + 1.0f) / 2.0f, (vetnor[1] + 1.0f) / 2.0f, (vetnor[2] + 1.0f) / 2.0f);
		glVertex3f(vetpos[0], vetpos[1], vetpos[2]);
	}
	glEnd();

	//for (int i = 0; i<g_pmesh->faces.size(); i++)
	//{
	//	int indtmp0 = g_pmesh->faces[i][0];
	//	int indtmp1 = g_pmesh->faces[i][1];
	//	int indtmp2 = g_pmesh->faces[i][2];
	//	
	//	point p0 = g_pmesh->vertices[indtmp0];
	//	point p1 = g_pmesh->vertices[indtmp1];
	//	point p2 = g_pmesh->vertices[indtmp2];
	//	GLubyte m_R, m_G, m_B;

	//	glBegin(GL_TRIANGLES);
	//	m_R = 255;
	//	m_G = 255;
	//	m_B = 255;
	//	glColor4ub(m_R, m_G, m_B, 255);
	//	glVertex3f(p0[0], p0[1], p0[2]);

	//	m_R = 255;
	//	m_G = 255;
	//	m_B = 255;
	//	glColor4ub(m_R, m_G, m_B, 255);
	//	glVertex3f(p1[0], p1[1], p1[2]);

	//	m_R = 255;
	//	m_G = 255;
	//	m_B = 255;
	//	glColor4ub(m_R, m_G, m_B, 255);
	//	glVertex3f(p2[0], p2[1], p2[2]);
	//	glEnd();
	//}

	
	GLubyte *poutbuffer = (GLubyte*)malloc(3 * fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height, GL_RGB, GL_UNSIGNED_BYTE, poutbuffer);

	IplImage *pimg = cvCreateImage(cvSize(fbo_width, fbo_height), 8, 3);

	for (int i = 0; i<fbo_height; i++)
	{
		for (int j = 0; j<fbo_width; j++)
		{
			int nind = i*fbo_width + j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width + j;
#else
			int nindgl = (fbo_height - i - 1)*fbo_width + j;
#endif

			pimg->imageData[nind * 3] = poutbuffer[nindgl * 3 + 2];
			pimg->imageData[nind * 3 + 1] = poutbuffer[nindgl * 3 + 1];
			pimg->imageData[nind * 3 + 2] = poutbuffer[nindgl * 3];

		}
	}

	//sprintf(filename,"%s_%d.png",g_prefixname,i);
	//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
	sprintf(filename, g_outputdir.c_str(), mview, imgindex);

	cvSaveImage(filename, pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);
	

	glPopAttrib();
	glPopMatrix();

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);



}


void render1rotate(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	//update_matrix_read(g_array_proj,g_array_mv,&g_projmat);
	//global_xf = xform(g_projmat.mv_matrix);
	//end of adding


	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
// 	glMultMatrixd(global_xf);
// 	glMultMatrixd(g_xforms);
	g_extrinsic = global_xf*g_xforms;
	//glMultMatrixf(g_extrinsic);

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}

	
	//GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	
	point scenectr = global_xf * global_bsph.center;
	

	xform mat_trans = xform::trans(-scenectr[0],-scenectr[1],-scenectr[2]);
	g_extrinsic = mat_trans*g_extrinsic ;
	//for volker canon rotation axle: 0.395764 -0.083446 -0.914553
	//point rot_axle(0.395764, -0.083446, -0.914553);//old
	//point rot_axle(-0.914553, -0.083446,0.395764);
	
	//point rot_axle(1, 0,0);//for canon seq
	point rot_axle(0, 1,0);//for gopro seq

	int mrotangle = -g_rotindex;	

	xform mat_rot = xform::rot(mrotangle/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	g_extrinsic = mat_rot*g_extrinsic ;	

	xform mat_negtrans = xform::trans(scenectr[0],scenectr[1],scenectr[2]);
	g_extrinsic = mat_negtrans*g_extrinsic ;

	



	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic		
		);

	
	
	//point initlit(-0.868953,0.467307,0.162925);
	//point initlit(-0.000985, - 0.000173, 1.000000);
	//point initlit(0.006493, 0.211984, -0.977252);
	point initlit(0.016347, 0.057697, -0.998200);
	point newlit = initlit;// mat_rot * initlit;

	g_projmat.light_position[1][0] = newlit[0];
	g_projmat.light_position[1][1] = newlit[1];
	g_projmat.light_position[1][2] = newlit[2];
	g_projmat.light_position[1][3] = 0.000000;

	 


	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[1][0], 
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2],
		g_projmat.light_position[1][3]
	);

	//fprintf(g_logfile,"%f %f %f %f\n",
	//	g_projmat.light_position[1][0],
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]);

    glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//glPolygonMode(	GL_FRONT,GL_LINE);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);



	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);

	glPopMatrix();
	glutSwapBuffers();

}


void render1(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_prog.program);

	//added by chenglei for iccv video
	//update_matrix_read(g_array_proj,g_array_mv,&g_projmat);
	//global_xf = xform(g_projmat.mv_matrix);
	//end of adding


	point m_virtualctr(0, 0, VIEW_DIST * global_bsph.r);
 	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	//g_camera.setupGL(m_virtualctr, global_bsph.r);


	glPushMatrix();
// 	glMultMatrixd(global_xf);
// 	glMultMatrixd(g_xforms);
	g_extrinsic = global_xf*g_xforms;
	//glMultMatrixf(g_extrinsic);

 	glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}

	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	glUniformMatrix4fv(
		g_prog.p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);

	glUniformMatrix4fv(
		g_prog.mv_matrix,
		1, GL_FALSE,
		g_extrinsic		
		);

	//float *lightdir = g_projmat.light_position[1];
	//GLfloat light0_position[] = { lightdir[0], lightdir[1], lightdir[2], 0 };
	//GLfloat light1_position[] = { -lightdir[0], -lightdir[1], -lightdir[2], 0 };
	//GLfloat light2_position[] = { lightdir[2], 0, -lightdir[0], 0 };
	//GLfloat light3_position[] = { -lightdir[2], 0, lightdir[0], 0 };
	//GLfloat light4_position[] = { 0, lightdir[2], -lightdir[1], 0 };
	//GLfloat light5_position[] = { 0, -lightdir[2], lightdir[1], 0 };

	//for (int i = 0; i < 6;i++)
	//	glUniform4f(
	//	g_prog.light_position[i],
	//	g_projmat.light_position[i][0],
	//	g_projmat.light_position[i][1],
	//	g_projmat.light_position[i][2],
	//	g_projmat.light_position[i][3]
	//);

	std::cout << g_projmat.light_position[0][0] << g_projmat.light_position[0][1] << g_projmat.light_position[0][2] << std::endl;

	for (int i = 0; i < 3; i++)
		glUniform4f(
		g_prog.light_position[i],
		g_projmat.light_position[i][0],
		g_projmat.light_position[i][1],
		g_projmat.light_position[i][2],
		g_projmat.light_position[i][3]
		);


	//fprintf(g_logfile,"%f %f %f %f\n",
	//	g_projmat.light_position[0][0],
	//	g_projmat.light_position[0][1],
	//	g_projmat.light_position[0][2],
	//	g_projmat.light_position[0][3]);

    glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);

	glEnableVertexAttribArray(g_prog.position);
	glEnableVertexAttribArray(g_prog.normal);
	glEnableVertexAttribArray(g_prog.diffuse);
	glEnableVertexAttribArray(g_prog.shininess);
	glEnableVertexAttribArray(g_prog.specular);
	glEnableVertexAttribArray(g_prog.uvcoord);

	//glPolygonMode(	GL_FRONT,GL_LINE);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		g_prog.position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		g_prog.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		g_prog.diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		g_prog.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		g_prog.specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);



	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	glDisableVertexAttribArray(g_prog.position);
	glDisableVertexAttribArray(g_prog.normal);
	glDisableVertexAttribArray(g_prog.diffuse);
	glDisableVertexAttribArray(g_prog.shininess);
	glDisableVertexAttribArray(g_prog.specular);
	glDisableVertexAttribArray(g_prog.uvcoord);


	//if (GlobalAppState::getInstance().s_b_renderHair){
	//	glLineWidth(2.0f);
	//	for (int i = 0; i < gvec_hairmap.size(); i++)
	//	{
	//		glBegin(GL_LINE_STRIP);
	//		for (int j = 0; j < gvec_hairmap[i].size(); j++){
	//		//for (int j = gvec_hairmap[i].size()-1; j > gvec_hairmap[i].size()-3; j--){
	//			glVertex3f(gvec_hairmap[i][j][0], gvec_hairmap[i][j][1], gvec_hairmap[i][j][2]);  // V0		
	//			glColor3ub(255, 0, 0);
	//		}
	//		glEnd();
	//	}
	//}

	glPopMatrix();
	glutSwapBuffers();

}

void OffScreenRenderICCV2(int mview,int imgindex)
{
	char filename[256];
	char norname[256];
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	glPolygonMode(GL_FRONT,GL_LINE);
	glLineWidth(10);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	mrtprogram->use();

	//g_prog.program = mrtprogram->getHandle();

	//GLint position, normal, diffuse, shininess, specular;
	//position = glGetAttribLocation(g_prog.program, "position");
	//normal = glGetAttribLocation(g_prog.program, "normal");
	//diffuse = glGetAttribLocation(g_prog.program, "diffuse");
	//shininess = glGetAttribLocation(g_prog.program, "shininess");
	//specular = glGetAttribLocation(g_prog.program, "specular");

	// 	mrtprogram->sendUniform("p_matrix",g_projmat.p_matrix,GL_FALSE,4);
	// 	mrtprogram->sendUniform("mv_matrix",g_projmat.mv_matrix,GL_FALSE,4);

	//GLint p_matrix,mv_matrix;
	//p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
	//mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");

	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);
	

	g_extrinsic = global_xf*g_xforms;

 	//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
 	glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);

	//vec litdir = g_camera.light();
	//std::cout<<litdir[0]<<" "<<litdir[1]<<" "<<litdir[2]<<std::endl;

	//update_matrix_read(g_array_proj,g_array_mv,&g_projmat);
	//for (int i=0;i<6;i++)
	//{
	//	glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	//}
	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);


	//for scan magzine
	//g_projmat.light_position[1][0] = 0.705069;
	//g_projmat.light_position[1][1] = -0.709007;
	//g_projmat.light_position[1][2] =  -0.013673 ;
	//g_projmat.light_position[1][3] = 0.000000;


	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;


	glUniformMatrix4fv(
		p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);
	glUniformMatrix4fv(
		mv_matrix,
		1, GL_FALSE,
		g_extrinsic
		);


	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[1][0], 
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2],
		g_projmat.light_position[1][3]
	);



	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);




	//glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClearColor(1.0f, 1.0f, 1.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer


	//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
	//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);
	//memcpy(g_projmat.mv_matrix,g_extrinsic,16*sizeof(float));




	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

	IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);

	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			int nind = i*fbo_width+j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}

	//sprintf(filename,"%s_%d.png",g_prefixname,i);
	//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
	sprintf(filename, g_outputdir.c_str(), mview, imgindex);
	sprintf(norname, resultdir.c_str(), filename);
	


	cvSaveImage(norname,pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);


	

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	mrtprogram->disable();


	


	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);



}

void OffScreenRenderROTATE(int mview,int imgindex)
{
	char filename[256];
	char norname[256];
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	mrtprogram->use();

	vec litdir = g_camera.light();
	//std::cout<<"litdir:"<<litdir[0]<<" "<<litdir[1]<<" "<<litdir[2]<<std::endl;
	//fprintf(g_logfile,"lit:%f %f %f\n",litdir[0],litdir[1],litdir[2]);
	litdir[0] = 0.000608;litdir[1] = 0.000794; litdir[2] = -1.000000;


	update_matrix_read(g_array_proj,g_array_mv,&g_projmat,g_width,g_height);
	global_xf = xform(g_projmat.mv_matrix);

	g_extrinsic = global_xf;
	xform mat_trans = xform::trans(-g_staticpt[0],-g_staticpt[1],-g_staticpt[2]);
	g_extrinsic = mat_trans*g_extrinsic ;
	//for volker canon rotation axle: 0.395764 -0.083446 -0.914553
	//point rot_axle(0.395764, -0.083446, -0.914553);//old
	//point rot_axle(-0.914553, -0.083446,0.395764);
	
	//point rot_axle(1, 0,0);//for canon seq
	point rot_axle(0, 1,0);//for gopro seq

	int mrotangle = g_rotindex;
	mrotangle = -45;

	xform mat_rot = xform::rot(mrotangle/180.0f*M_PI,rot_axle[0],rot_axle[1],rot_axle[2]);
	g_extrinsic = mat_rot*g_extrinsic ;

	//for debug
	////point rot_axle2(0, 1,0);//for gopro seq	
	//point rot_axle2(1,0,0);//for gopro seq	
	//mrotangle = 10;
	//mat_rot = xform::rot(mrotangle/180.0f*M_PI,rot_axle2[0],rot_axle2[1],rot_axle2[2]);
	//g_extrinsic = mat_rot*g_extrinsic ;
	//end of debug

	xform mat_negtrans = xform::trans(g_staticpt[0],g_staticpt[1],g_staticpt[2]);
	g_extrinsic = mat_negtrans*g_extrinsic ;


	//for (int i=0;i<6;i++)
	//{
	//	glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	//}
	//GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

	//fprintf(g_logfile,"lit:%f %f %f %f\n",g_projmat.light_position[1][0],g_projmat.light_position[1][1],g_projmat.light_position[1][2],g_projmat.light_position[1][3]);
	//g_projmat.light_position[1][0] = -0.000439;
	//g_projmat.light_position[1][1] = -0.000898;
	//g_projmat.light_position[1][2] = 1.000000;
	//g_projmat.light_position[1][3] = 0.000000;
	
	//g_projmat.light_position[1][0] = 0.464272;
	//g_projmat.light_position[1][1] = -0.169543;
	//g_projmat.light_position[1][2] = 0.869314;
	//g_projmat.light_position[1][3] = 0.000000;

	////canon sequence setting 0 for volker and miguel canon
	g_projmat.light_position[1][0] = 0.569416;
	g_projmat.light_position[1][1] = -0.068804;
	g_projmat.light_position[1][2] = 0.819165;
	g_projmat.light_position[1][3] = 0.000000;



	////canon sequence setting 1
	//g_projmat.light_position[1][0] = -0.000339;
	//g_projmat.light_position[1][1] = -0.000941;
	//g_projmat.light_position[1][2] = 1.000000;
	//g_projmat.light_position[1][3] =  0.000000;

	//canon sequence setting 2
	//g_projmat.light_position[1][0] = 0.104785 ;
	//g_projmat.light_position[1][1] = -0.017968;
	//g_projmat.light_position[1][2] = 0.994333;
	//g_projmat.light_position[1][3] = 0.000000;




	////gopro sequence 0
	//g_projmat.light_position[1][0] = -0.000195;
	//g_projmat.light_position[1][1] = -0.000981;
	//g_projmat.light_position[1][2] = 1.000000;
	//g_projmat.light_position[1][3] = 0.000000;   

	//gopro sequence 1
	//g_projmat.light_position[1][0] = -0.064254;
	//g_projmat.light_position[1][1] = -0.287380;
	//g_projmat.light_position[1][2] = 0.955659;
	//g_projmat.light_position[1][3] = 0.000000;  

	//gopro sequence 2
	//g_projmat.light_position[1][0] = 0.146512;
	//g_projmat.light_position[1][1] = 0.431115;
	//g_projmat.light_position[1][2] = 0.890322;
	//g_projmat.light_position[1][3] = 0.000000;  

	////gopro seq 3
	//g_projmat.light_position[1][0] = -0.353671;
	//g_projmat.light_position[1][1] = 0.844444;
	//g_projmat.light_position[1][2] = 0.402282;
	//g_projmat.light_position[1][3] = 0.000000;  

	//gopro seq 4
	//g_projmat.light_position[1][0] = -0.452749;
	//g_projmat.light_position[1][1] = -0.775398;
	//g_projmat.light_position[1][2] =  -0.440201;
	//g_projmat.light_position[1][3] = 0.000000;  


	//gopro seq 5 for miguel
	//g_projmat.light_position[1][0] = -0.019129;
	//g_projmat.light_position[1][1] = 0.714151;
	//g_projmat.light_position[1][2] =  0.699731 ;
	//g_projmat.light_position[1][3] = 0.000000;  

	//gopro seq 6 for miguel
	//g_projmat.light_position[1][0] = -0.180744;
	//g_projmat.light_position[1][1] = 0.826851;
	//g_projmat.light_position[1][2] = 0.532588 ;
	//g_projmat.light_position[1][3] = 0.000000; 

	//gopro seq 2
	//g_projmat.light_position[1][0] = 0.286401;
	//g_projmat.light_position[1][1] = 0.710698;
	//g_projmat.light_position[1][2] = 0.642559;
	//g_projmat.light_position[1][3] = 0.000000;



	   

	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;


	glUniformMatrix4fv(
		p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);
	glUniformMatrix4fv(
		mv_matrix,
		1, GL_FALSE,
		//g_projmat.mv_matrix
		g_extrinsic
		);


	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[1][0], 
		g_projmat.light_position[1][1],
		g_projmat.light_position[1][2],
		g_projmat.light_position[1][3]
	);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);



	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular,
		uvcoord = g_prog.uvcoord;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);
	glEnableVertexAttribArray(uvcoord);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);




	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer


	//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
	//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);
	//memcpy(g_projmat.mv_matrix,g_extrinsic,16*sizeof(float));




	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

#ifdef ROTATE_IMAGE
	IplImage *pimg = cvCreateImage(cvSize(fbo_height,fbo_width),8,3);
	
	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			//int nind = i*fbo_width+j;
			int nind = j*fbo_height+fbo_height-1-i;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];
		}
	}

#else
	IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);
	
	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			int nind = i*fbo_width+j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}

#endif

	//sprintf(filename,"%s_%d.png",g_prefixname,i);
	//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
	//sprintf(filename,"rendering//rendering_%d_%04d.png",mview,imgindex);
	sprintf(filename, g_outputdir.c_str(), mview, imgindex);
	sprintf(norname, resultdir.c_str(), filename);


	cvSaveImage(norname,pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);


	

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);
	glDisableVertexAttribArray(uvcoord);

	// disable the shader
	mrtprogram->disable();


	


	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);



}







void OffScreenRenderICCV(int mview,int imgindex)
{
	char filename[256];
	char norname[256];
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	//glPolygonMode(GL_FRONT,GL_LINE);
	//glLineWidth(1);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	mrtprogram->use();

	//g_prog.program = mrtprogram->getHandle();

	//GLint position, normal, diffuse, shininess, specular;
	//position = glGetAttribLocation(g_prog.program, "position");
	//normal = glGetAttribLocation(g_prog.program, "normal");
	//diffuse = glGetAttribLocation(g_prog.program, "diffuse");
	//shininess = glGetAttribLocation(g_prog.program, "shininess");
	//specular = glGetAttribLocation(g_prog.program, "specular");

	// 	mrtprogram->sendUniform("p_matrix",g_projmat.p_matrix,GL_FALSE,4);
	// 	mrtprogram->sendUniform("mv_matrix",g_projmat.mv_matrix,GL_FALSE,4);

	//GLint p_matrix,mv_matrix;
	//p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
	//mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");

	g_camera.setupGL(/*m_virtualctr*/ global_xf * global_bsph.center, global_bsph.r);

	vec litdir = g_camera.light();
	//std::cout<<"litdir:"<<litdir[0]<<" "<<litdir[1]<<" "<<litdir[2]<<std::endl;
	//fprintf(g_logfile,"lit:%f %f %f\n",litdir[0],litdir[1],litdir[2]);
	
	//litdir[0] = 0.000608;litdir[1] = 0.000794; litdir[2] = -1.000000;


	update_matrix_read(g_array_proj, g_array_mv, &g_projmat, g_width, g_height);
	for (int i=0;i<6;i++)
	{
		glGetLightfv(GL_LIGHT0+i, GL_POSITION, g_projmat.light_position[i]);
	}
	GLCamera::GL2CVprojmatrix(g_projmat.p_matrix,g_width,g_height,g_intrinsic);

 

	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;


	glUniformMatrix4fv(
		p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);
	glUniformMatrix4fv(
		mv_matrix,
		1, GL_FALSE,
		g_projmat.mv_matrix
		);


	vector<CVec4f> vec_lits(6);
	vec_lits[0] = CVec4f(lit_pos[0], lit_pos[1], lit_pos[2], 0);
	vec_lits[1] = CVec4f(-lit_pos[0], -lit_pos[1], -lit_pos[2], 0);
	vec_lits[2] = CVec4f(lit_pos[2], 0, -lit_pos[0], 0);
	vec_lits[3] = CVec4f(-lit_pos[2], 0, lit_pos[0], 0);
	vec_lits[4] = CVec4f(0, lit_pos[2], -lit_pos[1], 0);
	vec_lits[5] = CVec4f(0, -lit_pos[2], lit_pos[1], 0);

	for (int i = 0; i < 3; i++)
		glUniform4f(
		g_prog.light_position[i],
		vec_lits[i][0],
		vec_lits[i][1],
		vec_lits[i][2],
		vec_lits[i][3]
		);


	//glUniform4f(
	//	g_prog.light_position[0],
	//	g_projmat.light_position[1][0], 
	//	g_projmat.light_position[1][1],
	//	g_projmat.light_position[1][2],
	//	g_projmat.light_position[1][3]
	//);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, g_task_mesh.texture_buffer);
	glUniform1i(g_prog.texture, 0);



	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular,
		uvcoord = g_prog.uvcoord;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);
	glEnableVertexAttribArray(uvcoord);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);
	glVertexAttribPointer(
		g_prog.uvcoord,
		2, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, uvcoord)
		);




	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	//glClearColor(1.0f, 1.0f, 1.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer


	//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
	//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);
	//memcpy(g_projmat.mv_matrix,g_extrinsic,16*sizeof(float));




	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

	glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

#ifdef ROTATE_IMAGE
	IplImage *pimg = cvCreateImage(cvSize(fbo_height,fbo_width),8,3);
	
	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			//int nind = i*fbo_width+j;
			int nind = j*fbo_height+fbo_height-1-i;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];
		}
	}

#else
	IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);
	
	for (int i=0;i<fbo_height;i++)
	{
		for (int j=0;j<fbo_width;j++)
		{
			int nind = i*fbo_width+j;
#ifdef JUGEN_SIG
			int nindgl = (i)*fbo_width+j;
#else
			int nindgl = (fbo_height-i-1)*fbo_width+j;
#endif
			
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}

#endif

	sprintf(filename,g_outputdir.c_str(),mview,imgindex);
	


	cvSaveImage(filename, pimg);
	cvReleaseImage(&pimg);

	free(poutbuffer);


	

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);
	glDisableVertexAttribArray(uvcoord);

	// disable the shader
	mrtprogram->disable();


	


	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);



}



void OffScreenRender(void)
{
	char filename[256];
	char norname[256];
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,fbo_width,fbo_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	mrtprogram->use();

	//g_prog.program = mrtprogram->getHandle();

	//GLint position, normal, diffuse, shininess, specular;
	//position = glGetAttribLocation(g_prog.program, "position");
	//normal = glGetAttribLocation(g_prog.program, "normal");
	//diffuse = glGetAttribLocation(g_prog.program, "diffuse");
	//shininess = glGetAttribLocation(g_prog.program, "shininess");
	//specular = glGetAttribLocation(g_prog.program, "specular");

	// 	mrtprogram->sendUniform("p_matrix",g_projmat.p_matrix,GL_FALSE,4);
	// 	mrtprogram->sendUniform("mv_matrix",g_projmat.mv_matrix,GL_FALSE,4);

	//GLint p_matrix,mv_matrix;
	//p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
	//mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");


	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;

	glUniform4f(
		g_prog.light_position[0],
		g_projmat.light_position[0][0], 
		g_projmat.light_position[0][1],
		g_projmat.light_position[0][2],
		g_projmat.light_position[0][3]
	);


	
	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);



	for (int i=0;i<12;i++)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer


		sprintf(norname,"proj%02d.dat",i);
		sprintf(filename, resultdir.c_str(), norname);
		update_matrix_read (filename,&g_projmat);

		g_camera.setupGL(global_xf * global_bsph.center, global_bsph.r);
		g_extrinsic = global_xf*g_xforms;


		//glGetFloatv(GL_MODELVIEW_MATRIX,g_projmat.mv_matrix); 
		//glGetFloatv(GL_PROJECTION_MATRIX,g_projmat.p_matrix);
		//memcpy(g_projmat.mv_matrix,g_extrinsic,16*sizeof(float));

		glUniformMatrix4fv(
			p_matrix,
			1, GL_FALSE,
			g_projmat.p_matrix
			);
		glUniformMatrix4fv(
			mv_matrix,
			1, GL_FALSE,
			g_projmat.mv_matrix
			);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
		glDrawElements(
			GL_TRIANGLES,
			g_task_mesh.element_count,
			GL_UNSIGNED_INT,
			(void*)0
			);

		GLubyte *poutbuffer = (GLubyte*)malloc(3*fbo_width*fbo_height*sizeof(GLubyte));

		glReadPixels(0, 0, fbo_width, fbo_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

		IplImage *pimg = cvCreateImage(cvSize(fbo_width,fbo_height),8,3);

		for (int i=0;i<fbo_height;i++)
		{
			for (int j=0;j<fbo_width;j++)
			{
				int nind = i*fbo_width+j;
				int nindgl = (fbo_height-i-1)*fbo_width+j;
				pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
				pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
				pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

			}
		}

		//sprintf(filename,"%s_%d.png",g_prefixname,i);
		//sprintf(filename,"D://RenderTarget//Render//woman_%d_%04d_segmented.png",i,g_index);
		sprintf(filename,"angel_%d.png",i);
		

		cvSaveImage(filename,pimg);
		cvReleaseImage(&pimg);
	}

	

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	mrtprogram->disable();


	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}


void render(void)
{
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,g_width,g_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	mrtprogram->use();

	//g_prog.program = mrtprogram->getHandle();

	//GLint position, normal, diffuse, shininess, specular;
	//position = glGetAttribLocation(g_prog.program, "position");
	//normal = glGetAttribLocation(g_prog.program, "normal");
	//diffuse = glGetAttribLocation(g_prog.program, "diffuse");
	//shininess = glGetAttribLocation(g_prog.program, "shininess");
	//specular = glGetAttribLocation(g_prog.program, "specular");

	// 	mrtprogram->sendUniform("p_matrix",g_projmat.p_matrix,GL_FALSE,4);
	// 	mrtprogram->sendUniform("mv_matrix",g_projmat.mv_matrix,GL_FALSE,4);

	//GLint p_matrix,mv_matrix;
	//p_matrix = glGetUniformLocation(g_prog.program, "p_matrix");
	//mv_matrix = glGetUniformLocation(g_prog.program, "mv_matrix");


	GLint p_matrix, mv_matrix;
	p_matrix = g_prog.p_matrix;
	mv_matrix = g_prog.mv_matrix;

	glUniformMatrix4fv(
		p_matrix,
		1, GL_FALSE,
		g_projmat.p_matrix
		);
	glUniformMatrix4fv(
		mv_matrix,
		1, GL_FALSE,
		g_projmat.mv_matrix
		);

	GLint position = g_prog.position, 
		normal = g_prog.normal, 
		diffuse = g_prog.diffuse, 
		shininess = g_prog.shininess, 
		specular = g_prog.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, g_task_mesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(flag_vertex),
		(void*)offsetof(flag_vertex, specular)
		);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_task_mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		g_task_mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
		);

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	//glutSwapBuffers();

	// disable the shader
	mrtprogram->disable();


	GLubyte *poutbuffer = (GLubyte*)malloc(3*g_width*g_height*sizeof(GLubyte));

	glReadPixels(0, 0, g_width, g_height,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

	IplImage *pimg = cvCreateImage(cvSize(g_width,g_height),8,3);

	for (int i=0;i<g_height;i++)
	{
		for (int j=0;j<g_width;j++)
		{
			int nind = i*g_width+j;
			int nindgl = (g_height-i-1)*g_width+j;
			pimg->imageData[nind*3] = poutbuffer[nindgl*3+2];
			pimg->imageData[nind*3+1] = poutbuffer[nindgl*3+1];
			pimg->imageData[nind*3+2] = poutbuffer[nindgl*3];

		}
	}
	char filename[256];
	sprintf(filename,"D:\\outimg.bmp");
	cvSaveImage(filename,pimg);
	cvReleaseImage(&pimg);


	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);


	// Now we clear the default frame buffer we are going to render to
	glClearColor(0.0f, 0.0f, 0.2f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	// Now we get ready to draw the first cube with the first texture attached
	glLoadIdentity();

	// Now bind the first texture to use it
	glBindTexture(GL_TEXTURE_2D, img);
	glEnable(GL_TEXTURE_2D);

	glTranslatef(-1.2f,0.0f,-2.0f);
	glRotatef(-xrot,1.0f,0.0f,0.0f);
	glRotatef(-yrot,0.0f,1.0f,0.0f);

	glColor4f(1.0f,1.0f,1.0f,1.0f);

	GLfloat m_value = 0.5;

	// This time it's a textured spinning cube!
	// The texture being the scene we just rendered!
	glBegin(GL_QUADS);
	// Front Face
	glNormal3f( 0.0f, 0.0f, 1.0);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-m_value, -m_value,  m_value);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( m_value, -m_value,  m_value);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( m_value,  m_value,  m_value);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-m_value,  m_value,  m_value);
	// Back Face
	glNormal3f( 0.0f, 0.0f,-1.0);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-m_value, -m_value, -m_value);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-m_value,  m_value, -m_value);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( m_value,  m_value, -m_value);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( m_value, -m_value, -m_value);
	// Top Face
	glNormal3f( 0.0f, 1.0, 0.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-m_value,  m_value, -m_value);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-m_value,  m_value,  m_value);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( m_value,  m_value,  m_value);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( m_value,  m_value, -m_value);
	// Bottom Face
	glNormal3f( 0.0f,-1.0, 0.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-m_value, -m_value, -m_value);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( m_value, -m_value, -m_value);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( m_value, -m_value,  m_value);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-m_value, -m_value,  m_value);
	// Right face
	glNormal3f( 1.0, 0.0f, 0.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( m_value, -m_value, -m_value);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( m_value,  m_value, -m_value);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( m_value,  m_value,  m_value);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( m_value, -m_value,  m_value);
	// Left Face
	glNormal3f(-1.0, 0.0f, 0.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-m_value, -m_value, -m_value);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-m_value, -m_value,  m_value);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-m_value,  m_value,  m_value);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-m_value,  m_value, -m_value);
	glEnd();



	xrot+=xspeed;
	yrot+=yspeed;

	glutSwapBuffers ( );
	// Swap The Buffers To Not Be Left With A Clear Screen


	


}

void display(void)   
{
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,g_width,g_height);

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	glLoadIdentity();

	glTranslatef(0.0f,0.0f,-2.0f);
	glRotatef(xrot,1.0f,0.0f,0.0f);
	glRotatef(yrot,0.0f,1.0f,0.0f);

	// enable our MRT shader
	mrtprogram->use();
	{
		glBegin(GL_QUADS);
		// Front Face
		glColor4f(0.0f,1.0f,0.0f,1.0f);
		glVertex3f(-0.5f, -0.5,  0.5);
		glVertex3f( 0.5, -0.5,  0.5);
		glVertex3f( 0.5,  0.5,  0.5);
		glVertex3f(-0.5,  0.5,  0.5);
		// Back Face
		glColor4f(1.0f,0.0f,0.0f,1.0f);
		glVertex3f(-0.5, -0.5, -0.5);
		glVertex3f(-0.5,  0.5, -0.5);
		glVertex3f( 0.5,  0.5, -0.5);
		glVertex3f( 0.5, -0.5, -0.5);
		// Top Face
		glColor4f(0.0f,0.0f,1.0f,1.0f);
		glVertex3f(-0.5,  0.5, -0.5);
		glVertex3f(-0.5,  0.5,  0.5);
		glVertex3f( 0.5,  0.5,  0.5);
		glVertex3f( 0.5,  0.5, -0.5);
		// Bottom Face
		glColor4f(0.0f,1.0f,1.0f,1.0f);
		glVertex3f(-0.5, -0.5, -0.5);
		glVertex3f( 0.5, -0.5, -0.5);
		glVertex3f( 0.5, -0.5,  0.5);
		glVertex3f(-0.5, -0.5,  0.5);
		// Right face
		glColor4f(1.0f,1.0f,0.0f,1.0f);
		glVertex3f( 0.5, -0.5, -0.5);
		glVertex3f( 0.5,  0.5, -0.5);
		glVertex3f( 0.5,  0.5,  0.5);
		glVertex3f( 0.5, -0.5,  0.5);
		// Left Face
		glColor4f(1.0f,1.0f,1.0f,1.0f);
		glVertex3f(-0.5, -0.5, -0.5);
		glVertex3f(-0.5, -0.5,  0.5);
		glVertex3f(-0.5,  0.5,  0.5);
		glVertex3f(-0.5,  0.5, -0.5);
		glEnd();

	}
	// disable the shader
	mrtprogram->disable();

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

	// Now we clear the default frame buffer we are going to render to
	glClearColor(0.0f, 0.0f, 0.2f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	// Now we get ready to draw the first cube with the first texture attached
	glLoadIdentity();

	// Now bind the first texture to use it
	glBindTexture(GL_TEXTURE_2D, img);
	glEnable(GL_TEXTURE_2D);

	glTranslatef(-1.2f,0.0f,-2.0f);
	glRotatef(-xrot,1.0f,0.0f,0.0f);
	glRotatef(-yrot,0.0f,1.0f,0.0f);

	glColor4f(1.0f,1.0f,1.0f,1.0f);

	// This time it's a textured spinning cube!
	// The texture being the scene we just rendered!
	glBegin(GL_QUADS);
	// Front Face
	glNormal3f( 0.0f, 0.0f, 1.0);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5f, -0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5,  0.5,  0.5);
	// Back Face
	glNormal3f( 0.0f, 0.0f,-1.0);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5, -0.5);
	// Top Face
	glNormal3f( 0.0f, 1.0, 0.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5,  0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	// Bottom Face
	glNormal3f( 0.0f,-1.0, 0.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5, -0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5,  0.5);
	// Right face
	glNormal3f( 1.0, 0.0f, 0.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5,  0.5);
	// Left Face
	glNormal3f(-1.0, 0.0f, 0.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glEnd();

	// Now we move on to draw the cube again with the 2nd texture
	glLoadIdentity();

	// Now bind the first texture to use it
	glBindTexture(GL_TEXTURE_2D, img2);
	glEnable(GL_TEXTURE_2D);

	glTranslatef(1.2f,0.0f,-2.0f);
	glRotatef(-xrot,1.0f,0.0f,0.0f);
	glRotatef(-yrot,0.0f,1.0f,0.0f);

	glColor4f(1.0f,1.0f,1.0f,1.0f);

	// This time it's a textured spinning cube!
	// The texture being the scene we just rendered!
	glBegin(GL_QUADS);
	// Front Face
	glNormal3f( 0.0f, 0.0f, 1.0);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5f, -0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5,  0.5,  0.5);
	// Back Face
	glNormal3f( 0.0f, 0.0f,-1.0);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5, -0.5);
	// Top Face
	glNormal3f( 0.0f, 1.0, 0.0f);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5,  0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	// Bottom Face
	glNormal3f( 0.0f,-1.0, 0.0f);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5, -0.5, -0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5,  0.5);
	// Right face
	glNormal3f( 1.0, 0.0f, 0.0f);
	glTexCoord2f(1.0f, 0.0f); glVertex3f( 0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f( 0.5,  0.5, -0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f( 0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 0.0f); glVertex3f( 0.5, -0.5,  0.5);
	// Left Face
	glNormal3f(-1.0, 0.0f, 0.0f);
	glTexCoord2f(0.0f, 0.0f); glVertex3f(-0.5, -0.5, -0.5);
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5, -0.5,  0.5);
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5,  0.5,  0.5);
	glTexCoord2f(0.0f, 1.0f); glVertex3f(-0.5,  0.5, -0.5);
	glEnd();

	glDisable(GL_TEXTURE_2D);

	xrot+=xspeed;
	yrot+=yspeed;

	glutSwapBuffers ( );
	// Swap The Buffers To Not Be Left With A Clear Screen
}

// Handle mouse button and motion events
static unsigned buttonstate = 0;

void mousemotionfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};

	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else if (buttonstate & (1 << 30))
		b = Mouse::LIGHT;
	else
		b = physical_to_logical_map[buttonstate & 7];

	xform tmp_xf = global_xf * g_xforms;
	g_camera.mouse(x, y, b,
		tmp_xf * g_pmesh->bsphere.center,
		g_pmesh->bsphere.r,
		tmp_xf);
	g_xforms = inv(global_xf) * tmp_xf;
	update_bsph();

 	if (b != Mouse::NONE)
	{
		if(g_start_mouse_cap)
		{
			int mview = 0;
			OffScreenRenderICCV2(mview,g_mouseindx);
			g_mouseindx++;
		}
 		glutPostRedisplay();
	}
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	static timestamp last_click_time;
	static unsigned last_click_buttonstate = 0;
	static float doubleclick_threshold = 0.25f;

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		buttonstate |= (1 << 30);
	else
		buttonstate &= ~(1 << 30);

	if (state == GLUT_DOWN) {
		buttonstate |= (1 << button);
		if (buttonstate == last_click_buttonstate &&
			now() - last_click_time < doubleclick_threshold) {
				//doubleclick(button, x, y);
				last_click_buttonstate = 0;
		} else {
			last_click_time = now();
			last_click_buttonstate = buttonstate;
		}
	} else {
		buttonstate &= ~(1 << button);
	}

	mousemotionfunc(x, y);
}

void write_camera_para()
{
	char filename[256];
	sprintf(filename,"img%02d.ppm",g_imgindex);
	printf("\nSaving parameters of camera for image %s... ", filename);

	FILE *pfile;

	if (g_imgindex>0)
	{
		pfile = fopen(g_xffilenames.c_str(),"a");
	}
	else
	{
		pfile = fopen(g_xffilenames.c_str(),"w");
	}

	fprintf(pfile,"%d\n",g_imgindex);
	//for (int i=0;i<3;i++)
	{
		int i = 0;
		fprintf(pfile, "%f %f %f\n", g_intrinsic[i * 4] * viewportsclx, g_intrinsic[i * 4 + 1] * viewportsclx, g_intrinsic[i * 4 + 2] * viewportsclx);

		i = 1;
		fprintf(pfile, "%f %f %f\n", g_intrinsic[i * 4] * viewportscly, g_intrinsic[i * 4 + 1] * viewportscly, g_intrinsic[i * 4 + 2] * viewportscly);

		i = 2;
		fprintf(pfile, "%f %f %f\n", g_intrinsic[i * 4], g_intrinsic[i * 4 + 1], g_intrinsic[i * 4 + 2]);
	}
	fprintf(pfile,"0.0 0.0\n");


	//for (int i = 0; i<3; i++)
	//{
	//	fprintf(pfile, "%f %f %f %f\n", g_extrinsic[i], g_extrinsic[4 + i], g_extrinsic[8 + i], g_extrinsic[12 + i]);		
	//}

	
	CMtx4x4f extmat;
	extmat.MakeI();

	for (int i=0;i<3;i++)
	{
		if (i==0)
		{
			fprintf(pfile, "%f %f %f %f\n", g_extrinsic[i], g_extrinsic[4 + i], g_extrinsic[8 + i], g_extrinsic[12 + i]);			
			extmat(i, 0) = g_extrinsic[i];
			extmat(i, 1) = g_extrinsic[4 + i];
			extmat(i, 2) = g_extrinsic[8 + i];
			extmat(i, 3) = g_extrinsic[12 + i];
		}
		else
		{
			fprintf(pfile, "%f %f %f %f\n", -g_extrinsic[i], -g_extrinsic[4 + i], -g_extrinsic[8 + i], -g_extrinsic[12 + i]);			

			extmat(i, 0) = -g_extrinsic[i];
			extmat(i, 1) = -g_extrinsic[4 + i];
			extmat(i, 2) = -g_extrinsic[8 + i];
			extmat(i, 3) = -g_extrinsic[12 + i];
		}
	}

	CVec3f camctr = -extmat.Submat3x3().T()*CVec3f(extmat(0, 3), extmat(1, 3), extmat(2, 3));
	float tmpdist = (camctr - CVec3f(g_meshctr[0], g_meshctr[1], g_meshctr[2])).Magnitude();
	std::cout << "cam "<<g_imgindex<<" = "<< tmpdist << std::endl;


	//for (int i=0;i<3;i++)
	//{
	//	if (i==2)
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",-g_extrinsic[i],-g_extrinsic[4+i],-g_extrinsic[8+i],-g_extrinsic[12+i]);
	//	}
	//	else
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",g_extrinsic[i],g_extrinsic[4+i],g_extrinsic[8+i],g_extrinsic[12+i]);
	//	}
	//}

	fprintf(pfile,"\n");

	//for (int i=0;i<3;i++)
	//{
	//	if (i==2)
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",-g_extrinsic[i],-g_extrinsic[4+i],-g_extrinsic[8+i],-g_extrinsic[12+i]);
	//	}
	//	else
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",-g_extrinsic[i],-g_extrinsic[4+i],-g_extrinsic[8+i],-g_extrinsic[12+i]);
	//	}
	//}

	//for (int i=0;i<3;i++)
	//{
	//	if (i==2)
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",g_extrinsic[i],g_extrinsic[4+i],g_extrinsic[8+i],g_extrinsic[12+i]);
	//	}
	//	else
	//	{
	//		fprintf(pfile,"%f %f %f %f\n",-g_extrinsic[i],-g_extrinsic[4+i],-g_extrinsic[8+i],-g_extrinsic[12+i]);
	//	}
	//}
	//fprintf(pfile,"\n");

	fclose(pfile);



}
// Save the current image to a PPM file.
// Uses the next available filename matching filenamepattern
void dump_image()
{
	// Find first non-used filename
	const char filenamepattern[] = "img%d.ppm";
// 	int imgnum = 0;
// 	FILE *f;
// 	while (1) {
// 		char filename[1024];
// 		sprintf(filename, filenamepattern, imgnum++);
// 		f = fopen(filename, "rb");
// 		if (!f) {
// 			f = fopen(filename, "wb");
// 			printf("\n\nSaving image %s... ", filename);
// 			fflush(stdout);
// 			break;
// 		}
// 		fclose(f);
// 	}

	char filename[1024];
	sprintf(filename, filenamepattern, g_imgindex++);
	FILE *f = fopen(filename, "wb");
	printf("\n\nSaving image %s... ", filename);


	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	char *buf = new char[width*height*3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n#\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;

	printf("Done.\n\n");
}


void dump_imagePNG()
{
	// Find first non-used filename
	char filenamepattern[256];
	sprintf(filenamepattern, "view%04d.png", g_imgindex);

	printf("\n\nSaving image %s... ", filenamepattern);


	// Read pixels
	GLint V[4];
	glGetIntegerv(GL_VIEWPORT, V);
	GLint width = V[2], height = V[3];
	cv::Mat recordimg = cv::Mat(height, width, CV_8UC3);

	char *buf = new char[width*height * 3];
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height / 2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	memcpy(recordimg.data, buf, width*height * 3 * sizeof(BYTE));

	// Write out file
	cv::imwrite(filenamepattern, recordimg);
	delete[] buf;

	printf("Done.\n\n");
}


void loadDefaultCameraPos(std::string bspname)
{
	FILE *ptmpfile = fopen(bspname.c_str(), "r");
	fscanf(ptmpfile,"%f %f %f %f\n",&global_bsph.center[0],&global_bsph.center[1],&global_bsph.center[2], &global_bsph.r);
	fscanf(ptmpfile, "%f %f %f %f\n",
		&lit_pos[0],
		&lit_pos[1],
		&lit_pos[2],
		&lit_pos[3]);
	XForm<float> tmpxf;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			float tmpval;
			fscanf(ptmpfile, "%f", &tmpval);
			tmpxf[i + 4 * j] = tmpval;
		}
	}
	global_xf = tmpxf;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			float tmpval;
			fscanf(ptmpfile, "%f", &tmpval);
			tmpxf[i + 4 * j] = tmpval;
		}
	}
	g_xforms = tmpxf;

	fclose(ptmpfile);
}



#define Ctrl (1-'a')
void keyboardfunc(unsigned char key, int x, int y)
{
	FILE *pcamfile;

	XForm<float> tmpxf;

	switch (key) {
		case Ctrl+'r':
			g_imgindex = 0;
		
			break;
		case Ctrl+'b':
			g_start_mouse_cap = true;
		
			break;

			//setting the brdf according to the keyboard input
			//first specular albedo
		case Ctrl+'m':
			g_brdf_sp_scl += 0.1;
			if(g_brdf_sp_scl>255.0f)
				g_brdf_sp_scl = 255.0f;
			std::cout<<"specular alb:"<<g_brdf_sp_scl<<std::endl;
			break;					
		case Ctrl+'n':
			g_brdf_sp_scl /= 2.0;		
			std::cout<<"specular alb:"<<g_brdf_sp_scl<<std::endl;
			break;
		//then shinyness					
		case Ctrl+'u':
			g_brdf_shy_scl += 0.1;
			if(g_brdf_shy_scl>10.0f)
				g_brdf_shy_scl = 10.0f;
			std::cout<<"shinies:"<<g_brdf_shy_scl<<std::endl;
			break;					
		case Ctrl+'y':
			g_brdf_shy_scl /= 2.0;		
			std::cout<<"shinies:"<<g_brdf_shy_scl<<std::endl;
			break;
		//then diffuse albedo color
		case Ctrl+'j':
			g_brdf_diffuse_col += 0.05;
			if(g_brdf_diffuse_col>2*M_PI)
				g_brdf_diffuse_col = 0.0f;
			break;					
		case Ctrl+'h':
			g_brdf_diffuse_col -= 0.05;
			if(g_brdf_diffuse_col<0.0f)
				g_brdf_diffuse_col = 2*M_PI;		
			break;
		case 't':
			g_b_texture = !g_b_texture;		
			break;
			//end of the setting


		case 'I':
			printf("Writing %s\n", g_xffilenames.c_str());
			//g_extrinsic.write(g_xffilenames);
			write_camera_para();
			//dump_image(); 
			//dump_imagePNG();
			render1_writeout(g_viewind, g_imgindex++);
			break;
		case ' ':
			//if(g_index>=g_nummesh)
			//{
			//	printf("Finished!\n");
			//	break;
			//}

			g_index %= g_nummesh;
			g_pmesh = g_vec_mesh[g_index];
			update_mesh(g_vec_mesh[g_index],&g_task_mesh);

			//OffScreenRenderICCV(g_viewnum,g_index+g_initial_index);//the version for siggraph
			//OffScreenRenderICCV2(g_viewnum,g_index+g_initial_index);
			render1_writeout(g_viewind,g_index*g_frameinterval+g_initial_index);

			//glutPostRedisplay();			

			g_index++;
			break;

		case 'P':			
			pcamfile = fopen(GlobalAppState::getInstance().s_renderParams.c_str(), "w");
			fprintf(pcamfile,"%f %f %f %f\n",global_bsph.center[0],global_bsph.center[1],global_bsph.center[2], global_bsph.r);
			fprintf(pcamfile,"%f %f %f %f\n",
				g_projmat.light_position[0][0],
				g_projmat.light_position[0][1],
				g_projmat.light_position[0][2],
				g_projmat.light_position[0][3]);

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					fprintf(pcamfile, "%f ", global_xf[i + 4 * j]);
				}
				fprintf(pcamfile, "\n");
			}

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					fprintf(pcamfile, "%f ", g_xforms[i + 4 * j]);
				}
				fprintf(pcamfile, "\n");
			}
			fclose(pcamfile);

			

			break;


		case 'L':

			pcamfile = fopen(GlobalAppState::getInstance().s_renderParams.c_str(), "r");
			fscanf(pcamfile, "%f %f %f %f\n", &global_bsph.center[0], &global_bsph.center[1], &global_bsph.center[2], &global_bsph.r);
			fscanf(pcamfile, "%f %f %f %f\n",
				&lit_pos[0],
				&lit_pos[1],
				&lit_pos[2],
				&lit_pos[3]);			
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					float tmpval;
					fscanf(pcamfile, "%f", &tmpval);
					tmpxf[i + 4 * j] = tmpval;
				}
			}
			global_xf = tmpxf;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					float tmpval;
					fscanf(pcamfile, "%f", &tmpval);
					tmpxf[i + 4 * j] = tmpval;
				}
			}
			g_xforms = tmpxf;
			fclose(pcamfile);

			glutPostRedisplay();	

			break;

		case '\033': // Esc
		case Ctrl+'q':
		case 'Q':
		case 'q':
			exit(0);
	}
	glutPostRedisplay();

	
}

static void update_writeout(void)
{
	if(g_index>=g_nummesh)
	{
		std::cout<<"writing to images is done!"<<std::endl;
		exit(0);
		//return;
	}

	g_index %= g_nummesh;
	g_pmesh = g_vec_mesh[g_index];
	update_mesh(g_vec_mesh[g_index],&g_task_mesh);

	int newpos = g_index + g_initial_index - g_rotate_init+45;
	int m_tmpdenorm = newpos / 90;
	int m_tmpnorm = newpos % 90;
	if(m_tmpdenorm%2==0)
		g_rotindex = m_tmpnorm-45;
	else
		g_rotindex = 45 - m_tmpnorm;


	if (!g_pointcloud_render)
	{
		if (g_mode > 3)
			OffScreenRenderICCV(g_viewind, g_index + g_initial_index);
		else
			render1_writeout(g_viewind, g_index*g_frameinterval + g_initial_index);
	}
	else
	{
		render_PointCloud_writeout(g_viewind, g_index*g_frameinterval + g_initial_index);
	}

		
	g_index++;

	glutPostRedisplay();	

	//Sleep(40);


}


static void update_writeout_rotate(void)
{
	if (g_index >= GlobalAppState::getInstance().s_numberFrmRendered)
	{
		std::cout << "writing to images is done!" << std::endl;
		exit(-1);
		//return;
	}

	//g_index %= g_nummesh;
	g_pmesh = g_vec_mesh[g_index%g_nummesh];
	update_mesh(g_pmesh, &g_task_mesh);

	int newpos = g_index + g_initial_index - g_rotate_init + 45;
	int m_tmpdenorm = newpos / 90;
	int m_tmpnorm = newpos % 90;
	if (m_tmpdenorm % 2 == 0)
		g_rotindex = m_tmpnorm - 45;
	else
		g_rotindex = 45 - m_tmpnorm;

	if (!g_pointcloud_render)
		render1rotate_writeout(g_viewind, g_index);

	g_index++;

	glutPostRedisplay();

	Sleep(40);


}


static void update_writeout2(void)
{
	
	g_rotindex += 1;

	render1rotate_writeout(g_viewind,g_rotindex);		

	glutPostRedisplay();	

	//Sleep(10);


}

static void update(void)
{
	int milliseconds = glutGet(GLUT_ELAPSED_TIME);
	GLfloat seconds = (GLfloat)milliseconds * (1.0f/1000.0f);

	g_index %= g_nummesh;

	int newpos = g_index + g_initial_index - g_rotate_init+45;
	int m_tmpdenorm = newpos / 90;
	int m_tmpnorm = newpos % 90;
	if(m_tmpdenorm%2==0)
		g_rotindex = m_tmpnorm-45;
	else
		g_rotindex = 45 - m_tmpnorm;


	
	g_pmesh = g_vec_mesh[g_index];

	update_mesh(g_vec_mesh[g_index],&g_task_mesh);

	g_index++;                                 
	glutPostRedisplay();

	//Sleep(40);

}


static void update2(void)
{
	int milliseconds = glutGet(GLUT_ELAPSED_TIME);
	GLfloat seconds = (GLfloat)milliseconds * (1.0f/1000.0f);

	
	g_rotindex +=1;



	
	
	//g_pmesh = g_vec_mesh[g_index];
	//update_mesh(g_vec_mesh[g_index],&g_task_mesh);

	//g_index++;                                 
	glutPostRedisplay();

	//Sleep(8);

}

void init_geometry_buffer()
{
	//Generate the vertex buffer and element buffer
	glGenBuffers(1, &g_task_mesh.vertex_buffer);
	glGenBuffers(1, &g_task_mesh.element_buffer);
	
#ifdef TEXTURE_MAPPING
	glGenTextures(1, &g_task_mesh.texture_buffer); 
#endif
	

}

void visibility_call()
{
	char filename[256];
	sprintf(filename,"C:\\PerfCap\\tatjana_seq1s1\\modelDensified\\model0006.off");
	sprintf(g_prefixname, resultdir.c_str(), filename);
	g_pmesh = TriMesh::read(filename);

	int numvertices = g_pmesh->vertices.size();

	g_pmesh->need_bsphere();
	g_pmesh->need_normals();

	if (g_pmesh->colors.size()!=numvertices)
	{
		g_pmesh->colors.resize(numvertices);
	}

	for (int i=0;i<numvertices;++i)
	{
		Color vec_color(1.0,1.0,1.0);
		g_pmesh->colors[i] = vec_color;
	}

	update_mesh(g_pmesh,&g_task_mesh);

	g_xforms = xform();
	g_xffilenames = xfname("camera.ply");
	resetview();


	RenderVisibilityTriangle();
	//RenderVisibilityVertex();
	//OffScreenRender();

}


void read_in_calibration(std::string calibfilename, int viewind)
{
	CCameraArray cameras = CCameraArray(calibfilename.c_str());
	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			g_array_proj[i][j] = cameras.cameras[viewind].A(i, j);
		}
	}

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			g_array_mv[j * 4 + i] = cameras.cameras[viewind].Rt(i, j);
		}
	}
}

int read_in_calibration(char *filename)
{
	FILE *infile = fopen(filename, "r");
	if (infile == NULL)
	{
		fprintf(stderr, "Can't open '%s'! Exiting...\n", filename);
		exit(-1);
	}

	float dump_proj[3][4];

	fscanf(infile, "P:\n");
	fscanf(infile, "%f %f %f %f\n", &dump_proj[0][0], &dump_proj[0][1], &dump_proj[0][2], &dump_proj[0][3]);
	fscanf(infile, "%f %f %f %f\n", &dump_proj[1][0], &dump_proj[1][1], &dump_proj[1][2], &dump_proj[1][3]);
	fscanf(infile, "%f %f %f %f\n", &dump_proj[2][0], &dump_proj[2][1], &dump_proj[2][2], &dump_proj[2][3]);

// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[0][0], dump_proj[0][1], dump_proj[0][2], dump_proj[0][3]);
// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[1][0], dump_proj[1][1], dump_proj[1][2], dump_proj[1][3]);
// 	fprintf(g_pcalfile,"%f %f %f %f\n", dump_proj[2][0], dump_proj[2][1], dump_proj[2][2], dump_proj[2][3]);
// 	fprintf(g_pcalfile,"\n");



	fscanf(infile, "\nK:\n");
	fscanf(infile, "%f %f %f\n", &g_array_proj[0][0], &g_array_proj[0][1], &g_array_proj[0][2]);
	fscanf(infile, "%f %f %f\n", &g_array_proj[1][0], &g_array_proj[1][1], &g_array_proj[1][2]);
	fscanf(infile, "%f %f %f\n", &g_array_proj[2][0], &g_array_proj[2][1], &g_array_proj[2][2]);

	
// 	for (int j=0;j<3;j++)
// 	{
// 		array_proj[0][j] = array_proj[0][j]*fbo_width/g_width;
// 		array_proj[1][j] = array_proj[1][j]*fbo_height/g_height;
// 
// 	}

	fscanf(infile, "\nM:\n");
	fscanf(infile, "%f %f %f %f\n", &g_array_mv[0], &g_array_mv[4], &g_array_mv[8], &g_array_mv[12]);
	fscanf(infile, "%f %f %f %f\n", &g_array_mv[1], &g_array_mv[5], &g_array_mv[9], &g_array_mv[13]);
	fscanf(infile, "%f %f %f %f\n", &g_array_mv[2], &g_array_mv[6], &g_array_mv[10], &g_array_mv[14]);
	fscanf(infile, "%f %f %f %f\n", &g_array_mv[3], &g_array_mv[7], &g_array_mv[11], &g_array_mv[15]);


	fclose(infile);
}

void smoothAnimation(float strength,vector<float> &vec_weight)
{
	
	for (unsigned int f=1; f<g_nummesh-1; f++)
	{
		for (unsigned int v=0; v<g_vec_mesh[f]->vertices.size(); v++)
		{
			float m_tmpstr = strength * vec_weight[v];

			vec p0 = g_vec_mesh[f-1]->vertices[v];
			vec p1 = g_vec_mesh[f]->vertices[v];
			vec p2 = g_vec_mesh[f+1]->vertices[v];
			vec d;
			d[0]= p0[0] - p1[0] * 2.0 + p2[0];
			d[1]= p0[1] - p1[1] * 2.0 + p2[1];
			d[2]= p0[2] - p1[2] * 2.0 + p2[2];

			g_vec_mesh[f]->vertices[v][0] += m_tmpstr * d[0];
			g_vec_mesh[f]->vertices[v][1] += m_tmpstr * d[1];
			g_vec_mesh[f]->vertices[v][2] += m_tmpstr * d[2];			
		}
	}
}


//int build_uv_coord_from_projections(TriMesh *pmesh,vector<CvPoint2D32f> &vec_uvcoord)
//{
//	pmesh->need_normals();
//	if(vec_uvcoord.size()!=pmesh->vertices.size())
//		vec_uvcoord.resize(pmesh->vertices.size());
//
//	char filename[256];
//	sprintf(filename,resultdir,"uv.txt");
//	FILE *pfile = fopen(filename,"w");
//	for(int i=0;i<pmesh->vertices.size();i++)
//	{
//		CVec3f tmppt (pmesh->vertices[i][0],
//			pmesh->vertices[i][1],
//			pmesh->vertices[i][2]);
//		CVec2f imgpos;
//		m_cameras->WorldToImgCoords(tmppt,imgpos,0);
//		CVec3f tmpctr = m_cameras->cameras[0].Center;
//		CVec3f tmpdir = (tmpctr - tmppt).Unit();
//		CVec3f tmpnor ( pmesh->normals[i][0],
//			pmesh->normals[i][1],
//			pmesh->normals[i][2]);
//		CvPoint2D32f uv;
//		//if(tmpnor*tmpdir<=0)
//		//{
//		//	uv.x = 1/(float)g_width;
//		//	uv.y = 1/(float)g_height;		
//		//}else
//		{
//			uv.x = imgpos.x/(float)g_width;
//			uv.y = imgpos.y/(float)g_height;
//		}
//
//		vec_uvcoord[i] = uv;
//		fprintf(pfile,"%f %f\n",uv.x,uv.y);
//	}
//	fclose(pfile);
//
//	return 1;
//}


int test()
{
	TriMesh *pmesh = TriMesh::read("newmichijawwide.ply");
	for(int i=0;i<3;i++)
		subdiv(pmesh);

	//pmesh->need_neighbors();
	//int numvertices = pmesh->vertices.size();
	//std::cout<<numvertices<<std::endl;
	//vector<float> vec_colors(numvertices);
	//for(int time=0;time<2;time++)
	//{
	//	for(int i=0;i<numvertices;i++)
	//	{
	//		vector<float> vectmp_val;
	//		for(int j=0;j<pmesh->neighbors[i].size();j++)
	//		{
	//			int tmpind = pmesh->neighbors[i][j];
	//			vectmp_val.push_back(pmesh->colors[tmpind][0]);		
	//		}
	//		sort(vectmp_val.begin(),vectmp_val.end());

	//		int length = vectmp_val.size();
	//		vec_colors[i] = vectmp_val[length/2];		
	//	}
	//	for(int i=0;i<numvertices;i++)
	//	{
	//		pmesh->colors[i][0] = vec_colors[i];
	//		pmesh->colors[i][1] = vec_colors[i];
	//		pmesh->colors[i][2] = vec_colors[i];	
	//	}
	//}

	pmesh->write("newmichijawwidedense.ply");

	return 1;



}


void generateChessBoardImage()
{
	int grid = 100;
	int m_tmpwidth = grid * 16;
	int m_tmpheight = grid * 12;
	IplImage *pimg = cvCreateImage(cvSize(m_tmpwidth, m_tmpheight), 8, 1);

	
	int m_tmpsizegrid = grid;
	BYTE *pixels = (BYTE*)malloc( m_tmpwidth*m_tmpheight * sizeof(BYTE));
	for (int i = 0; i<m_tmpwidth; i++)
	{
		for (int j = 0; j<m_tmpheight; j++)
		{
			int m_tmpind = j*m_tmpwidth + i;
			BYTE tmpval = 0;

			int m_x = i / m_tmpsizegrid;
			int m_y = j / m_tmpsizegrid;

			bool b_x = true;
			bool b_y = true;
			if (m_x % 2 == 0)
				b_x = false;
			if (m_y % 2 == 0)
				b_y = false;

			int tmpadd = m_x / 2 + m_y / 2;

			int val = 0;
			if (tmpadd % 2 == 0)
				val = 255;

			if (b_x&&b_y)
			{
				pixels[m_tmpind] = val;
				
			}
			else
			{
				pixels[m_tmpind] = val;
			}

			//pixels[m_tmpind*3] = 255;
			//pixels[m_tmpind*3+1] = 255;
			//pixels[m_tmpind*3+2] = 255;
		}
	}

	for (int i = 0; i < (m_tmpwidth * m_tmpheight); i++)
	{
		pimg->imageData[i] = pixels[i];
	}

	


	cvSaveImage("chessboard.png", pimg);
	cvReleaseImage(&pimg);



}


void printcmd()
{
	FILE *pfile = fopen("run.bat", "w");
	fprintf(pfile, "python camera_intrinsics.py --nsquaresx 7 --nsquaresy 5 --squaresize 30 --chessboard --width 4096 --height 3072 ..//..//output//80006// ");
	int numfr = 899;

	//char stname[256] = " ..//..//calib2//80006//80006_ch_calib2_%03d.png";
	//char stname[256] = " 80006_ch_calib2_%03d.png";
	char stname[256] = " %03d.png";
	
	for( int i = 0; i < numfr; i+=2)
	{

		fprintf(pfile, stname, i);
	
	}
	fclose(pfile);

}


void printcmd2()
{
	FILE *pfile = fopen("run.bat", "w");
	fprintf(pfile, "python camera_extrinsics.py");
	int numfr = 899;

	//char stname[256] = " ..//..//calib2//80006//80006_ch_calib2_%03d.png";
	//char stname[256] = " 80006_ch_calib2_%03d.png";

	//fprintf(pfile, " --intrin1 ..//output//80006//");
	//fprintf(pfile, " --intrin2 ..//output//80007//");
	
	int interval = 10;

	fprintf(pfile, " --img1");
	for (int i = 0; i < numfr; i += interval)
	{

		char stname[256] = " 80006//%03d.png";

		fprintf(pfile, stname, i);

	}

	fprintf(pfile, " --img2");
	for (int i = 0; i < numfr; i += interval)
	{

		char stname[256] = " 80007//%03d.png";

		fprintf(pfile, stname, i);

	}

	fprintf(pfile, " --nsquaresx 7 --nsquaresy 5 --squaresize 30");

	fprintf(pfile, "  ..//output// ");



	fclose(pfile);

}


void contaganentMeshes()
{


}

void createMeshBoard()
{
	TriMesh *pmesh = new TriMesh();
	pmesh->vertices.resize(4);

	pmesh->vertices[0] = point(0, 0, 0);
	pmesh->vertices[1] = point(0, 500, 0);
	pmesh->vertices[2] = point(500, 0, 0);
	pmesh->vertices[3] = point(500, 500, 0);

	pmesh->faces.resize(2);
	pmesh->faces[0].v[0] = 0;
	pmesh->faces[0].v[1] = 1;
	pmesh->faces[0].v[2] = 3;

	pmesh->faces[1].v[0] = 0;
	pmesh->faces[1].v[1] = 3;
	pmesh->faces[1].v[2] = 2;

	pmesh->write("board.ply");

}

void combineMeshes(TriMesh *ptotmesh, const TriMesh *pbdmesh)
{

	int oldnumvertices = ptotmesh->vertices.size();
	
	
	for (int j = 0; j < pbdmesh->vertices.size(); j++)
	{
		ptotmesh->vertices.push_back(pbdmesh->vertices[j]);
	}
	for (int j = 0; j < pbdmesh->faces.size(); j++)
	{
		TriMesh::Face tmpface = pbdmesh->faces[j];
		tmpface.v[0] += oldnumvertices;
		tmpface.v[1] += oldnumvertices;
		tmpface.v[2] += oldnumvertices;
		ptotmesh->faces.push_back(tmpface);
	}
}

void synthesizeHairMesh(const vector<CVec3f> &hairline, TriMesh &hairmesh, float hairwidth, int samplerate)
{
	
	int numpts = hairline.size();
	assert(numpts >= 2);

	vector<CVec3f> vec_refdir(numpts);
	vector<CVec3f> vec_othdir(numpts);

	vector<CVec3f> vec_dir(numpts - 1);
	for (int i = 0; i < (numpts - 1); i++)
	{
		CVec3f tmpdir = (hairline[i + 1] - hairline[i]).Unit();
		vec_dir[i] = tmpdir;
	}

	vec_othdir[0] = vec_dir[0];
	vec_othdir[numpts - 1] = vec_dir[numpts - 2];
	for (int i = 1; i < (numpts-1); i++)
	{
		vec_othdir[i] = (vec_dir[i - 1] + vec_dir[i]).Unit();
	}

	vector<CVec3f> vec_tmprefdir(numpts - 2);
	for (int i = 0; i < (numpts - 2); i++)
	{
		CVec3f tmpdir = vec_dir[i].Cross(vec_dir[(i + 1)%(numpts-1)]);
		int next = (i + 2);
		while (tmpdir.Magnitude() == 0 && next<2*numpts)
		{
			tmpdir = vec_dir[i].Cross(vec_dir[next%(numpts-1)]);
			next++;
		}
		if (tmpdir.Magnitude()>0)
			vec_tmprefdir[i] = tmpdir.Unit();
		else{
			vec_tmprefdir[i] = CVec3f(-vec_dir[i][1]- vec_dir[i][2],vec_dir[i][0]-vec_dir[i][2],vec_dir[i][0]+vec_dir[i][1]).Unit();
		}
	}

	//for (int i = 1; i < (numpts - 2); i++)
	//{
	//	if (vec_tmprefdir[0] * vec_tmprefdir[i] < 0)
	//		vec_tmprefdir[i] = -vec_tmprefdir[i];
	//}	
	//vec_refdir[0] = vec_tmprefdir[0];
	//vec_refdir[numpts - 1] = vec_tmprefdir[numpts - 3];
	//for (int i = 1; i < (numpts - 1); i++)
	//{
	//	vec_refdir[i] = vec_tmprefdir[i - 1];
	//}

	CVec3f refdir = CVec3f(-vec_dir[0][1] - vec_dir[0][2], vec_dir[0][0] - vec_dir[0][2], vec_dir[0][0] + vec_dir[0][1]).Unit();
	if (vec_tmprefdir.size() > 0)
		refdir = vec_tmprefdir[0];

	for (int i = 0; i < (numpts); i++)
	{
		vec_refdir[i] = refdir;
	}

	
	hairmesh.vertices.clear();

	for (int i = 0; i < numpts; i++)
	{
		CVec3f tmpxdir = (vec_othdir[i].Cross(vec_refdir[i])).Unit();
		CVec3f tmpydir = (vec_othdir[i].Cross(tmpxdir)).Unit();

		float thikness = hairwidth * (numpts - i - 1) / (numpts - 1);

		for (int j = 0; j < samplerate; j++)
		{
			float x = thikness * cos(j * 2 * M_PI / samplerate);
			float y = thikness * sin(j * 2 * M_PI / samplerate);

			CVec3f tmppos = tmpxdir * x + tmpydir * y + hairline[i];

			hairmesh.vertices.push_back(point(tmppos.x,tmppos.y,tmppos.z));
		}
	}

	hairmesh.faces.clear();
	for (int i = 0; i < (numpts - 1); i++)
	{
		for (int j = 0; j < samplerate; j++)
		{
			TriMesh::Face tmpface;
			tmpface[0] = i*samplerate + j;
			tmpface[1] = i*samplerate + (j+1)%samplerate;
			tmpface[2] = (i + 1)*samplerate + j;
			hairmesh.faces.push_back(tmpface);

			tmpface[0] = i*samplerate + (j + 1) % samplerate;
			tmpface[1] = (i+1)*samplerate + (j + 1) % samplerate;
			tmpface[2] = (i + 1)*samplerate + j;
			hairmesh.faces.push_back(tmpface);
		}
	}


	

}

void growHairCylinderGeometry(TriMesh *pmesh, const vector<vector<CVec3f> > &vec_hairlines, float hairwidth, int samplerate)
{

	int numsegs = vec_hairlines.size();
	for (int i = 0; i < numsegs; i++)
	{
		TriMesh hairmesh;
		synthesizeHairMesh(vec_hairlines[i], hairmesh, hairwidth, samplerate);
		combineMeshes(pmesh, &hairmesh);
	}
}

int main(int argc, char* argv[])
{	
	//createMeshBoard();
	//return 0;

	if (argc < 2)
	{
		std::cout << "configuration file missing!\n" << std::endl;
		return -1;
	}

	const std::string &configfile = argv[1];
	ParameterFile parameterFileGlobalApp(configfile);
	GlobalAppState::getInstance().readMembers(parameterFileGlobalApp);
	GlobalAppState::getInstance().print();


	g_mode = GlobalAppState::getInstance().s_mode;

	g_rotate_init = 30;
	g_initial_index = GlobalAppState::getInstance().s_dynSeqFrameSt;
	g_end_index = GlobalAppState::getInstance().s_dynSeqFrameEnd;
	g_frameinterval = GlobalAppState::getInstance().s_dynSeqFrameInterval;
	g_viewind = GlobalAppState::getInstance().s_renderViewPoint;
	
	if (g_mode == 3 || g_mode ==5){
		g_nummesh = g_end_index - g_initial_index + 1;
		g_nummesh /= g_frameinterval;
	}
	else
		g_nummesh = 1;

	resultdir = GlobalAppState::getInstance().s_outputResultsDirectory;
	g_outputdir = GlobalAppState::getInstance().s_outputImagePath;
	std::string meshfilenamepath = GlobalAppState::getInstance().s_inputMeshFileDirectory;

	int winsizex = GlobalAppState::getInstance().s_renderWindowSizeX;
	int winsizey = GlobalAppState::getInstance().s_renderWindowSizeY;

	g_width = GlobalAppState::getInstance().s_renderImageSizeX;
	g_height = GlobalAppState::getInstance().s_renderImageSizeY;

	std::string uvtextureimg = GlobalAppState::getInstance().s_UVtexturePath;

	int magfactor = GlobalAppState::getInstance().s_offrenderMagFactor;

	fbo_width = g_width*magfactor;
	fbo_height = g_height * magfactor;
	
	//vector<CvPoint2D32f> vec_uvcoord;

	camfilename = GlobalAppState::getInstance().s_inputCalibFileDirectory;

	if (GlobalAppState::getInstance().s_b_renderHair)
	{
		gvec_hairmap.clear();
		ifstream is(GlobalAppState::getInstance().s_hairMapPath);
		vector<CVec3f> p;
		for (int n; is >> n; gvec_hairmap.push_back(p)) {
			p.resize(n);
			for (auto& x : p) {
				is >> x[0] >> x[1] >> x[2];
			}
			string str;
			getline(is, str); // skip optional numbers
		}
	}

	
	g_logfile = fopen("log.txt","w");


	glutInit(&argc, (char**)argv);
	glutInitDisplayMode ( GLUT_RGB | GLUT_DEPTH| GLUT_DOUBLE );		// Display Mode
	glutInitWindowSize(winsizex, winsizey);	
	glutCreateWindow( "Mesh Rendering" );
	if (g_mode>1)
		glutHideWindow();

	glewInit();
	if (!GLEW_VERSION_2_0) {
		fprintf(stderr, "OpenGL 2.0 not available\n");
		return 1;
	}

	init_gl_state();
	make_resources();

	init_geometry_buffer();



	g_vec_mesh.resize(g_nummesh);

	int initial_index = g_initial_index;
	//int viewnum = g_viewnum;
	
	if (g_mode == 3 || g_mode==5)
	{
		for (int fr = 0; fr < g_nummesh; fr++)
		{
			int tmpindex = fr*g_frameinterval + initial_index;

			char filename[256];
			char prefixnametmp[256];
			
			sprintf(filename, meshfilenamepath.c_str(), tmpindex);

			TriMesh *ptmpmesh = TriMesh::read(filename);
			
			//TriMesh *ptmpmesh = TriMesh::read("C:\\tempsyv1\\YaserICT\\renderMaterial\\00_Final.ply");

			int numvertices = ptmpmesh->vertices.size();

			ptmpmesh->need_bsphere();

			g_meshctr = ptmpmesh->bsphere.center;

			//if(ptmpmesh->normals.size() > 0)
			//	vector<vec> ().swap(ptmpmesh->normals);

			ptmpmesh->need_normals();

			std::cout << ptmpmesh->uvcoord.size() << " vs " << numvertices << std::endl;
			if (ptmpmesh->uvcoord.size() != numvertices)
				ptmpmesh->uvcoord.resize(numvertices, Vec<2, float>(0, 0));

			//if (ptmpmesh->colors.size() != numvertices)
			{
				ptmpmesh->colors.resize(numvertices);
				for (int i = 0; i < numvertices; ++i)
				{
					Color vec_color(149,150,208);
					//Color vec_color(123, 128, 164);
					//Color vec_color(200,0,0);

					float m_scl = 1.0;///255.0f;
					ptmpmesh->colors[i][0] = vec_color[0] * m_scl;
					ptmpmesh->colors[i][1] = vec_color[1] * m_scl;
					ptmpmesh->colors[i][2] = vec_color[2] * m_scl;

					//debug here
					//ptmpmesh->colors[i] = pcheckmesh->colors[i];
				}
			}
			
			g_vec_mesh[fr] = ptmpmesh;
		}
	}
	else
	{
		TriMesh *ptmpmesh = TriMesh::read(meshfilenamepath.c_str());

		if (ptmpmesh->faces.size() == 0)
			g_pointcloud_render = true;
		


		int oldnumvertices = ptmpmesh->vertices.size();
		if (GlobalAppState::getInstance().s_boardMeshDir != "" && !g_pointcloud_render)
		{
			TriMesh *pbdmesh = TriMesh::read(GlobalAppState::getInstance().s_boardMeshDir.c_str());			
			combineMeshes(ptmpmesh, pbdmesh);			
		}


		//grow hair geometry
		if (GlobalAppState::getInstance().s_b_renderHair)
		{
			growHairCylinderGeometry(ptmpmesh, gvec_hairmap,GlobalAppState::getInstance().s_hairThickness,10 );
			ptmpmesh->normals.clear();
			ptmpmesh->colors.clear();
			ptmpmesh->write("meshwithHair.ply");
		}

		int numvertices = ptmpmesh->vertices.size();
		ptmpmesh->need_bsphere();
		g_meshctr = ptmpmesh->bsphere.center;

		//if(ptmpmesh->normals.size() > 0)
		//	vector<vec> ().swap(ptmpmesh->normals);

		ptmpmesh->need_normals();

		std::cout << ptmpmesh->uvcoord.size() << " vs " << numvertices << std::endl;
		if (ptmpmesh->uvcoord.size() != numvertices)
			ptmpmesh->uvcoord.resize(numvertices, Vec<2, float>(0, 0));

		
		if (ptmpmesh->colors.size() != numvertices)
		{
			ptmpmesh->colors.resize(numvertices);
			for (int i = 0; i < oldnumvertices; ++i)
			{
				//Color vec_color(149,150,208);
				//Color vec_color(123, 128, 164);
				Color vec_color(0.7f, 0.7f, 0.7f);
				//Color vec_color(200,0,0);

				float m_scl = 1.0;///255.0f;
				ptmpmesh->colors[i][0] = vec_color[0] * m_scl;
				ptmpmesh->colors[i][1] = vec_color[1] * m_scl;
				ptmpmesh->colors[i][2] = vec_color[2] * m_scl;

				//debug here
				//ptmpmesh->colors[i] = pcheckmesh->colors[i];
			}			
			for (int i = oldnumvertices; i < numvertices; ++i)
			{
				//Color vec_color(149,150,208);
				//Color vec_color(123, 128, 164);
				Color vec_color(0, 0, 0);
				//Color vec_color(200,0,0);

				float m_scl = 1.0;///255.0f;
				ptmpmesh->colors[i][0] = vec_color[0] * m_scl;
				ptmpmesh->colors[i][1] = vec_color[1] * m_scl;
				ptmpmesh->colors[i][2] = vec_color[2] * m_scl;

				//debug here
				//ptmpmesh->colors[i] = pcheckmesh->colors[i];
			}
		}


		g_vec_mesh[0] = ptmpmesh;
	}

	g_pmesh = g_vec_mesh[0];

	int numvertices = g_pmesh->vertices.size();	
	read_in_calibration(camfilename, g_viewind);

	g_staticpt = g_pmesh->bsphere.center;
	//g_refpt[0] = g_array_mv[0]*g_staticpt[0]+g_array_mv[4]*g_staticpt[1]+g_array_mv[8]*g_staticpt[2]+g_array_mv[12];
	//g_refpt[1] = g_array_mv[1]*g_staticpt[0]+g_array_mv[5]*g_staticpt[1]+g_array_mv[9]*g_staticpt[2]+g_array_mv[13];
	//g_refpt[2] = g_array_mv[2]*g_staticpt[0]+g_array_mv[6]*g_staticpt[1]+g_array_mv[10]*g_staticpt[2]+g_array_mv[14];

	
	if (!g_pointcloud_render){

		if (GlobalAppState::getInstance().s_b_renderWithUV)
			update_mesh(g_pmesh, &g_task_mesh, uvtextureimg);
		else
			update_mesh(g_pmesh, &g_task_mesh);
	}

	//OffScreenRender();
	printf("Done!");

	// Loading done!!!!!!!!!!!!!!!!!!!!!!!!
	g_xforms = xform();
	g_xffilenames = xfname("camera.ply");
	resetview();

	std::string bspfilename = GlobalAppState::getInstance().s_renderParams;
	if (g_mode>1)
		loadDefaultCameraPos(bspfilename);

	if (g_mode == 0)
		glutDisplayFunc(&render1);//the up to date rendering route

	if (g_mode == 1)
		glutDisplayFunc(&renderICCV);//with the same camera view



	glutMouseFunc(mousebuttonfunc);
	glutMotionFunc(mousemotionfunc);
	glutKeyboardFunc(keyboardfunc);

	
	//glutPostRedisplay();	

	if (g_mode <2)
		glutIdleFunc(idle);
	else
	{
		if (g_mode == 6)
			glutIdleFunc(update_writeout_rotate);
		else
			glutIdleFunc(update_writeout);
	}

	


	//   	glutReshapeFunc     ( reshape );
	//   	glutKeyboardFunc    ( keyboard );
	//   	glutIdleFunc		( idle );
	glutMainLoop();			// Run the main GLUT loop for rendering


	
	fclose(g_logfile);

	return 0;

}

