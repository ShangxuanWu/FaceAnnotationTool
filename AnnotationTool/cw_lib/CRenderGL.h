#ifndef CRENDERGL_H
#define CRENDERGL_H

#define USE_OPENGL

#ifdef USE_OPENGL

#include <stdlib.h>//should before glew.h
#include <math.h>
#include <stdio.h>

#ifdef WIN32
#include <conio.h>
#endif

#include "GL/glew.h"
#include "glut.h"


#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>

#include "GLSL/GLSLProgram.hpp"

#include "CCameraArray.h"
#include "PlyMeshIO.h" 


struct render_mesh {
	GLuint vertex_buffer, element_buffer;
	GLsizei element_count;
	GLuint texture;
};

struct render_vertex {
	GLfloat position[4];
	GLfloat normal[4];
	GLfloat diffuse[4];
	GLfloat shininess;
	GLubyte specular[4];
};


struct render_uniform_varible
{
	GLfloat p_matrix[16], mv_matrix[16];
	GLfloat light_position[7][4];
};

struct render_program
{
	GLuint vertex_shader, fragment_shader, program;

	GLint position, normal, diffuse, shininess, specular;


	GLint p_matrix,mv_matrix;
	GLint light_position[7];
};


namespace wRender{
	enum Render_Mode { NONE, OFFSCREEN, STATIC, DYNAMIC};
}

class CRenderGL
{
public:

	CRenderGL(int argc, char**argv, wRender::Render_Mode emode = wRender::OFFSCREEN, int width = 1, int height = 1, int scl = 1)
		:m_eMode(emode),m_imgwidth(width),m_imgheight(height),m_scl(scl)
	{
		m_width = width*scl;
		m_height = height*scl;
		switch(emode)
		{
		case wRender::OFFSCREEN:
			initGLUT_hidden(argc,argv);
			break;
		case wRender::STATIC:
		case wRender::DYNAMIC:
			initGLUT(argc,argv);
			break;
		default:
			break;
		}
	}

	~CRenderGL()
	{
		if (m_iTaskMesh.vertex_buffer != 0)
		{
			glBindBuffer(1, m_iTaskMesh.vertex_buffer);
			glDeleteBuffers(1, &m_iTaskMesh.vertex_buffer);
			m_iTaskMesh.vertex_buffer = 0;
		}

		if (m_iTaskMesh.element_buffer != 0)
		{
			glBindBuffer(1, m_iTaskMesh.element_buffer);
			glDeleteBuffers(1, &m_iTaskMesh.element_buffer);
			m_iTaskMesh.element_buffer = 0;
		}

		if (m_nTexture != 0){
			glDeleteTextures(1, &m_nTexture);
			m_nTexture = 0;
		}
		if (m_nTexture2 != 0){
			glDeleteTextures(1, &m_nTexture2);
			m_nTexture2 = 0;
		}
		if (m_nDepthBuffer != 0){
			glDeleteRenderbuffers(1, &m_nDepthBuffer);
			m_nDepthBuffer = 0;
		}
		if (m_nFbo != 0)
		{
			glDeleteFramebuffers(1, &m_nFbo);
			m_nFbo = 0;
		}


		if(m_pMrtProgram)
			delete m_pMrtProgram;

		if(m_cvCameras)
			delete m_cvCameras;	

	};

	void RenderingPrepartory(const std::string &vetShdName, const std::string &fragShdName);
	
	void SendModelForRender(PlyMeshIO &pmesh, int option = 0, CVec3f m_overlay = CVec3f(1.0, 0.0, 0.0))
	{
		m_pMesh = &pmesh;
		pmesh.computeNormals();
				

		if (m_pMesh->m_colors.size()!=m_pMesh->m_vertices.size())
		{
			m_pMesh->m_colors.resize(m_pMesh->m_vertices.size());
			option = 1;
		}

		if (option==1)
		{
			for (unsigned int i=0;i<m_pMesh->m_vertices.size();++i)
			{
				m_pMesh->m_colors[i] = m_overlay;
			}

			m_color_overlay = m_overlay;
		}	

		if (option==2)
		{
			for (unsigned int i=0;i<m_pMesh->m_vertices.size();++i)
			{
				m_pMesh->m_colors[i] = m_overlay - m_pMesh->m_colors[i];
			}

			m_color_overlay = m_overlay;
		}

		update_mesh(m_pMesh, &m_iTaskMesh);
	}
	

	void texturingMesh(const std::string &imgpath, int numview, float downsizefactor, const vector<vector<int> > &vec_visicams, PlyMeshIO &mesh);


	void OffScreenRender(const std::string &imgpath, const std::string &resultdir, int numview, float downsizefactor)
	{
		char filename[256];
		
		vector<cv::Mat> RGBimgs(numview);
		for (int i = 0; i < numview; i++)
		{
			sprintf(filename, imgpath.c_str(), i);
			cv::Mat img = cv::imread(filename, 1);

			float downsizefactorx = downsizefactor;
			float downsizefactory = downsizefactor;

			cv::Size ndsize1(img.size().width * downsizefactorx, img.size().height*downsizefactory);
			cv::Mat dsimg;
			cv::resize(img, dsimg, ndsize1, 0, 0, CV_INTER_CUBIC);

			RGBimgs[i] = dsimg;
		}

		render_to_image_overlay(RGBimgs, resultdir);
	}




	void renderDepthMultiView(PlyMeshIO &mesh, vector<cv::Mat> &retdepths, bool bLocalCoordinate = false);
	int computeVetexVisibilityCamera(PlyMeshIO &mesh, vector<vector<int> > &vec_visicams, int occulusionbd=5, float depthqn = 1.0f);

	int computeVetexVisibilityCameraCuda(PlyMeshIO &mesh, vector<vector<int> > &vec_visicams, int occulusionbd=5, float depthqn = 1.0f);
	

	void setupRenderCam(const char *pszfile )
	{
		m_cvCameras = new CCameraArray(pszfile);
		update_matrix_read (m_cvCameras);
	}

	void setupRenderCam(const CCameraArray *incameras )
	{
		m_cvCameras = new CCameraArray(incameras);
		update_matrix_read (m_cvCameras);	
	}


	CCameraArray *GetCvCameraParam(){return m_cvCameras;}
		
	void render_to_mask(const std::string &resultdir, CVec3i refcol);

	void render_to_depth_map(std::vector<cv::Mat>& depths);

	CVec3f m_color_overlay;
protected:
	void initGLUT_hidden(int argc, char** argv);
	void initGLUT(int argc, char** argv);

	void init_gl_state();

	void init_FBO();

	void init_geometry_buffer();

	void update_mesh(PlyMeshIO * pmesh, render_mesh *pout_mesh);

	int make_resources(const std::string &vetShdName, const std::string &fragShdName);

	void  update_matrix_read (const char * filename,render_uniform_varible * resource);

	void  update_matrix_read (CCamera &m_cameras,render_uniform_varible * resource);

	void  update_matrix_read (CCameraArray *m_glCameras);
	
		
	void render_to_image_overlay(const vector<cv::Mat> &RGBimg, const std::string &resultdir);
	
	
			

private:
	wRender::Render_Mode m_eMode;
	
	int m_width;		// The width of the texture we'll be rendering to
	int m_height;		// The  hight of the texture we'll be rendering to

	int m_imgwidth;
	int m_imgheight;

	int m_scl;
		
	CCameraArray *m_cvCameras;
	
	PlyMeshIO *m_pMesh;	

private:
	GLuint m_nFbo;					// Our handle to the FBO
	GLuint m_nDepthBuffer;			// Our handle to the depth render buffer
	GLuint m_nTexture,m_nTexture2;			// Our handle to a textures
	
	GLSLProgram * m_pMrtProgram;	// The shader we will be using for MRT rendering output
	render_mesh m_iTaskMesh;
	render_uniform_varible m_iInputUniform;
	vector<render_uniform_varible> vec_glrenderCam;//a sequence cam to render

	render_program m_iRenderProg;

	//for save the image of the rendering result
	int m_nImgindex;
	int m_nIndex;
		
};


#endif

#endif
