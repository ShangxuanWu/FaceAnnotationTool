#include "CRenderGL.h"


//#include <cuda_runtime.h>



#ifdef USE_OPENGL

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif



CVec3f obtainColorFromImage(const cv::Mat &img, float x, float y)
{
	int width = img.cols;
	int height = img.rows;

	const int xm = floor(x);
	const float  a = x - xm;
	const int ym = floor(y);
	const float b = y - ym;

	if (xm <= 0 || xm > (width - 2) || ym <= 0 || ym > (height - 2))
		return CVec3f(0,0,0);

	cv::Vec3b val0 = img.at<cv::Vec3b>(ym, xm);
	cv::Vec3b val1 = img.at<cv::Vec3b>(ym, xm+1);
	cv::Vec3b val2 = img.at<cv::Vec3b>(ym+1, xm);
	cv::Vec3b val3 = img.at<cv::Vec3b>(ym+1, xm+1);
	
	CVec3f blendval0,blendval1;

	for (int i = 0; i < 3; i++){
		blendval0[i] = val0[i] * (1.0f - a) + val1[i] * a;
		blendval1[i] = val2[i] * (1.0f - a) + val3[i] * a;
	}

	CVec3f retval = (blendval0*(1.0f - b) + blendval1 * b);
	return retval;
	

}
void CRenderGL::initGLUT_hidden(int argc, char** argv) 
{
	// init GLUT
	int dummy = 1;
	glutInit (&argc, argv);

	// Generate a hidden dummy window
	glutInitWindowSize (1, 1); 
	glutCreateWindow (""); 
	glutHideWindow();



}

void CRenderGL::initGLUT(int argc, char** argv)
{
	glutInit(&argc, (char**)argv);
	glutInitDisplayMode ( GLUT_RGB | GLUT_DEPTH| GLUT_DOUBLE );		// Display Mode
	glutInitWindowSize(m_width,m_height);
	glutCreateWindow( "FrameBuffer Object Example 2 - Press ESC to exit" );



}

void CRenderGL::init_geometry_buffer()
{
	//Generate the vertex buffer and element buffer
	glGenBuffers(1, &m_iTaskMesh.vertex_buffer);
	glGenBuffers(1, &m_iTaskMesh.element_buffer);

}


void CRenderGL::update_mesh(PlyMeshIO * pmesh, render_mesh *pout_mesh)
{	
	int numvertices = pmesh->m_vertices.size();

	render_vertex *vertex_data
		= (render_vertex*) malloc(numvertices * sizeof(render_vertex));

	for (int i=0;i<numvertices;i++)
	{
		vertex_data[i].position[0] = pmesh->m_vertices[i][0];
		vertex_data[i].position[1] = pmesh->m_vertices[i][1];
		vertex_data[i].position[2] = pmesh->m_vertices[i][2];
		vertex_data[i].position[3] = 1.0f;

		vertex_data[i].normal[0] = pmesh->m_normals[i][0];
		vertex_data[i].normal[1] = pmesh->m_normals[i][1];
		vertex_data[i].normal[2] = pmesh->m_normals[i][2];
		vertex_data[i].normal[3] = 1.0f;

		vertex_data[i].diffuse[0] = pmesh->m_colors[i][0];
		vertex_data[i].diffuse[1] = pmesh->m_colors[i][1];
		vertex_data[i].diffuse[2] = pmesh->m_colors[i][2];
		vertex_data[i].diffuse[3] = 1.0f;

		vertex_data[i].shininess   = 0.0f;
		vertex_data[i].specular[0] = 0;
		vertex_data[i].specular[1] = 0;
		vertex_data[i].specular[2] = 0;
		vertex_data[i].specular[3] = 0;
	}



	GLsizei numfaces = pmesh->m_faces.size();

	GLsizei element_count = 3 * numfaces;

	GLuint *element_data
		= (GLuint*) malloc(element_count * sizeof(GLuint));
	GLuint index;

	for (int i=0;i<numfaces;i++)
	{
		element_data[i*3] = pmesh->m_faces[i][0];
		element_data[i*3+1] = pmesh->m_faces[i][1];
		element_data[i*3+2] = pmesh->m_faces[i][2];
	}


	pout_mesh->element_count = element_count;

	glBindBuffer(GL_ARRAY_BUFFER, pout_mesh->vertex_buffer);
	glBufferData(
		GL_ARRAY_BUFFER,
		numvertices * sizeof(render_vertex),
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
	free((void*)vertex_data);

}

void CRenderGL::init_gl_state()     
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	// Setup our FBO
	glGenFramebuffersEXT(1, &m_nFbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_nFbo);

	// Create the render buffer for depth	
	glGenRenderbuffersEXT(1, &m_nDepthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, m_nDepthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, m_width, m_height);

	// Now setup the first texture to render to
	glGenTextures(1, &m_nTexture);
	glBindTexture(GL_TEXTURE_2D, m_nTexture);
	//glTexImage2D(GL_TEXTURE_2D, 0, 3,  m_width, m_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F,  m_width, m_height, 0, GL_RGB, GL_FLOAT, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_nTexture, 0);

	// Attach the depth render buffer to the FBO as it's depth attachment
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_nDepthBuffer);

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(status != GL_FRAMEBUFFER_COMPLETE_EXT)
		exit(1);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	// Unbind the FBO for now

}


void CRenderGL::init_FBO()     
{
	glShadeModel(GL_SMOOTH);
	glClearColor(0.0f, 0.0f, 0.2f, 0.5f);
	glClearDepth(1.0f);					
	glEnable(GL_DEPTH_TEST);			
	glDepthFunc(GL_LEQUAL);				
	glViewport(0,0,800,600);

	// Setup our FBO
	glGenFramebuffersEXT(1, &m_nFbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_nFbo);

	// Create the render buffer for depth	
	glGenRenderbuffersEXT(1, &m_nDepthBuffer);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, m_nDepthBuffer);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, m_width, m_height);

	// Now setup the first texture to render to
	glGenTextures(1, &m_nTexture);
	glBindTexture(GL_TEXTURE_2D, m_nTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,  m_width, m_height, 0, GL_RGBA, GL_FLOAT/*GL_UNSIGNED_BYTE*/, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_nTexture, 0);

	// after that setup the second texture to render to
	glGenTextures(1, &m_nTexture2);
	glBindTexture(GL_TEXTURE_2D, m_nTexture2);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,  m_width, m_height, 0, GL_RGBA, GL_FLOAT/*GL_UNSIGNED_BYTE*/, NULL);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// And attach it to the FBO so we can render to it
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, m_nTexture2, 0);

	// Attach the depth render buffer to the FBO as it's depth attachment
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_nDepthBuffer);

	// Define an array which contains the targets we wish to render to...
	GLenum mrt[] = {GL_COLOR_ATTACHMENT0_EXT, GL_COLOR_ATTACHMENT1_EXT};
	// ... then inform OpenGL that we wish to render to these two targets
	glDrawBuffers(2,mrt);

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(status != GL_FRAMEBUFFER_COMPLETE_EXT)
		exit(1);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);	// Unbind the FBO for now

	m_pMrtProgram = new GLSLProgram("perf.v.glsl","perf.f.glsl");
}





int CRenderGL::make_resources(const std::string &vetShdName, const std::string &fragShdName)
{
	m_pMrtProgram = new GLSLProgram(vetShdName,fragShdName);

	m_iRenderProg.program = m_pMrtProgram->getHandle();

	m_iRenderProg.p_matrix = glGetUniformLocation(m_iRenderProg.program, "p_matrix");
	m_iRenderProg.mv_matrix = glGetUniformLocation(m_iRenderProg.program, "mv_matrix");
	m_iRenderProg.light_position[0] = glGetUniformLocation(m_iRenderProg.program,"light0_position");

	m_iRenderProg.position = glGetAttribLocation(m_iRenderProg.program, "position");
	m_iRenderProg.normal = glGetAttribLocation(m_iRenderProg.program, "normal");
	m_iRenderProg.diffuse = glGetAttribLocation(m_iRenderProg.program, "diffuse");
	m_iRenderProg.shininess = glGetAttribLocation(m_iRenderProg.program, "shininess");
	m_iRenderProg.specular = glGetAttribLocation(m_iRenderProg.program, "specular");


	const float lit_pos[4] = { 0.408248, -0.816497, 0.408248, 0.0 };
 	for (int i=0;i<4;i++)
 	{
 		m_iInputUniform.light_position [0][i] = lit_pos[i];
 	}


	return 1;
}

void  CRenderGL::update_matrix_read (CCameraArray *array_glCameras)
{
	vec_glrenderCam.resize(array_glCameras->cameras.size());
	for(int ncam =0;ncam<array_glCameras->cameras.size();ncam++)
	{
		CCamera &m_glCameras = array_glCameras->cameras[ncam];
		render_uniform_varible * resource = &vec_glrenderCam[ncam];


		float array_proj[3][4];
		for (int i=0;i<3;i++)
		{
			for (int j=0;j<4;j++)
			{
				array_proj[i][j] = m_glCameras.A.El(i,j)*m_scl;
			}
		}
		array_proj[2][2] = 1.0f;

		float array_mv[16];
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++)
			{
				array_mv[j*4+i] = m_glCameras.Rt.El(i,j);
			}
		}


		float nearplane = 0.01;
		float farplane = 10000.0f;


		CV2GLprojmatrix(array_proj, m_width, m_height,nearplane, farplane, resource->p_matrix );

		for (int i=0;i<16;i++)
		{
			resource->mv_matrix[i] = array_mv[i];
		}

	}

}



void  CRenderGL::update_matrix_read (CCamera &m_glCameras,render_uniform_varible * resource)
{

	float array_proj[3][4];
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<4;j++)
		{
			array_proj[i][j] = m_glCameras.A.El(i,j)*m_scl;
		}
	}
	array_proj[2][2] = 1.0f;


	float array_mv[16];
	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			array_mv[j*4+i] = m_glCameras.Rt.El(i,j);
		}
	}


	float m_scale = 3.0;
	float nearplane = 0.01f;// dist_cam2near - m_scale*m_pMesh->bsphere.r;
	float farplane = 5000.0f;// dist_cam2near + m_scale*m_pMesh->bsphere.r;

	CV2GLprojmatrix(array_proj, m_width, m_height,nearplane, farplane, resource->p_matrix );

	for (int i=0;i<16;i++)
	{
		resource->mv_matrix[i] = array_mv[i];
	}

}


void  CRenderGL::update_matrix_read (const char * filename,render_uniform_varible * resource)
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



	float array_proj[3][4];
	fscanf(infile, "\nK:\n");
	fscanf(infile, "%f %f %f\n", &array_proj[0][0], &array_proj[0][1], &array_proj[0][2]);
	fscanf(infile, "%f %f %f\n", &array_proj[1][0], &array_proj[1][1], &array_proj[1][2]);
	fscanf(infile, "%f %f %f\n", &array_proj[2][0], &array_proj[2][1], &array_proj[2][2]);

	//ajust the intrinsic according to the rendering resolution
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<4;j++)
		{
			array_proj[i][j] *= m_scl;
		}
	}
	array_proj[2][2] = 1.0f;


	float array_mv[16];
	fscanf(infile, "\nM:\n");
	fscanf(infile, "%f %f %f %f\n", &array_mv[0], &array_mv[4], &array_mv[8], &array_mv[12]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[1], &array_mv[5], &array_mv[9], &array_mv[13]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[2], &array_mv[6], &array_mv[10], &array_mv[14]);
	fscanf(infile, "%f %f %f %f\n", &array_mv[3], &array_mv[7], &array_mv[11], &array_mv[15]);

	fclose(infile);

	float m_scale = 3.0;
	float nearplane = 0.1f;// dist_cam2near - m_scale*m_pMesh->bsphere.r;
	float farplane = 5000.0f;// dist_cam2near + m_scale*m_pMesh->bsphere.r;

	CV2GLprojmatrix(array_proj, m_width, m_height,nearplane, farplane, resource->p_matrix );
	
	for (int i=0;i<16;i++)
	{
		resource->mv_matrix[i] = array_mv[i];
	}
}





void CRenderGL::RenderingPrepartory(const std::string &vetShdName, const std::string &fragShdName)
{

	glewInit();
	if (!GLEW_VERSION_2_0) {
		fprintf(stderr, "OpenGL 2.0 not available\n");
	}

	init_gl_state();
	make_resources(vetShdName,fragShdName);
	init_geometry_buffer();
}

int CRenderGL::computeVetexVisibilityCamera(PlyMeshIO &mesh, vector<vector<int> > &vec_visicams, int occulusionbd, float depthqn)
{
	//int occulusionbd = 5;
	//float depthqn = 1.0f;
	//float grazethre = 0.3f;

	vector<cv::Mat> depthimgs;
	renderDepthMultiView(mesh, depthimgs);

	int numviews = depthimgs.size();
	int numvertices = mesh.m_vertices.size();

	if (vec_visicams.size() != numvertices)
		vec_visicams.resize(numvertices);

//#pragma omp parallel for
	for (int v = 0; v < numvertices; v++)
	{
		const CVec3f &tmpvet = mesh.m_vertices[v];

		for (int c = 0; c < numviews; c++)
		{
			const CCamera &cam = m_cvCameras->cameras[c];
			CVec2f tmpimgpos;

			float distval = (cam.Center - tmpvet).Magnitude();

			cam.WorldToImgCoords(tmpvet, tmpimgpos);

			bool bOcculd = false;
			for (int i = -occulusionbd; i <= occulusionbd; i++)
			{
				for (int j = -occulusionbd; j <= occulusionbd; j++)
				{
					int tmpx = tmpimgpos.x + i + 0.5f;
					int tmpy = tmpimgpos.y + j + 0.5f;

					if (tmpx<0 || tmpx>(m_imgwidth - 1) || tmpy<0 || tmpy>(m_imgheight - 1))
					{
						bOcculd = true;
						break;
					}
					else
					{
						cv::Vec3f tmpdepthpt = depthimgs[c].at<cv::Vec3f>(tmpy, tmpx);
						CVec3f tmpsam(tmpdepthpt[0], tmpdepthpt[1], tmpdepthpt[2]);
						float tmpcompval = (tmpsam - cam.Center).Magnitude();

						if (fabs(tmpcompval - distval) > depthqn)
						{
							bOcculd = true;
							break;
						}

					}
				}
			}

			if (!bOcculd)
				vec_visicams[v].push_back(c);
		}

		CVec3f tmpnor = mesh.m_normals[v];
		int tmpnum = vec_visicams[v].size();
		vector<pair<float, int> > vectmp_pair(tmpnum);
		for (int c = 0; c < tmpnum; c++)
		{
			CVec3f tmpdir = m_cvCameras->cameras[c].Center - tmpvet;
			float tmpdot = tmpnor * tmpdir.Unit();
			vectmp_pair[c] = pair<float, int>(tmpdot, vec_visicams[v][c]);
		}
		std::sort(vectmp_pair.begin(), vectmp_pair.end());

		for (int c = 0; c < tmpnum; c++)
		{
			vec_visicams[v][c] = vectmp_pair[tmpnum - 1-c].second;
			//std::cout << vectmp_pair[tmpnum-1-c].second << " ";
		}
		//std::cout << std::endl;
		//getch();
	}


	int maxnum = 0;
	for (int v = 0; v < numvertices; v++)
	{
		if (vec_visicams[v].size()>maxnum)
			maxnum = vec_visicams[v].size();		
	}
	std::cout << "Max number of visible camera " << maxnum << std::endl;

	return maxnum;
}




int write2PLYwithFloatImageInput(const cv::Mat &meshimg, const std::string &plyname)
{
	FILE *file = fopen(plyname.c_str(), "w");
	if (file == NULL)
	{
		std::cout << "CANNOT create file " << plyname << std::endl;
		return 0;
	}
	int numtt = 0;
	int numfaces = 0;

	int width = meshimg.cols;
	int height = meshimg.rows;

	int bufferDimOutput = width * height;
	vector<bool> p_mask(bufferDimOutput, true);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			cv::Vec3f tmpval = meshimg.at<cv::Vec3f>(i, j);
			if (tmpval[0] == 0 && tmpval[1] == 0 && tmpval[2] == 0){				
				p_mask[i*width + j] = false;
			}
			else{
				p_mask[i*width + j] = true;
				numtt += 1;
			}
		}
	}	
	std::cout << "valid vertex " << numtt << std::endl;

	for (int j = 0; j < (height - 1); j++)
	{
		for (int i = 0; i < (width - 1); i++)
		{
			int tmpind0 = j * width + i;
			int tmpind1 = j * width + i + 1;
			int tmpind2 = (j + 1) * width + i;

			if ((p_mask[tmpind0] ) && (p_mask[tmpind1] ) && (p_mask[tmpind2] ))
				numfaces += 1;

			//fprintf(file, "3 %d %d %d\n", tmpind0, tmpind2, tmpind1);

		}
	}

	for (int j = 1; j < height; j++)
	{
		for (int i = 1; i < width; i++)
		{
			int tmpind0 = j * width + i;
			int tmpind1 = j * width + i - 1;
			int tmpind2 = (j - 1) * width + i;

			if ((p_mask[tmpind0] ) && (p_mask[tmpind1] ) && (p_mask[tmpind2] ))
				numfaces += 1;

			//fprintf(file, "3 %d %d %d\n", tmpind0, tmpind2, tmpind1);

		}
	}
	std::cout << "valid faces " << numfaces << std::endl;

	//vector<vector<int> > vec_mapping(height, vector<int>(width, -1));
	vector<int> vec_mapping(bufferDimOutput, -1);


	fprintf(file, "ply\n");
	fprintf(file, "format ascii 1.0\n");
	fprintf(file, "element vertex %d\n", numtt);
	fprintf(file, "property float x\n");
	fprintf(file, "property float y\n");
	fprintf(file, "property float z\n");	
	fprintf(file, "element face %d\n", numfaces);
	fprintf(file, "property list uchar int vertex_indices\n");
	fprintf(file, "end_header\n");

	int ptval = 0;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			int tmpind = j*width + i;
			cv::Vec3f tmpval = meshimg.at<cv::Vec3f>(j, i);

			if (p_mask[tmpind] ){
				fprintf(file, "%f %f %f\n", tmpval[0], tmpval[1], tmpval[2]);
				//vec_mapping[j][i] = ptval;
				vec_mapping[tmpind] = ptval;
				ptval += 1;
			}
		}
	}

	//for (int j = 0; j < height; j++)
	//{
	//	for (int i = 0; i < width; i++)
	//	{
	//		int tmpind = j*width + i;
	//		if (p_mask[tmpind] > MASK_THRESHOLD_VALUE){
	//			std::cout << vec_mapping[tmpind] << std::endl;
	//			getch();
	//		}
	//	}
	//}

	for (int j = 0; j < (height - 1); j++)
	{
		for (int i = 0; i < (width - 1); i++)
		{
			int tmpind0 = j * width + i;
			int tmpind1 = j * width + i + 1;
			int tmpind2 = (j + 1) * width + i;

			cv::Vec3f tmpval = meshimg.at<cv::Vec3f>(j, i);

			if ((p_mask[tmpind0]) && (p_mask[tmpind1] ) && (p_mask[tmpind2]))			
			{
				
				fprintf(file, "3 %d %d %d\n", vec_mapping[tmpind0], vec_mapping[tmpind2], vec_mapping[tmpind1]);

			}
		}
	}

	for (int j = 1; j < height; j++)
	{
		for (int i = 1; i < width; i++)
		{
			int tmpind0 = j * width + i;
			int tmpind1 = j * width + i - 1;
			int tmpind2 = (j - 1) * width + i;

			if ((p_mask[tmpind0] ) && (p_mask[tmpind1] ) && (p_mask[tmpind2] ))
			{
				fprintf(file, "3 %d %d %d\n", vec_mapping[tmpind0], vec_mapping[tmpind2], vec_mapping[tmpind1]);
			}
		}
	}
	fclose(file);

	return 1;

}


void CRenderGL::renderDepthMultiView(PlyMeshIO &mesh, vector<cv::Mat> &retdepths, bool bLocalCoordinate)
{
	std::cout << "renderDepthMultiView...";

	int numviews = m_cvCameras->cameras.size();
	if (retdepths.size() != numviews)
		retdepths.resize(numviews);

	int numvertices = mesh.m_vertices.size();

	vector<CVec2f> vec_vet_range(3, CVec2f(10000, -10000));	

	for (int i = 0; i < numvertices; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (mesh.m_vertices[i][j] < vec_vet_range[j][0])
				vec_vet_range[j][0] = mesh.m_vertices[i][j];
			if (mesh.m_vertices[i][j] > vec_vet_range[j][1])
				vec_vet_range[j][1] = mesh.m_vertices[i][j];
		}
	}

	for (int j = 0; j < 3; j++)
	{
		vec_vet_range[j][0] -= 0.01f;
		vec_vet_range[j][1] += 0.01f;
		//std::cout << vec_vet_range[j][0] << " " << vec_vet_range[j][1] << std::endl;
	}

	
	mesh.m_colors.resize(numvertices);
	for (int i = 0; i < numvertices; i++)
	{

		mesh.m_colors[i][0] = (mesh.m_vertices[i][0] - vec_vet_range[0][0]) / (vec_vet_range[0][1] - vec_vet_range[0][0]);
		mesh.m_colors[i][1] = (mesh.m_vertices[i][1] - vec_vet_range[1][0]) / (vec_vet_range[1][1] - vec_vet_range[1][0]);
		mesh.m_colors[i][2] = (mesh.m_vertices[i][2] - vec_vet_range[2][0]) / (vec_vet_range[2][1] - vec_vet_range[2][0]);
	}
	update_mesh(&mesh, &m_iTaskMesh);


	char filename[256];
	char norname[256];
	
	//numofcam = vec_glrenderCam.size();
	
	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_nFbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, m_imgwidth, m_imgheight);//**control the rendering resolution

	
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	m_pMrtProgram->use();

	GLint p_matrix, mv_matrix;
	p_matrix = m_iRenderProg.p_matrix;
	mv_matrix = m_iRenderProg.mv_matrix;

	glUniform4f(
		m_iRenderProg.light_position[0],
		m_iInputUniform.light_position[0][0],
		m_iInputUniform.light_position[0][1],
		m_iInputUniform.light_position[0][2],
		m_iInputUniform.light_position[0][3]
		);

	GLint position = m_iRenderProg.position,
		normal = m_iRenderProg.normal,
		diffuse = m_iRenderProg.diffuse,
		shininess = m_iRenderProg.shininess,
		specular = m_iRenderProg.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_iTaskMesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, specular)
		);

	for (int view = 0; view<numviews; view++)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

		update_matrix_read(m_cvCameras->cameras[view], &m_iInputUniform);

		glUniformMatrix4fv(
			p_matrix,
			1, GL_FALSE,
			m_iInputUniform.p_matrix			
			);
		glUniformMatrix4fv(
			mv_matrix,
			1, GL_FALSE,
			m_iInputUniform.mv_matrix			
			);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iTaskMesh.element_buffer);
		glDrawElements(
			GL_TRIANGLES,
			m_iTaskMesh.element_count,
			GL_UNSIGNED_INT,
			(void*)0
			);
				
		GLfloat *poutbuffer = new GLfloat[3 * m_imgwidth*m_imgheight];

		glReadPixels(0, 0, m_imgwidth, m_imgheight, GL_RGB, GL_FLOAT, poutbuffer);
		
		cv::Mat img = cv::Mat(m_imgheight, m_imgwidth, CV_32FC3);
		
#pragma omp parallel for
		for (int i = 0; i<m_imgheight; i++)
		{
			for (int j = 0; j<m_imgwidth; j++)
			{
				int nind = i*m_imgwidth + j;
				int nindgl = (m_imgheight - i - 1)*m_imgwidth + j;

				float tmpr = poutbuffer[nindgl * 3 + 0];
				float tmpg = poutbuffer[nindgl * 3 + 1];
				float tmpb = poutbuffer[nindgl * 3 + 2];

				if (tmpr == 0 && tmpg == 0 && tmpb == 0)
				{
					img.at<cv::Vec3f>(i, j)[0] = tmpr;
					img.at<cv::Vec3f>(i, j)[1] = tmpg;
					img.at<cv::Vec3f>(i, j)[2] = tmpb;

				}
				else
				{
					CVec3f tmpvet;
					tmpvet[0] = tmpr * (vec_vet_range[0][1] - vec_vet_range[0][0]) + vec_vet_range[0][0];
					tmpvet[1] = tmpg * (vec_vet_range[1][1] - vec_vet_range[1][0]) + vec_vet_range[1][0];
					tmpvet[2] = tmpb * (vec_vet_range[2][1] - vec_vet_range[2][0]) + vec_vet_range[2][0];

					if (bLocalCoordinate)
					{
						CVec3f tmpnewvet = m_cvCameras->cameras[view].Rt.transform_point(tmpvet);
						img.at<cv::Vec3f>(i, j)[0] = tmpnewvet[0];
						img.at<cv::Vec3f>(i, j)[1] = tmpnewvet[1];
						img.at<cv::Vec3f>(i, j)[2] = tmpnewvet[2];
					}
					else
					{
						img.at<cv::Vec3f>(i, j)[0] = tmpvet[0];
						img.at<cv::Vec3f>(i, j)[1] = tmpvet[1];
						img.at<cv::Vec3f>(i, j)[2] = tmpvet[2];
					}					
				}
			}
		}

		retdepths[view] = img;

		//write2PLYwithFloatImageInput(img, "heloo.ply");
		//getch();

		//sprintf(filename, "overlay%04d.exr", view);
		//cv::imwrite(std::string(filename), img);

		delete[] poutbuffer;


	}



	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	m_pMrtProgram->disable();

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	
	std::cout << "Done!"<<std::endl;

}





void CRenderGL::texturingMesh(const std::string &imgpath, int numview, float downsizefactor, const vector<vector<int> > &vec_visicams, PlyMeshIO &mesh)
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

	int numvertices = mesh.m_vertices.size();
	if (mesh.m_colors.size() != numvertices)
		mesh.m_colors.resize(numvertices);

	for (int v = 0; v < numvertices; v++)
	{
		const vector<int> &vec_tmpind = vec_visicams[v];
		const CVec3f &vetpos = mesh.m_vertices[v];

		if (vec_tmpind.size() == 0)
			mesh.m_colors[v] = CVec3f(1.0f, 0, 0);
		else
		{
			CVec3f tmpsum(0, 0, 0);
			float sumwt = 0;
			for (int c = 0; c < vec_tmpind.size(); c++)
			{
				int camind = vec_tmpind[c];
				const CCamera &cam = m_cvCameras->cameras[camind];
				CVec2f imgpos;
				cam.WorldToImgCoords(vetpos, imgpos);

				CVec3f tmpcolor = obtainColorFromImage(RGBimgs[camind], imgpos.x, imgpos.y);
				float tmpwt = 1.0f;
				sumwt += tmpwt;
				tmpsum += tmpcolor*tmpwt;
			}			
			mesh.m_colors[v] = tmpsum / sumwt/255.0f;
		}
	}
}

void CRenderGL::render_to_mask(const std::string &resultdir, CVec3i refcol)
{
	char filename[256];
	char norname[256];

	float m_overratio = 0.5f;


	int numofcam = 1;
	numofcam = vec_glrenderCam.size();


	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_nFbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0, 0, m_imgwidth, m_imgheight);//**control the rendering resolution

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	m_pMrtProgram->use();

	GLint p_matrix, mv_matrix;
	p_matrix = m_iRenderProg.p_matrix;
	mv_matrix = m_iRenderProg.mv_matrix;

	glUniform4f(
		m_iRenderProg.light_position[0],
		m_iInputUniform.light_position[0][0],
		m_iInputUniform.light_position[0][1],
		m_iInputUniform.light_position[0][2],
		m_iInputUniform.light_position[0][3]
		);



	GLint position = m_iRenderProg.position,
		normal = m_iRenderProg.normal,
		diffuse = m_iRenderProg.diffuse,
		shininess = m_iRenderProg.shininess,
		specular = m_iRenderProg.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_iTaskMesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, specular)
		);

	for (int view = 0; view<numofcam; view++)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

		update_matrix_read(m_cvCameras->cameras[view], &m_iInputUniform);

		glUniformMatrix4fv(
			p_matrix,
			1, GL_FALSE,
			m_iInputUniform.p_matrix
			//vec_glrenderCam[view].p_matrix
			);
		glUniformMatrix4fv(
			mv_matrix,
			1, GL_FALSE,
			m_iInputUniform.mv_matrix
			//vec_glrenderCam[view].mv_matrix
			);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iTaskMesh.element_buffer);
		glDrawElements(
			GL_TRIANGLES,
			m_iTaskMesh.element_count,
			GL_UNSIGNED_INT,
			(void*)0
			);

		GLubyte *poutbuffer = new GLubyte[3 * m_imgwidth*m_imgheight*sizeof(GLubyte)];

		glReadPixels(0, 0, m_imgwidth, m_imgheight, GL_RGB, GL_UNSIGNED_BYTE, poutbuffer);

		cv::Mat tmpimg = cv::Mat(m_imgheight, m_imgwidth, CV_8UC1);

#ifdef USEOPENMP
#pragma omp parallel for
#endif	
		for (int i = 0; i<m_imgheight; i++)
		{
			for (int j = 0; j<m_imgwidth; j++)
			{
				int nind = i*m_imgwidth + j;
				int nindgl = (m_imgheight - i - 1)*m_imgwidth + j;

				float tmprat = 0.0f;
				unsigned char nB = poutbuffer[nindgl * 3 + 2];
				unsigned char nG = poutbuffer[nindgl * 3 + 1];
				unsigned char nR = poutbuffer[nindgl * 3];

				//if (nR || nB || nG)
				if (CVec3i(nR, nG, nB) == refcol)
				{
					tmpimg.at<unsigned char>(i, j) = 255;
				}
				else
					tmpimg.at<unsigned char>(i, j) = 0;
			}
		}

		sprintf(filename, "overlay%04d.png", view);
		cv::imwrite(resultdir + "/" + filename, tmpimg);

		delete[] poutbuffer;
	}

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	m_pMrtProgram->disable();

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void CRenderGL::render_to_image_overlay(const vector<cv::Mat> &RGBimg,const std::string &resultdir)
{
	char filename[256];
	char norname[256];

	float m_overratio = 0.5f;


	int numofcam = 1;
	numofcam = vec_glrenderCam.size();


	// First we bind the FBO so we can render to it
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_nFbo);

	// Save the view port and set it to the size of the texture
	glPushAttrib(GL_VIEWPORT_BIT);
	glViewport(0,0,m_imgwidth,m_imgheight);//**control the rendering resolution

	// Then render as normal
	// Today's scene is a wonderful multi-coloured spinning cube ;)
	// The destination of the data is controlled by the fragment shader we are using
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	m_pMrtProgram->use();

	GLint p_matrix, mv_matrix;
	p_matrix = m_iRenderProg.p_matrix;
	mv_matrix = m_iRenderProg.mv_matrix;

	glUniform4f(
		m_iRenderProg.light_position[0],
		m_iInputUniform.light_position[0][0], 
		m_iInputUniform.light_position[0][1],
		m_iInputUniform.light_position[0][2],
		m_iInputUniform.light_position[0][3]
	);


	
	GLint position = m_iRenderProg.position, 
		normal = m_iRenderProg.normal, 
		diffuse = m_iRenderProg.diffuse, 
		shininess = m_iRenderProg.shininess, 
		specular = m_iRenderProg.specular;

	glEnableVertexAttribArray(position);
	glEnableVertexAttribArray(normal);
	glEnableVertexAttribArray(diffuse);
	glEnableVertexAttribArray(shininess);
	glEnableVertexAttribArray(specular);

	//set the format for the vertex buffer
	glBindBuffer(GL_ARRAY_BUFFER, m_iTaskMesh.vertex_buffer);
	glVertexAttribPointer(
		position,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, position)
		);
	glVertexAttribPointer(
		normal,
		3, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, normal)
		);
	glVertexAttribPointer(
		diffuse,
		4, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, diffuse)
		);
	glVertexAttribPointer(
		shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, shininess)
		);
	glVertexAttribPointer(
		specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(render_vertex),
		(void*)offsetof(render_vertex, specular)
		);

	for (int view=0;view<numofcam;view++)
	{
		glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

		update_matrix_read (m_cvCameras->cameras[view],&m_iInputUniform);
		
		glUniformMatrix4fv(
			p_matrix,
			1, GL_FALSE,
			m_iInputUniform.p_matrix
			//vec_glrenderCam[view].p_matrix
			);
		glUniformMatrix4fv(
			mv_matrix,
			1, GL_FALSE,
			m_iInputUniform.mv_matrix
			//vec_glrenderCam[view].mv_matrix
			);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iTaskMesh.element_buffer);
		glDrawElements(
			GL_TRIANGLES,
			m_iTaskMesh.element_count,
			GL_UNSIGNED_INT,
			(void*)0
			);

		GLubyte *poutbuffer = new GLubyte[3*m_imgwidth*m_imgheight*sizeof(GLubyte)];

		glReadPixels(0, 0, m_imgwidth, m_imgheight,GL_RGB,GL_UNSIGNED_BYTE,poutbuffer);

		cv::Mat tmpimg = cv::Mat(m_imgheight, m_imgwidth, CV_8UC3);

#ifdef USEOPENMP
#pragma omp parallel for
#endif	
		for (int i=0;i<m_imgheight;i++)
		{
			for (int j=0;j<m_imgwidth;j++)
			{
				int nind = i*m_imgwidth+j;
				int nindgl = (m_imgheight-i-1)*m_imgwidth+j;

				unsigned char oldR = RGBimg[view].at<cv::Vec3b>(i, j)[2];
				unsigned char oldG = RGBimg[view].at<cv::Vec3b>(i, j)[1];
				unsigned char oldB = RGBimg[view].at<cv::Vec3b>(i, j)[0];

				float tmprat = 0.0f;
				int tmpsum = poutbuffer[nindgl*3+2]+poutbuffer[nindgl*3+1]+poutbuffer[nindgl*3];

				if(tmpsum)
				{
					tmprat = m_overratio;
					//pimg->imageData[nind*3]   = (unsigned char)(MIN(10*cvRound(fabs(m_oldB*(1.0f-m_tmprat) - poutbuffer[nindgl*3+2]*m_tmprat)),255));
					//pimg->imageData[nind*3+1] = (unsigned char)(MIN(10*cvRound(fabs(m_oldG*(1.0f-m_tmprat) - poutbuffer[nindgl*3+1]*m_tmprat)),255));
					//pimg->imageData[nind*3+2] = (unsigned char)(MIN(10*cvRound(fabs(m_oldR*(1.0f-m_tmprat) - poutbuffer[nindgl*3]*m_tmprat)),255));


				}else
				{					
					//pimg->imageData[nind*3]   = (unsigned char)cvRound(fabs(m_oldB*(1.0f-m_tmprat) - poutbuffer[nindgl*3+2]*m_tmprat));
					//pimg->imageData[nind*3+1] = (unsigned char)cvRound(fabs(m_oldG*(1.0f-m_tmprat) - poutbuffer[nindgl*3+1]*m_tmprat));
					//pimg->imageData[nind*3+2] = (unsigned char)cvRound(fabs(m_oldR*(1.0f-m_tmprat) - poutbuffer[nindgl*3]*m_tmprat));
				}

				
				tmpimg.at<cv::Vec3b>(i, j)[0] = (unsigned char)(poutbuffer[nindgl * 3 + 2] * tmprat + oldB*(1.0f - tmprat));
				tmpimg.at<cv::Vec3b>(i, j)[1] = (unsigned char)(poutbuffer[nindgl * 3 + 1] * tmprat + oldG*(1.0f - tmprat));
				tmpimg.at<cv::Vec3b>(i, j)[2] = (unsigned char)(poutbuffer[nindgl * 3] * tmprat + oldR*(1.0f - tmprat));

			}
		}

		sprintf(filename,"overlay%04d.png",view);

		cv::imwrite(resultdir + "/" + std::string(filename), tmpimg);
		
		delete[] poutbuffer;
	}

	glDisableVertexAttribArray(position);
	glDisableVertexAttribArray(normal);
	glDisableVertexAttribArray(diffuse);
	glDisableVertexAttribArray(shininess);
	glDisableVertexAttribArray(specular);

	// disable the shader
	m_pMrtProgram->disable();

	// Restore old view port and set rendering back to default frame buffer
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}


#endif
