
#include <iostream>
#include <fstream>
#include <sstream>

#include "CCameraArray.h"
#include "CVec.cpp"

using namespace std;

CCameraArray::CCameraArray(void){};

CCameraArray::~CCameraArray(void){};

CCameraArray::CCameraArray(const char * pszFile){
	m_distType = 5;
	LoadCalibrationData(pszFile);
};

CCameraArray::CCameraArray(const char * pszFile, int distType)
{
	m_distType = distType;
	LoadCalibrationData(pszFile);
};

CCameraArray::CCameraArray(const CCameraArray * input){
	if(!input)
		std::cout<<"invalid input CameraArray!\n";

	for(int i=0;i<input->cameras.size();i++)
	{
		CCamera camera_no;
		camera_no=CCamera(input->cameras[i].A,input->cameras[i].Rt,input->cameras[i].k0,input->cameras[i].k1,i);
		cameras.push_back(camera_no);
	}
	
};


CMtx3x3f CCameraArray::GetFundamentalMatrix(int ind_left, int ind_right)
{
	CMtx3x3f mat_r0, mat_r1;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j<3; j++)
		{
			mat_r0(i, j) = this->cameras[ind_left].Rt(i, j);
			mat_r1(i, j) = this->cameras[ind_right].Rt(i, j);
		}
	}

	CMtx3x3f mat_refr = mat_r1*mat_r0.T();

	CVec3f vec_reft = this->cameras[ind_right].Center - this->cameras[ind_left].Center;
	vec_reft = mat_r0*vec_reft;//pay attention to this, this should be rotated to left camera

	CMtx3x3f res, tmpmat;
	tmpmat = (this->cameras[ind_right].A).Inv();
	res = tmpmat.T()*mat_refr;
	tmpmat = res * (vec_reft.CrossMatrix());
	res = tmpmat * ((this->cameras[ind_left].A).Inv());


	return res;
}


CMtx3x3f CCameraArray::GetFundamentalMatrix2(int ind_left, int ind_right)
{
	CMtx3x3f mat_r0, mat_r1;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j<3; j++)
		{
			mat_r0(i, j) = this->cameras[ind_left].Rt(i, j);
			mat_r1(i, j) = this->cameras[ind_right].Rt(i, j);
		}
	}

	CMtx3x3f mat_refr = mat_r1*mat_r0.T();

	CVec3f vec_reft = this->cameras[ind_right].Center - this->cameras[ind_left].Center;
	//vec_reft = mat_r0*vec_reft;//pay attention to this, this should be rotated to left camera

	CMtx3x3f res, tmpmat;
	tmpmat = (this->cameras[ind_right].A).Inv();
	res = tmpmat.T()*mat_r1;
	tmpmat = res * (vec_reft.CrossMatrix());
	tmpmat = tmpmat * mat_r0.T();
	res = tmpmat * ((this->cameras[ind_left].A).Inv());


	return res;
}


CMtx3x3f CCameraArray::GetFundamentalMatrix2(int ind_left, int ind_right, float downsizefactor)
{
	CMtx3x3f mat_r0, mat_r1;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j<3; j++)
		{
			mat_r0(i, j) = this->cameras[ind_left].Rt(i, j);
			mat_r1(i, j) = this->cameras[ind_right].Rt(i, j);
		}
	}

	CMtx3x3f mat_refr = mat_r1*mat_r0.T();

	CMtx3x3f leftIntMat = this->cameras[ind_left].A;
	CMtx3x3f rightIntMat = this->cameras[ind_right].A;

	//downsize the intrinsics due to downsizing the image
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			leftIntMat(i, j) *= downsizefactor;
			rightIntMat(i, j) *= downsizefactor;
		}
	}
	


	CVec3f vec_reft = this->cameras[ind_right].Center - this->cameras[ind_left].Center;
	//vec_reft = mat_r0*vec_reft;//pay attention to this, this should be rotated to left camera

	CMtx3x3f res, tmpmat;
	tmpmat = (rightIntMat).Inv();
	res = tmpmat.T()*mat_r1;
	tmpmat = res * (vec_reft.CrossMatrix());
	tmpmat = tmpmat * mat_r0.T();
	res = tmpmat * ((leftIntMat).Inv());


	return res;
}

void CCameraArray::DownSizeCameras(float downsizefactor)
{
	for (int i = 0; i < cameras.size(); i++)
	{
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 3; k++)
				cameras[i].A(j, k) *= downsizefactor;
	}
}


void CCameraArray::printout(const std::string &outname)
{
	FILE *pfile = fopen(outname.c_str(), "w");
	for (int c = 0; c < cameras.size(); c++)
	{
		const CCamera &cam = cameras[c];
		fprintf(pfile, "%d\n", cam.Id);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				fprintf(pfile, "%f ", cam.A(i, j));
			}
			fprintf(pfile, "\n");
		}

		fprintf(pfile, "%f %f %f %f %f\n", cam.k0, cam.k1, cam.p0, cam.p1, cam.k2);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				fprintf(pfile, "%f ", cam.Rt(i, j));
			}
			fprintf(pfile, "\n");
		}

		fprintf(pfile, "\n");
	}
	fclose(pfile);

}

void CCameraArray::DownSizeCameras(float downsizefactorx, float downsizefactory)
{
	for (int i = 0; i < cameras.size(); i++)
	{		
		for (int k = 0; k < 3; k++)
			cameras[i].A(0, k) *= downsizefactorx;
		for (int k = 0; k < 3; k++)
			cameras[i].A(1, k) *= downsizefactory;
	}
}

CMtx3x3f CCameraArray::GetFundamentalMatrix2(int ind_left, int ind_right, float downsizefactorx, float downsizefactory)
{
	CMtx3x3f mat_r0, mat_r1;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j<3; j++)
		{
			mat_r0(i, j) = this->cameras[ind_left].Rt(i, j);
			mat_r1(i, j) = this->cameras[ind_right].Rt(i, j);
		}
	}

	CMtx3x3f mat_refr = mat_r1*mat_r0.T();

	CMtx3x3f leftIntMat = this->cameras[ind_left].A;
	CMtx3x3f rightIntMat = this->cameras[ind_right].A;

	//downsize the intrinsics due to downsizing the image
	
	for (int j = 0; j < 3; j++)
	{
		leftIntMat(0, j) *= downsizefactorx;
		rightIntMat(0, j) *= downsizefactorx;
	}
	for (int j = 0; j < 3; j++)
	{
		leftIntMat(1, j) *= downsizefactory;
		rightIntMat(1, j) *= downsizefactory;
	}




	CVec3f vec_reft = this->cameras[ind_right].Center - this->cameras[ind_left].Center;
	//vec_reft = mat_r0*vec_reft;//pay attention to this, this should be rotated to left camera

	CMtx3x3f res, tmpmat;
	tmpmat = (rightIntMat).Inv();
	res = tmpmat.T()*mat_r1;
	tmpmat = res * (vec_reft.CrossMatrix());
	tmpmat = tmpmat * mat_r0.T();
	res = tmpmat * ((leftIntMat).Inv());


	return res;
}

CCamera CCameraArray::ReCenterCamera(int refind, int curind)
{
	CMtx4x4f invRt = cameras[refind].Rt.Inv();

	CMtx4x4f newRt = cameras[curind].Rt * invRt;

	CCamera retcam = cameras[curind];
	retcam.Rt = newRt;
	return retcam;
}

vector<int> CCameraArray::ViewSelection(int refind, float min, float max)
{
	float minval = min / 180.0f*M_PI;
	float maxval = max / 180.0f*M_PI;

	//to be implemented
	vector<CVec3f> vec_dirs(cameras.size(), CVec3f(0, 0, 1));
	for (int i = 0; i < vec_dirs.size(); i++)
	{
		CMtx3x3f invR;
		for (int k = 0; k < 3; k++)
		{
			for (int l = 0; l < 3; l++)
				invR(k, l) = cameras[i].Rt(l, k);
		}
		vec_dirs[i] = (invR * CVec3f(0, 0, 1.0f)).Unit();
		
	}
	CVec3f refdir = vec_dirs[refind];
	vector<int> vec_ret;

	CVec3f ctrpos = cameras[refind].Center;
	vector<pair<float, int> > vectmp_pair;
	for (int i = 0; i < cameras.size(); i++)
	{
		if (i == refind)
			continue;
		CVec3f tmppos = cameras[i].Center;
		float tmpbaseline = (tmppos - ctrpos).Magnitude();
		float tmpdot = refdir * vec_dirs[i];
		float tmpangle = acos(tmpdot);
		if (tmpangle > minval && tmpangle < maxval)
		{
#ifdef SELECT_CAMERA_ACCORDING_ANGLE
			vectmp_pair.push_back(pair<float, int>(tmpangle, i));
#else
			vectmp_pair.push_back(pair<float, int>(tmpbaseline, i));
#endif
		}
	}

	std::sort(vectmp_pair.begin(), vectmp_pair.end());
	float numneiview = vectmp_pair.size();
	for (int c = 0; c < numneiview; c++)
	{
		vec_ret.push_back( vectmp_pair[c].second);
		//std::cout << vectmp_pair[tmpnum-1-c].second << " ";
	}

	
	//to show the baseline
	std::cout << "baseline: " << std::endl;
	for (int i = 0; i < vec_ret.size(); i++)
	{
		int tmpind = vec_ret[i];
		CVec3f tmppos = cameras[tmpind].Center;
		std::cout << i << " Camera Center " << tmppos.x << "  " << tmppos.y << " " << tmppos.z << std::endl;
		float tmpbaseline = (tmppos - ctrpos).Magnitude();
		std::cout << refind << " to " << tmpind << ": " << tmpbaseline << std::endl;
	}

	return vec_ret;
}


void CCameraArray::scaleCameraCoordinates(float scl)
{
	for (int i = 0; i < cameras.size(); i++)
	{
		cameras[i].scaleCamera(scl);
	}
}




bool CCameraArray::LoadCalibrationData(const char * pszFile){

	ifstream fp_in;	
	unsigned int	camId;		// camera Id number
	CMtx3x3f		A;			// intrinsics
	CMtx4x4f		Rt;			// extrinsics
	CMtxf			P;			// projection matrix
	float			k0,k1,k2;		// radial distortion coefficients
	float			p0, p1;

	// Open input file
	fp_in.open(pszFile, ios::in);   
	if (fp_in.fail()){
		cameras.clear();
		camNums.clear();
		cout << "Failed!\n";
		return false;
	}

	// Initialize camera params
	A.MakeI();
	Rt.MakeI();
	P=0;

	char calibfile[]="ProjMatrix.txt";
	FILE *cfp=fopen(calibfile,"w");

	int index=0;

	std::cout << "distortion Type " << m_distType << std::endl;

	// Load calibration data
	// Try to read camNum first, will get eof if no data left.
	while(fp_in>>camId && !fp_in.eof()){
			fp_in >> A(0,0);
			fp_in >> A(0,1);
			fp_in >> A(0,2);
			fp_in >> A(1,0);
			fp_in >> A(1,1);
			fp_in >> A(1,2);
			fp_in >> A(2,0);
			fp_in >> A(2,1);
			fp_in >> A(2,2);	

			std::string discoeff;
			getline(fp_in, discoeff);
			getline(fp_in, discoeff);
			std::cout << discoeff << std::endl;
			switch (m_distType)
			{
			case 0:
				break;
			case 2:
				//fp_in >> k0;
				//fp_in >> k1;
				sscanf(discoeff.c_str(), "%f %f", &k0, &k1);
				break;
			case 4:
				//fp_in >> k0;
				//fp_in >> k1;
				//fp_in >> p0;
				//fp_in >> p1;
				sscanf(discoeff.c_str(), "%f %f %f %f", &k0, &k1,&p0,&p1);
				break;
			case 5:
				//fp_in >> k0;
				//fp_in >> k1;
				//fp_in >> p0;
				//fp_in >> p1;
				//fp_in >> k2;
				sscanf(discoeff.c_str(), "%f %f %f %f %f", &k0, &k1, &p0, &p1,&k2);
				break;
			default:
				std::cout << "CANNOT find distortation mode " << std::endl;
			}			
			fp_in >> Rt(0,0);
			fp_in >> Rt(0,1);
			fp_in >> Rt(0,2);
			fp_in >> Rt(0,3);
			fp_in >> Rt(1,0);
			fp_in >> Rt(1,1);
			fp_in >> Rt(1,2);
			fp_in >> Rt(1,3);
			fp_in >> Rt(2,0);
			fp_in >> Rt(2,1);
			fp_in >> Rt(2,2);
			fp_in >> Rt(2,3);
			CCamera camera_no;


			switch (m_distType)
			{
			case 0:
				camera_no = CCamera(A, Rt, camId);
				break;
			case 2:
				camera_no = CCamera(A, Rt, k0, k1, camId);
				break;
			case 4:
				camera_no = CCamera(A, Rt, k0, k1, p0, p1, camId);
				break;
			case 5:
				camera_no = CCamera(A, Rt, k0, k1, p0, p1, k2, camId);
				break;
			default:
				std::cout << "CANNOT find distortation mode!!! " << std::endl;
			}
			
			cameras.push_back(camera_no);
			camNums.push_back(camId);
	
			fprintf(cfp,"%d\n",index);
			for(int r=0;r<3;r++){
				for(int c=0;c<4;c++){
					if(c==3){fprintf(cfp,"%f",camera_no.P(r,c));}
					else if(c==2)
					{
						fprintf(cfp,"%f\t",camera_no.P(r,c));
					}
					else
					{
						fprintf(cfp,"%f\t",camera_no.P(r,c));
					}
				}
				fprintf(cfp,"\n");
			}
			index++;
	}
	fclose(cfp);
	fp_in.close();
	return true;
};



void CV2GLprojmatrix(float cparam[3][4], int width, int height, float
	gnear, float gfar, float m[16])
{
	float   p[3][3], q[4][4];


	/* Camera parameters are converted openGL representation. */
	/* Camera parameter represents the transformation from camera coordinates
	to screen coordinates[pixel]. OpenGL projection matrix represents the
	transformation from the camera coordinates to normalized view volume. */
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			p[i][j] = cparam[i][j];// / cparam[2][2];
		}
	}
	q[0][0] = (2.0 * p[0][0] / width);
	q[0][1] = (2.0 * p[0][1] / width);
	q[0][2] = ((2.0 * p[0][2] / width) - 1.0);
	q[0][3] = 0.0;

	q[1][0] = 0.0;
	q[1][1] = -(2.0 * p[1][1] / height);
	q[1][2] = -((2.0 * p[1][2] / height) - 1.0);
	q[1][3] = 0.0;

	q[2][0] = 0.0;
	q[2][1] = 0.0;
	q[2][2] = (gfar + gnear) / (gfar - gnear);
	q[2][3] = -2.0 * gfar * gnear / (gfar - gnear);

	q[3][0] = 0.0;
	q[3][1] = 0.0;
	q[3][2] = 1.0;
	q[3][3] = 0.0;

	for (int i = 0; i<4; i++)
	{
		for (int j = 0; j<4; j++)
		{
			m[j * 4 + i] = q[i][j];
		}
	}
}