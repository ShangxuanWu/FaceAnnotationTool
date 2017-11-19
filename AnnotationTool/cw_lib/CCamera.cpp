

#include "CCamera.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>

#include "CVec.cpp"

CCamera::CCamera(void){
	A.MakeI();
	Rt.MakeI();
	P.MakeI();
	this->k0 = 0;
	this->k1 = 0;
	this->p0 = 0;
	this->p1 = 0;
	this->k2 = 0;
};

CCamera::CCamera(CMtx3x3f A, CMtx4x4f Rt, unsigned int camId)
{
	this->k0 = 0;
	this->k1 = 0;
	this->p0 = 0;
	this->p1 = 0;
	this->k2 = 0;

	Init(A, Rt, camId);
}


CCamera::CCamera(CMtx3x3f A, CMtx4x4f Rt, 
				 float k0, float k1, unsigned int camId)
{
	
	this->k0 = k0;
	this->k1 = k1;
	this->p0 = 0;
	this->p1 = 0;
	this->k2 = 0;
	Init(A,Rt,camId);
};

CCamera::CCamera(CMtx3x3f A, CMtx4x4f Rt,
	float k0, float k1, float p0, float p1, unsigned int camId)
{
	Init(A, Rt, camId);
	this->k0 = k0;
	this->k1 = k1;
	this->p0 = p0;
	this->p1 = p1;
	this->k2 = 0;
};

CCamera::CCamera(CMtx3x3f A, CMtx4x4f Rt,
	float k0, float k1, float p0, float p1, float k2, unsigned int camId)
{
	Init(A, Rt, camId);
	this->k0 = k0;
	this->k1 = k1;
	this->p0 = p0;
	this->p1 = p1;
	this->k2 = k2;
};

CCamera::~CCamera(void){};

void CCamera::Init(CMtx3x3f A, CMtx4x4f Rt, 
				   unsigned int camId){
					   
	CMtx4x4f	tmpA;
	CMtx3x3f	R;
	CVec3f t;

	this->A = A;
	this->Rt = Rt;	
	this->Id = camId;
	this->invA = A.Inv();


	tmpA.MakeI();
	for(int r=0;r<3;r++){
		for(int c=0;c<3;c++){
			tmpA(r,c) = A(r,c);
		}
	}
	P = tmpA*Rt;


	// Compute camera center
	for(int r=0;r<3;r++){
		for(int c=0;c<3;c++){
			R(r,c) = Rt(r,c);
			P3x3trans(r,c) = P(r,c);
		}
	}

	P3x3trans = P3x3trans.T();
	t = CVec3f(Rt(0,3),Rt(1,3),Rt(2,3));
	Center = -R.T()*t;
}

void CCamera::ComputeExtraInfo()
{
	CMtx4x4f	tmpA;
	CMtx3x3f	R;
	CVec3f t;

	tmpA.MakeI();
	for (int r = 0; r<3; r++){
		for (int c = 0; c<3; c++){
			tmpA(r, c) = A(r, c);
		}
	}
	P = tmpA*Rt;


	// Compute camera center
	for (int r = 0; r<3; r++){
		for (int c = 0; c<3; c++){
			R(r, c) = Rt(r, c);
			P3x3trans(r, c) = P(r, c);
		}
	}

	P3x3trans = P3x3trans.T();
	t = CVec3f(Rt(0, 3), Rt(1, 3), Rt(2, 3));
	Center = -R.T()*t;


}

void CCamera::printout(const char *filename) const
{
	FILE *pfile = fopen(filename, "w");
	fprintf(pfile, "0\n");
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++)
			fprintf(pfile, "%f ", A(i, j));
		fprintf(pfile, "\n");
	}
	fprintf(pfile, "%f %f %f %f %f\n", k0, k1, p0, p1,k2);
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 4; j++)
			fprintf(pfile, "%f ", Rt(i, j));
		fprintf(pfile, "\n");
	}
	fclose(pfile);

	

}


CCamera CCamera::reCenterCamera(const CCamera &curindcam)
{
	CMtx4x4f invRt = this->Rt.Inv();

	CMtx4x4f newRt = curindcam.Rt * invRt;

	CCamera retcam = curindcam;
	retcam.Rt = newRt;

	retcam.ComputeExtraInfo();
	return retcam;
}


void CCamera::Dump(FILE *fp){
	fprintf(fp,"Cam ID: %d\n",Id);
	fprintf(fp,"A = ");
	//A.Dump(fp);
	fprintf(fp,"Rt = ");
	//Rt.Dump(fp);
	fprintf(fp,"P = ");
	//P.Dump(fp);
	fprintf(fp,"k0 = %f \nk1= %f\n\n", k0, k1);
};


// Convert from world 3D coords to cam 3D coords 
CVec3f CCamera::WorldToEyeCoords(CVec3f p){
	CVec4f pt(p.x,p.y,p.z,1.0f);
	CVec4f eyePt = Rt*pt;
	eyePt = eyePt/eyePt.w;
	CVec3f result(eyePt.x,eyePt.y,eyePt.z);
	return result;
};


// Convert from world 3D coords to cam 3D coords 
void CCamera::WorldToImgCoords(const CVec3f &wp, CVec2f &ip) const{
	CVec4f pt(wp.x,wp.y,wp.z,1.0f);
	CVec4f imgPt;

	imgPt = P*pt;
	assert(imgPt.z!=0);
	ip.x = imgPt.x/imgPt.z;
	ip.y = imgPt.y/imgPt.z;
};

// Convert from world 3D coords to cam 3D coords
// accounting for radial distortion
void CCamera::WorldToImgCoordsRad(const CVec3f &wp, CVec2f &ip){
	CVec4f pt(wp.x,wp.y,wp.z,1.0f);
	CVec4f imgPt;
	float rSquared, rToTheFourth;

	// Project into normalized camera coordinates
	imgPt = Rt*pt;	
	assert(imgPt.z!=0);
	imgPt.x = imgPt.x/imgPt.z;
	imgPt.y = imgPt.y/imgPt.z;
	// Account for radial distortion
	rSquared = imgPt.x*imgPt.x + imgPt.y*imgPt.y;
	rToTheFourth = rSquared*rSquared;
	imgPt.x = imgPt.x * (1 + k0*rSquared + k1*rToTheFourth);
	imgPt.y = imgPt.y * (1 + k0*rSquared + k1*rToTheFourth);
	// Convert to image coordinates
	ip.x = A(0,2) + A(0,0)*imgPt.x + A(0,1)*imgPt.y;
	ip.y = A(1,2) + A(1,1)*imgPt.y;
};


// Fill OpenGL modelview matrix
void CCamera::GetModelViewMatrix(float *m){
	unsigned int n,r,c;

	n=0;
	// OpenGL matrices are column order
	for(c=0;c<4;c++){
		for(r=0;r<4;r++){
			m[n++]=Rt(r,c);
		}
	}
};

// Computes projection matrix for specified image dimensions and z clip planes
void CCamera::GetProjectionMatrix(float *p, unsigned int w, unsigned int h, 
								  float nearClip, float farClip){

	// Need xform from right-handed eye coords to left-handed NDC coords.
	// No change to x. fx*x + sx*y + cx*1, then div by W/2 and sub 1.0
	// Y needs to be inverted then offset by (h-Cy)
	// Z stays the same! 
	// Viewing frustum (last step): 
	// Must map Znear to -1, Zfar to +1
	// Remember, normalized device coordinates (NDC):
	// - after divide by w, get image coords by ORTHOGRAPHIC projection,
	//   i.e. you just use the x and y coordinates as is.
	// - the z coord is used just for clipping to (-1,1).
	// - NDC's are a right handed coordinate system again.
	// (OpenGL camera coordinates are not, but that doesn't matter
	// because I set the modelview and projection matrices myself.)
	//
	p[0] = 2.0f*A(0,0)/(float)w;
	p[4] = 2.0f*A(0,1)/(float)w;
	p[8] = 2.0f*A(0,2)/(float)w - 1.0f;
	p[12] = 0.0f;
	p[1] = 0.0f;
	p[5] = -2.0f*A(1,1)/(float)h;
	p[9] = 2.0f*(h-A(1,2))/(float)h - 1.0f;
	p[13] = 0.0f;
	p[2] = 0.0f;
	p[6] = 0.0f;
	p[10] = (farClip + nearClip)/(farClip - nearClip);
	p[14] = -2.0f * farClip * nearClip/
		(farClip - nearClip);
	p[3] = 0.0f;
	p[7] = 0.0f;
	p[11] = 1.0f;
	p[15] = 0.0f;
};

// This does what it's expected to do. 
// Verified starting from p, projection matrix. 
// Double checked in Matlab starting from line and P.
void CCamera::ImgLineToWorldNormal(const CVec3f &line, CVec3f &normal){
	normal = P3x3trans * line;
	normal = normal.Unit();
};


void CCamera::downsize(float downsizefactorx, float downsizefactory){
	for (int k = 0; k < 3; k++)
		A(0, k) *= downsizefactorx;
	for (int k = 0; k < 3; k++)
		A(1, k) *= downsizefactory;

	ComputeExtraInfo();

}

void CCamera::scaleCamera(float scl)
{
	for (int i = 0; i < 3; i++)
		this->Rt(i, 3) *= scl;

	ComputeExtraInfo();

}



