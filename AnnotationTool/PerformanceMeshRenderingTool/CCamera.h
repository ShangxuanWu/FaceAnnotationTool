//////////////////////////////////////////////////////////////
//
// Camera.h 
//

#pragma once


#include <assert.h>
#include "CVec.h"

#ifndef _H_CAMERA
#define _H_CAMERA

class CCamera {

public:
	CCamera(void);
	CCamera(CMtx3x3f A, CMtx4x4f Rt, 
		float k0, float k1, unsigned int camId);
	CCamera(CMtx3x3f A, CMtx4x4f Rt,unsigned int camId);
	CCamera(CMtx3x3f A, CMtx4x4f Rt,
		float k0, float k1, float p0, float p1, unsigned int camId);
	CCamera(CMtx3x3f A, CMtx4x4f Rt,
		float k0, float k1, float p0, float p1, float k2, unsigned int camId);

	void Init(CMtx3x3f A, CMtx4x4f Rt, unsigned int camId);
	
	~CCamera(void);
	void Dump(FILE *fp = stdout);			// Dumps camera params
	CMtx4x4f GetRt() {return Rt;};
	CMtx4x4f GetA() {return A;};
	CVec3f WorldToEyeCoords(CVec3f p);
	void WorldToImgCoords(const CVec3f &wp, CVec2f &ip) const;
	void WorldToImgCoordsRad(const CVec3f &wp, CVec2f &ip);
	void GetModelViewMatrix(float *MVM);
	void GetProjectionMatrix(float *PM, unsigned int w, unsigned int h,
		float nearClip, float farClip);
	CVec3f GetCenter(void){return Center;}; 
	void ImgLineToWorldNormal(const CVec3f &line, CVec3f &normal);

	void printout(const char *filename) const;

	void downsize(float dx, float dy);

	void scaleCamera(float scl);

	CCamera reCenterCamera(const CCamera &curindcam);
	
	void ComputeExtraInfo(); // Computes Center, P, P3x3trans from A & Rt

public:
	CMtx3x3f	A;		// intrinsics
	CMtx3x3f	invA;		// inverse of intrinsics
	CMtx4x4f	Rt;		// extrinsics
	CMtx4x4f	P;		// projection matrix
	float			k0,k1,k2;	// radial distortion coefficients
	float			p0, p1;
	unsigned int	Id;		// Camera Id number
	CVec3f		Center; // Camera center in world coords = -inv(R)*t
	CMtx3x3f	P3x3trans; // P(1:3,1:3)', use for mapping lines to normals 
//private:
	
};




#endif

