#pragma once


#include <stdlib.h>
#include <vector>
#include "CCamera.h"
#include "CVec.h"

#ifndef _H_CAMERAARRAY
#define _H_CAMERAARRAY

class CCameraArray{

public:
	CCameraArray(void);
	CCameraArray(const char * pszFile);
	CCameraArray(const char * pszFile, int distType);
	CCameraArray(const CCameraArray * input);
	~CCameraArray(void);
	bool LoadCalibrationData(const char * pszFile);
	CCamera GetNthCam(unsigned int n){return cameras[n];};
	CVec3f GetNthCamCenter(unsigned int n){return((cameras[n]).GetCenter());};
	unsigned int ArrayNumCams(){return (unsigned int)(cameras.size());};
	void WorldToImgCoords(const CVec3f &wp, CVec2f &ip, unsigned int n)
		{return (cameras[n]).WorldToImgCoords(wp,ip);};
	void ImgLineNthCamToWorldNormal(const CVec3f &line, CVec3f &normal, 
		unsigned int cam){(cameras[cam]).ImgLineToWorldNormal(line, normal);};
	//	CCamera CCameraArray::GetCamById(unsigned int camID);

	CMtx3x3f GetFundamentalMatrix(int camind0, int camind1);
	CMtx3x3f GetFundamentalMatrix2(int camind0, int camind1);
	CMtx3x3f GetFundamentalMatrix2(int ind_left, int ind_right, float downsizefactor);
	CMtx3x3f GetFundamentalMatrix2(int ind_left, int ind_right, float downsizefactorx, float downsizefactory);

	void scaleCameraCoordinates(float scl);

	CCamera ReCenterCamera(int refind, int curind);
	vector<int> ViewSelection(int refind, float min, float max);

	void DownSizeCameras(float downsizefactor);
	void DownSizeCameras(float downsizefactorx, float downsizefactory);

	void printout(const std::string &outname);
	
public:
	std::vector<CCamera>		cameras;
	std::vector<unsigned int>	camNums;

	int m_distType;
};


void CV2GLprojmatrix(float cparam[3][4], int width, int height, float gnear, float gfar, float m[16]);



#endif
