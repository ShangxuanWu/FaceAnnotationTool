#pragma once

#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
//#include "CVec.h"

#include "parameterFile.h"

#define X_GLOBAL_APP_STATE_FIELDS \
	X(int, s_mode) \
	X(int, s_renderWindowSizeX) \
	X(int, s_renderWindowSizeY) \
	X(int, s_renderImageSizeX) \
	X(int, s_renderImageSizeY) \
	X(int, s_offrenderMagFactor) \
	X(std::string, s_inputMeshFileDirectory) \
	X(std::string, s_inputCalibFileDirectory) \
	X(std::string, s_outputImagePath) \
	X(std::string, s_outputResultsDirectory) \
	X(int, s_renderViewPoint) \
	X(bool, s_b_renderWithUV) \
	X(std::string, s_UVtexturePath) \
	X(std::string, s_glShaderPath) \
	X(std::string, s_renderParams) \
	X(int, s_dynSeqFrameSt) \
	X(int, s_dynSeqFrameEnd) \
	X(int, s_dynSeqFrameInterval) \
	X(std::string, s_boardMeshDir) \
	X(int, s_numberFrmRendered) \
	X(std::string, s_hairMapPath) \
	X(bool, s_b_renderHair) \
	X(float, s_hairThickness)

	







//#define checkSizeArray(a, d)( (((sizeof a)/(sizeof a[0])) >= d))

class GlobalAppState
{
public:

#define X(type, name) type name;
	X_GLOBAL_APP_STATE_FIELDS
#undef X

		//! sets the parameter file and reads
	void readMembers(const ParameterFile& parameterFile) {
		s_ParameterFile = parameterFile;
		readMembers();
	}

	//! reads all the members from the given parameter file (could be called for reloading)
	void readMembers() {
#define X(type, name) \
	if (!s_ParameterFile.readParameter(std::string(#name), name)) {std::cout<<(std::string(#name).append(" ").append("uninitialized"))<<std::endl;	name = type();}
		X_GLOBAL_APP_STATE_FIELDS
#undef X
	

		m_bIsInitialized = true;
	}

	void print() const {
#define X(type, name) \
	std::cout << #name " = " << name << std::endl;
		X_GLOBAL_APP_STATE_FIELDS
#undef X
	}

	static GlobalAppState& getInstance() {
		static GlobalAppState s;
		return s;
	}
	static GlobalAppState& get() {
		return getInstance();
	}

	//! constructor
	GlobalAppState() {
		m_bIsInitialized = false;		
	}

	//! destructor
	~GlobalAppState() {
	}	
		

private:
	bool m_bIsInitialized;
	ParameterFile s_ParameterFile;	
};
