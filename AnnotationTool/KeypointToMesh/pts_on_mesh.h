#pragma once
// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu

#ifndef __PTS_ON_MESH_H_INCLUDED__
#define __PTS_ON_MESH_H_INCLUDED__


#include "../camera_propagation/myheader.h"
#include "../camera_propagation/pts_3d_conf.h"

class pts_on_mesh : public pts_3d_conf {
public:
	int vertice_id;						// vertice id

	pts_on_mesh(int pts_id, double tx, double ty, double tz, double tconf)			: pts_3d_conf(tx, ty, tz, tconf), vertice_id(pts_id) {}
	pts_on_mesh(int pts_id, double tx, double ty, double tz)						: pts_3d_conf(tx, ty, tz), vertice_id(pts_id) {}
	pts_on_mesh(int pts_id, float tx, float ty, float tz)							: pts_3d_conf(tx, ty, tz), vertice_id(pts_id){}
	pts_on_mesh(int pts_id, float tx, float ty, float tz, double tconf)				: pts_3d_conf(tx, ty, tz, tconf), vertice_id(pts_id) {}
	pts_on_mesh(int pts_id, int tx, int ty, int tz)									: pts_3d_conf(tx, ty, tz), vertice_id(pts_id) {}
	pts_on_mesh(int pts_id, int tx, int ty, int tz, double tconf)					: pts_3d_conf(tx, ty, tz, tconf), vertice_id(pts_id) {}

	cv::Point3d convert_to_point3d();
	std::vector<double> convert_to_pts_vec();
	std::vector<double> convert_to_pts_vec_conf();
	void print(int prec = default_precision);
};


#endif