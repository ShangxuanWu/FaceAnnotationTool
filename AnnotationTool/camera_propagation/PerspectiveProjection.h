#ifndef _PERSPECTIVE_PROJECTION_H_
#define _PERSPECTIVE_PROJECTION_H_ 

#include <Eigen/Eigen>

//////////////////////////////////////////////////////////////////////////////////////////
// basics

// global 3d pos -> camera-local 3d pos
inline Eigen::Vector3d TransformGlobalPosIntoCameraCoordSystem(
	const Eigen::Vector3d& global_pos,
	const Eigen::MatrixXd& cam_extrinsics)
{
	return cam_extrinsics * global_pos.homogeneous();
}

// camera-local 3d pos -> normalized image coords
inline Eigen::Vector3d TransformLocalPosToNormalizedImageCoords(
	const Eigen::Vector3d& local_pos
	)
{
	double pos_z = local_pos[2];
	return local_pos / pos_z;
}

// normalized image coords -> pixel coords
inline Eigen::Vector3d TransformNormalizedImageCoordsToPixelCoords(
	const Eigen::Vector3d& normalized_image_coords,
	const Eigen::Matrix3d& cam_intrinsics
	)
{
	return cam_intrinsics * normalized_image_coords;
}

//////////////////////////////////////////////////////////////////////////////////////////
// wrapper

// global 3d pos -> normalized image coords
inline Eigen::Vector3d TransformGlobalPos3DToNormalizedImageCoords(
	const Eigen::Vector3d& global_pos,
	const Eigen::MatrixXd& cam_extrinsics)
{
	return
		TransformLocalPosToNormalizedImageCoords(TransformGlobalPosIntoCameraCoordSystem(global_pos, cam_extrinsics));
}

// global 3d pos -> pixel coords
inline Eigen::Vector3d TransformGlobalPos3DToPixelCoords(
	const Eigen::Vector3d& global_pos,
	const Eigen::MatrixXd& cam_extrinsics,
	const Eigen::Matrix3d& cam_intrinsics)
{
	return 
		TransformNormalizedImageCoordsToPixelCoords(
			TransformLocalPosToNormalizedImageCoords(TransformGlobalPosIntoCameraCoordSystem(global_pos, cam_extrinsics)),
			cam_intrinsics);
}



#endif  // _PERSPECTIVE_PROJECTION_H_