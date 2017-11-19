#include "EpipolarConstraint.h"

int main()
{
	std::string root_folder = "C:\\Users\\shangxuanu\\Desktop\\image_original";

	EpipolarConstraint epipolar_constraint(root_folder);
	epipolar_constraint.getClickedPoint();
	epipolar_constraint.calculateEpipolarLine();
	epipolar_constraint.showEpipolarLine();
	
	return 0;
}