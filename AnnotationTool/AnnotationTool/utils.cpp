#include "utils.h"

// get the last 2 char of the folder name as ID
std::string getIDFromFilename(std::string str)
{
	return str.substr(str.length() - 16, 2);
}

std::string getFolderFromFilename(std::string full_path)
{
	return full_path.substr(full_path.length() - 20, 6);
}

std::string getIDFromFoldername(std::string folder_name)
{
	return folder_name.substr(folder_name.length() - 2, 2);
}

std::string getFolderFromID(std::string ID)
{
	return "3300"+ID;
}

std::string getFullPathFromFolder(std::string root, std::string folder_name)
{
	return root + "/cam" + folder_name + "/image0000.png";
}

mycamera getMycameraFromFullPath(QString imagePaths, std::vector<mycamera> camera_cluster)
{
	std::string folder_name = getFolderFromFilename(imagePaths.toStdString());
	for (int i = 0; i < camera_cluster.size(); i++)
	{
		if (camera_cluster[i].name == folder_name) 
		{
			return camera_cluster[i];
		}
	}
}

void rot90(cv::Mat &matImage, int rotflag) {
	//1=CW, 2=CCW, 3=180
	if (rotflag == RIGHT) {
		transpose(matImage, matImage);
		flip(matImage, matImage, 1); //transpose+flip(1)=CW
	}
	else if (rotflag == LEFT) {
		transpose(matImage, matImage);
		flip(matImage, matImage, 0); //transpose+flip(0)=CCW     
	}
	//else if (rotflag == 3) {
	//	flip(matImage, matImage, -1);    //flip(-1)=180          
	//}
	//else if (rotflag != 0) { //if not 0,1,2,3:
	else
	{
		qDebug() << "Unknown rotation flag(" << rotflag << ")" << endl;
	}
}

void helper_map_to_vec_tool_to_conf(std::map<std::string, std::vector<pts_2d_conf>> pts_src)
{
	;
}
