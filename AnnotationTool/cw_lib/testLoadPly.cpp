#include "PlyMeshIO.h"

#include <string.h>
#include <iostream>
#include <stdio.h>

using namespace std;

PlyMeshIO sourceMesh;

//std::sprintf(info, "%05d", mdb.startFrame);
//if (fromNeutral) { meshFile = neutralTakeDir + "/scans/" + string(info) + ".ply"; }
//else { meshFile = takeDir + "/scans/" + string(info) + ".ply"; }

string meshFile = "D:\oculus\mesh\sahana-expressions\67200\mesh.ply";
sourceMesh.loadPlyMesh(meshFile);
//cout << "Loading " << meshFile << endl;
//if (sourceMesh.m_vertices.size() == 0) { std::cout << "...Error!" << std::endl; }

//std::sprintf(info, "%05d", mdb.nextFrame);
//meshFile = takeDir + "/scans/" + info + ".ply";
//targetMesh.loadPlyMesh(meshFile);
//std::cout << "Loading " << meshFile << endl;
//if (targetMesh.m_vertices.size() == 0) { std::cout << "...Error!" << std::endl; }
