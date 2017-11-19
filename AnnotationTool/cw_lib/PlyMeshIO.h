//*******************************************
// PlyMeshIO.h
//
// generated by Chenglei Wu, ORP
//
// copyright reserved
//***************************************

#pragma once


#include <cstdlib>
#include <stdio.h>
#include "CVec.h"



class PlyMeshIO
{
	public:

		PlyMeshIO(){ m_bNeighborReOdered = false; };
		~PlyMeshIO(){};


		void loadPlyMesh(const std::string &plyname);
		void writePlyMesh(const std::string &plyname, bool binary = false, bool write_norm = false, bool float_color = false);

		
		void computeNormals();
		void computeAdjacentfaces();
		void reoderNeighors();
		void computeNeighbors();
		int obtainMaxNumofNeighbors();
		void smoothNormals(int Iters);

		void triangleOrderOpt();

		void NaNRemove();

		void remap_verts(const std::vector<int> &remap_table);
		void remove_vertices(const vector<bool> &toremove);
		void remove_unused_vertices();
		void remove_faces(const vector<bool> &toremove);

		bool connected(int f1, int f2, bool conn_vert);
		void find_connected(vector<int> &comps, vector<int> &compsizes, int f, int whichcomponent, bool conn_vert);
		void find_comps(vector<int> &comps, vector<int> &compsizes, bool conn_vert /* = false */);
		void select_big_comps(const vector<int> &comps, const vector<int> &compsizes, int min_size, int total_largest /* = std::numeric_limits<int>::max() */);
		void select_big_comps_auto(const vector<int> &comps, const vector<int> &compsizes, float ratio);

		void select_biggest_comps();

		void reOrderViaGraphGrouping(int upboundnum);

		int find_vetex_connected(vector<int> &comps, vector<int> &compsizes, int v, int whichcomponent, int upboundnum);
		void find_vertex_comps(vector<int> &comps, vector<int> &compsizes, int upboundnum);


		bool read_ply(FILE *f);

		bool read_faces_bin(FILE *f, bool need_swap,
			int nfaces, int face_len, int face_count, int face_idx);
		bool read_faces_asc(FILE *f, int nfaces,
			int face_len, int face_count, int face_idx, bool read_to_eol = false);


		bool read_verts_bin(FILE *f, bool &need_swap,
			int nverts, int vert_len, int vert_pos, int vert_norm,
			int vert_color, bool float_color, int vert_conf, int vert_uv);

		bool read_verts_asc(FILE *f,
			int nverts, int vert_len, int vert_pos, int vert_norm,
			int vert_color, bool float_color, int vert_conf, int vert_uv);


		void write_verts_bin(FILE *f, bool need_swap,
			bool write_norm, bool write_color, bool float_color,
			bool write_conf);
		void write_verts_asc(FILE *f,
			const char *before_vert,
			const char *before_norm,
			const char *before_color,
			bool float_color,
			const char *before_conf,
			const char *after_line);

		void write_faces_bin(FILE *f, bool need_swap,
			int before_face_len, const char *before_face,
			int after_face_len, const char *after_face);
		void write_faces_asc(FILE *f, const char *before_face, const char *after_line);
		
		void write_ply_header(FILE *f, const char *format, bool write_norm,  bool float_color);
		
		void write_ply_ascii(FILE *f, bool write_norm,  bool float_color);
		void write_ply_binary(FILE *f, bool need_swap, bool write_norm, bool float_color);

		void remove_duplicated_vertices();
		
		double computeLaplacianSum();
		double computeCotLaplacianSum();

		void computeCotangentWeights();
		void computeEdgeMappedVertices();

			
		vector<CVec3f> m_vertices;
		vector<CVec3f> m_normals;
		vector<CVec3f> m_colors;
		vector<float> m_confidences;

		vector<vector<int> > m_neighbors;
		vector<vector<int> > m_adjacentfaces;
		bool m_bNeighborReOdered;
		vector<vector<int> > m_edgemappedvetindices;

		vector<vector<float> > m_cotweights;


		vector<CVec3i> m_faces;

		vector<CVec3f> m_facenormals;

		vector<CVec2f> m_uvcoord;


		//vector<int> m_grid;
		//int m_grid_width;
		//int m_grid_height;

	private:
		
};

