#include "scene/mesh.hpp"
#include <map>

namespace _462 {

	static void compute_normals(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles)
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			vertices[i].normal = Vector3::Zero;
		}

		for (size_t i = 0; i < triangles.size(); i++)
		{
			unsigned int index0 = triangles[i].vertices[0];
			unsigned int index1 = triangles[i].vertices[1];
			unsigned int index2 = triangles[i].vertices[2];

			Vector3 surface_normal = normalize(cross(vertices[index1].position - vertices[index0].position, vertices[index2].position - vertices[index0].position));
			vertices[index0].normal += surface_normal;
			vertices[index1].normal += surface_normal;
			vertices[index2].normal += surface_normal;
		}
		//average all the surface normals
		for (size_t i = 0; i < vertices.size(); i++)
		{
			vertices[i].normal = normalize(vertices[i].normal);
		}
	}

	static void add_new_dots(unsigned int new_vertex_index, Mesh::MeshTriangleList& new_triangles, unsigned int old_triangle_index, unsigned int a_vertex_tri_index, unsigned int b_vertex_tri_index)
	{
		unsigned int new_triangle_index[] = {old_triangle_index*4, old_triangle_index*4, old_triangle_index*4 + 3};
		unsigned int new_triangle_vertex_slot[3] = {};
		if (a_vertex_tri_index > b_vertex_tri_index)
		{
			std::swap(a_vertex_tri_index, b_vertex_tri_index);
		}

		switch (a_vertex_tri_index)
		{
		case 0:
		{
			switch (b_vertex_tri_index)
			{
			case 1:
				new_triangle_index[1] += 1;
				new_triangle_vertex_slot[0] += 1;
				new_triangle_vertex_slot[1] += 2;
				break;
			case 2:
				new_triangle_index[1] += 2;
				new_triangle_vertex_slot[0] += 2;
				new_triangle_vertex_slot[1] += 1;
				new_triangle_vertex_slot[2] += 2;
				break;
			}
			break;
		}
		case 1:
		{
			// the b_vertex_tri_index then musti be  2
			new_triangle_index[0] += 1;
			new_triangle_index[1] += 2;
			new_triangle_vertex_slot[0] += 1;
			new_triangle_vertex_slot[1] += 2;
			new_triangle_vertex_slot[2] += 1;
			break;
		}
		}

		for (int i = 0; i < 3; i++)
		{
			MeshTriangle& new_triangle = new_triangles[new_triangle_index[i]];
			new_triangle.vertices[new_triangle_vertex_slot[i]] = new_vertex_index;
		}
	}

static void first_pass(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles, Mesh::MeshEdgeList& edges, Mesh::MeshTriangleList& new_triangles)
{	
	const float a_weight_interior = 3.0f / 8.0f;
	const float b_weight_interior = a_weight_interior;
	const float c_weight = 1.0f / 8.0f;
	const float d_weight = c_weight;
	const float a_weight_boundary = 0.5f;
	const float b_weight_boundary = 0.5f;
	unsigned int a_vertex_tri_index;
	unsigned int b_vertex_tri_index;
	unsigned int a_vertex_tri_1_index;
	unsigned int b_vertex_tri_1_index;

	for (int i = 0; i < edges.size(); i++)
	{
		MeshVertex v;
		const MeshEdge& e = edges[i];
		unsigned int a_index = e.vertices[0];
		unsigned int b_index = e.vertices[1];
		const MeshVertex& a = vertices[a_index];
		const MeshVertex& b = vertices[b_index];

		const MeshTriangle& triangle = triangles[e.triangles[0]];
		for (int j = 0; j < 3; j++)
		{
			if (triangle.vertices[j] == a_index)
			{
				a_vertex_tri_index = j;
			}
			else if (triangle.vertices[j] == b_index)
			{
				b_vertex_tri_index = j;
			}
		}

		if (e.triangle_size == 1)
		{
			//boundary case
			v.position = a_weight_boundary * a.position + b_weight_boundary * b.position;
			add_new_dots(vertices.size(), new_triangles, e.triangles[0], a_vertex_tri_index, b_vertex_tri_index);
		}
		else
		{
			//interior case
			unsigned int c_index, d_index;
			const MeshTriangle& tri_0 = triangles[e.triangles[0]];
			const MeshTriangle& tri_1 = triangles[e.triangles[1]];
			for (int j = 0; j < 3; j++)
			{
				if (tri_0.vertices[j] != a_index && tri_0.vertices[j] != b_index)
				{
					c_index = tri_0.vertices[j];
				}

				if (tri_1.vertices[j] != a_index && tri_1.vertices[j] != b_index)
				{
					d_index = tri_1.vertices[j];
				}
			}
			const MeshVertex& c = vertices[c_index];
			const MeshVertex& d = vertices[d_index];

			v.position = a_weight_interior * a.position + b_weight_interior * b.position + c_weight * c.position + d_weight * d.position;	

			const MeshTriangle& triangle1 = triangles[e.triangles[1]];
			for (int j = 0; j < 3; j++)
			{
				if (triangle1.vertices[j] == a_index)
				{
					a_vertex_tri_1_index = j;
				}
				else if (triangle1.vertices[j] == b_index)
				{
					b_vertex_tri_1_index = j;
				}
			}

			add_new_dots(vertices.size(), new_triangles, e.triangles[0], a_vertex_tri_index, b_vertex_tri_index);
			add_new_dots(vertices.size(), new_triangles, e.triangles[1], a_vertex_tri_1_index, b_vertex_tri_1_index);
		}
		vertices.push_back(v);
	}
}

static void second_pass(Mesh::MeshVertexList& vertices, Mesh::MeshEdgeList& edges, unsigned int old_vertices_size)
{
	const float v_weight_boundary = 0.75f;
	const float a_weight_boundary = 1.0f / 8.0f;
	const float b_weight_boundary = a_weight_boundary;	

	for (int i = 0; i < old_vertices_size; i++)
	{
		MeshVertex& v = vertices[i];
		if (v.is_boundary || v.edges.size() == 2)
		{
			//boundary case
			unsigned int a_index, b_index;
			bool found_a = false;

			for (int j = 0; j < v.edges.size(); j++)
			{
				const MeshEdge& e = edges[v.edges[j]];
				if (e.triangle_size == 1)
				{
					if (!found_a)
					{
						a_index = e.vertices[0] != i ? e.vertices[0] : e.vertices[1];
						found_a = true;
					}
					else
					{
						b_index = e.vertices[0] != i ? e.vertices[0] : e.vertices[1];
						break;
					}
				}
			}
			MeshVertex& a = vertices[a_index];
			MeshVertex& b = vertices[b_index];
			v.position = v_weight_boundary * v.position + a_weight_boundary * a.position + b_weight_boundary * b.position;
		}
		else
		{
			//interiors case
			const int N = v.edges.size();
			float beta = (0.625f - pow(0.375f + 0.25f*cos(2 * PI / N), 2)) / N;
			Vector3 sum = Vector3::Zero;

			for (int j = 0; j < N; j++)
			{
				const MeshEdge& e = edges[v.edges[j]];
				unsigned int u_index = e.vertices[0] != i ? e.vertices[0] : e.vertices[1];
				const MeshVertex& u = vertices[u_index];
				sum += u.position;
			}

			v.position = ((1 - beta * N) * v.position) + beta * (sum);
		}
	}
}

static void assemble_new_triangles(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles, Mesh::MeshTriangleList& new_triangles)
{
	for (int i = 0; i < triangles.size(); i++)
	{
		const MeshTriangle& old_triangle = triangles[i];
		for (int i = 0; i < 3; i++)
		{
			MeshTriangle t;
			t.vertices[0] = old_triangle.vertices[i];
			new_triangles.push_back(t);
		}
		//vertices for triangle index 3 all consist of new vertices from first pass
		new_triangles.push_back(MeshTriangle());
	}
}

static void recreate_new_edges(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles, Mesh::MeshEdgeList& edges, unsigned int old_vertices_size)
{
	typedef std::map<std::pair<int, int>, unsigned int> EdgeMap;
	EdgeMap edge_map;	
	unsigned int edge_idx_counter = 0; // current edge index, for creating new edges

	for (int i = 0; i < triangles.size(); i++)
	{
		MeshTriangle& tri = triangles[i];

		for (int j = 0; j < 3; j++)
		{
			int next = (j + 1) >= 3 ? 0 : j + 1;
			unsigned int cur_v_idx = tri.vertices[j];
			unsigned int next_v_idx = tri.vertices[next];
			if (cur_v_idx > next_v_idx)
			{
				std::swap(cur_v_idx, next_v_idx);
			}
			std::pair< EdgeMap::iterator, bool> edgemap_check_1 = edge_map.insert(std::make_pair(std::make_pair(cur_v_idx, next_v_idx), edge_idx_counter));
			if (edgemap_check_1.second) //  if edge that is connected by two vertices doesnt exist yet
			{
				MeshEdge e;
				e.vertices[0] = cur_v_idx;
				e.vertices[1] = next_v_idx;
				e.triangles[0] = i;
				e.triangle_size = 1;
				edges.push_back(e);
				tri.edges[j] = edgemap_check_1.first->second;
				edge_idx_counter++;
			}
			else
			{
				unsigned int edge_index;
				edge_index = edgemap_check_1.first->second;				
				//update the triangle
				MeshEdge& e = edges[edge_index];
				e.triangles[1] = i;
				e.triangle_size++;
				tri.edges[j] = edge_index;
			}
		}
	}

	//clear all the old edges on the old vertices
	for (int i = 0; i < old_vertices_size; i++)
	{
		vertices[i].edges.clear();
	}

	//marking boundary vertices
	for (int i = 0; i < edges.size(); i++)
	{
		const MeshEdge& e = edges[i];
		for (int j = 0; j < 2; j++)
		{
			MeshVertex& v = vertices[e.vertices[j]];
			v.edges.push_back(i);
			if (e.triangle_size == 1)
			{
				//mark as boundary vertex
				v.is_boundary = true;
			}
		}
	}
}

bool Mesh::subdivide()
{
	unsigned int old_vertices_size = vertices.size();
	MeshTriangleList new_triangles;
	MeshEdgeList new_edges;
	new_triangles.reserve(triangles.size() * 4);
	new_edges.reserve(new_triangles.size() * 2);
	assemble_new_triangles(vertices, triangles, new_triangles);

	first_pass(vertices, triangles, edges, new_triangles);
	second_pass(vertices, edges, old_vertices_size);
	//
	triangles = new_triangles;
	recreate_new_edges(vertices, triangles, new_edges, old_vertices_size);
	edges = new_edges;
	compute_normals(vertices, triangles);
	create_gl_data(); // recreate the gl data
    return true;
}

} /* _462 */
