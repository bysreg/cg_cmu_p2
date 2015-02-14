#include "scene/mesh.hpp"

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

static void first_pass(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles, Mesh::MeshEdgeList& edges)
{
	Mesh::MeshVertexList oddVertices;
	const float a_weight_interior = 3.0f / 8.0f;
	const float b_weight_interior = a_weight_interior;
	const float c_weight = 1.0f / 1.0f;
	const float d_weight = c_weight;
	const float a_weight_boundary = 0.5f;
	const float b_weight_boundary = 0.5f;

	for (int i = 0; i < edges.size(); i++)
	{
		MeshVertex v;
		const MeshEdge& e = edges[i];
		unsigned int a_index = e.vertices[0];
		unsigned int b_index = e.vertices[1];
		const MeshVertex& a = vertices[a_index];
		const MeshVertex& b = vertices[b_index];
		if (e.triangle_size == 1)
		{
			//boundary case
			v.position = a_weight_boundary * a.position + b_weight_boundary * b.position;
		}
		else
		{
			//interior case
			unsigned int c_index, d_index;
			const MeshTriangle& tri_0 = triangles[e.triangles[0]];
			const MeshTriangle& tri_1 = triangles[e.triangles[1]];
			for (int i = 0; i < 3; i++)
			{
				if (tri_0.vertices[i] != a_index && tri_0.vertices[i] != b_index)
				{
					c_index = tri_0.vertices[i];
				}

				if (tri_1.vertices[i] != a_index && tri_1.vertices[i] != b_index)
				{
					d_index = tri_1.vertices[i];
				}
			}
			const MeshVertex& c = vertices[c_index];
			const MeshVertex& d = vertices[d_index];

			v.position = a_weight_interior * a.position + b_weight_interior * b.position + c_weight * c.position + d_weight * d.position;	
		}
		oddVertices.push_back(v);
	}
}

static void second_pass(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles, Mesh::MeshEdgeList& edges)
{
	const float v_weight_boundary = 0.75f;
	const float a_weight_boundary = 1.0f / 0.8f;
	const float b_weight_boundary = a_weight_boundary;	

	for (int i = 0; i < vertices.size(); i++)
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
			//interios case
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

static void assemble_new_triangles(Mesh::MeshVertexList& vertices, Mesh::MeshTriangleList& triangles)
{

}

bool Mesh::subdivide()
{
	first_pass(vertices, triangles, edges);
	second_pass(vertices, triangles, edges);
	assemble_new_triangles(vertices, triangles);
	compute_normals(vertices, triangles);
	create_gl_data(); // recreate the gl data
    return true;
}

} /* _462 */
