#include "inmost.h"
#include <stdio.h>

double MESH_CENTER_X = 3.679;
double MESH_CENTER_Y = 2.746;

using namespace INMOST;
using namespace std;

void make_vec_to_center_tag(Mesh *m)
{
	Tag tagVectortoCenter = m->CreateTag("Vector_to_Center", DATA_REAL, CELL, NONE, 3);
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		double *coords = new double[3];
		cell.Barycenter(coords);
		cell.RealArray(tagVectortoCenter)[0] = MESH_CENTER_X - coords[0];
		cell.RealArray(tagVectortoCenter)[1] = MESH_CENTER_Y - coords[1];
		cell.RealArray(tagVectortoCenter)[2] = 0;
	}
}


double dist_between_nodes(Node first_node, Node second_node) {
	double first_comp = first_node.Coords()[0] - second_node.Coords()[0], sec_comp = first_node.Coords()[1] - second_node.Coords()[1];
	return sqrt(first_comp * first_comp + sec_comp * sec_comp);
}


double mesh_diam(Mesh* m)
{
	double max_diam = 0;
	for (Mesh::iteratorCell icell = m->BeginCell(); icell != m->EndCell(); icell++) {
		Cell cell = icell->getAsCell();
		ElementArray<Node> nodes = cell.getNodes();

		double f_dist = dist_between_nodes(nodes[0], nodes[1]);
		double s_dist = dist_between_nodes(nodes[0], nodes[2]);
		double t_dist = dist_between_nodes(nodes[1], nodes[2]);

		double cell_diam = max(max(f_dist, s_dist), t_dist);
		if (cell_diam > max_diam) {
			max_diam = cell_diam;
		}
	}
	return max_diam;
}

void make_cells_count_tag(Mesh* m)
{
	Tag tagCoord = m->CreateTag("Cells count", DATA_INTEGER, NODE, NONE, 1);
	for (Mesh::iteratorNode inode = m->BeginNode(); inode != m->EndNode(); inode++) {
		Node node = inode->getAsNode();
		node.Integer(tagCoord) = node.getCells().size();
	}
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		cout << "Usage: mesh <mesh.vtk>" << endl;
		return 1;
	}
    Mesh *m = new Mesh;
	m->Load(argv[1]);
	make_cells_count_tag(m);
	m->Save("res.vtk");
	delete m;
	cout << "Success! \n";
	return 0;
}
