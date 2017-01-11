#include <iostream>
#include <igl/facet_components.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOFF.h>
#include <Eigen/Core>
#include <vector>
#include <map>
#include <set>

void extractComponent(const std::string &s, const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F, const std::vector<int> &facets)
{
	std::set<int> usedVerts;
	for(int i=0; i<facets.size(); i++)
	{
		for(int j=0; j<3; j++)
		{
			usedVerts.insert(F(facets[i],j));
		}
	}
	int nverts = usedVerts.size();

	Eigen::MatrixX3d cV(nverts, 3);
	int idx=0;
	std::map<int, int> invmap;
	for(std::set<int>::iterator it = usedVerts.begin(); it != usedVerts.end(); ++it)
	{
		cV.row(idx) = V.row(*it);
		invmap[*it] = idx;
		idx++;
	}
	Eigen::MatrixX3i cF(facets.size(), 3);
	for(int i=0; i<facets.size(); i++)
	{
		for(int j=0; j<3; j++)
			cF(i,j) = invmap[F(facets[i], j)];
	}
	igl::writeOFF(s, cV, cF);
}

int main(int argc, char *argv[])
{
	if(argc != 2)
		return -1;
	Eigen::MatrixX3d V;
	Eigen::MatrixX3i F;
	igl::read_triangle_mesh(argv[1], V, F);
	Eigen::VectorXd C(F.rows());
	igl::facet_components(F, C);
	std::map<int, std::vector<int> > components;
	for(int i=0; i<F.rows(); i++)
		components[C[i]].push_back(i);

	int largestid = -1;
	int largestsize = 0;	
	for(std::map<int, std::vector<int> >::iterator it = components.begin(); it != components.end(); ++it)
	{
		if(it->second.size() > largestsize)
		{
			largestid = it->first;
			largestsize = it->second.size();
		}
	}	
	
	// extract largest component
	extractComponent("largest.off", V, F, components[largestid]);
	std::vector<int> rest;
	for(std::map<int, std::vector<int> >::iterator it = components.begin(); it != components.end(); ++it)
	{
		if(it->first == largestid)
			continue;
		for(int i=0; i<it->second.size(); i++)
			rest.push_back(it->second[i]);
	}
	extractComponent("rest.off", V, F, rest);
}
