#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <string>
// Domain 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
// Polyhedron type
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
using namespace std;
int main(int argc, char *argv[])
{
  // Load a polyhedron
	//std::cerr << "One" << endl;
  Polyhedron poly;
  std::ifstream input(argv[1]);
  input >> poly;

	//std::cerr << "Two" << endl;
  Polyhedron rest;
  std::ifstream rests(argv[2]);
  rests >> rest;

	//std::cerr << "Three" << endl;
	string output = string(argv[3]);
  // Create a vector with only one element: the pointer to the polyhedron.
  std::vector<Polyhedron*> poly_ptrs_vector;
	//std::cerr << "Four" << endl;
  poly_ptrs_vector.push_back(&poly);
	//std::cerr << "Five" << endl;
  poly_ptrs_vector.push_back(&rest);
	//std::cerr << "Six" << endl;
  // Create a polyhedral domain, with only one polyhedron,
  // and no "bounding polyhedron", so the volumetric part of the domain will be
  // empty.
  Mesh_domain domain(poly_ptrs_vector.begin(), poly_ptrs_vector.end());
	//std::cerr << "Seven" << endl;
  
  // Get sharp features
 // domain.detect_features(); //includes detection of borders
  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.025,
                         facet_angle = 25,
                         facet_size = 0.1,
                         facet_distance = 0.01);
	//std::cerr << "Eight" << endl;

  //Mesh_criteria criteria(edge_size = 0.1);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());
	//std::cerr << "Nine" << endl;
  // Output the facets of the c3t3 to an OFF file. The facets will not be
  // oriented.
  std::ofstream off_file((output + "out.off").c_str());
	//std::cerr << "Ten" << endl;
  c3t3.output_boundary_to_off(off_file);
	//std::cerr << "Eleven" << endl;
  int val = off_file.bad() ? EXIT_FAILURE : EXIT_SUCCESS;
	//std::cerr << "End Val: " << val << endl;
	return val;
}

