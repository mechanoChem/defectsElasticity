/**
*base class for dislocation representing by force dipole 
*derived from IGA
*/


#ifndef IGA_dislocationClass_h
#define IGA_dislocationClass_h

#include "IGA.h"
#include "model/model_dislocation.h"
using namespace std;
template <int dim>
class IGA_dislocation: public IGA<dim>
{
public:
	/**
	*IGA_dislocation constructor
	*/
  IGA_dislocation (NURBSMesh<dim>& _mesh, parametersClass& _params);
	
	/**
	*IGA_dislocation destructor
	*/
  ~IGA_dislocation();
	
	/**
	*entrance of run IGA_dislocation
	*User may use it directly
	*/
  void run();
	
	/**
	*apply dirchlet boundary condition
	*User may change it accordingly and add neumman boundary condition if needed
	*for applying neumman boundary condition please see example of structure problem  
	*/
  void apply_boundary_conditions();
	
	/**
	*Essential to dislocation represented by force dipole
	*please see mark_plane.cc for specific implementation
	*Mark special plane and defects flags
	*We have 7 slots in planeFlags. 0-5 denotes 6 planes. 6: denote plane inside cell or just cell surface. 
	*planeFlags[i]=1:i:normal vector of plane;  1->negtive half plane; =2->postive plane;i=0-5
	*planeFlags[6]=1 mark plane inside of cell. planeFlags[6]=0 (default) mark plane on cell surface.
	*using Sp_planeQuad_Point to specify the position of plane inside of cell.
	*NOTE: planed marked should be compatible with mesh for accuracy:
	* e.g. the upper half plane and lower half plane may be separated by element edge
	*We have 6 slots in defectFlags. 0-2 for direction of burger vector. 3:TBD 4:edge, 5:screw
	*defectFlags[i]=n: i:direction of burger vector, n:number of half planes overlaped. i<=2;
	*/
  void mark_plane();
	
	/**
	*this function generates vtk file
	*User may use it directly
	*/
  void output (unsigned int _cycle);
	
	/**
	*this function assembles jacobian and righ hand side 
	*use "dislocationModel" defined in "run"
	*User may use it directly
	*/
  void assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end);
  
	/**
	*include model type
	*/
	model_dislocation<Sacado::Fad::DFad<double>,dim>* dislocationModel; 
};

#endif