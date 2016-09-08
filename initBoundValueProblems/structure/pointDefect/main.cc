/*
*main file
*parameters for numerical simulation (e.g. solver parameters) and material parameters are defined in one class
*User could defined two separated parameter class for clearity
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <ctime>
#include <vector>

//method: IGA
#include "initBoundValueProb/IGA_structure.h"

//Constant
#define DIMS 3
#define NUM_QUAD_POINTS 5 //NUM_QUAD_POINTS<=5 implemented
#define NUM_THREADS 16



using namespace std;

int main(int argc, char *argv[]){
  std::clock_t start=std::clock();
  dealii::MultithreadInfo::set_thread_limit	(NUM_THREADS);

  //set material parameters genric 
  parametersClass params;
  params.setDouble("lambda", 1.0);
  params.setDouble("mu", 1.0);
  params.setDouble("muSG", 1.0);
  params.setDouble("l", 0.0);
  params.setDouble("load", 3e-6);
  params.setDouble("Gamma", 0.0);
  params.setDouble("C", 5.0);
	
	params.setInt("DOF",0);
  //params.setString("bcType", "line");
  params.setString("bcType", "pointDefect");
  params.setString("order", "Quadratic");
  params.setBool("enforceWeakBC", true);
	params.setBool("finiteStrain", true);
  params.setInt("knots",5);

	//sover parameters
	params.setDouble("res",1);
	params.setDouble("tol",1.0e-12);
	params.setDouble("abs_tol",1.0e-16);
	params.setDouble("initial_norm",0);
	params.setDouble("current_norm",0);
	params.setDouble("max_iteration",3);
	
	//output path default:current folder
	params.setString("output_path", "output0");
	
  printf("reading environmental variables...\n");
  //NURBS file prefix
  char fileName[100];
  std::sprintf (fileName, "%uD%s%u", DIMS, params.getString("order").c_str(), params.getInt("knots"));
  std::string filePrefix(fileName);
		

  char meshFile[100];
	std::sprintf (meshFile, "/home/wzhenlin/workspace/meshes/defectsGradientElasicity/lineDefect/C1/IGAMesh%s.h5", filePrefix.c_str());
  //readHDF5<DIMS>(meshFile, geometry); 

  //Read NURBS geometry  
  nativeNURBSStructure<DIMS> geometry(meshFile);
  //Generate IGA mesh and data structures (control nodes, knot spans, basis functions, etc)
	
  NURBSMesh<DIMS> mesh(geometry, DIMS, NUM_QUAD_POINTS);
	
	//problem & Method
  IGA_structure<DIMS> problem(mesh,params);
  problem.run();
	
  //Stats
  printf ("\nTime taken:%10.2e sec\n", (std::clock()-start)/((double)CLOCKS_PER_SEC));
}

