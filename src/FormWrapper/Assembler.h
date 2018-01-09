/*
 * Assembler.h
 *
 *  Created on: Oct 10, 2014
 *      Author: lurker
 */

#ifndef ASSEMBLER_PRIVATE_ASSEMBLER_H_
#define ASSEMBLER_PRIVATE_ASSEMBLER_H_

#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <mexplus.h>

#include "common.h"

using namespace std;
using namespace mexplus;


#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

class Assembler {
public:
	Assembler();
	virtual ~Assembler();

	/*
	 * public methods
	 */
	/*
	 * Reference Matrix
	 */
	void Reference(M_Ptr&, M_Ptr&, M_Ptr&,M_Ptr,M_Ptr);
	void Reference(M_Ptr& DX, M_Ptr& DY, M_Ptr Points);
	/*
	 * Mass Matrix
	 */
	void AssembleMass(double*&, double*&, double*&, M_Ptr, M_Ptr,
			M_Ptr, M_Ptr, M_Ptr, M_Ptr);
	/*
	 * Stiffness Matrix
	 */
	// if kernel is a function
	void AssembleStiff(double*&, double*&, double*&,M_Ptr, M_Ptr,
			M_Ptr, M_Ptr, M_Ptr, M_Ptr, M_Ptr Arealist);

	// with a symmetric/diagonal matrix kernel
	void AssembleStiff(double*&, double*&, double*&,M_Ptr, M_Ptr,
			M_Ptr, M_Ptr, M_Ptr, M_Ptr, M_Ptr, M_Ptr);
	/*
	 * Load Vector from Robin Boundary
	 */
	void AssembleLoad(double*& pLoad, M_Ptr Nodes,
			M_Ptr QNodes, M_Ptr Elems,M_Ptr Ref,
			M_Ptr Weights, M_Ptr Fcn, M_Ptr Arealist);

	/*
	 * Integral for point source < Fcn , delta(x)> = Fcn(x) in infinite dimensional space
	 * if x is not a node, it will be some weighted expression(localized) function.
	 */
	void AssembleLoad(double*& pLoad, M_Ptr _point, M_Ptr Fcn);

	// 1d Integral
	void Reference(M_Ptr&, M_Ptr&, M_Ptr, M_Ptr);

	void AssembleBC(double*& pNeumann,  M_Ptr Nodes,
			M_Ptr QNodes, M_Ptr eNeumann,
			M_Ptr Ref, M_Ptr Weights,  M_Ptr Fcn);
	void AssembleBC(double* &pI, double* &pJ, double* &pV,
			M_Ptr Nodes, M_Ptr eRobin,
			M_Ptr Ref, M_Ptr Weights, M_Ptr Fcn);

	// Auxiliary
	void Qnodes2D(double*& Coords, M_Ptr Nodes, M_Ptr QNodes, M_Ptr Elems);
	void Qnodes1D(double*& Coords, M_Ptr Nodes, M_Ptr QNodes, M_Ptr Edges);


	// modules
	/*
	 * integral on int \phi_i_x \phi_j
	 */
	void AssembleGradXFunc(double* &pI, double* &pJ, double* &pV,
			M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX,M_Ptr RefY,
			M_Ptr Weights, M_Ptr Fcn);
	/*
	 * integral on int \phi_i_y \phi_j
	 */
	void AssembleGradYFunc(double* &pI, double* &pJ, double* &pV,
			M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX, M_Ptr RefY,
			M_Ptr Weights, M_Ptr Fcn);
	/*
	 * combine them all
	 */
	void AssembleGradXYFunc(double* &pI, double* &pJ, double* &pV,double* &pW,
			M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX,M_Ptr RefY,
			M_Ptr Weights, M_Ptr Fcn_X, M_Ptr Fcn_Y);



	/*
	 * load matrix
	 */
	void AssembleLoadMatrix(double*& pI, double*& pJ, double*& pV, M_Ptr Nodes,
			 M_Ptr Elems,M_Ptr Ref,
			M_Ptr Weights, M_Ptr Fcn);


	/*
	 * calculate elementwise inner product
	 *
	 * w_K = u' M_{K} v
	 */
	void AssembleOverElement(double*& w, M_Ptr Nodes, M_Ptr Elems,
			M_Ptr Ref, M_Ptr RefX,
			M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn_s, M_Ptr Fcn_a,
			M_Ptr u, M_Ptr v);


	void AssembleOverNode(double*& w, M_Ptr Nodes, M_Ptr Elems,
			M_Ptr Ref, M_Ptr RefX,
			M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn_s, M_Ptr Fcn_a,
			M_Ptr u, M_Ptr v);
	/*
	 * Calculate gradient at each node. It is known that the gradient will suffer from non-uniqueness.
	 * but the error estimate is bounded when high order FEM is used.
	 *
	 * To make things compact, we put grad here in form builder. But it is not a form yet.
	 */

	void CalculateGrad(double*& wx, double *& wy, M_Ptr u, M_Ptr Nodes,
			M_Ptr Elems, M_Ptr DX, M_Ptr DY);

};

#endif /* ASSEMBLER_PRIVATE_ASSEMBLER_C_ */
