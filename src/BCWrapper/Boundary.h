/*
 * Boundary.h
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#ifndef BOUNDARY_PRIVATE_BOUNDARY_H_
#define BOUNDARY_PRIVATE_BOUNDARY_H_

#include <cstdlib>
#include <vector>
#include <cmath>

#include <iostream>
#include <iterator>
#include <string.h>

#include <unordered_set>

#include <mexplus.h>
#include "exprtk.hpp"

#include "common.h"

using namespace std;
using namespace mexplus;


enum Boundary_t {DIRICHLET = 0, NEUMANN, ROBIN};

class Boundary {
public:
	Boundary(M_Ptr);
	/*
	 * Default constructor, placeholder
	 */
	Boundary();
	virtual ~Boundary();

	/*
	 * For the case with placeholder, delayed construction.
	 */
	virtual void setDirichlet(M_Ptr);

	unordered_set<int32_t> b_edge_set;
	vector<std::string> b_expr;
	// apply boundary condition on LHS and RHS

};

class DirichletBC:public Boundary {
public:
	explicit DirichletBC(M_Ptr _edges) : Boundary(_edges) {}
	/*
	 * placeholder, for later use
	 */
	explicit DirichletBC() : Boundary() {}
};

#endif /* BOUNDARY_PRIVATE_BOUNDARY_H_ */
