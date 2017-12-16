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

/*
 * boundary is all lazy evaluated.
 */
class Boundary {
public:
	Boundary(M_Ptr);
	virtual ~Boundary();
	std::string id;
	vector<std::string> expr;
};

#endif /* BOUNDARY_PRIVATE_BOUNDARY_H_ */
