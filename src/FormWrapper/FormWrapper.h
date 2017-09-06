/*
 * FormWrapper.h
 *
 *  Created on: Sep 5, 2017
 *      Author: lurker
 */

#ifndef FormWrapper_H
#define FormWrapper_H

#include <cassert>
#include "mexplus.h"
#include "common.h"


class FormWrapper {
public:
	FormWrapper();
	virtual ~FormWrapper();

	int qdeg;

	/*
	 * 2D and 1D references,
	 */
	void Reference(M_Ptr &R, M_Ptr &Rx, M_Ptr &Ry, M_Ptr points, M_Ptr qpoints);
	void Reference(M_Ptr &R, M_Ptr points, M_Ptr qpoints);

};

#endif /* FormWrapper.h_H */
