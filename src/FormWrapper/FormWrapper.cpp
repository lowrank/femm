/*
 * FormWrapper.cpp
 *
 *  Created on: Sep 5, 2017
 *      Author: lurker
 */

#include "FormWrapper.h"

FormWrapper::FormWrapper() {
	// TODO Auto-generated constructor stub

}

FormWrapper::~FormWrapper() {
	// TODO Auto-generated destructor stub
}

void FormWrapper::Reference(M_Ptr &R, M_Ptr &Rx, M_Ptr &Ry, M_Ptr points, M_Ptr qpoints){

	auto np = mxGetN(points);
	auto nq = mxGetN(qpoints);

	auto T = dMat(np, np);
	auto F = dMat(np, nq);
	auto Fx = dMat(np, nq);
	auto Fy = dMat(np, nq);

	int deg = round((sqrt(8 * np + 1) - 3)/2);

	assert((deg + 1) * (deg + 2) == np);

	auto Tp = mxGetPr(T);
	auto Fp = mxGetPr(F);
	auto Fxp = mxGetPr(Fx);
	auto Fyp = mxGetPr(Fy);

	auto pp = mxGetPr(points);
	auto qp = mxGetPr(qpoints);

	/*
	 * build Vandermonte matrix, x^(i-j) * y^j
	 */
	for (auto col = 0; col < np;++col) {
		for (auto i = 0; i < deg + 1; ++i) {
			for (auto j = 0; j < i + 1; ++j) {
				*Tp++ = pow(pp[2 * col], i - j) * pow(pp[2 * col + 1], j);
			}
		}
	}

	for (auto col = 0; col < nq;++col) {
		for (auto i = 0; i < deg + 1; ++i) {
			for (auto j = 0; j < i + 1; ++j) {
				*Fp++ = pow(qp[2 * col], i - j) * pow(qp[2 * col + 1], j);
			}
		}
	}

	for (auto col = 0; col < nq;++col) {
		for (auto i = 0; i < deg + 1; ++i) {
			for (auto j = 0; j < i + 1; ++j) {
				if (j == i) {
					*Fxp++ = 0.;
				}
				else {
					*Fxp++ = pow(qp[2 * col], i - j - 1) * pow(qp[2 * col + 1], j);
				}
			}
		}
	}


	for (auto col = 0; col < nq;++col) {
		for (auto i = 0; i < deg + 1; ++i) {
			for (auto j = 0; j < i + 1; ++j) {
				if (j == 0) {
					*Fyp++ = 0.;
				}
				else {
					*Fyp++ = pow(qp[2 * col], i - j) * pow(qp[2 * col + 1], j - 1);
				}
			}
		}
	}

	M_Ptr RHS_f[] = {T, F};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	M_Ptr RHS_x[] = {T, Fx};
	mexCallMATLAB(1, &Rx, 2, RHS_x, "mldivide");

	M_Ptr RHS_y[] = {T, Fy};
	mexCallMATLAB(1, &Ry, 2, RHS_y, "mldivide");


	mxDestroyArray(T);
	mxDestroyArray(F);
	mxDestroyArray(Fx);
	mxDestroyArray(Fy);
}


using namespace mexplus;
template class mexplus::Session<FormWrapper>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<FormWrapper>::create(new FormWrapper()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<FormWrapper>::destroy(input.get(0));
}

MEX_DEFINE(ref)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);

	auto form = Session<FormWrapper>::get(input.get(0));

	auto np = mxGetN(prhs[1]);
	auto nq = mxGetN(prhs[2]);

	plhs[0] = dMat(np, nq);
	plhs[1] = dMat(np, nq);
	plhs[2] = dMat(np, nq);

	form->Reference(plhs[0], plhs[1], plhs[2], C_CAST(prhs[1]), C_CAST(prhs[2]));
}

}

MEX_DISPATCH
