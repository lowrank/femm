/*
 * Assembler.cc
 *
 *  Created on: Oct 9, 2014
 *      Author: lurker
 */
#include "Assembler.h"

Assembler::Assembler(){


}

Assembler::~Assembler() {
#ifdef DEBUG
	mexPrintf("Assembler detached\n");
#endif
}


/*
 * Extract information of basis on Reference Triangle
 */
void Assembler::Reference(M_Ptr &F, M_Ptr &DX, M_Ptr &DY,
		M_Ptr Points, M_Ptr QPoints){

	auto _numberofpoints = mxGetN(Points);
	auto _numberofqpoints = mxGetN(QPoints);


	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	auto Vander = mxCreateNumericMatrix(_numberofpoints,_numberofpoints,mxDOUBLE_CLASS, mxREAL);
	auto VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderY = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	int deg = round((sqrt(8*_numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != _numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	// extract all nodes from promoted nodes.
	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderF_ptr = mxGetPr(VanderF);
	auto VanderX_ptr = mxGetPr(VanderX);
	auto VanderY_ptr = mxGetPr(VanderY);

	auto nodes_ptr   = mxGetPr(Points);
	auto qnodes_ptr  = mxGetPr(QPoints);

	// basis on all nodes

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*Vander_ptr++ = pow(nodes_ptr[2*col], i - j)*pow(nodes_ptr[2*col + 1], j);
			}
		}
	}



	// basis on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*VanderF_ptr++ = pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j);
			}
		}
	}

	// partial x on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == i){
					*VanderX_ptr++ = 0.;
				}
				else{
					*VanderX_ptr++ = (i - j) * pow(qnodes_ptr[2*col], i - j - 1)*pow(qnodes_ptr[2*col + 1], j);
				}
			}
		}
	}

	// partial y on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == 0){
					*VanderY_ptr++ = 0.;
				}
				else{
					*VanderY_ptr++ = j* pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j - 1);
				}
			}
		}
	}



	mxArray* RHS_f[] = {Vander, VanderF};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	mxArray* RHS_y[] = {Vander, VanderY};
	mexCallMATLAB(1, &DY, 2, RHS_y, "mldivide");


	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);
	mxDestroyArray(VanderY);


}

/*
 * 1D reference
 */
void Assembler::Reference(M_Ptr& F, M_Ptr& DX, M_Ptr Degree, M_Ptr QPoints){
	// Reference segment [-1 , 1],
	// with (Degree - 1) equal-spaced points in between.

	auto _numberofpoints = static_cast<int32_t>(*(mxGetPr(Degree)) + 1);

	vector<double> Points(_numberofpoints);
	Points[0] = -1.;
	Points[1] = 1. ;
	for (size_t i = 2; i < _numberofpoints; i++) {
		Points[i] = (-1.) + 2.0*(i - 1)/static_cast<double>(_numberofpoints - 1);
	}

	// 1D vector-row-majored
	auto _numberofqpoints = mxGetN(QPoints);
	auto qnodes_ptr       = mxGetPr(QPoints);

	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	auto Vander = mxCreateNumericMatrix(_numberofpoints, _numberofpoints ,mxDOUBLE_CLASS, mxREAL);
	auto VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	// extract all nodes from promoted nodes.
	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderF_ptr = mxGetPr(VanderF);
	auto VanderX_ptr = mxGetPr(VanderX);

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < _numberofpoints; i++){
			*Vander_ptr++ = pow(Points[col], i);
		}
	}

	for (size_t col = 0; col < _numberofqpoints; col++) {
		for (size_t i = 0; i < _numberofpoints; i++) {
			*VanderF_ptr++ = pow(qnodes_ptr[col], i);
		}
	}

	for (size_t col = 0; col < _numberofqpoints; col++) {
		for (size_t i = 0 ; i < _numberofpoints ; i++) {
			if (i == 0) {
				*VanderX_ptr++ = 0.;
			}
			else{
				*VanderX_ptr++ = i * pow(qnodes_ptr[col], i - 1);
			}
		}
	}

	mxArray* RHS_f[] = {Vander, VanderF};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	Points.clear();
	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);

}

void Assembler::AssembleBC(double* &pI, double* &pJ, double* &pV,
		M_Ptr Nodes,M_Ptr eRobin,
		M_Ptr Ref, M_Ptr Weights, M_Ptr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(eRobin);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofedges          = mxGetN(eRobin);
	auto numberofnodesperedge   = mxGetM(eRobin);
	auto numberofqnodes         = mxGetN(Ref);

	mwSize vertex_1, vertex_2;
	double length;

	/*
	 * More codes but there is only one judgment at beginning.
	 * Performance does not drop.
	 */
	if (mxGetNumberOfElements(Fcn) == numberofedges*numberofqnodes){

		for (size_t i = 0; i < numberofedges; i++) {


			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;


			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				for (size_t k = 0; k < numberofnodesperedge; k++) {
					*pI++ = pedge_ptr[i*numberofnodesperedge + j];
					*pJ++ = pedge_ptr[i*numberofnodesperedge + k];
					*pV = 0.;

					for (size_t l = 0; l < numberofqnodes; l++) {
						*pV = *pV + Interp[i*numberofqnodes + l] *
								reference[j + l*numberofnodesperedge] *
								reference[k + l*numberofnodesperedge] *
								weights[l];
					}
					*pV *= length/2.0; pV++;
				}
			}
		}
	}
	else {
		for (size_t i = 0; i < numberofedges; i++) {


			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;


			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				for (size_t k = 0; k < numberofnodesperedge; k++) {
					*pI++ = pedge_ptr[i*numberofnodesperedge + j];
					*pJ++ = pedge_ptr[i*numberofnodesperedge + k];
					*pV = 0.;

					for (size_t l = 0; l < numberofqnodes; l++) {
						*pV = *pV +
								reference[j + l*numberofnodesperedge] *
								reference[k + l*numberofnodesperedge] *
								weights[l];
					}
					*pV *= *(Interp)*length/2.0; pV++;
				}
			}
		}
	}
}


void Assembler::AssembleBC(double*& pNeumann, M_Ptr Nodes,
		M_Ptr QNodes, M_Ptr eNeumann,
		M_Ptr Ref, M_Ptr Weights,  M_Ptr Fcn) {
	// Fcn needs to be same size as eNeumann * qnodes

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(eNeumann);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofedges          = mxGetN(eNeumann);
	auto numberofnodesperedge   = mxGetM(eNeumann);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2;
	double length, tmp;

	if (mxGetNumberOfElements(Fcn)  == numberofedges*numberofqnodes) {
		/*
		 * Fcn is a matrix.
		 */
		for (size_t i = 0; i < numberofedges; i++) {
			// integral over each interval
			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;

			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				// each basis
				tmp = 0.;
				for (size_t l = 0 ; l < numberofqnodes; l++) {
					tmp += Interp[i*numberofqnodes + l] * reference[j + l*numberofnodesperedge] * weights[l];
				}
				pNeumann[pedge_ptr[i*(numberofnodesperedge) + j] - 1] += tmp * length/2.;
			}
		}
	}
	else {
		/*
		 * Fcn is a number
		 */
		for (size_t i = 0; i < numberofedges; i++) {
			// integral over each interval
			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;

			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				// each basis
				tmp = 0.;
				for (size_t l = 0 ; l < numberofqnodes; l++) {
					tmp += reference[j + l*numberofnodesperedge] * weights[l];
				}
				pNeumann[pedge_ptr[i*(numberofnodesperedge) + j] - 1] += *(Interp)*tmp * length/2.;
			}
		}
	}
}

void Assembler::AssembleMass(double* &pI, double* &pJ, double* &pV,
		M_Ptr Nodes, M_Ptr Elems,M_Ptr Ref,
		M_Ptr Weights, M_Ptr Fcn, M_Ptr A){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);

	auto area_ptr               = mxGetPr(A);

	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;


	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){
		for (size_t i =0; i < numberofelem; i++){

			area = area_ptr[i];

			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l]*
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = (*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
	else if (mxGetNumberOfElements(Fcn) == numberofelem) {
		for (size_t i = 0; i < numberofelem; i++) {
			area = area_ptr[i];

			auto r = Interp[i];
			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = r*(*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){

			area = area_ptr[i];
			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = *(Interp)*(*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
}

void Assembler::Qnodes2D(double*& Coords, M_Ptr Nodes, M_Ptr QNodes, M_Ptr Elems){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(QNodes);


	mwSize vertex_1, vertex_2 , vertex_3;

	for (size_t i = 0; i < numberofelem; i++) {

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		for (size_t l = 0; l < numberofqnodes; l++) {
			Coords[2*(i*numberofqnodes + l)] =
					pnodes_ptr[2 * vertex_1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3]*qnodes_ptr[2*l + 1];
			Coords[2*(i*numberofqnodes + l) + 1] =
					pnodes_ptr[2 * vertex_1 + 1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2 + 1]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3 + 1]*qnodes_ptr[2*l + 1];
		}
	}
}


void Assembler::Qnodes1D(double*& Coords, M_Ptr Nodes, M_Ptr QNodes, M_Ptr Edges){
	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(Edges);

	auto numberofedges          = mxGetN(Edges);
	auto numberofnodesperedge   = mxGetM(Edges);
	auto numberofqnodes         = mxGetN(QNodes);

	mwSize vertex_1, vertex_2;

	if (mxGetM(Nodes) == 1) {
		// 1D problem
		for (size_t i = 0; i < numberofedges; i++) {
			vertex_1 = pedge_ptr[numberofnodesperedge*i] - 1;
			vertex_2 = pedge_ptr[numberofnodesperedge*i + 1] - 1;

			for (size_t l = 0; l < numberofqnodes; l++) {
				Coords[(i*numberofqnodes + l)] =
						pnodes_ptr[vertex_1] + pnodes_ptr[vertex_2] +
						(pnodes_ptr[vertex_2] - pnodes_ptr[vertex_1])*qnodes_ptr[l];

			}
		}
	}
	else {
		// 2D problem
		for (size_t i = 0; i < numberofedges; i++) {

			vertex_1 = pedge_ptr[numberofnodesperedge*i] - 1;
			vertex_2 = pedge_ptr[numberofnodesperedge*i + 1] - 1;

			for (size_t l = 0; l < numberofqnodes; l++) {
				Coords[2*(i*numberofqnodes + l)] =
						(pnodes_ptr[2*vertex_1] + pnodes_ptr[2*vertex_2] +
						(pnodes_ptr[2*vertex_2] - pnodes_ptr[2*vertex_1])*qnodes_ptr[l])/2.0;
				Coords[2*(i*numberofqnodes + l) + 1] =
						(pnodes_ptr[2*vertex_1 + 1] + pnodes_ptr[2*vertex_2 + 1] +
						(pnodes_ptr[2*vertex_2 + 1] - pnodes_ptr[2*vertex_1 + 1])*qnodes_ptr[l])/2.0;
			}
		}
	}

}



// calculate integral on boundary
void Assembler::AssembleLoad(double*& pLoad, M_Ptr Nodes,
		M_Ptr QNodes, M_Ptr Elems,M_Ptr Ref,
		M_Ptr Weights, M_Ptr Fcn, M_Ptr A) {

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);

	auto area_ptr               = mxGetPr(A);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area, tmp;


	// Fcn cannot be a function handle, too slow

	/* @Revised: Fcn can be a function handle for one pass.
	 * Which means extra space to store all the values. However,
	 * we did not do it because Matlab can handle this easily.
	 */
	auto Fcn_ptr = mxGetPr(Fcn);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		area = area_ptr[i];

		if(mxGetNumberOfElements(Fcn) == numberofqnodes * numberofelem) {
			// Fcn has numberofqnodes * numberofelem elements
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Fcn_ptr[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp*area;
			}
		}
		else if (mxGetNumberOfElements(Fcn) == 1) {
			// Fcn is constant
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp += reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += *(Fcn_ptr)*tmp*area;
			}
		}
		else {
			// Other cases does not match
			mexErrMsgTxt("Error:Assembler:AssembleLoad::Failed with unexpected Fcn.\n");
		}
	}
}


void Assembler::AssembleStiff(double* &pI, double* &pJ, double*&pV,
		M_Ptr Nodes, M_Ptr Elems, M_Ptr RefX,
		M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);


	mwSize vertex_1, vertex_2, vertex_3;
	double det, area;
	double Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn)  == numberofelem*numberofqnodes) {
		// Fcn is a matrix
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else if (mxGetNumberOfElements(Fcn)  == numberofelem){
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			auto r = Interp[i];
			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = r*(*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = *(Interp)*(*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
}


/*
 * Stiff Kernel as two functions
 */

void Assembler::AssembleStiff(double* &pI, double* &pJ, double*&pV,
		M_Ptr Nodes, M_Ptr Elems, M_Ptr RefX,
		M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn_X, M_Ptr Fcn_Y) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp_X             = mxGetPr(Fcn_X);
	auto  Interp_Y             = mxGetPr(Fcn_Y);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);


	mwSize vertex_1, vertex_2, vertex_3;
	double det, area;
	double Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn_X)  == numberofelem * numberofqnodes &&
			mxGetNumberOfElements(Fcn_Y) == numberofelem * numberofqnodes) {
		// Fcn is a matrix
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								Interp_X[i*numberofqnodes + l] *
								(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								)
								+
								Interp_Y[i*numberofqnodes + l] *
								(
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								(
								*(Interp_X)*
								(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								)
								+
								*(Interp_Y)*
								(
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
}


void Assembler::AssembleGradXFunc(double* &pI, double* &pJ, double* &pV,
		M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX, M_Ptr RefY,
		M_Ptr Weights, M_Ptr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;
	double Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];
					}
					*pV  = *(Interp)*(*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}
}

void Assembler::AssembleGradYFunc(double* &pI, double* &pJ, double* &pV,
		M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX,M_Ptr RefY,
		M_Ptr Weights, M_Ptr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;
	double Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = *(Interp)*(*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}
}


void Assembler::AssembleGradXYFunc(double* &pI, double* &pJ, double* &pV,double* &pW,
		M_Ptr Nodes, M_Ptr Elems, M_Ptr Ref, M_Ptr RefX,M_Ptr RefY,
		M_Ptr Weights, M_Ptr Fcn_X, M_Ptr Fcn_Y){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp_X             = mxGetPr(Fcn_X);
	auto  Interp_Y             = mxGetPr(Fcn_Y);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area;
	double Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn_X) == numberofelem*numberofqnodes &&
			mxGetNumberOfElements(Fcn_Y) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					*pW = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){


						*pV = *pV + Interp_X[i*numberofqnodes + l] *
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];

						*pW = *pW + Interp_Y[i*numberofqnodes + l] *
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];

					}
					*pV  = (*pV)/2.0;
					*pW  = (*pW)/2.0;
					pI++; pJ++; pV++;pW++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					*pW = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){


						*pV = *pV +
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];

						*pW = *pW +
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];

					}
					*pV  = *(Interp_X)*(*pV)/2.0;
					*pW  = *(Interp_Y)*(*pW)/2.0;
					pI++; pJ++; pV++;pW++;
				}
			}
		}//end for
	}
}//end all

// build matrix form of load vector.
void Assembler::AssembleLoadMatrix(double*& pI, double*& pJ, double*& pV, M_Ptr Nodes,
		M_Ptr Elems,M_Ptr Ref,
		M_Ptr Weights, M_Ptr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
//	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	double det, area, tmp;


	auto Fcn_ptr = mxGetPr(Fcn);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		if(mxGetNumberOfElements(Fcn) == numberofqnodes * numberofelem) {

			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Fcn_ptr[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem]*weights[l];
				}

				// which node
				*pI = pelem_ptr[i*numberofnodesperelem + j];
				// which element
				*pJ = (i + 1);
				// value
				*pV = tmp * area;

				pI++;pJ++;pV++;

			}
		}
		else {
			mexErrMsgTxt("AssemblerExtension::AssembleLoadMatrix::Dimension does not match.\n");
		}
	}
}


void Assembler::AssembleOverElement(double*& w, M_Ptr Nodes, M_Ptr Elems,
		M_Ptr Ref, M_Ptr RefX,
		M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn_S, M_Ptr Fcn_A, M_Ptr u,
		M_Ptr v){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);

	auto  Interp_S             = mxGetPr(Fcn_S);
	auto  Interp_A             = mxGetPr(Fcn_A);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);
	auto ptru                   = mxGetPr(u);
	auto ptrv                   = mxGetPr(v);


	mwSize vertex_1, vertex_2, vertex_3;
	double det, area;
	double Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn_S)  == numberofelem &&
			mxGetNumberOfElements(Fcn_A) == numberofelem) {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);


			auto a = Interp_A[i];
			auto s = Interp_S[i];

			int32_t I, J;
			double K1, K2, K;


			w[i] = 0;
			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					I = pelem_ptr[i*numberofnodesperelem + j] - 1;
					J = pelem_ptr[i*numberofnodesperelem + k] - 1;
					K1 = 0.;
					K2 = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						K1 = K1 + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
						K2 = K2 +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					K1 = s*(K1)/4.0/area;
					K2 = a*(K2) * area;
					K  = K1 + K2;

					w[i] += K * ptru[I] * ptrv[J];
					if (j != k){
						w[i] += K * ptru[J] * ptrv[I];
					}
				}
			}
		}
	}
}


void Assembler::AssembleOverNode(double*& w, M_Ptr Nodes, M_Ptr Elems,
		M_Ptr Ref, M_Ptr RefX,
		M_Ptr RefY, M_Ptr Weights, M_Ptr Fcn_S, M_Ptr Fcn_A,
		M_Ptr u, M_Ptr v) {

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);

	auto  Interp_S             = mxGetPr(Fcn_S);
	auto  Interp_A             = mxGetPr(Fcn_A);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofnodes          = mxGetN(Nodes);
	auto numberofqnodes         = mxGetN(RefX);
	auto ptru                   = mxGetPr(u);
	auto ptrv                   = mxGetPr(v);


	mwSize vertex_1, vertex_2, vertex_3;
	double det, area;
	double Jacobian[2][2];

	// for performance purpose, not to use fill
	memset(w, 0,sizeof(double) * numberofnodes);

	if (mxGetNumberOfElements(Fcn_S)  == numberofnodes &&
			mxGetNumberOfElements(Fcn_A) == numberofnodes) {

		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			int32_t I, J, K;
			double K1, K2, K3;

			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					I = pelem_ptr[i*numberofnodesperelem + j] - 1;
					J = pelem_ptr[i*numberofnodesperelem + k] - 1;

					/*
					 * calculate (I,J) element
					 */

					/*
					 * (I,J) element uses all qnodes.
					 */
					for (size_t l = 0; l < numberofqnodes; l++) {
						/*
						 * l-th component
						 */
						K1 =(
							(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
							 Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] +
							 Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
							+
							(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
							 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] +
							 Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
							)*weights[l];
						K2 = reference[j+ l*numberofnodesperelem]*
							 reference[k+ l*numberofnodesperelem]*
							 weights[l];
						for (size_t m = 0; m < numberofnodesperelem; m++) {
							// m th node in this element
							K = pelem_ptr[numberofnodesperelem * i + m] - 1;
							K3  = reference[m + l * numberofnodesperelem] *
									(Interp_S[K] * (K1/4.0/area) + Interp_A[K] * K2 * area);
							w[K] += K3 * ptru[I] * ptrv[J];
							if (j != k) {
								w[K] += K3 * ptru[J] * ptrv[I];
							}
						}// end of m
					}// end of l
				}// end of k
			}//end of j
		}//end of i
	}
}


template class mexplus::Session<Assembler>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Assembler>::create(new Assembler()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Assembler>::destroy(input.get(0));
}

MEX_DEFINE(reference2D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofpoints = mxGetN(prhs[1]);
	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);

	assembler->Reference(plhs[0], plhs[1], plhs[2], C_CAST(prhs[1]), C_CAST(prhs[2]));


}

MEX_DEFINE(reference1D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 2);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofpoints = static_cast<size_t>(*mxGetPr(C_CAST(prhs[1])) + 1);

	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	assembler->Reference(plhs[0], plhs[1], C_CAST(prhs[1]), C_CAST(prhs[2]));

}

MEX_DEFINE(assems)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler->AssembleStiff(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]));
}

// advanced interface with matrix kernel
MEX_DEFINE(assemsm)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler->AssembleStiff(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]));
}

MEX_DEFINE(assema)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler->AssembleMass(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),  C_CAST(prhs[6]));


}

MEX_DEFINE(asseml) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	double* pLoad = mxGetPr(plhs[0]);
	assembler->AssembleLoad(pLoad, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]));

}

MEX_DEFINE(qnodes2D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofelem           = mxGetN(prhs[3]);
	size_t numberofqnodes         = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(2, numberofelem * numberofqnodes,  mxDOUBLE_CLASS, mxREAL);
	double* Coords = mxGetPr(plhs[0]);
	assembler->Qnodes2D(Coords,C_CAST(prhs[1]), C_CAST(prhs[2]),
			C_CAST(prhs[3]));
}

MEX_DEFINE(qnodes1D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofedges          = mxGetN(prhs[3]);
	size_t numberofqnodes         = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[1]), numberofedges * numberofqnodes,  mxDOUBLE_CLASS, mxREAL);
	double* Coords = mxGetPr(plhs[0]);
	assembler->Qnodes1D(Coords,C_CAST(prhs[1]), C_CAST(prhs[2]),
			C_CAST(prhs[3]));
}

MEX_DEFINE(assemrbc) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	double* pNeumann = mxGetPr(plhs[0]);
	assembler->AssembleBC(pNeumann, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]));

}

MEX_DEFINE(assemlbc) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes         = mxGetN(prhs[1]);
	size_t numberofedges         = mxGetN(prhs[2]);
	size_t numberofnodesperedges = mxGetM(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);
	double* pJ = mxGetPr(plhs[1]);
	double* pV = mxGetPr(plhs[2]);
	assembler->AssembleBC(pI, pJ, pV,  C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]));

}

MEX_DEFINE(assemex_gradfunc_x) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradXFunc(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]));
}

MEX_DEFINE(assemex_gradfunc_y)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradYFunc(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]));
}

MEX_DEFINE(assemex_gradfunc_xy)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 9);
	OutputArguments output(nlhs, plhs, 4);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	plhs[3] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pW = mxGetPr(plhs[3]);

	assembler_ex->AssembleGradXYFunc(pI, pJ, pV, pW, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]),C_CAST(prhs[8]));
}

MEX_DEFINE(assemex_lm) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	auto numberofnodesperelem = mxGetM(prhs[2]);
	auto numberofelem    =  mxGetN(prhs[2]);


//	M_Ptr Nodes,
//			M_Ptr QNodes, M_Ptr Elems,M_Ptr Ref,
//			M_Ptr Weights, M_Ptr Fcn

	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem , 1,  mxDOUBLE_CLASS, mxREAL);
	double* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem , 1,  mxDOUBLE_CLASS, mxREAL);
	double* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleLoadMatrix(pI, pJ, pV, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]));

}

MEX_DEFINE(assemloe)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 11);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes          = mxGetN(prhs[1]);
	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	double* w = mxGetPr(plhs[0]);

	assembler->AssembleOverElement(w, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]), C_CAST(prhs[8]),
			C_CAST(prhs[9]), C_CAST(prhs[10]));

}


MEX_DEFINE(assemlon)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 11);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes          = mxGetN(prhs[1]);
	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, 1, mxDOUBLE_CLASS, mxREAL);
	double* w = mxGetPr(plhs[0]);

	assembler->AssembleOverNode(w, C_CAST(prhs[1]),
			C_CAST(prhs[2]), C_CAST(prhs[3]),
			C_CAST(prhs[4]), C_CAST(prhs[5]),
			C_CAST(prhs[6]), C_CAST(prhs[7]), C_CAST(prhs[8]),
			C_CAST(prhs[9]), C_CAST(prhs[10]));

}

}


MEX_DISPATCH


