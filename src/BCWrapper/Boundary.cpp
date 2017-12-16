/*
 * Boundary.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: lurker
 */

#include "Boundary.h"

Boundary::Boundary(M_Ptr name) {
	id = mxArrayToString(name);
}

Boundary::~Boundary() {
	expr.clear();

#ifdef DEBUG
	mexPrintf("Boundary detached\n");
#endif
}
template class mexplus::Session<Boundary>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);

	output.set(0,
			Session<Boundary>::create(new Boundary(
			C_CAST(input.get(0))
	)));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs,
		const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Boundary>::destroy(input.get(0));
}

MEX_DEFINE(push_expr) (int nlhs, mxArray* plhs[], int nrhs,
		const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 2);
	OutputArguments output(nlhs, plhs, 0);
	auto bc = Session<Boundary>::get(input.get(0));

	std::string expr_string(mxArrayToString(prhs[1]));
	bc->expr.push_back(expr_string);
}

MEX_DEFINE(filter_segments) (int nlhs, mxArray* plhs[], int nrhs,
		const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	auto bc = Session<Boundary>::get(input.get(0));

	OutputArguments output(nlhs, plhs, (int) (bc->expr.size()));

	auto edges_ptr = (int*) mxGetPr(C_CAST(prhs[1]));
	auto nodes_ptr = mxGetPr(C_CAST(prhs[2]));

	auto _num_edge = mxGetN(prhs[1]);
	auto _num_node = mxGetM(prhs[1]);

	auto _num_expr = bc->expr.size();

	if (_num_node < 2) {
		mexErrMsgTxt(
				"Error: Boundary:filter_segments::Each edge should contain at least 2 nodes.\n");
	}

	if (_num_expr < 1) {
		mexErrMsgTxt(
				"Error: Boundary:filter_segments::boundary filter cannot be null");
	}

	/*
	 * stores edges for each boundary
	 */
	vector<vector<int32_t>> _global;
	_global.resize(_num_expr);

	exprtk::symbol_table<double> symbol_table;

	double x = 0., y = 0.;

	symbol_table.add_variable("x", x);
	symbol_table.add_variable("y", y);
	symbol_table.add_constants();

	vector<exprtk::expression<double>> expressions(_num_expr);
	vector<exprtk::parser<double>> parsers(_num_expr);

	for (size_t i = 0; i < _num_expr; i++) {
		expressions[i].register_symbol_table(symbol_table);
		parsers[i].compile(bc->expr[i], expressions[i]);
	}

	double _x_1, _y_1, _x_2, _y_2, result_1, result_2;

	for (size_t i = 0; i < _num_edge; i++) {

		_x_1 = nodes_ptr[2 * (edges_ptr[i * _num_node] - 1)];
		_y_1 = nodes_ptr[2 * (edges_ptr[i * _num_node] - 1) + 1];
		_x_2 = nodes_ptr[2 * (edges_ptr[i * _num_node + 1] - 1)];
		_y_2 = nodes_ptr[2 * (edges_ptr[i * _num_node + 1] - 1) + 1];

		for (size_t j = 0; j < _num_expr; j++) {
			x = _x_1;
			y = _y_1;
			result_1 = fabs(expressions[j].value());
			x = _x_2;
			y = _y_2;
			result_2 = fabs(expressions[j].value());

			if (result_1 < MEX_EPS && result_2 < MEX_EPS) {
				for (size_t k = 0; k < _num_node; k++) {
					_global[j].push_back(edges_ptr[i * _num_node + k]);
				}
			} else {
				continue;
			}
		}
	}
	for (size_t i = 0; i < _num_expr; i++) {
		plhs[i] = mxCreateNumericMatrix(_num_node,
				_global[i].size() / _num_node, mxINT32_CLASS, mxREAL);
		// implicit convert pointer
		memcpy(mxGetPr(plhs[i]), &_global[i][0],
				sizeof(int32_t) * _global[i].size());
	}
	_global.clear();
	expressions.clear();
	parsers.clear();
}


}



MEX_DISPATCH
