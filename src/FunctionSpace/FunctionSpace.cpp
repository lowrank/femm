#include <cassert>
#include "FunctionSpace.h"

void FunctionSpace::build(double* p_ptr, int* s_ptr, int* t_ptr, int* e_ptr, int* n_ptr, int np, int ns, int nt, int ne, int _deg){
    // assertion on degree.
    assert(deg > 0 && deg < 20);
    // initialize.
    nodePtr = p_ptr;
    segPtr  = s_ptr;
    triPtr  = t_ptr;
    edgePtr = e_ptr;
    neighPtr= n_ptr;

    deg = _deg;

    this->clear();
    int  dof = (deg + 1) * (deg + 2)/2;
    this->resize(2 * (np  +  (deg - 1) * ne + (dof - 3 * deg) * nt),
           dof * nt,
           (deg + 1) * ns,
           ne);
    // fill neighbors.
    this->neighbors.resize((unsigned long) (3 * nt));

    double P1_Coord_X, P1_Coord_Y, P2_Coord_X,P2_Coord_Y;
    double Delta_X, Delta_Y;

    // copy neighbors. todo: use memcopy.
    for (auto i = 0; i < 3 * nt; ++i) {
        neighbors[i] = neighPtr[i];
    }

    for (auto i = 0; i < np; ++i) {
        nodes[2 * i]     = _Point_X(i);
        nodes[2 * i + 1] = _Point_Y(i);
    }

    for (auto i = 0; i < ne; ++i) {
        P1_Coord_X = _Point_X(_Edge_L(i));
        P1_Coord_Y = _Point_Y(_Edge_L(i));
        P2_Coord_X = _Point_X(_Edge_R(i));
        P2_Coord_Y = _Point_Y(_Edge_R(i));

        Delta_X = (P2_Coord_X - P1_Coord_X)/deg;
        Delta_Y = (P2_Coord_Y - P1_Coord_Y)/deg;

        std::vector<int> EdgePoints((unsigned long) (deg - 1));

        for (auto j = 0; j < deg - 1; ++j) {
            EdgePoints[j] = np + i * (deg - 1) + j;

            nodes[2 * EdgePoints[j]    ] = P1_Coord_X + Delta_X * (j + 1);
            nodes[2 * EdgePoints[j] + 1] = P1_Coord_Y + Delta_Y * (j + 1);
        }

        edges.insert(make_pair(
                to_string(_Edge_L(i)) + "-" + to_string(_Edge_R(i)), EdgePoints
        ));
    }

    unordered_set<string> boundary_set;
    boundary_set.reserve((unsigned long) ns);
    unordered_map<string, int> boundary_elems;

    for (auto index = 0; index < ns; ++index) {
        auto bound_l = to_string(_Seg_R(index));
        auto bound_r = to_string(_Seg_L(index));

        boundary_set.insert(bound_l + "-" + bound_r);
    }

    int counter = 2 * np+ 2 * (deg - 1) * ne;

    for (auto index = 0; index < nt; ++index) {
        elems[dof * index    ] = _Tri_U(index);
        elems[dof * index + 1] = _Tri_V(index);
        elems[dof * index + 2] = _Tri_W(index);

        auto tri_u = to_string(_Tri_U(index));
        auto tri_v = to_string(_Tri_V(index));
        auto tri_w = to_string(_Tri_W(index));

        // edge 0 - 1.
        auto EdgeIterator = edges.find(tri_u + "-" + tri_v);
        if (EdgeIterator != edges.end()) {
            copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                 elems.begin() + dof * index + 3);
        }
        else {
            EdgeIterator = edges.find(tri_v + "-" + tri_u);
            reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                         elems.begin() + dof * index + 3);
        }

        // edge 1 - 2.

        EdgeIterator = edges.find(tri_v + "-" + tri_w);
        if (EdgeIterator != edges.end()) {
            copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                 elems.begin() + dof * index + 3 + (deg - 1));
        }
        else {
            EdgeIterator = edges.find(tri_w + "-" + tri_v);
            reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                         elems.begin() + dof * index + 3 + (deg - 1));
        }

        // edge 2 - 0.
        EdgeIterator = edges.find(tri_w + "-" + tri_u);
        if (EdgeIterator != edges.end()) {
            copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                 elems.begin() + dof * index + 3 + 2 * (deg - 1));
        }
        else {
            EdgeIterator = edges.find(tri_u + "-" + tri_w);
            reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(),
                         elems.begin() + dof * index + 3 + 2 * (deg - 1));
        }

        auto BoundaryIterator = boundary_set.find(tri_u + "-" + tri_v);
        if (BoundaryIterator != boundary_set.end()) {
            boundary_elems.insert(make_pair(tri_u + "-" + tri_v, index));
        }

        BoundaryIterator = boundary_set.find(tri_v + "-" + tri_w);
        if (BoundaryIterator != boundary_set.end()) {
            boundary_elems.insert(make_pair(tri_v + "-" + tri_w, index));
        }

        BoundaryIterator = boundary_set.find(tri_w + "-" + tri_u);
        if (BoundaryIterator != boundary_set.end()) {
            boundary_elems.insert(make_pair(tri_w + "-" + tri_u, index));
        }


        int internal_counter = 0;
        for (auto i = 1; i < deg - 1; ++i) {
            for (auto j = 1; j < deg - i; ++j) {
                int k = deg - i - j;

                elems[dof * index + 3 * deg + internal_counter] = counter / 2;
                internal_counter += 1;

                nodes[counter] =
                        (double)k/(double)deg * _Point_X(_Tri_U(index)) +
                        (double)i/(double)deg * _Point_X(_Tri_V(index)) +
                        (double)j/(double)deg * _Point_X(_Tri_W(index));

                nodes[counter + 1] =
                        (double)k/(double)deg * _Point_Y(_Tri_U(index)) +
                        (double)i/(double)deg * _Point_Y(_Tri_V(index)) +
                        (double)j/(double)deg * _Point_Y(_Tri_W(index));
                counter += 2;
            }
        }
    }

    for (auto index = 0; index < ns; index++) {
        // orientation reversed.
        boundary[(deg + 1) * index + 1] = _Seg_L(index);
        boundary[(deg + 1) * index    ] = _Seg_R(index);

        auto bound_l = to_string(_Seg_R(index));
        auto bound_r = to_string(_Seg_L(index));

        auto EdgeIterator = edges.find(bound_l + "-" + bound_r);
        auto BoundaryElemIterator = boundary_elems.find(bound_l + "-" + bound_r);

        if (EdgeIterator != edges.end()) {
            copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), boundary.begin() + (deg + 1) * index + 2);
        }
        else {
            EdgeIterator = edges.find(bound_r + "-" + bound_l);
            reverse_copy(EdgeIterator->second.begin(), EdgeIterator->second.end(), boundary.begin() + (deg + 1) * index + 2);
        }


        if (BoundaryElemIterator != boundary_elems.end()) {
            boundary_index.push_back(BoundaryElemIterator->second);
        }
        else {
            BoundaryElemIterator = boundary_elems.find(bound_r + "-" + bound_l);
            if (BoundaryElemIterator != boundary_elems.end()) {
                boundary_index.push_back(BoundaryElemIterator->second);
            }
        }

    }
}

void FunctionSpace::clear() {
    nodes.clear(); edges.clear();
    elems.clear(); boundary.clear();
    boundary_index.clear(); neighbors.clear();
}

void FunctionSpace::resize(int nsize, int esize, int bsize, int edsize) {
    nodes.resize((unsigned long) nsize);
    elems.resize((unsigned long) esize);
    boundary.resize((unsigned long) bsize);
    edges.reserve((unsigned long) edsize);
}

FunctionSpace::FunctionSpace() {}

FunctionSpace::~FunctionSpace() {

}

using namespace mexplus;

template class Session<FunctionSpace>;

namespace {
    MEX_DEFINE(new)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 0);
        OutputArguments output(nlhs, plhs, 1);
        output.set(0, Session<FunctionSpace>::create(new FunctionSpace()));
    }

    MEX_DEFINE(delete)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        Session<FunctionSpace>::destroy(input.get(0));
    }

    MEX_DEFINE(build)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 7);
        OutputArguments output(nlhs, plhs, 0);
        auto funcspace = Session<FunctionSpace>::get(input.get(0));
        auto deg       = (int)(*M_Cast<double>(C_CAST(input.get(6))));

        /*
         * build a tTriangleInfo struct.
         */

        auto p_ptr = M_Cast<double>(C_CAST(prhs[1])); auto np = mxGetN(prhs[1]);
        auto s_ptr = M_Cast<int>(C_CAST(prhs[2]));    auto ns = mxGetN(prhs[2]);
        auto t_ptr = M_Cast<int>(C_CAST(prhs[3]));    auto nt = mxGetN(prhs[3]);
        auto e_ptr = M_Cast<int>(C_CAST(prhs[4]));    auto ne = mxGetN(prhs[4]);
        auto n_ptr = M_Cast<int>(C_CAST(prhs[5]));

        funcspace->build(p_ptr, s_ptr, t_ptr, e_ptr, n_ptr, np, ns, nt, ne, deg);

    }

    MEX_DEFINE(getData)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 5);
        auto funcspace = Session<FunctionSpace>::get(input.get(0));

        int deg = funcspace->deg;

        plhs[0] = mxCreateNumericMatrix(2, funcspace->nodes.size()/2, mxDOUBLE_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[0]), &funcspace->nodes[0], funcspace->nodes.size() * sizeof(double) );

        plhs[1] = mxCreateNumericMatrix((deg + 1) * (deg + 2)/2, funcspace->elems.size()/((deg + 1) * (deg + 2)/2), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[1]), &funcspace->elems[0], funcspace->elems.size() * sizeof(int));
        M_Ptr temp_1[] = {plhs[1],mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[1], 2, temp_1, "plus");

        plhs[2] = mxCreateNumericMatrix((deg + 1), funcspace->boundary.size()/(deg + 1), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[2]), &funcspace->boundary[0], funcspace->boundary.size() * sizeof(int));
        M_Ptr temp_2[] = {plhs[2],mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[2], 2, temp_2, "plus");

        plhs[3] = mxCreateNumericMatrix(1, funcspace->boundary_index.size(), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[3]), &funcspace->boundary_index[0], funcspace->boundary_index.size() * sizeof(int));
        M_Ptr temp_3[] = {plhs[3],mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[3], 2, temp_3, "plus");

        plhs[4] = mxCreateNumericMatrix(3, funcspace->neighbors.size()/3, mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[4]), &funcspace->neighbors[0], funcspace->neighbors.size() * sizeof(int));
        M_Ptr temp_4[] = {plhs[4], mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[4], 2, temp_4, "plus");

    }



}

MEX_DISPATCH

