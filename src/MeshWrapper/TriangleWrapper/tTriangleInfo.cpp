//
// Created by lurker on 4/14/17.
//


#include "tTriangleInfo.h"
#include "mexplus.h"
#include "common.h"

tTriangleInfo::tTriangleInfo() {
    /*
     * initialize all information on meta.
     *
     *  primitive types.
     */
    _meta.numberoftriangles = 0;
    _meta.numberofpointattributes = 0;
    _meta.numberofsegments = 0;
    _meta.numberoftriangleattributes = 0;
    _meta.numberofcorners = 0;
    _meta.numberofregions = 0;
    _meta.numberofedges = 0;
    _meta.numberofholes = 0;
    _meta.numberofpoints = 0;


    /*
     * initialize all pointers.
     */
    _meta.pointlist = nullptr;
    _meta.pointattributelist = nullptr;
    _meta.pointmarkerlist = nullptr;

    _meta.trianglelist = nullptr;
    _meta.triangleattributelist = nullptr;
    _meta.trianglearealist = nullptr;

    _meta.segmentlist = nullptr;
    _meta.segmentmarkerlist = nullptr;

    _meta.edgelist = nullptr;
    _meta.edgemarkerlist = nullptr;

    _meta.regionlist = nullptr;
    _meta.holelist = nullptr;
    _meta.normlist = nullptr;
    _meta.neighborlist = nullptr;
}

tTriangleInfo::~tTriangleInfo() {
    /*
     * release all information of the pointers.
     *
     * todo: use unique_ptr to describe the mesh. [minor]
     *
     */
    trifree((void*)_meta.pointlist);
    trifree((void*)_meta.pointattributelist);
    trifree((void*)_meta.pointmarkerlist);

    trifree((void*)_meta.trianglelist);
    trifree((void*)_meta.triangleattributelist);

    trifree((void*)_meta.segmentlist);
    trifree((void*)_meta.segmentmarkerlist);;

    trifree((void*)_meta.edgelist);
    trifree((void*)_meta.edgemarkerlist);

    trifree((void*)_meta.regionlist);
    trifree((void*)_meta.holelist);
    trifree((void*)_meta.normlist);
    trifree((void*)_meta.neighborlist);
}

void tTriangleInfo::set_points(Array<double> &_points) {
    auto len = _points.get_size();
    assert(len % 2 == 0);
    /*
     * memory is handled by parent.
     */
    _meta.pointlist = (double*) malloc(len * sizeof(double));
    assert(_meta.pointlist != nullptr);
    _meta.numberofpoints = (int) (len / 2);
    memcpy(_meta.pointlist, _points.get_data(), len * sizeof(double));

}

void tTriangleInfo::set_facets(Array<int> &_facets) {
    auto len = _facets.get_size();
    assert(len % 2 == 0);
    /*
     * memory is handled by parent.
     */
    _meta.segmentlist = (int*) malloc(len * sizeof(int));
    assert(_meta.segmentlist != nullptr);
    _meta.numberofsegments = (int) (len / 2);
    memcpy(_meta.segmentlist, _facets.get_data(), len * sizeof(int));
}

void tTriangleInfo::build(std::string switches, tTriangleInfo& out) {
    tTriangleInfo vor;
    triangulate((char*)switches.c_str(), &_meta, &out._meta, &vor._meta);

}


void tTriangleInfo::refine(std::string switches, tTriangleInfo &out) {
    tTriangleInfo vor;
    triangulate((char*)switches.c_str(), &_meta, &out._meta, &vor._meta);
}

using namespace mexplus;

template class Session<tTriangleInfo>;

namespace {
    MEX_DEFINE(new)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 0);
        OutputArguments output(nlhs, plhs, 1);

        output.set(0, Session<tTriangleInfo>::create(new tTriangleInfo()));
    }

    MEX_DEFINE(delete)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        Session<tTriangleInfo>::destroy(input.get(0));
    }

    MEX_DEFINE(verionInfo_tri) (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 0);
        OutputArguments output(nlhs, plhs, 0);

        mexPrintf("tTriangleInfo : Wrapper of triangle library.\n");
        mexPrintf("Version: 0.1\n");
    }

    MEX_DEFINE(getInfo_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        auto mesh = Session<tTriangleInfo>::get(input.get(0));

        mexPrintf("tTriangleInfo: \n");
        mexPrintf("\t number of points: %d\n", mesh->_meta.numberofpoints);
        mexPrintf("\t number of segments: %d\n", mesh->_meta.numberofsegments);
        mexPrintf("\t number of triangles: %d\n", mesh->_meta.numberoftriangles);
        mexPrintf("\t number of edges: %d\n", mesh->_meta.numberofedges);

    }

    MEX_DEFINE(getData_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 5);

        auto mesh = Session<tTriangleInfo>::get(input.get(0));

        plhs[0] = mxCreateNumericMatrix(2, (int32_t)(mesh->_meta.numberofpoints), mxDOUBLE_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[0]), mesh->_meta.pointlist,2 *  (int32_t)(mesh->_meta.numberofpoints) * sizeof(double));


        plhs[1] = mxCreateNumericMatrix(2, (int32_t)(mesh->_meta.numberofsegments), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[1]), mesh->_meta.segmentlist, 2 * (int32_t)(mesh->_meta.numberofsegments) * sizeof(int32_t));


        M_Ptr temp_1[] = {plhs[1], mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[1], 2, temp_1, "plus");

        plhs[2] = mxCreateNumericMatrix(3, (int32_t)(mesh->_meta.numberoftriangles), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[2]), mesh->_meta.trianglelist, 3 * (int32_t)(mesh->_meta.numberoftriangles) * sizeof(int32_t));

        M_Ptr temp_2[] = {plhs[2], mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[2], 2, temp_2, "plus");

        plhs[3] = mxCreateNumericMatrix(2, (int32_t)(mesh->_meta.numberofedges), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[3]), mesh->_meta.edgelist, 2 * (int32_t)(mesh->_meta.numberofedges) * sizeof(int32_t));

        M_Ptr temp_3[] = {plhs[3], mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[3], 2, temp_3, "plus");

        plhs[4] = mxCreateNumericMatrix(3, (int32_t)(mesh->_meta.numberoftriangles), mxINT32_CLASS, mxREAL);
        memcpy(mxGetPr(plhs[4]), mesh->_meta.neighborlist, 3 * (int32_t)(mesh->_meta.numberoftriangles) * sizeof(int32_t));

        M_Ptr temp_4[] = {plhs[4], mxCreateDoubleScalar(1.0)};
        mexCallMATLAB(1, &plhs[4], 2, temp_4, "plus");

    }

    MEX_DEFINE(set_points_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 2);
        OutputArguments output(nlhs, plhs, 0);

        auto mesh = Session<tTriangleInfo>::get(input.get(0));

        /*
         * when array is also double, we do not have to convert.
         */
        Array<double> points(M_Cast<double>(C_CAST(prhs[1])), mxGetM(prhs[1]), false);
        mesh->set_points(points);
    }

    MEX_DEFINE(set_facets_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 2);
        OutputArguments output(nlhs, plhs, 0);

        auto mesh = Session<tTriangleInfo>::get(input.get(0));

        /*
         * when array is not double, we need to convert.
         */
        Array<int> facets; facets.convert(M_Cast<double>(C_CAST(prhs[1])), mxGetM(prhs[1]));
        mesh->set_facets(facets);
    }


    MEX_DEFINE(build_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 2);
        OutputArguments output(nlhs, plhs, 0);
        auto mesh_in = Session<tTriangleInfo>::get(input.get(0));
        auto mesh_out = Session<tTriangleInfo>::get(input.get(1));

        mesh_in->build("pcz", *mesh_out);
    }


    MEX_DEFINE(refine_tri)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 3);
        OutputArguments output(nlhs, plhs, 0);

        auto mesh_in = Session<tTriangleInfo>::get(input.get(0));
        auto mesh_out = Session<tTriangleInfo>::get(input.get(1));
        auto switches = mxArrayToString(input.get(2));

        mesh_in->refine((("pzcen") + std::string(switches)).c_str(), *mesh_out);
    }

    MEX_DEFINE(connectivity)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 3);


        auto mesh_in = Session<tTriangleInfo>::get(input.get(0));

        auto numberofelem = mesh_in->_meta.numberoftriangles;
        auto numberofsegs = mesh_in->_meta.numberofsegments;
        auto numberofedge = mesh_in->_meta.numberofedges;

        plhs[0] = mxCreateNumericMatrix((numberofedge - numberofsegs) * 2, 1, mxDOUBLE_CLASS, mxREAL);
        plhs[1] = mxCreateNumericMatrix((numberofedge - numberofsegs) * 2, 1, mxDOUBLE_CLASS, mxREAL);
        plhs[2] = mxCreateNumericMatrix((numberofedge - numberofsegs) * 2, 1, mxDOUBLE_CLASS, mxREAL);


        double* pI = mxGetPr(plhs[0]);
        double* pJ = mxGetPr(plhs[1]);
        double* pV = mxGetPr(plhs[2]);

        for (int elem_id = 0; elem_id < mesh_in->_meta.numberoftriangles; ++elem_id) {
        	for (int edge_id = 0; edge_id < 3; ++edge_id) {
        		*pI = elem_id + 1;
        		*pJ = mesh_in->_meta.neighborlist[3 * elem_id + edge_id] + 1;
        		if (*pJ != 0) {
        			pI++;pJ++; *pV = 1; pV++;
        		}
        	}
        }
    }
}

MEX_DISPATCH



