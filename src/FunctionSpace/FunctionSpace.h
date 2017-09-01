//
// FunctionSpace initialize space dofs, based on geometry information.
//

#ifndef FEM_FUNCTIONSPACE_H
#define FEM_FUNCTIONSPACE_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstring>
#include <algorithm>

#include "mexplus.h"
#include "common.h"

using namespace std;

class FunctionSpace {
public:
    FunctionSpace(const FunctionSpace&) = delete;
    FunctionSpace();
    ~FunctionSpace();

    void build(double* p_ptr, int* s_ptr, int* t_ptr, int* e_ptr, int* n_ptr, int np, int ns, int nt, int ne, int deg);

    // function space only contains points to be evaluated.
    // for integration, quadrature order will be adaptive.
    // this part stems from original Mesh Promote.
    // CG is now the only choice.
    double  *nodePtr;
    int     *segPtr, *triPtr, *edgePtr, *neighPtr;
    int deg;

    std::vector<double> nodes;
    std::unordered_map<std::string,std::vector<int>> edges;
    std::vector<int> elems;
    std::vector<int> boundary;
    std::vector<int> boundary_index;
    std::vector<int> neighbors;

    void clear();
    void resize(int nsize, int esize, int bsize, int edsize);

private:
    inline double _Point_X(int _index){
        return nodePtr[2*_index ];
    }

    inline double _Point_Y(int _index){
        return nodePtr[2*_index + 1];
    }

    inline int _Edge_L(int _index){
        return edgePtr[2*_index];
    }

    inline int _Edge_R(int _index){
        return edgePtr[2* _index + 1];
    }

    inline int _Tri_U(int _index){
        return triPtr[3 * _index];
    }

    inline int _Tri_V(int _index){
        return triPtr[3 * _index + 1];
    }

    inline int _Tri_W(int _index){
        return triPtr[3 * _index + 2];
    }

    inline int _Seg_L(int _index){
        return segPtr[2 * _index ];
    }

    inline int _Seg_R(int _index){
        return segPtr[2 * _index + 1];
    }



};


#endif //FEM_FUNCTIONSPACE_H
