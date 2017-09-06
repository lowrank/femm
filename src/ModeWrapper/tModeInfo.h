//
// Created by lurker on 4/18/17.
//

#ifndef MESHWRAPPER_TMODEINFO_H
#define MESHWRAPPER_TMODEINFO_H


#include <vector>
#include <cstdlib>
#include <quadmath.h>

#include <cstring>
#include "mexplus.h"
#include "common.h"


using namespace mexplus;


using std::vector;

class QRule_tet {
public:
    vector<double > points_x;
    vector<double > points_y;
    vector<double > points_z;
    vector<double > weights;
    size_t          degree;

    void resize(size_t n) {
        points_x.resize(n);
        points_y.resize(n);
        points_z.resize(n);
        weights.resize(n);
    }
};

class tModeInfo_tet {
public:
    QRule_tet table;

    tModeInfo_tet();
    ~tModeInfo_tet();

    void get_vr_data(size_t deg);
};

class QRule_tri {
public:
    vector<__float128 > points_x;
    vector<__float128 > points_y;
    vector<__float128 > weights;
    size_t          degree;

    void resize(size_t n) {
        points_x.resize(n);
        points_y.resize(n);
        weights.resize(n);
    }
};

class tModeInfo_tri {
public:
    QRule_tri table;

    tModeInfo_tri();
    ~tModeInfo_tri();

    void get_vr_data(size_t deg);
};


class QRule_lin {
public:
	vector<double> points_x;
	vector<double> weights;

	void resize(size_t n) {
		points_x.resize(n);
		weights.resize(n);
	}
};

class tModeInfo_lin {
public:
	QRule_lin table;
	tModeInfo_lin();
	~tModeInfo_lin();

	void get_vr_data(size_t deg);
};


#endif //MESHWRAPPER_TMODEINFO_H
