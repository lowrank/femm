//
// Created by lurker on 4/14/17.
//

#ifndef TMESHINFO_H
#define TMESHINFO_H

#include <iostream>
#include <cassert>
#include <cstring>
#include <unordered_map>

extern "C" {
#include "triangle.h"
}

template <typename T>
class Array {
protected:
    T* data;
    size_t size;
    bool owner;
public:
    Array() {
        owner = false;
    }
    Array(size_t _n) {
        data = (T*)malloc(_n * sizeof(T));
        memset(data, 0, _n * sizeof(T));
        size = _n;
        owner = true;
    }

    Array(size_t _n, bool _o) {
        data = (T*)malloc(_n * sizeof(T));
        size = _n;
        owner = _o;
    }

    Array(T* _data, size_t _n, bool _o) {
        if (_o) {
            memcpy(data, _data, _n * sizeof(T));
            size = _n;
            owner = _o;
        }
        else {
            data = _data;
            size = _n;
            owner = _o;
        }
    }

    ~Array() {
        if (owner) {
            free(data);
            data = nullptr;
        }
    }
    size_t get_size() {
        return size;
    }

    void convert(double* _data, size_t _n) {
        data = (T*)malloc(_n * sizeof(T));
        for (int i = 0; i < _n; ++i) {
            data[i] = static_cast<T>(_data[i]);
        }
        size = _n;
        owner = true;
    }

    bool get_owner() {
        return owner;
    }

    T* get_data() {
        return data;
    }

    T& operator[](size_t i) {
        assert(i < size);
        return data[i];
    }

    const T& operator[](size_t i) const{
        assert(i < size);
        return data[i];
    }
};

class tTriangleInfo {
public:
    tTriangleInfo(const tTriangleInfo&) = delete;
    tTriangleInfo();
    ~tTriangleInfo();

    void set_points(Array<double>& _points);
    void set_facets(Array<int>& _facets);
    void build(std::string, tTriangleInfo& out);
    void refine(std::string, tTriangleInfo& out);

    struct triangulateio _meta;
};


#endif //TMESHINFO_H
