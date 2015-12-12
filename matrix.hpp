//
//  matrix.hpp
//  Simplex
//
//  Created by Nic Mon on 11/6/15.
//  Copyright © 2015 Nic Mon. All rights reserved.
//


#include <stdio.h>

#ifndef __QS_MATRIX_H
#define __QS_MATRIX_H

#include <vector>
#include <iomanip>
using namespace std;

template <typename T> class QSMatrix {
private:
    std::vector<std::vector<T> > mat;
    unsigned rows;
    unsigned cols;
    
public:
    QSMatrix(unsigned _rows, unsigned _cols, const T& _initial);
    QSMatrix(std::vector<std::vector<T> > matrix);
    QSMatrix(const QSMatrix<T>& rhs);
    virtual ~QSMatrix();
    
    // Operator overloading, for "standard" mathematical matrix operations
    QSMatrix<T>& operator=(const QSMatrix<T>& rhs);
    
    // Matrix mathematical operations
    QSMatrix<T> operator+(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator+=(const QSMatrix<T>& rhs);
    QSMatrix<T> operator-(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator-=(const QSMatrix<T>& rhs);
    QSMatrix<T> operator*(const QSMatrix<T>& rhs);
    QSMatrix<T>& operator*=(const QSMatrix<T>& rhs);
    QSMatrix<T> transpose();
    
    // Matrix/scalar operations
    QSMatrix<T> operator+(const T& rhs);
    QSMatrix<T> operator-(const T& rhs);
    QSMatrix<T> operator*(const T& rhs);
    QSMatrix<T> operator/(const T& rhs);
    
    // Matrix/vector operations
    std::vector<T> operator*(const std::vector<T>& rhs);
    std::vector<T> diag_vec();
    
    // Access the individual elements
    T& operator()(const unsigned& row, const unsigned& col);
    const T& operator()(const unsigned& row, const unsigned& col) const;
    
    // Access the row and column sizes
    unsigned get_rows() const;
    unsigned get_cols() const;
    
    //my added stuff below
    //row operations
    std::vector<T> row(int row_num);
    void row_scalar_mult(int row_num, double scalar);
    
    //add or subtract scalar*row 1 to or from row 2 (row 1 unchanged)
    void add_row(int row1, int row2, double scalar);
    
    //forms tableau
    void form_tableau( std::vector<T> b, std::vector<T> c);
    
    // returns -1 if optimal, index of most neg value of obj function
    int piv_col_index();
    
    int piv_row_index(int piv_col);
    
    void pivot(int row_index, int col_index);
    
    void print();
    
    void simplex();
};


class formatted_output
{
private:
    int width;
    ostream& stream_obj;
    
public:
    formatted_output(ostream& obj, int w): width(w), stream_obj(obj) {}
    
    template<typename T>
    formatted_output& operator<<(const T& output)
    {
        //stream_obj.width(10);
        stream_obj << output;
        
        return *this;
    }
    
    formatted_output& operator<<(ostream& (*func)(ostream&))
    {
        func(stream_obj);
        return *this;
    }
};

#include "matrix.cpp"

#endif
