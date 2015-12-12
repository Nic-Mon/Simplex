


#ifndef __QS_MATRIX_CPP
#define __QS_MATRIX_CPP

#include "matrix.hpp"
#include <iostream>



// Parameter Constructor
template<typename T>
QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _initial) {
    mat.resize(_rows);
    for (unsigned i=0; i<mat.size(); i++) {
        mat[i].resize(_cols, _initial);
    }
    rows = _rows;
    cols = _cols;
}

template<typename T>
QSMatrix<T>::QSMatrix(std::vector< std::vector<T> > matrix) {
    mat = matrix;
    rows = matrix.size();
    cols = matrix[0].size();
}

// Copy Constructor
template<typename T>
QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) {
    mat = rhs.mat;
    rows = rhs.get_rows();
    cols = rhs.get_cols();
}

// (Virtual) Destructor
template<typename T>
QSMatrix<T>::~QSMatrix() {}

// Assignment Operator
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) {
    if (&rhs == this)
        return *this;
    
    unsigned new_rows = rhs.get_rows();
    unsigned new_cols = rhs.get_cols();
    
    mat.resize(new_rows);
    for (unsigned i=0; i<mat.size(); i++) {
        mat[i].resize(new_cols);
    }
    
    for (unsigned i=0; i<new_rows; i++) {
        for (unsigned j=0; j<new_cols; j++) {
            mat[i][j] = rhs(i, j);
        }
    }
    rows = new_rows;
    cols = new_cols;
    
    return *this;
}

// Addition of two matrices
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const QSMatrix<T>& rhs) {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] + rhs(i,j);
        }
    }
    
    return result;
}

// Cumulative addition of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            this->mat[i][j] += rhs(i,j);
        }
    }
    
    return *this;
}

// Subtraction of this matrix and another
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] - rhs(i,j);
        }
    }
    
    return result;
}

// Cumulative subtraction of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator-=(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            this->mat[i][j] -= rhs(i,j);
        }
    }
    
    return *this;
}

// Left multiplication of this matrix and another
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            for (unsigned k=0; k<rows; k++) {
                result(i,j) += this->mat[i][k] * rhs(k,j);
            }
        }
    }
    
    return result;
}

// Cumulative left multiplication of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix<T>& rhs) {
    QSMatrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}

// Calculate a transpose of this matrix
template<typename T>
QSMatrix<T> QSMatrix<T>::transpose() {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[j][i];
        }
    }
    
    return result;
}

// Matrix/scalar addition
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] + rhs;
        }
    }
    
    return result;
}

// Matrix/scalar subtraction
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] - rhs;
        }
    }
    
    return result;
}

// Matrix/scalar multiplication
template<typename T>
QSMatrix<T> QSMatrix<T>::operator*(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] * rhs;
        }
    }
    
    return result;
}

// Matrix/scalar division
template<typename T>
QSMatrix<T> QSMatrix<T>::operator/(const T& rhs) {
    QSMatrix result(rows, cols, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result(i,j) = this->mat[i][j] / rhs;
        }
    }
    
    return result;
}

// Multiply a matrix with a vector
template<typename T>
std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs) {
    std::vector<T> result(rhs.size(), 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        for (unsigned j=0; j<cols; j++) {
            result[i] = this->mat[i][j] * rhs[j];
        }
    }
    
    return result;
}

// Obtain a vector of the diagonal elements
template<typename T>
std::vector<T> QSMatrix<T>::diag_vec() {
    std::vector<T> result(rows, 0.0);
    
    for (unsigned i=0; i<rows; i++) {
        result[i] = this->mat[i][i];
    }
    
    return result;
}

// Access the individual elements
template<typename T>
T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) {
    return this->mat[row][col];
}

// Access the individual elements (const)
template<typename T>
const T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) const {
    return this->mat[row][col];
}

// Get the number of rows of the matrix
template<typename T>
unsigned QSMatrix<T>::get_rows() const {
    return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned QSMatrix<T>::get_cols() const {
    return this->cols;
}

//my added stuff below
//row operations
template<typename T>
std::vector<T> QSMatrix<T>::row(int row_num) {
    return mat[row_num + 1];
}

template<typename T>
void QSMatrix<T>::row_scalar_mult(int row_num, double scalar){
    for (auto& x : mat[row_num-1]) {    x *= scalar;    }
}

//add row 1 to row 2 (row 1 unchanged)
template<typename T>
void QSMatrix<T>::add_row(int row1, int row2, double scalar){
    for(int i=0; i<get_cols(); i++){
        mat[row2-1][i] += ( scalar * mat[row1-1][i] );
    }
}

//forms tableau
template<typename T>
void QSMatrix<T>::form_tableau(std::vector<T> b, std::vector<T> c){
    
    for(int i=0; i<rows; i++){
        mat[i].push_back(b[i]);
    }
    
    c.push_back(0);
    mat.push_back(c);
    
    rows++; cols++;
}

template <typename T>
int QSMatrix<T>::piv_col_index() {
    double smallest = 0;
    
    for(auto x : mat[rows-1]) {     if(x<smallest) smallest = x;    }
    
    if (smallest == 0) return -1;
    else {
        return find(mat[rows-1].begin(), mat[rows-1].end(), smallest) - mat[rows-1].begin();
    }
}

template<typename T>
int QSMatrix<T>::piv_row_index(int piv_col_index) {
    int b = cols-1;
    std::vector<double> v;
    for(int i=0; i<rows-1; i++) {
        v.push_back(mat[i][b] / mat[i][piv_col_index]);
    }
    
    double smallest = v[0];
    int smallest_index = 0;
    for(int i=1; i<rows-1; i++) {
        if (v[i] < smallest) {
            smallest = v[i];
            smallest_index = i;
        }
    }
    
    return smallest_index;
}

template<typename T>
void QSMatrix<T>::pivot(int row_index, int col_index) {
    row_scalar_mult(row_index+1, ( 1. / mat[row_index][col_index]));
    
    for(int i=0; i<rows; i++){
        if ( i != row_index && mat[i][col_index] != 0) {
            add_row(row_index+1, i+1, -mat[i][col_index]);
        }
    }
}

template<typename T>
void QSMatrix<T>::print() {
    for (int i=0; i<rows; i++) {
        std::cout << "[";
        
        for (int j=0; j<cols-1; j++) {
            std::cout << std::setw(4);
            std::cout << mat[i][j] << ", ";
        }
        std::cout << std::setw(5);
        std::cout << mat[i][cols-1];
        std::cout << "]";
        std::cout << std::endl;
    }
}

template<typename T>
void QSMatrix<T>::simplex() {
    
    int step = 1;
    
    while(piv_col_index() != -1) {
        pivot(piv_row_index(piv_col_index()), piv_col_index());
        
        std::cout << "STEP " << step << ": " << std::endl;
        step++;
        print();
        std::cout << std::endl;
    }
    
}
#endif