//
//  main.cpp
//  Simplex
//
//  Created by Nic Mon on 11/6/15.
//  Copyright © 2015 Nic Mon. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include "matrix.hpp"

void read_func(vector<string>& v, vector<double> c) {
    string func;
    cin >> func;
    
    istringstream in(func);
    
    double d;
    string s;
    
    while( in.peek() != '\n' ) {
        in >> d;
        c.push_back(d);
        in >> s;
        v.push_back(s);
    }
}

void print(double i) { std::cout << i <<" "; return; }


int main(int argc, const char * argv[]) {
    
    cout << "           LINEAR OPTIMZER             " << endl;
    cout << "---------------------------------------" << endl << endl;
    
    /*
    vector<string> vars;
    vector<double> c;
    
    cout << "What equation to maximize?" << endl;
    read_func(vars, c);
    
    for(auto x: vars) { cout << x << " "; }
    cout << endl;
    for(auto x: c) { cout << x << " "; }
    */
    
    
    
    std::vector<double> v1 = { 2, 1, 1, 1, 0, 0 };
    std::vector<double> v2 = { 1, 3, 2, 0, 1, 0 };
    std::vector<double> v3 = { 2, 1, 2, 0, 0, 1 };
    
    std::vector<double> b = { 180, 300, 240 };
    std::vector<double> e = { -6, -5, -4, 0, 0, 0 };
    
    std::vector< std::vector<double> > a = { v1, v2, v3};
    
    QSMatrix<double> matrix = QSMatrix<double>(a);
    matrix.form_tableau(b, e);
    
    matrix.print();
    std::cout << std::endl;
    
    matrix.simplex();
    
    
    return 0;
}
