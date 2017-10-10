#pragma once
#include <iostream>
#include "Set.h"

class Matrix {
private:
	Matrix();
	int n, m;
	double **matrix;
public:
	Matrix(const Matrix&);
	Matrix(int, int);
	Matrix(const Matrix&, int, int);
	Matrix(const Matrix&, const Matrix&, const Set&);
	double& get(int, int);
	Matrix operator+(const Matrix&) const;
	friend Matrix operator*(double, const Matrix&);
	friend std::istream& operator>>(std::istream&, const Matrix&);
	friend std::ostream& operator<<(std::ostream&, const Matrix&) ;
	Matrix operator-(const Matrix&) const;
	Matrix operator*(const Matrix&) const;
	Matrix& operator=(const Matrix&);
	~Matrix();
	int getN() const;
	int getM() const;
	Matrix inverse() const;
	Matrix transpose() const;
	double determinant() const;
	Matrix part(const Set&) const;
	Matrix part(const Set&, const Set&) const;
};