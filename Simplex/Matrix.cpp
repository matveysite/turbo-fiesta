#include "Matrix.h"
#include <algorithm>

Matrix::Matrix() {
}

Matrix::Matrix(const Matrix& m) : Matrix(m, 0, 0) {
}

Matrix::Matrix(int n, int m) {
	this->n = n;
	this->m = m;
	matrix = new double*[n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[m]();
	}
}

Matrix::Matrix(const Matrix& m, int dn, int dm) :Matrix(m.n + dn, m.m + dm) {
	for (int i = 0; i < std::min(m.n, m.n + dn); i++) {
		for (int j = 0; j < std::min(m.m, m.m + dm); j++) {
			matrix[i][j] = m.matrix[i][j];
		}
	}
}

Matrix::Matrix(const Matrix& mb, const Matrix& mn, const Set& base) :Matrix(mb.n + mn.n, 1) {
	int mbi=0, mni = 0;
	for (int i = 0; i < n; i++) {
		if (base.check(i)) {
			matrix[i][0] = mb.matrix[mbi++][0];
		}
		else {
			matrix[i][0] = mn.matrix[mni++][0];
		}
	}
}

double& Matrix::get(int i, int j) {
	return matrix[i][j];
}

Matrix Matrix::operator*(const Matrix& m2) const {
	Matrix result(this->n, m2.m);
	for (int i = 0; i < this->n; i++) {
		for (int j = 0; j < m2.m; j++) {
			for (int k = 0; k < this->m; k++) {
				result.matrix[i][j] += matrix[i][k]* m2.matrix[k][j];
			}
		}
	}
	return result;
}

Matrix & Matrix::operator=(const Matrix& m) {
	if ((n != m.n) || (this->m != m.m)){
		delete[] matrix;
		n = m.n;
		this->m = m.m;
		matrix = new double*[n];
	}
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[this->m];
		for (int j = 0; j < this->m; j++) {
			matrix[i][j] = m.matrix[i][j];
		}
	}
	return *this;
}

Matrix Matrix::operator+(const Matrix& m2) const {
	Matrix result(this->n, this->m);
	for (int i = 0; i < this->n; i++) {
		for (int j = 0; j < this->m; j++) {
			result.matrix[i][j] = this->matrix[i][j] + m2.matrix[i][j];
		}
	}
	return result;
}

Matrix operator*(double t, const Matrix& m) {
	Matrix result(m.n, m.m);
	for (int i = 0; i < result.n; i++) {
		for (int j = 0; j < result.m; j++) {
			result.matrix[i][j] = m.matrix[i][j] * t;
		}
	}
	return result;
}

Matrix Matrix::operator-(const Matrix& m2) const {
	return (*this) + (-1)*m2;
}

Matrix::~Matrix() {
	for (int i = 0; i < n; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

std::istream& operator>>(std::istream& in, const Matrix& m) {
	for (int i = 0; i < m.n; i++) {
		for (int j = 0; j < m.m; j++) {
			in >> m.matrix[i][j];
		}
	}
	return in;
}

std::ostream& operator<<(std::ostream& out, const Matrix& m) {
	for (int i = 0; i < m.n; i++) {
		for (int j = 0; j < m.m; j++) {
			out << m.matrix[i][j] << " ";
		}
		if (i != m.n - 1) {
			out << std::endl;
		}
	}
	return out;
}

int Matrix::getN() const {
	return n;
}

int Matrix::getM() const {
	return m;
}

Matrix Matrix::inverse() const {
	Matrix adj(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			Set rows(n), cols(m);
			rows.add(i), cols.add(j);
			rows.inverse(), cols.inverse();
			double sign = ((i + j) % 2) ? -1 : 1;
			adj.matrix[i][j] = sign*part(rows, cols).determinant();
		}
	}
	return (1 / determinant())*adj.transpose();
}

Matrix Matrix::transpose() const {
	Matrix result(m, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result.matrix[j][i] = matrix[i][j];
		}
	}
	return result;
}

double Matrix::determinant() const {
	double result = 0;
	if (n > 1) {
		int sign = 1;
		for (int j = 0; j < m; j++) {
			Set rows(n), cols(m);
			rows.add(0);
			rows.inverse();
			cols.add(j);
			cols.inverse();
			result += sign*part(rows, cols).determinant()*matrix[0][j];
			sign *= -1;
		}
	}
	else result = matrix[0][0];
	return result;
}

Matrix Matrix::part(const Set& cols) const {
	Set rows(n);
	rows.inverse();
	return part(rows, cols);
}

Matrix Matrix::part(const Set& rows, const Set& cols) const {
	Matrix result(rows.getCount(), cols.getCount());
	int ii = 0, jj = 0;
	for (int i = 0; i < n; i++) {
		if (rows.check(i)) {
			jj = 0;
			for (int j = 0; j < m; j++) {
				if (cols.check(j)) {
					result.matrix[ii][jj] = matrix[i][j];
					jj++;
				}
			}
			ii++;
		}
	}
	return result;
}
