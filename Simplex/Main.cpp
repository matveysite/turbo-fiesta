#include <iostream>
#include <math.h>
#include "LP_algorithms.h"
#include "Matrix.h"
#include "Set.h"

using namespace std;

int main() {
	setlocale(LC_ALL, "Russian");
	int n, m;
	bool success, needSecondPhase = true;

	cout << "n: ";
	cin >> n;
	cout << "m: ";
	cin >> m;

	Matrix initialA(m, n), B(m, 1), initialC(n, 1), initialDl(n, 1), initialDh(n, 1);
	cout << "������� A: " << endl;
	cin >> initialA;
	cout << "������ B: " << endl;
	cin >> B;
	cout << "������ C: " << endl;
	cin >> initialC;
	cout << "������ D_*: " << endl;
	cin >> initialDl;
	cout << "������ D^*: " << endl;
	cin >> initialDh;

	Matrix X = initialDl, W = B - initialA*X;
	cout << "������ ����, x: " << X.transpose() << endl;
	cout << "������ ����, w: " << W.transpose() << endl;
	int lastNullInW = 0;
	while ((lastNullInW < m) && (W.get(lastNullInW, 0) == 0)) {
		lastNullInW++;
	}
	Matrix simplexA(initialA, 0, m), simplexDl(initialDl, m, 0), simplexDh(initialDh, m, 0), simplexC(n + m, 1);
	Set base(n + m);
	X = Matrix(X, m, 0);
	if (lastNullInW != m) {
		for (int i = 0; i < m; i++) {
			double c = sign(W.get(i, 0));
			if (c == 0) {
				c = 1;
			}
			simplexA.get(i, n + i) = c;
			simplexDh.get(n + i, 0) = abs(W.get(i, 0));
			simplexC.get(n + i, 0) = -1;
			base.add(n + i);
			X.get(n + i, 0) = abs(W.get(i, 0));
		}
		needSecondPhase = true;
	}
	else {
		simplexC = Matrix(initialC, m, 0);
		for (int i = n; i < n + m; i++) {
			base.add(i);
			simplexA.get(i - n, i) = 1;
		}
		cout << "������ ���� �� �����, ������ �������� ������." << endl;
		needSecondPhase = false;
	}
	cout << "������ ����, A: " << endl << simplexA << endl;
	simplexMethod(simplexA, B, simplexC, simplexDl, simplexDh, X, base);

	if (needSecondPhase) {
		int i = n;
		while ((i < n + m) && (X.get(i, 0) == 0)) {
			i++;
		}
		if ((i != (n + m))) {
			cout << "������� ���, �� ��� ������������� ���������� ����������!" << endl;
			success = false;
		}
		else {
			cout << endl << "������ ����: " << endl;
			int artInBase = 0;
			for (int i = n; i < n + m; i++) {
				if (base.check(i)) {
					artInBase++;
				}
			}
			X = Matrix(X, -m + artInBase, 0);
			simplexA = Matrix(initialA, 0, artInBase);
			simplexC = Matrix(initialC, artInBase, 0);
			simplexDl = Matrix(initialDl, artInBase, 0);
			simplexDh = Matrix(initialDh, artInBase, 0);
			Set secondPhaseBase(n + artInBase);
			secondPhaseBase.insert(base, n);
			if (artInBase) {
				cout << "� ������ �������� ������������� ����������, ������ �������� ������!" << endl;
				int localI = 0;
				for (int i = n; i < n + m; i++) {
					if (base.check(i)) {
						simplexA.get(i - n, n + localI) = 1;
						secondPhaseBase.add(n + localI);
						localI++;
					}
				}
			}
			base = secondPhaseBase;
			simplexMethod(simplexA, B, simplexC, simplexDl, simplexDh, X, base);
			success = true;
		}
	} else {
		success = true;
	}

//-----------------------------------------------------------------------------------------
	cout << "������������ ������: " << endl;
	cout << "( " << B.transpose() << (-1) * initialDl.transpose() << initialDh.transpose() << ") * lambda -> min" << endl;
	Matrix dualA(initialA.transpose(), 0, 2*n);
	for(int i = 0; i < n; i++){
		dualA.get(i, m + i) = -1;
		dualA.get(i, m + n + i) = 1;
	}
	cout << endl;
	cout << dualA << " * lambda = ( " << initialC.transpose() << ")^T" << endl << endl;
	cout << "v >= 0" << endl << "w >= 0" << endl << endl;
	cout << "������������ ������������ ����: " << endl;
	Matrix w(n, 1), v(n, 1), y(n, 1);
	for(int i = 0; i < n; i++) {
		double c = initialC.get(i, 0);
		if(c < 0) {
			w.get(i, 0) = 0;
			v.get(i, 0) = -c;
		} else {
			w.get(i, 0) = c;
			v.get(i, 0) = 0;
		}
	}
	cout << "lambda = ( " << y.transpose() << v.transpose() << w.transpose() << ")" << endl << endl;
	if (success) {
		int fictInBase = 0;
		for(int i = n; i < base.getN(); i++) {
			if(base.check(i)) {
				fictInBase++;
			}
		}
		if(!fictInBase) {
			cout << "����������� �������� ������������ ����(�� ������ � �� ������ ��������, ���� �� ��������): " << endl;
			Matrix dualA_base = initialA.part(base), c_base = initialC.transpose().part(base).transpose();
			Matrix dualA_bti = dualA_base.transpose().inverse();
			Matrix dualU = dualA_bti*c_base;
			cout << "u: " << dualU.transpose() << endl;
			Matrix dualDelta = initialC - initialA.transpose()*dualU;
			cout << "delta: " << dualDelta.transpose() << endl;
			for (int i = 0; i < n; i++) {
				double d = dualDelta.get(i, 0);
				if (d < 0) {
					w.get(i, 0) = 0;
					v.get(i, 0) = -d;
				}
				else {
					w.get(i, 0) = d;
					v.get(i, 0) = 0;
				}
			}
			cout << "lambda = ( " << dualU.transpose() << v.transpose() << w.transpose() << ")" << endl << endl;
		} else {
			cout << "� ������ �������� ��������� ����������, ����������� �������� ������������ ���� �� ��������� ������ ���!" << endl;
		}
	}

	//-----------------------------------------------------------------------------------------
	Matrix dsA(initialA, 0, m), dsC(initialC, m, 0), dsDl(initialDl, m, 0), dsDh(initialDh, m, 0);
	Set dBase(n + m);
	for (int i = n; i < n + m; i++) {
		dBase.add(i);
		dsA.get(i - n, i) = 1;
	}
	cout << endl << endl << endl << endl << "������������ �������� �����: " << endl;
	dualSimplexMethod(dsA, B, dsC, dsDl, dsDh, dBase);
	system("pause");
}