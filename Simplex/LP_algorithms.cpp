#include <iostream>
#include "Matrix.h"
#include "Set.h"
#include "LP_algorithms.h"

using namespace std;

double sign(double a) {
	return (a > 0) - (a < 0);
}

void dualSimplexMethod(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Set& base) {
	int n = A.getM(), m = A.getN();
	bool crit = false, solvable = true;
	Matrix x(n, 1);
	int iteration = 0;
	while (!crit && solvable) {
		cout << "----------------------------------------------------" << endl;
		cout << ++iteration << " итерация" << endl;
		cout << "----------------------------------------------------" << endl;
		cout << "Базис: " << base << endl;
		Matrix A_base = A.part(base), c_base = c.transpose().part(base).transpose();
		base.inverse();
		Matrix A_n = A.part(base);
		base.inverse();
		Matrix A_bti = A_base.transpose().inverse();
		Matrix u = A_bti*c_base;
		cout << "u: " << u.transpose() << endl;
		Matrix delta = c - A.transpose()*u;
		cout << "delta: " << delta.transpose() << endl;

		Matrix kappa_n(n - m, 1);
		int cc = 0;
		for (int i = 0; i < n; i++) {
			if (!base.check(i)) {
				kappa_n.get(cc++, 0) = (delta.get(i, 0) <= 0) ? dl.get(i, 0) : dh.get(i, 0);
			}
		}
		Matrix kappa_b = A_base.inverse()*(b - A_n*kappa_n);
		Matrix kappa(kappa_b, kappa_n, base);
		cout << "kappa: " << kappa.transpose() << endl;
		crit = true;
		int jz = -1;
		double xn = 0;
		int jzl = m;
		for (int i = n - 1; (i >= 0) && crit; i--) {
			if (base.check(i)) {
				if (kappa.get(i, 0) < dl.get(i, 0)) {
					crit = false;
					jz = i;
					xn = dl.get(i, 0);
				}
				else if (kappa.get(i, 0) > dh.get(i, 0)) {
					crit = false;
					jz = i;
					xn = dh.get(i, 0);
				}
				jzl--;
			}
		}
		if (!crit) {
			cout << "j*: " << jz + 1 << endl;
			Matrix ejz(m, 1);
			ejz.get(jzl, 0) = -1;
			Matrix pu = sign(kappa.get(jz, 0) - xn)*A_bti*ejz;
			cout << "pu: " << pu.transpose() << endl;
			Matrix pdelta = -1 * A_n.transpose()*pu;
			cout << "p_delta: " << pdelta.transpose() << endl;
			double minSigma = -1;
			int j0 = -1;
			int c = 0;
			for (int i = 0; i < n; i++) {
				if (!base.check(i)) {
					if (pdelta.get(c, 0)*delta.get(i, 0) < 0) {
						double sigma = -delta.get(i, 0) / pdelta.get(c, 0);
						cout << "Sigma_" << i + 1 << ": " << sigma << endl;
						if ((sigma < minSigma) || (minSigma == -1)) {
							minSigma = sigma;
							j0 = i;
						}
					}
					else {
						cout << "Sigma_" << i + 1 << ": +oo" << endl;
					}
					c++;
				}
			}
			cout << "min Sigma: ";
			if (minSigma == -1) {
				cout << "+oo" << endl;
				solvable = false;
			}
			else {
				cout << minSigma << endl;
				cout << "j0: " << j0 + 1 << endl;
				base.remove(jz);
				base.add(j0);
			}
		}
		else {
			x = kappa;
		}
	}
	if (solvable) {
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
		cout << "Критерий выполнен!" << endl;
		cout << "Базис: " << base << endl;
		cout << "X: " << x.transpose() << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}
	else {
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
		cout << "Решения нет (двойственная ЦФ бесконечно убывает), или задача двойственно вырождена!" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}
}

void simplexMethod(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Matrix& x, Set& base) {
	int n = A.getM(), m = A.getN();
	bool crit = false;
	int iteration = 0;
	while (!crit) {
		cout << "----------------------------------------------------" << endl;
		cout << ++iteration << " итерация" << endl;
		cout << "----------------------------------------------------" << endl;
		cout << "Базис: " << base << endl;
		cout << "X: " << x.transpose() << endl;
		Matrix A_base = A.part(base), c_base = c.transpose().part(base).transpose();
		Matrix u = A_base.transpose().inverse()*c_base;
		cout << "u: " << u.transpose() << endl;
		Matrix delta = c - A.transpose()*u;
		cout << "delta: " << delta.transpose() << endl;
		crit = true;
		int j0 = -1;
		for (int i = 0; i < n; i++) {
			if (!base.check(i)) {
				double dj = delta.get(i, 0);
				double xj = x.get(i, 0);
				double dlj = dl.get(i, 0);
				double dhj = dh.get(i, 0);
				if (!((dj <= 0) && (xj == dlj) || (dj >= 0) && (xj == dhj))) {
					crit = false;
					j0 = i;
					break;
				}
			}
		}
		if (!crit) {
			cout << "j0: " << j0 + 1 << endl;
			int j0n = j0;
			for (int i = 0; i < j0; i++) {
				if (base.check(i)) {
					j0n--;
				}
			}
			Matrix ln(n - m, 1);
			double dj0 = delta.get(j0, 0);
			ln.get(j0n, 0) = sign(dj0);
			cout << "l_н: " << ln.transpose() << endl;
			Set Aj0(n);
			Aj0.add(j0);
			Matrix lb = (-sign(dj0)) * A_base.inverse()*A.part(Aj0);
			cout << "l_б: " << lb.transpose() << endl;
			Matrix l(lb, ln, base);
			double minT = -1;
			int jn = -1;
			for (int j = 0; j < n; j++) {
				double lj = l.get(j, 0);
				if (base.check(j) || j == j0) {
					if (lj != 0) {
						double theta = (((lj < 0) ? dl.get(j, 0) : dh.get(j, 0)) - x.get(j, 0)) / lj;
						cout << "Theta_" << j + 1 << ": " << theta << endl;
						if ((theta < minT) || (minT == -1)) {
							minT = theta;
							jn = j;
						}
					}
					else {
						cout << "Theta_" << j + 1 << ": +oo" << endl;
					}
				}
			}
			cout << "min Theta: " << minT << endl;
			cout << "j*: " << jn + 1 << endl;
			base.add(j0);
			base.remove(jn);
			x = x + minT*l;
		}
	}
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "Критерий выполнен!" << endl;
	cout << "Базис: " << base << endl;
	cout << "X: " << x.transpose() << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
}

void problemSensitivity(Matrix& A, Matrix& b, Matrix& c, Matrix& dl, Matrix& dh, Set& optimalBase) {
}
