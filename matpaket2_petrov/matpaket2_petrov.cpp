#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void printMatrix(const vector<vector<double>>& A, const vector<double>& b) {
	int n = b.size();
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << setw(9) << fixed << setprecision(4) << A[i][j] << " ";
		}
		cout << "| " << b[i] << endl;
	}
}

void printMultiplierMatrix(const vector<vector<double>>& multipliers) {
	int n = multipliers.size();
	cout << "\nMatrica M:\n\n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j < i) {
				cout << setw(9) << fixed << setprecision(4) << multipliers[i][j] << " ";
			}
			else if (j == i) {
				cout << setw(9) << fixed << setprecision(4) << 1.0 << " ";
			}
			else {
				cout << setw(9) << fixed << setprecision(4) << 0.0 << " ";
			}
		}
		cout << endl;
	}
}

void GaussElimination(vector<vector<double>>& A, vector<double>& b) {
	int n = b.size();

	vector<vector<double>> P(n, vector<double>(n, 0.0));
	for (int i = 0; i < n; i++) {
		P[i][i] = 1.0;
	}

	for (int i = 0; i < (n - 1); ++i) {

		int maxRow = i;
		for (int k = i + 1; k < n; ++k) {
			if (abs(A[k][i]) > abs(A[maxRow][i])) {
				maxRow = k;
			}
		}

		swap(A[i], A[maxRow]);
		swap(b[i], b[maxRow]);

		swap(P[i], P[maxRow]);

		vector<vector<double>> M(n, vector<double>(n, 0.0));
		for (int k = i + 1; k < n; ++k) {
			double factor = A[k][i] / A[i][i];
			M[k][i] = factor;
			for (int j = i; j < n; ++j) {
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}

		printMultiplierMatrix(M);

		cout << "\nMatrica P:\n\n";
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				cout << setw(9) << fixed << setprecision(4) << P[i][j] << " ";
			}
			cout << "\n";
		}

		cout << "\nSLAU (Matrica A) posle " << (i + 1) << " shaga perestanovki: \n\n";
		printMatrix(A, b);
	}

	vector<double> x(n);
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}

	cout << "\nReshenie:\n\n";
	for (int i = 0; i < n; ++i) {
		cout << setw(5) << "X[" << i << "]  =  " << x[i] << endl;
	}
}

int main() {

	cout << "----------------------------------zadanie 1-------------------------------\n";

	int n = 5;
	double pi = 3.141592653589793238462643383279502884;

	vector<vector<double>> A(n, vector<double>(n, 0));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j == 0 || i == j + 1 || (i == 0 && j == n - 1)) {
				A[i][j] = 1;
			}
			else if (i == j) {
				A[i][j] = 1 + (j) * 0.01;
			}
			else {
				A[i][j] = 0;
			}
		}
	}

	cout << "Matrica A:\n\n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << setw(9) << fixed << setprecision(4) << A[i][j] << " ";
		}
		cout << "\n";
	}

	vector<double> b(n, 0);
	for (int i = 0; i < n; ++i) {
		b[i] = sin(pi / (n - i));
	}

	cout << "\nVector b:\n\n"; for (int i = 0; i < n; ++i) {
		cout << setw(9) << b[i] << "\n";
	}

	cout << "\n----------------------------------------------------------\n\n";
	cout << "SLAU:\n\n";
	printMatrix(A, b);
	GaussElimination(A, b);


	cout << "\n-------------------------zadanie 4-------------------------\n\n";

	int np = 10;

	vector<vector<double>> Ap(np, vector<double>(np, 0));

	for (int i = 0; i < np; ++i) {
		for (int j = 0; j < np; ++j) {
			if (i == j) {
				Ap[i][j] = 3.1 * (i + 1);
			}
			if (i == (j - 1)) {
				Ap[i][j] = -2 * (i + 1);
			}
			if (i == (j + 1)) {
				Ap[i][j] = i + 1;
			}
		}
	}

	vector<double> d(np, 1.0);
	for (int i = 0; i < np; ++i) {
		double k = i + 1;
		d[i] = (2.1 * k * k + 7.2 * k + 2) / (k * k + 3 * k + 2);
	}

	cout << "Matrix:\n\n";
	for (int i = 0; i < np; ++i) {
		for (int j = 0; j < np; ++j) {
			cout << setw(9) << fixed << setprecision(4) << Ap[i][j] << " ";
		}
		cout << "| " << d[i] << endl;
	}

	vector<double> gamma(np), beta(np), xp(np);

	gamma[0] = Ap[0][0];
	beta[0] = d[0] / gamma[0];
	for (int i = 1; i < np; ++i) {
		gamma[i] = Ap[i][i] - (Ap[i][i - 1] * Ap[i - 1][i]) / gamma[i - 1];
		beta[i] = (d[i] - Ap[i][i - 1] * beta[i - 1]) / gamma[i];
	}

	xp[np - 1] = beta[np - 1];
	for (int i = np - 2; i >= 0; --i) {
		xp[i] = beta[i] - (Ap[i][i + 1] * xp[i + 1]) / gamma[i];
	}

	cout << "\nSolution:\n\n";
	for (int i = 0; i < np; ++i) {
		cout << setw(5) << "X[" << i << "]  =  " << xp[i] << endl;
	}


	return 0;
}