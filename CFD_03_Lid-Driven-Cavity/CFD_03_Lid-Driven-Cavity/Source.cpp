//0.8 0.25 0.01
//0.8 0.05 

#include <iostream>
#include <math.h>
#include "write.h"
#include "math.h"

#define PI 3.1415926
#define N 81
#define Re 100
#define SCHEME 1
#define SOR_s 0.8
#define SOR_p 0.8
#define tol 1E-7
#define P 0
#define U 1
#define V 2

using namespace std;

//int p->0; u->1; v->2
void setBC(int);
void pseudo(int);
void scheme(int, int, int, int);
void corr();
double residual();
void result();

double u[N + 1][N] = { 0.0 };
double v[N][N + 1] = { 0.0 };
double p[N + 1][N + 1] = { 0.0 }, p_corr[N + 1][N] = { 0.0 };
double nu = (double)1 / (Re - 1);
double h = (double)1 / (N - 1);
double dm = 0;
double u_final[N][N], v_final[N][N], p_final[N][N], velocity_final[N][N];
char name1[] = "velocity.csv", name2[] = "u.csv", name3[] = "v.csv", name4[] = "p.csv";

int main() {
	int i, iteration = 0;
	setBC(U);
	setBC(V);
	setBC(P);

	do {
		iteration += 1;
		pseudo(U);
		setBC(U);
		pseudo(V);
		setBC(V);
		corr();
		if (iteration % 100 == 0) {
			cout << "iteration:	" << iteration << "	residual:	" << dm << endl;
			result();

		}
	} while (residual() > tol);//
	result();
	return 0;
}


/*Functions definition*/
void setBC(int k) {
	int j;
	switch (k) {
	case 0:
		//set boundary for p
		for (j = 0; j < N + 1; j++) {
			p[0][j] = p[1][j];
			p[N][j] = p[N - 1][j];
			p[j][0] = p[j][1];
			p[j][N] = p[j][N - 1];
		}
		break;
	case 1:
		//set boundary for u
		for (j = 0; j < N; j++) {
			u[0][j] = 2 - u[1][j];
			u[N][j] = -u[N - 1][j];
			u[j][0] = 0;
			u[j][N - 1] = 0;
		}
		break;
	case 2:
		//set boundary for v
		for (j = 0; j < N; j++) {
			v[0][j] = 0;
			v[N - 1][j] = 0;
			v[j][0] = -v[j][1];
			v[j][N] = -v[j][N - 1];
		}
		break;
	}
};

double d_u[N + 1][N] = { 0.0 };
double d_v[N][N + 1] = { 0.0 };
double a_u[N + 1][N] = { 0.0 };
double a_v[N][N + 1] = { 0.0 };
double u_E, u_W, v_N, v_S;
double a_E, a_W, a_N, a_S, a_P;
double temp = 0;
void pseudo(int k) {
	int i, j;
	switch (k) {
	case 1:
		//Pre-cal for u
		for (i = 1; i < N; i++) {
			for (j = 1; j < N - 1; j++) {
				scheme(U, SCHEME, i, j);
				a_E = -u_E * h / 2.0 + nu;
				a_W = u_W * h / 2.0 + nu;
				a_N = -v_N * h / 2.0 + nu;
				a_S = v_S * h / 2.0 + nu;
				a_u[i][j] = (u_E - u_W + v_N - v_S) * h / 2.0 + 4 * nu;
				d_u[i][j] = h / a_u[i][j];
				temp = u[i][j];
				u[i][j] = (a_E * u[i][j + 1] + a_W * u[i][j - 1] + a_N * u[i - 1][j] + a_S * u[i + 1][j] - h * (p[i][j + 1] - p[i][j])) / a_u[i][j];
				u[i][j] = temp + SOR_s * (u[i][j] - temp);
			}
		}
		break;

	case 2:
		//Pre-cal for v
		for (i = 1; i < N - 1; i++) {
			for (j = 1; j < N; j++) {
				scheme(V, SCHEME, i, j);
				a_E = -u_E * h * 0.5 + nu;
				a_W = u_W * h * 0.5 + nu;
				a_N = -v_N * h * 0.5 + nu;
				a_S = v_S * h * 0.5 + nu;
				a_v[i][j] = (u_E - u_W + v_N - v_S) * h * 0.5 + 4 * nu;
				d_v[i][j] = h / a_v[i][j];
				temp = v[i][j];
				v[i][j] = (a_E * v[i][j + 1] + a_W * v[i][j - 1] + a_N * v[i - 1][j] + a_S * v[i + 1][j] - h * (p[i][j] - p[i + 1][j])) / a_v[i][j];//
				v[i][j] = temp + SOR_s * (v[i][j] - temp);
			}
		}
		break;
	}

};

void corr() {
	int i, j;
	double u_corr, v_corr;
	double b;
	for (i = 1; i < N; i++) {
		for (j = 1; j < N; j++) {
			p_corr[i][j] = { 0.0 };
		}
	}
	//correct p
	for (i = 1; i < N; i++) {
		for (j = 1; j < N; j++) {
			b = h * (u[i][j] - u[i][j - 1] + v[i - 1][j] - v[i][j]);
			a_P = h * d_u[i][j - 1] + h * d_u[i][j] + h * d_v[i - 1][j] + h * d_v[i][j];
			p_corr[i][j] = (p_corr[i][j - 1] * h * d_u[i][j - 1] + p_corr[i][j + 1] * h * d_u[i][j] + p_corr[i - 1][j] * h * d_v[i - 1][j] + p_corr[i + 1][j] * h * d_v[i][j] - b) / a_P;
			p[i][j] = p[i][j] + SOR_p * p_corr[i][j];

		}
	}
	setBC(P);
	//correct u
	for (i = 1; i < N; i++) {
		for (j = 1; j < N - 1; j++) {
			u_corr = d_u[i][j] * (p_corr[i][j] - p_corr[i][j + 1]);
			temp = u[i][j];
			u[i][j] += u_corr;
			temp = u[i][j];
		}
	}
	setBC(U);
	//correct v
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N; j++) {
			v_corr = d_v[i][j] * (p_corr[i + 1][j] - p_corr[i][j]);
			v[i][j] += v_corr;
		}
	}
	setBC(V);
};

double residual() {
	int i, j;
	dm = 0;
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			dm += abs(u[i][j - 1] - u[i][j] + v[i][j] - v[i - 1][j]);
		}
	}
	return dm;
}

double r_e, r_w, r_n, r_s;
double phi_e, phi_w, phi_n, phi_s;
void scheme(int k, int m, int i, int j) {
	if (k == U) {
		switch (m) {
		case 1:

			r_e = (u[i][j + 1] - u[i][j]) / (u[i][j] - u[i][j - 1]);
			r_w = (u[i][j] - u[i][j - 1]) / (u[i][j - 1] - u[i][j - 2]);
			r_n = (v[i - 1][j + 1] - v[i - 1][j]) / (v[i - 1][j] - v[i - 1][j - 1]);
			r_s = (v[i][j + 1] - v[i][j]) / (v[i][j] - v[i][j - 1]);
			//phi
			phi_e = r_e;
			phi_n = r_n;
			phi_s = r_s;
			phi_w = r_w;
			//
			u_E = u[i][j] + 0.5 * phi_e * (u[i][j] - u[i][j - 1]);
			u_W = u[i][j - 1] + 0.5 * phi_w * (u[i][j - 1] - u[i][j - 2]);
			v_N = v[i - 1][j] + 0.5 * phi_n * (v[i - 1][j] - v[i - 1][j - 1]);
			v_S = v[i][j] + 0.5 * phi_s * (v[i][j] - v[i][j - 1]);

			if (u[i][j] == u[i][j - 1])
				u_E = 0;
			if ((u[i][j - 1] == u[i][j - 2]))
				u_W = 0;
			if (v[i - 1][j] == v[i - 1][j - 1])
				v_N = 0;
			if ((v[i][j] == v[i][j - 1]))
				v_S = 0;

		case 2:
			r_e = (u[i][j + 1] - u[i][j]) / (u[i][j] - u[i][j - 1]);
			r_w = (u[i][j] - u[i][j - 1]) / (u[i][j - 1] - u[i][j - 2]);
			r_n = (v[i - 1][j + 1] - v[i - 1][j]) / (v[i - 1][j] - v[i - 1][j - 1]);
			r_s = (v[i][j + 1] - v[i][j]) / (v[i][j] - v[i][j - 1]);
			//phi
			phi_e = max(0.0, min(2 * r_e, min(2.0, 2 * (r_e + abs(r_e)) / (r_e + 3))));
			phi_n = max(0.0, min(2 * r_n, min(2.0, 2 * (r_n + abs(r_n)) / (r_n + 3))));
			phi_s = max(0.0, min(2 * r_s, min(2.0, 2 * (r_s + abs(r_s)) / (r_s + 3))));
			phi_w = max(0.0, min(2 * r_w, min(2.0, 2 * (r_w + abs(r_w)) / (r_w + 3))));
			//
			u_E = u[i][j] + 0.5 * phi_e * (u[i][j] - u[i][j - 1]);
			u_W = u[i][j - 1] + 0.5 * phi_w * (u[i][j - 1] - u[i][j - 2]);
			v_N = v[i - 1][j] + 0.5 * phi_n * (v[i - 1][j] - v[i - 1][j - 1]);
			v_S = v[i][j] + 0.5 * phi_s * (v[i][j] - v[i][j - 1]);
			break;
		case 3:
			r_e = (u[i][j + 1] - u[i][j]) / (u[i][j] - u[i][j - 1]);
			r_w = (u[i][j] - u[i][j - 1]) / (u[i][j - 1] - u[i][j - 2]);
			r_n = (v[i - 1][j + 1] - v[i - 1][j]) / (v[i - 1][j] - v[i - 1][j - 1]);
			r_s = (v[i][j + 1] - v[i][j]) / (v[i][j] - v[i][j - 1]);
			//phi
			phi_e = max(0.0, min(2 * r_e, min(2.0, (r_e + 1) / 2)));
			phi_n = max(0.0, min(2 * r_n, min(2.0, (r_n + 1) / 2)));
			phi_s = max(0.0, min(2 * r_s, min(2.0, (r_s + 1) / 2)));
			phi_w = max(0.0, min(2 * r_w, min(2.0, (r_w + 1) / 2)));
			//
			u_E = u[i][j] + 0.5 * phi_e * (u[i][j] - u[i][j - 1]);
			u_W = u[i][j - 1] + 0.5 * phi_w * (u[i][j - 1] - u[i][j - 2]);
			v_N = v[i - 1][j] + 0.5 * phi_n * (v[i - 1][j] - v[i - 1][j - 1]);
			v_S = v[i][j] + 0.5 * phi_s * (v[i][j] - v[i][j - 1]);
			break;

		}
	}
	else if (k == V) {
		switch (m) {
		case 1:
			u_E = (u[i][j] + u[i + 1][j]) * 0.5;
			u_W = (u[i][j - 1] + u[i + 1][j - 1]) * 0.5;
			v_N = (v[i - 1][j] + v[i][j]) * 0.5;
			v_S = (v[i][j] + v[i + 1][j]) * 0.5;
			break;
		case 2:
			r_e = (u[i][j] - u[i + 1][j]) / (u[i + 1][j] - u[i + 2][j]);
			r_w = (u[i][j - 1] - u[i + 1][j - 1]) / (u[i + 1][j - 1] - u[i + 2][j - 1]);
			r_n = (v[i - 1][j] - v[i][j]) / (v[i][j] - v[i + 1][j]);
			r_s = (v[i][j] - v[i + 1][j]) / (v[i + 1][j] - v[i + 2][j]);
			//phi
			phi_e = max(0.0, min(2 * r_e, min(2.0, 2 * (r_e + abs(r_e)) / (r_e + 3))));
			phi_n = max(0.0, min(2 * r_n, min(2.0, 2 * (r_n + abs(r_n)) / (r_n + 3))));
			phi_s = max(0.0, min(2 * r_s, min(2.0, 2 * (r_s + abs(r_s)) / (r_s + 3))));
			phi_w = max(0.0, min(2 * r_w, min(2.0, 2 * (r_w + abs(r_w)) / (r_w + 3))));
			//
			u_E = u[i + 1][j] + 0.5 * phi_e * (u[i + 1][j] - u[i + 2][j]);
			u_W = u[i + 1][j - 1] + 0.5 * phi_w * (u[i + 1][j - 1] - u[i + 2][j - 1]);
			v_N = v[i][j] + 0.5 * phi_n * (v[i][j] - v[i + 1][j]);
			v_S = v[i + 1][j] + 0.5 * phi_s * (v[i + 1][j] - v[i + 2][j]);
			break;
		case 3:
			r_e = (u[i][j] - u[i + 1][j]) / (u[i + 1][j] - u[i + 2][j]);
			r_w = (u[i][j - 1] - u[i + 1][j - 1]) / (u[i + 1][j - 1] - u[i + 2][j - 1]);
			r_n = (v[i - 1][j] - v[i][j]) / (v[i][j] - v[i + 1][j]);
			r_s = (v[i][j] - v[i + 1][j]) / (v[i + 1][j] - v[i + 2][j]);
			//phi
			phi_e = max(0.0, min(2 * r_e, min(2.0, (r_e + 1) / 2)));
			phi_n = max(0.0, min(2 * r_n, min(2.0, (r_n + 1) / 2)));
			phi_s = max(0.0, min(2 * r_s, min(2.0, (r_s + 1) / 2)));
			phi_w = max(0.0, min(2 * r_w, min(2.0, (r_w + 1) / 2)));
			//
			u_E = u[i + 1][j] + 0.5 * phi_e * (u[i + 1][j] - u[i + 2][j]);
			u_W = u[i + 1][j - 1] + 0.5 * phi_w * (u[i + 1][j - 1] - u[i + 2][j - 1]);
			v_N = v[i][j] + 0.5 * phi_n * (v[i][j] - v[i + 1][j]);
			v_S = v[i + 1][j] + 0.5 * phi_s * (v[i + 1][j] - v[i + 2][j]);
			break;
		}
	}

}



void result() {
	int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			u_final[i][j] = (u[i][j] + u[i + 1][j]) * 0.5;
			v_final[i][j] = (v[i][j] + v[i][j + 1]) * 0.5;
			p_final[i][j] = (p[i][j] + p[i][j + 1] + p[i + 1][j] + p[i + 1][j + 1]) * 0.25;
			velocity_final[i][j] = sqrt(u_final[i][j] * u_final[i][j] + v_final[i][j] * v_final[i][j]);
		}
	}
	write(&velocity_final[0][0], N, N, name1);
	write(&u_final[0][0], N, N, name2);
	write(&v_final[0][0], N, N, name3);
	write(&p_final[0][0], N, N, name4);
}

