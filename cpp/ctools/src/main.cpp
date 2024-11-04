#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include<ctime> // time

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


double Hellinger2(
	std::vector<double> &u, std::vector<double> &v, int &ii, int &jj
) {
	int N = u.size(), k=0;
	double p = 0, pu = 0., pv = 0., maxp=0, s=0.;
	std::vector<bool> possible;
	std::vector<double> ppu(N, 0.), ppv(N, 0.), peff(N, 0.);

	for (int j=0; j < N; j++){
		pu += u[j];
		pv += v[j];

		if (j == ii || j == jj) continue;
		ppu[k] = u[j];
		ppv[k] = v[j];
		k++;
	}
	ppu[N-2] = u[ii];
	ppu[N-1] = u[jj];
	ppv[N-2] = v[jj];
	ppv[N-1] = v[ii];

	if (pu == 0 || pv == 0) return  0.;

	for (int j=0; j < N; j++) {
		if (ppu[j] > 0 && ppv[j] > 0) {
			peff[j] = 0.5 * (log(ppu[j]) + log(ppv[j]) - log(pu) -log(pv));
			possible.push_back(true);
		}	else {
			possible.push_back(false);
		}
		
	}
	
	for (int j=0; j < N; j++) {
		if (possible[j])  {
			if (maxp > peff[j]) maxp = peff[j];
		}
	}
	if (maxp == 0) return 0.;

	for (int j=0; j < N; j++) {
		if (possible[j]) {
			s += exp(peff[j]-maxp);
		}
	}
	return abs(s * exp(maxp));
}

double cosine_similarity(
	std::vector<double> &u, std::vector<double> &v, int &ii, int &jj
) {
	int N = u.size();
	double uv=0., uu=0., vv=0.;
	for (int i=0; i < N; i++) {
		uu += u[i] * u[i];
		vv += v[i] * v[i];
		if (i == ii | i == jj) continue;
		uv += u[i] * v[i];
	}

	if (ii < N && jj < N) {
		uv += u[jj] * v[ii];
		uv += u[ii] * v[jj]; 
	}

	return uv / (sqrt(uu * vv));
}

PYBIND11_MODULE(ctools, m) {

  m.doc() = "Creates fast random networks";

	m.def(
    "Hellinger2",
    &Hellinger2,
    py::return_value_policy::reference_internal
  );

	m.def(
    "cosine_similarity",
    &cosine_similarity,
    py::return_value_policy::reference_internal
  );
}