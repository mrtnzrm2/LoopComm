#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include<ctime> // time

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

const int INF = std::numeric_limits<int>::max();
const double INF2 = std::numeric_limits<double>::max();


std::vector<std::vector<int> > create_id_matrix(
	std::vector<std::vector<bool> >& matrix,
	const int& M, const int& N
) {
	std::vector<std::vector<int> > id_matrix(N, std::vector<int>(M, 0));
	int id = 1;
	for (int i=0; i < M; i++){
		for (int j=0; j < N; j++){
			if (matrix[i][j] > 0){
				id_matrix[i][j] = id;
				id++;
			}
		}
	}
	return id_matrix;
}

class simquest {
  private:
    std::vector<double> linksim_matrix;
    std::vector<std::vector<double> > source_matrix;
    std::vector<std::vector<double> > target_matrix;
  public:
    simquest(
      std::vector<std::vector<bool> > BA,
			std::vector<std::vector<double> > D,
      std::vector<std::vector<double> > AKI,
      std::vector<std::vector<double> > AIK,
			std::vector<std::vector<bool> > BAKI,
  		std::vector<std::vector<bool> > BAIK,
      const int N,
      const int leaves,
      const int topology,
      const int index
    );
    ~simquest(){};
    std::vector<double> calculate_linksim_matrix(
      std::vector<std::vector<bool> >& bmatrix, const int& N, const int& leaves, const int& index
    );

    double similarity_index(
      std::vector<double> &u, std::vector<double> &v, std::vector<bool> &bu, std::vector<bool> &bv, double &D, int &ii, int &jj, const int &index
    );


    std::vector<std::vector<double> > calculate_nodesim_matrix(
      std::vector<std::vector<double> >& matrix, std::vector<std::vector<bool> >& bmatrix, std::vector<std::vector<double> > D, const int& N, const int& index
    );
		
		std::vector<double> get_linksim_matrix();
		std::vector<std::vector<double> > get_source_matrix();
		std::vector<std::vector<double> > get_target_matrix();

		double dist_sim(double &d);
		double cosine_similarity(std::vector<double> &u, std::vector<double> &v, int &ii, int &jj);
		double Hellinger2(std::vector<double> &u, std::vector<double> &v, int &ii, int&jj);
};

simquest::simquest(
	std::vector<std::vector<bool> > BA,
	std::vector<std::vector<double>> D,
  std::vector<std::vector<double> > AKI,
  std::vector<std::vector<double> > AIK,
	std::vector<std::vector<bool> > BAKI,
  std::vector<std::vector<bool> > BAIK,
	const int N,
	const int leaves,
  const int topology,
  const int index
){
	// MIX topology
	if (topology == 0) {
		source_matrix = calculate_nodesim_matrix(AIK, BAIK, D, N, index);
		target_matrix = calculate_nodesim_matrix(AKI, BAKI, D, N, index);
	}
	// SOURCE topology
	else if (topology == 1) {
		source_matrix = calculate_nodesim_matrix(AIK, BAIK, D, N, index);
		target_matrix = source_matrix;
	}
	// TARGET topology
	else if (topology == 2) {
		source_matrix = calculate_nodesim_matrix(AKI, BAKI, D, N, index);
		target_matrix = source_matrix;
	}
	linksim_matrix = calculate_linksim_matrix(BA, N, leaves, index);
}

std::vector<std::vector<double> > simquest::calculate_nodesim_matrix(
	 std::vector<std::vector<double> >& matrix, std::vector<std::vector<bool> >& bmatrix, std::vector<std::vector<double> > D, const int& N, const int& index
) {
	std::vector<std::vector<double> > node_sim_matrix(N, std::vector<double>(N, 0.));

	for (int i=0; i < N; i++) {
		for (int j=i; j < N; j++) {
			if (i == j) continue;
			node_sim_matrix[i][j] = similarity_index(matrix[i], matrix[j], bmatrix[i], bmatrix[j], D[i][j], i, j, index);
			node_sim_matrix[j][i] = node_sim_matrix[i][j];
		}
	}
	return node_sim_matrix;
}

std::vector<double> simquest::calculate_linksim_matrix(
	std::vector<std::vector<bool> >& matrix, const int& N, const int& leaves, const int& index
) {
	int t;
	std::vector<std::vector<int> > id_matrix = create_id_matrix(matrix, N, N);

	t = (int) ((leaves - 1.) * leaves / 2.);
	std::vector<double> link_similarity_matrix(t, 0.);
	int col_id, row_id;
	for (int i =0; i < N; i++) {
		for (int j=0; j < N; j++) {
			row_id = id_matrix[i][j];
			if (row_id == 0) continue;
			for (int k=j; k < N; k++) {
				col_id = id_matrix[i][k];
				if (k == j || col_id == 0) continue;
				t = leaves * (row_id - 1) + col_id - 1 - 2 * (row_id -1) - 1;
				t -= (int) ((row_id - 1.) * (row_id - 2.) / 2);
				link_similarity_matrix[t] = target_matrix[j][k];
			}
			for (int k=i; k < N; k++) {
				col_id = id_matrix[k][j];
				if (k == i || col_id == 0) continue;
				t = leaves * (row_id - 1) + col_id - 1 - 2 * (row_id -1) - 1;
				t -= (int) ((row_id - 1.) * (row_id - 2.) / 2);
				link_similarity_matrix[t] = source_matrix[i][k];
			}
		}
	}
	return link_similarity_matrix;
}

double simquest::dist_sim(
	double &d
) {
	return 1 / (1 + d);
}

double simquest::cosine_similarity(
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

double simquest::Hellinger2(
	std::vector<double> &u, std::vector<double> &v, int &ii, int &jj
) {
	int N = u.size(), k=0;
	double p = 0, pu = 0., pv = 0.,  maxp=0, s=0.;
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
		}	else{
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
	s =  abs(s * exp(maxp));
	if (s < 0)
		s = 0;
	else if (s > 1)
		s = 1;
	else return s;
}

double simquest::similarity_index(std::vector<double> &u, std::vector<double> &v, std::vector<bool> &bu, std::vector<bool> &bv, double &D, int &ii, int &jj, const int &index) {
	// Jaccard probability index
  if (index == 0) {
		return cosine_similarity(u, v, ii, jj);
	}
	else if (index == 1) {
		return dist_sim(D);
	}
	else if (index == 2) {
		return Hellinger2(u, v, ii, jj);
	}
  else {
    std::range_error("Similarity index must be a integer from 0 to 5\n");
  }
}

std::vector<double> simquest::get_linksim_matrix() {
	return linksim_matrix;
}

std::vector<std::vector<double> > simquest::get_source_matrix() {
	return source_matrix;
}

std::vector<std::vector<double> > simquest::get_target_matrix() {
	return target_matrix;
}

PYBIND11_MODULE(simquest, m) {
    py::class_<simquest>(m, "simquest")
        .def(
          py::init<
           std::vector<std::vector<bool> >,
						std::vector<std::vector<double> >,
						std::vector<std::vector<double> >,
						std::vector<std::vector<double> >,
						std::vector<std::vector<bool> >,
						std::vector<std::vector<bool> >,
						const int,
						const int,
						const int,
						const int
          >()
        )
        .def("get_linksim_matrix", &simquest::get_linksim_matrix)
        .def("get_source_matrix", &simquest::get_source_matrix)
				.def("get_target_matrix", &simquest::get_target_matrix);
}