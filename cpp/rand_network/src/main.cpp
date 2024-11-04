#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <ctime> // time
#include <map>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

template<typename T>
void display_progress(T &ini, T &total, int &carrier, int &sp) {
  int true_percent;
  double percent = static_cast<double>(ini) / static_cast<double>(total);
  percent *= 100;
  percent = floor(percent);
  true_percent = (int) percent;
  if (true_percent % sp == 0 && carrier != true_percent) {
    std::cout << true_percent << "%   ";
    carrier = true_percent;
  }
}

void network_density(
  std::vector<std::vector<int> > &A, int &nodes, double &den
) {
  // Create variables ----
  den = 0;
  for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      if (A[i][j] > 0)
        den++;
    }
  }
  den /= static_cast<double>(nodes * (nodes - 1));
}

void network_M(
  std::vector<std::vector<int> > &A, int &nrows, int&ncols, int &m
) {
  // Create variables ----
  m = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (A[i][j] > 0)
        m++;
    }
  }
}

void network_count(
  std::vector<std::vector<int> > &A, int &nodes, int &counter
) {
  // Create variables ----
  counter = 0;
  for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      counter += A[i][j];
    }
  }
}

void network_count_M(
  std::vector<std::vector<int> > &A, int &nrows, int &ncols, int &counter
) {
  // Create variables ----
  counter = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      counter += A[i][j];
    }
  }
}

void network_links(
  std::vector<std::vector<int> > &A, int &nodes, double &counter
) {
  // Create variables ----
  counter = 0;
  for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      if (A[i][j] > 0)
        counter++;
    }
  }
}

template<typename T>
void unique(std::vector<T> &v) {
  std::vector<int>::iterator ip;
  std::sort(v.begin(), v.end());
  ip = std::unique(
    v.begin(),
    v.begin() + v.size()
  );
  v.resize(
    std::distance(v.begin(), ip)
  );
}

// distbase***

void  get_subnets(
  std::vector<std::vector<double> > &net,
  std::vector<double> &bins,
  int &leaves, int& nbins,
  std::vector<std::vector<std::vector<double> > > &subnet
) {
  // Declare variables ----
  std::vector<std::vector<double> > sub;
  // Start loops ----
  for (int i = 0; i < nbins - 1; i++) {
    for (int j = 0; j < leaves; j++) {
      if (i < nbins - 2) {
        if (bins[i] <= net[j][2] && bins[i + 1] > net[j][2])
          sub.push_back(net[j]);
      } else {
        if (bins[i] <= net[j][2] && bins[i + 1] + 0.1 > net[j][2])
          sub.push_back(net[j]);
      }
      
    }
    subnet.push_back(sub);
    sub.clear();
  }
}

std::vector<std::vector<double> > find_bin(double &d, std::vector<std::vector<double> >  &net) {
  std::vector<std::vector<double> > subnet;
  for (int i = 0; i < net.size(); i++) {
    if (d >= net[i][2] - 1e-2 && d < net[i][2] + 1e-2)
      subnet.push_back(net[i]);  
  }
  return subnet;
}

int find_bin(double &d, std::vector<double> &bins, int & nbins) {
  int c = -1;
  for (int i = 0; i < nbins - 1; i++) {
    if (i < nbins - 2) {
      if (bins[i] <= d && bins[i + 1] > d) {
        c = i;
        return c;
      }
    } else {
      if (bins[i] <= d && bins[i + 1] + 0.1 > d) {
        c = i; 
        return c;
      }
    }
   
  }
  return c;
}

// distbase: EDR network generators ***

std::vector<std::vector<int> > sample_distbase(
  std::vector<std::vector<double> > &net,
  std::vector<double> &bins, int &nbins,
  int &nodes, int &ncols, int &leaves,
  double &rho, double &lb, double& loc
) {
  // Change random seed ----
  srand(time(0));
  // Declare variables ----
  int c, r, A, B, t = 0, sp = 5;
  double Rho = 0.0, d;
  // Declare network ----
  std::vector<std::vector<int> > NET(
    nodes, std::vector<int>(nodes, 0)
  );
  // Categorize distances ----
  std::vector<std::vector<std::vector<double> > > subnet;
  get_subnets(
    net, bins, leaves, nbins, subnet
  );
  std::cout << "Matching network density:\n";
  while (Rho < rho) {
    display_progress(Rho, rho, t, sp);
    // Generate random distance
    d =  (rand() % 1000) / 1000.;
    d = - log(1 - d) / lb + loc;
    if (d > bins[nbins - 1]) continue;
    c = find_bin(d, bins, nbins);
    if (c == -1) continue;
    if (subnet[c].size() == 0) continue;
    r = rand() % subnet[c].size();
    A = rand() % 2;
    if (A == 0) B = 1;
    else B = 0;
    NET[subnet[c][r][A]][subnet[c][r][B]]++;
    network_density(NET, ncols, Rho);
  }
  if (Rho > rho)
     NET[subnet[c][r][A]][subnet[c][r][B]]--;
  std::cout << "\nDone!\n";
  return NET;
}

std::vector<std::vector<int> > sample_distbase_trunc(
  std::vector<std::vector<double> > &net,
  std::vector<double> &bins, int &nbins,
  int &nodes, int &ncols, int &leaves,
  double &rho, double &lb, double &xmin, double &xmax
) {
  // Change random seed ----
  srand(time(0));
  // Declare variables ----
  int c, r, A, B, t = 0, sp = 5;
  double Rho = 0.0, d;
  // Declare network ----
  std::vector<std::vector<int> > NET(
    nodes, std::vector<int>(nodes, 0)
  );
  // Categorize distances ----
  std::vector<std::vector<std::vector<double> > > subnet;
  get_subnets(
    net, bins, leaves, nbins, subnet
  );
  std::cout << "Matching network density:\n";
  while (Rho < rho) {
    display_progress(Rho, rho, t, sp);
    // Generate random distance
    d =  (rand() % 100000) / 100000.;
    d = xmin - log(1. - d * (1. - exp(-lb * (xmax - xmin)))) / lb ;
    if (d > xmax) continue;
    c = find_bin(d, bins, nbins);
    if (c == -1) continue;
    if (subnet[c].size() == 0) continue;
    r = rand() % subnet[c].size();
    A = rand() % 2;
    if (A == 0) B = 1;
    else B = 0;
    NET[subnet[c][r][A]][subnet[c][r][B]]++;
    network_density(NET, ncols, Rho);
  }
  if (Rho > rho)
     NET[subnet[c][r][A]][subnet[c][r][B]]--;
  std::cout << "\nDone!\n";
  return NET;
}

std::vector<std::vector<int> > sample_distbase_M(
  std::vector<std::vector<double> > &net,
  std::vector<double> &bins, int &nbins,
  int &nrows, int &ncols, int &Drows,
  int &M, double &lb, double &loc
) {
  // Change random seed ----
  srand(time(0));
  // Declare variables ----
  int c, r, A, B, m = 0, t = 0, sp = 5;
  double d;
  // Declare network ----
  std::vector<std::vector<int> > NET(
    nrows, std::vector<int>(nrows, 0)
  );
  // Categorize distances ----
  std::vector<std::vector<std::vector<double> > > subnet;
  get_subnets(
    net, bins, Drows, nbins, subnet
  );
  std::cout << "Matching network density:\n";
  while (m <= M) {
    display_progress(m, M, t, sp);
    // Generate random distance
    d = (rand() % 1000000) / 1000000.;
    d = - log(1 - d) / lb + loc;
    if (d > bins[nbins - 1]) continue;
    c = find_bin(d, bins, nbins);
    if (c == -1) continue;
    if (subnet[c].size() == 0) continue;
    r = rand() % subnet[c].size();
    A = rand() % 2;
    if (A == 0) B = 1;
    else B = 0;
    NET[subnet[c][r][A]][subnet[c][r][B]]++;
    network_M(NET, nrows, ncols, m);
  }
  std::cout << "\nDone!\n";
  return NET;
}

std::vector<std::vector<int> > sample_distbase_trunc_M(
  std::vector<std::vector<double> > &net,
  std::vector<double> &bins, int &nbins,
  int &nrows, int &ncols, int &Drows,
  int &M, double &lb, double &loc, double &b
) {
  // Change random seed ----
  srand(time(0));
  // Declare variables ----
  int c, r, A, B, m = 0, t = 0, sp = 5;
  double d;
  // Declare network ----
  std::vector<std::vector<int> > NET(
    nrows, std::vector<int>(nrows, 0)
  );
  // Categorize distances ----
  std::vector<std::vector<std::vector<double> > > subnet;
  get_subnets(
    net, bins, Drows, nbins, subnet
  );
  std::cout << "Matching network density:\n";
  while (m <= M) {
    display_progress(m, M, t, sp);
    // Generate random distance
    d = (rand() % 1000000) / 1000000.;
    d = - log(1 - d * (1 - exp(-lb * (b - loc)))) / lb + loc;
    if (d > b) continue;
    c = find_bin(d, bins, nbins);
    if (c == -1) continue;
    if (subnet[c].size() == 0) continue;
    r = rand() % subnet[c].size();
    A = rand() % 2;
    if (A == 0) B = 1;
    else B = 0;
    NET[subnet[c][r][A]][subnet[c][r][B]]++;
    network_M(NET, nrows, ncols, m);
  }
  std::cout << "\nDone!\n";
  return NET;
}

// swapnet: generator for configuration model ***

std::vector<std::vector<double> > swap_one_k(
  std::vector<std::vector<double> > &G,
  const int &rows, const int &cols, int &swaps
) {
  // Change random seed ----
  srand(time(0));
  // Declare the shuffle nodes ----
  int A, B, C, D, keep = 0, t = 0, sp = 5;
  int luck;
  double m;
  // Copy network ----
  std::vector<std::vector<double> > GG = G;
  // Start loop!!!
  printf("Swapping edges %i times:\n", swaps);
  while (t < swaps) {
    display_progress(t, swaps, keep, sp);
    A = rand() % rows;
    B = rand() % rows;
    C = rand() % cols;
    D = rand() % cols;
    if (A == B || C == D) continue;
    if (A == D || B == C) continue;
    if (A == C || B == D) continue;
    if (GG[A][C] == 0 || GG[B][D] == 0) continue;
    if (GG[A][D] > 0 || GG[B][C] > 0) continue;
    // Create luck ----
    luck = rand() % 2;
    switch (luck) {
    // Horizontal cases ----
    case 0:
      GG[A][D] = GG[A][C];
      GG[A][C] *= 0.0;
      GG[B][C] = GG[B][D];
      GG[B][D] *= 0.0;
      break;
    // Vertical case ----
    case 1:
      GG[B][C] = GG[A][C];
      GG[A][C] *= 0.0;
      GG[A][D] = GG[B][D];
      GG[B][D] *= 0.0;
      break;
    default:
      break;
    }
    t++;
  }
  std::cout << "\nDone!\n";
  return GG;
}

std::vector<std::vector<double> > swap_dir_weights(
  std::vector<std::vector<double> > &G,
  const int &rows, const int &cols, int &swaps
) {
  // Change random seed ----
  srand(time(0));
  // Declare the shuffle nodes ----
  int A, B, C, keep = 0, t = 0, sp = 5;
  int luck;
  double m;
  // Copy network ----
  std::vector<std::vector<double> > GG = G;
  // Start loop!!!
  printf("Swapping weights %i times:\n", swaps);
  while (t < swaps) {
    display_progress(t, swaps, keep, sp);
    // Create luck ----
    luck = rand() % 2;
    switch (luck) {
    // Horizontal cases ----
    case 0:
      A = rand() % rows;
      B = rand() % rows;
      C = rand() % cols;
      if (A == B) continue;
      if (A == C || B == C) continue;
      if (GG[A][C] == 0 || GG[B][C] == 0) continue;
      m = GG[A][C];
      GG[A][C] = GG[B][C];
      GG[B][C] = m;
      break;
    // Vertical case ----
    case 1:
      A = rand() % rows;
      B = rand() % cols;
      C = rand() % cols;
      if (A == B) continue;
      if (A == C || B == C) continue;
      if (GG[A][B] == 0 || GG[A][C] == 0) continue;
      m = GG[A][B];
      GG[A][B] = GG[A][C];
      GG[A][C] = m;
      break;
    default:
      break;
    }
    t++;
  }
  std::cout << "\nDone!\n";
  return GG;
}


std::vector<std::vector<std::vector<double> > > swap_one_k_TWOMX(
  std::vector<std::vector<double> > &G, std::vector<std::vector<double> > &H,
  const int &rows, const int &cols, int &swaps
) {
  // Change random seed ----
  srand(time(0));
  // Declare the shuffle nodes ----
  int A, B, C, D, keep = 0, t = 0, sp = 5;
  int luck;
  double m;
  // Copy network ----
  std::vector<std::vector<double> > GG = G;
  std::vector<std::vector<double> > HH = H;
  // Start loop!!!
  printf("Swapping edges %i times:\n", swaps);
  while (t < swaps) {
    display_progress(t, swaps, keep, sp);
    A = rand() % rows;
    B = rand() % rows;
    C = rand() % cols;
    D = rand() % cols;
    if (A == B || C == D) continue;
    if (A == D || B == C) continue;
    if (A == C || B == D) continue;
    if (GG[A][C] == 0 || GG[B][D] == 0) continue;
    if (GG[A][D] > 0 || GG[B][C] > 0) continue;
    // Create luck ----
    luck = rand() % 2;
    switch (luck) {
    // Horizontal cases ----
    case 0:
      GG[A][D] = GG[A][C];
      GG[A][C] *= 0.0;
      GG[B][C] = GG[B][D];
      GG[B][D] *= 0.0;

      HH[A][D] = HH[A][C];
      HH[A][C] *= 0.0;
      HH[B][C] = HH[B][D];
      HH[B][D] *= 0.0;

      break;
    // Vertical case ----
    case 1:
      GG[B][C] = GG[A][C];
      GG[A][C] *= 0.0;
      GG[A][D] = GG[B][D];
      GG[B][D] *= 0.0;

      HH[B][C] = HH[A][C];
      HH[A][C] *= 0.0;
      HH[A][D] = HH[B][D];
      HH[B][D] *= 0.0;

      break;
    default:
      break;
    }
    t++;
  }
  
  std::vector<std::vector<std::vector<double> > > GH;
  GH.push_back(GG);
  GH.push_back(HH);

  std::cout << "\nDone!\n";
  return GH;
}

PYBIND11_MODULE(rand_network, m) {

  m.doc() = "Creates fast random networks";

  m.def(
    "sample_distbase",
    &sample_distbase,
    py::return_value_policy::reference_internal
  );

  m.def(
    "sample_distbase_trunc",
    &sample_distbase_trunc,
    py::return_value_policy::reference_internal
  );

  m.def(
    "sample_distbase_M",
    &sample_distbase_M,
    py::return_value_policy::reference_internal
  );

  m.def(
    "sample_distbase_trunc_M",
    &sample_distbase_trunc_M,
    py::return_value_policy::reference_internal
  );

  m.def(
    "swap_one_k",
    &swap_one_k,
    py::return_value_policy::reference_internal
  );

  m.def(
    "swap_dir_weights",
    &swap_dir_weights,
    py::return_value_policy::reference_internal
  );

  m.def(
    "swap_one_k_TWOMX",
    &swap_one_k_TWOMX,
    py::return_value_policy::reference_internal
  );

}