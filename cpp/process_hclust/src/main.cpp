#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric> 
#include <cmath>

#include "hclust-cpp/fastcluster.h"

#include "vite2.cpp"

#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

struct lcprops {
  int m = 0;
  int n = 0;
};

double sum(std::vector<double> &v) {
  if (v.size() == 0) {
    printf("Warning: mean of vector with zero size\n");
    return 0;
  }
  double mv = 0;
  for (int i = 0; i < v.size(); i++)
    mv += v[i];
  return mv;
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

template<typename T>
void intersection(
  std::vector<T> &a, std::vector<T> &b, std::vector<T> &result
) {
  std::set_intersection(
    a.begin(), a.end(),
    b.begin(), b.end(),
    back_inserter(result)
  );
}

bool search_key(std::map<int, std::vector<int> > &a, const int &key) {
  for (std::map<int, std::vector<int> >::iterator f=a.begin(); f != a.end(); ++f) {
    if (f->first == key) return true;
  }
  return false;
}

bool search_key(std::map<int, lcprops > &a, const int &key) {
  for (std::map<int, lcprops>::iterator f=a.begin(); f != a.end(); ++f) {
    if (f->first == key) return true;
  }
  return false;
}

bool cmp( std::pair<int, std::vector<int> > &a, std::pair <int, std::vector<int> > &b ){
   return a.second[0] > b.second[0];
}

std::map<int, std::vector<int>> sort_map(std::map<int, std::vector<int>> givenMap){
   std::vector<std::pair<int, std::vector<int> > > pairVec;
   std::map<int, std::vector<int>> newMap;
   for ( auto& it : givenMap ) {
      pairVec.push_back( it );
   }
   sort( pairVec.begin(), pairVec.end(), cmp); 
   for ( auto& it : pairVec ) {
      newMap.insert( { it.first, it.second } );
   }
   return newMap;
}

std::vector<double> simplify_height_to_k_end(
  int &n,
  double* height,
  std::vector<double>& sim_height,
  int  &size
) {
  double h = height[0];
  std::vector<double> sim_k;
  for (int i = 0; i < n - 1; i++) {
    if (i < n - 2) {
      if (height[i + 1] != h) {
        sim_k.push_back(n - i);
        sim_height.push_back(h);
        h = height[i + 1];
        ++(size);
      }
    } else {
      if (height[i] != height[i - 1]) {
        sim_k.push_back(n - i);
        h = height[i];
        sim_height.push_back(h);
        ++(size);
      }
    }
    
  }
  return sim_k;
}

std::vector<double> simplify_height_to_k_start(
  int &n,
  double* height,
  std::vector<double>& sim_height,
  int &size
) {
  double h = height[0];
  std::vector<double> sim_k;
  for (int i = 0; i < n - 1; i++) {
    if (i == 0) {
      sim_k.push_back(n - 1);
      sim_height.push_back(h);
      ++size;
    }
    if (height[i] != h && i != 0) {
      h = height[i];
      sim_k.push_back(n - i - 1);
      sim_height.push_back(h);
      ++size;
    }
  }
  return sim_k;
}

std::vector<double> complete_height_to_k(
  int &n,
  double* height,
  std::vector<double>& sim_height,
  int &size
) {
  std::vector<double> sim_k;
  for (int i = 0; i < n - 1; i++) {
    sim_k.push_back(n - i - 1);
    sim_height.push_back(height[i]);
    size++;
}
  return sim_k;
}

double Dc(int &m, int &n, bool &undirected) {
  double dc;
  if (!undirected)
    dc = (m - n + 1.) / pow(n - 1., 2.);
  else
    dc = (m - n + 1.) / ((n * (n - 1.) / 2.) - n + 1.);
  if (dc > 0) return dc;
  else return 0;
}

double Sc(int &m, int &n, int &M, int& N) {
  double pc;
  pc = (m - n + 1.) / (M - N + 1.);
  if (pc > 0) return -pc * log(pc);
  else return 0;
}

double Xsus(std::map<int, lcprops> &v, int &N, int &order) {
  double  x = 0; // percolation suceptability
  for (std::map<int, lcprops >::iterator it = v.begin(); it != v.end(); ++it) {
    if (it->second.m <= 1 || it->second.n <= 2) continue;
    if (it->second.m != order)
      x += pow(it->second.m, 2.);
  }
  return x / N;
}

double order_parameter(std::vector<int> &v, int &M) {
  return v[0] / (M * 1.);
}

double Xm(std::map<int, lcprops> &v) {
  double n = v.size();
  double xm2 = 0, xm = 0;
  std::map<int, int> v_count;
  for (std::map<int, lcprops>::iterator it = v.begin(); it != v.end(); ++it) {
    if (!search_key(v, it->second.m)) v_count[it->second.m] = 1;
    else v_count[it->second.m]++;
  }
  for (std::map<int, int >::iterator it = v_count.begin(); it != v_count.end(); ++it) {
    xm2 += pow(it->first, 2.0) * it->second;
    xm += it->first * it->second * 1.;
  }
  if (xm > 0)
    return xm2 / xm;
  else return 0;
}

class ph {
  private:

    std::vector<int> K;
    std::vector<double> Height;
    std::vector<int> NEC;
    std::vector<double> MU;
    std::vector<double> D;
    std::vector<int> ntrees;
    std::vector<double> X;
    std::vector<double> OrP;
    std::vector<double> XM;
    std::vector<double> S;
    std::vector<double> CC;
    std::vector<double> elinks;

    int max_level=0;

    int number_of_elements;
    std::vector<double> distane_matrix;
    std::vector<int> source_vertices;
    std::vector<int> target_vertices;
    int total_nodes;
    int Linkage;
    bool CUT;
    bool undirected;

  public:
    ph(
      const int n,
      std::vector<double> distmat,
      std::vector<int> source,
      std::vector<int> target,
      const int nodes,
      const int linkage,
      const bool cut,
      const bool uni
    );
    ~ph(){};

    void vite2();
    // std::vector<std::vector<int> > &link_communities, std::vector<double> & H, 
    template <typename T>
    void expand_vector(std::vector<T>& v, const int& N);

    std::vector<int> get_K();
    std::vector<double> get_Height();
    std::vector<int> get_NEC();
    std::vector<double> get_MU();
    std::vector<double> get_D();
    std::vector<int> get_ntrees();
    std::vector<double> get_X();
    std::vector<double> get_OrP();
    std::vector<double> get_XM();
    std::vector<double> get_S();
    std::vector<double> get_CC();
    std::vector<double> get_elinks();
    int get_max_level();
    void get_sizes(
      std::map<int, lcprops> &info_sizes,
      int* labels, std::vector<int> &lcsize,
      std::vector<int>& unique_labels,
      std::vector<int>& source,
      std::vector<int>& target,
      int& n
    );
};

ph::ph(
  const int n,
  std::vector<double> distmat,
  std::vector<int> source,
  std::vector<int> target,
  const int nodes,
  const int linkage,
  const bool cut,
  const bool uni
) {
  number_of_elements = n;
  distane_matrix = distmat;
  source_vertices = source;
  target_vertices = target;
  total_nodes = nodes;
  Linkage = linkage;
  CUT = cut;
  undirected = uni;
}

template <typename T>
void ph::expand_vector(std::vector<T>& v, const int& N) {
  v = std::vector<T>(N, 0);
}

std::vector<int> ph::get_K() {
  return K;
}
std::vector<double> ph::get_Height() {
  return Height;
}
std::vector<int> ph::get_NEC() {
  return NEC;
}
std::vector<double> ph::get_MU() {
  return MU;
}
std::vector<double> ph::get_D() {
  return D;
}
std::vector<int> ph::get_ntrees() {
  return ntrees;
}
std::vector<double> ph::get_X() {
  return X;
}
std::vector<double> ph::get_OrP() {
  return OrP;
}
std::vector<double> ph::get_XM() {
  return XM;
}

std::vector<double> ph::get_S() {
  return S;
}

std::vector<double> ph::get_CC() {
  return CC;
}

std::vector<double> ph::get_elinks() {
  return elinks;
}


int ph::get_max_level(){
  return max_level;
}

void ph::get_sizes(
  std::map<int, lcprops> &info_sizes,
  int* labels, std::vector<int> &lcsize,
  std::vector<int>& unique_labels,
  std::vector<int>& source,
  std::vector<int>& target,
  int& n
) {
  std::vector<std::set<int> > node_buffer(unique_labels.size());
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < unique_labels.size(); j++) {
      if (labels[i] == unique_labels[j]) {
        info_sizes[unique_labels[j]].m++;
        lcsize[j]++;
        node_buffer[j].insert(source[i]);
        node_buffer[j].insert(target[i]);
        break;
      }
    }
  }
  for (int j = 0; j < unique_labels.size(); j++)
    info_sizes[unique_labels[j]].n = node_buffer[j].size();
  sort(lcsize.begin(), lcsize.end(), std::greater<int>());
}

void ph::vite2() {
  // Various variables ----
  int nt, nec, length = 0, m, n, el;
  double dc, mtree;
  std::vector<double> sim_k, sim_height, dcv, scv;
  std::vector<int> lcsizes;
  std::map<int, link_groups> groups;
  // Condense distance matrix ----
  double* tri_distmat = new double[(number_of_elements * (number_of_elements - 1)) / 2];
  // hclust arrays ----
  int* merge = new int[2 * (number_of_elements - 1)];
  double* height = new double[number_of_elements-1];
  int* labels = new int[number_of_elements];
  // Get condense matrix ----
  for (int i = 0; i < distane_matrix.size(); i++) {
    tri_distmat[i] = distane_matrix[i];
  }
  // Get hierarchy!! ----
  hclust_fast(
    number_of_elements,
    tri_distmat,
    Linkage, // linkage method
    merge,
    height
  );
  if (CUT) {
    // Delete duplicated heights preserving the first K and height ----
    sim_k = simplify_height_to_k_start(number_of_elements, height, sim_height, length);
  } else {
    // Keep the all the steps ----
    sim_k = complete_height_to_k(number_of_elements, height, sim_height, length);
  }
  expand_vector(K, length);
  expand_vector(Height, length);
  expand_vector(NEC, length);
  expand_vector(D, length);
  expand_vector(S, length);

  expand_vector(CC, length);
  expand_vector(elinks, length);
  // THE GAME STARTS
  for (int i=0; i < length; i++) {
    K[i] = sim_k[i];
    Height[i] = sim_height[i];
    // Cut tree at given sim_k and get memberships ----
    cutree_k(
      number_of_elements,
      merge,
      K[i],
      labels
    );
    // Get number of links and number of nodes in link communities in order
    groups = group_linkcommunities(labels, source_vertices, target_vertices, number_of_elements);
    
    mtree = 0.;
    nec = 0;
    dcv = std::vector<double>(groups.size(), 0.);
    scv = std::vector<double>(groups.size(), 0.);

    el = 0;

    for (std::map<int, link_groups>::iterator it=groups.begin(); it != groups.end(); ++it) {
      m = it->second.edgelist.size();
      n = it->second.link_nodes.size();
      if (m > 1 && n > 2) {
        dcv[nec] = Dc(m, n, undirected) * m / number_of_elements;
        scv[nec] = Sc(m, n, number_of_elements, total_nodes);

        CC[i] += global_clustering_coefficient(it->second.link_nodes, it->second.edgelist);
        mtree += m - n + 1;
        nec++;
        el += m;
      }
    }

    S[i] = sum(scv);

    // CC[i]/= nec;
    elinks[i] = el;

    mtree = (number_of_elements - total_nodes + 1.) - mtree;
    mtree = mtree / (number_of_elements  - total_nodes + 1.);
    if (mtree > 0) S[i] += -mtree * log(mtree);
    NEC[i] = nec;
    D[i] = sum(dcv);
    // CC[i] = lc_mean_cc(groups);
    // std::cout << CC[i] << " ";
    groups.clear();
  }
  // std::cout<<"\n";
  // Delete pointers
  delete[] labels;
  delete[] merge;
  delete[] height;
  delete[] tri_distmat;
}

PYBIND11_MODULE(process_hclust, m) {
    py::class_<ph>(m, "ph")
        .def(
          py::init<
            const int,
            std::vector<double>,
            std::vector<int>,
            std::vector<int>,
            const int,
            const int,
            const bool,
            const bool
          >()
        )
        .def("vite2", &ph::vite2)
        .def("get_K", &ph::get_K)
				.def("get_Height", &ph::get_Height)
				.def("get_NEC", &ph::get_NEC)
        .def("get_MU", &ph::get_MU)
				.def("get_D", &ph::get_D)
			  .def("get_ntrees", &ph::get_ntrees)
        .def("get_X", &ph::get_X)
        .def("get_OrP", &ph::get_OrP)
			  .def("get_XM", &ph::get_XM)
        .def("get_S", &ph::get_S)
        .def("get_CC", &ph::get_CC)
        .def("get_elinks", &ph::get_elinks)
        .def("get_max_level", &ph::get_max_level);
}