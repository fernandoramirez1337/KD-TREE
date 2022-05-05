#include <cstdarg>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "KDTree.hpp"

using namespace std;

template <size_t N, typename IteratorType>
Point<N> point_from_range(IteratorType begin, IteratorType end) {
  Point<N> result;
  std::copy(begin, end, result.begin());
  return result;
}

int main() {

  const double data_points[][4] = {
      {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 0, 1, 1},
      {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 1},
      {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 0}, {1, 0, 1, 1},
      {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 1, 1, 1},
  };

  KDTree<4, size_t> test01;

  for (size_t i = 0; i < 16; ++i) {
    test01.insert(point_from_range<4>(data_points[i], data_points[i] + 4), i);
  }

  cout<< "Elementos insertados" <<endl;

  Point<4> testpoint01;
  testpoint01[0] = 0.0;
  testpoint01[1] = 0.0;
  testpoint01[2] = 0.0;
  testpoint01[3] = 0.0;


  if(test01.contains(testpoint01)) {
    cout<<"Contiene 0,0,0,0"<<endl;
  } else {
    cout<<"No contiene 0,0,0,0"<<endl;
  }

  
  return 0;
}
