// Control: semantics of the CAPD interval comparisons used by both checkers.
// '<' and '>' between intervals must mean CERTAINLY (every element), so overlapping or touching intervals compare false.
#include <iostream>
#include "capd/capdlib.h"
using namespace capd; using namespace std;
int main() {
  interval a(0., 2.), b(1., 3.), c(2.5, 4.);
  struct T { const char* what; bool got, must; } t[] = {
    {"[0,2] < [1,3]   (overlap)  ", a < b, false},
    {"[0,2] < [2.5,4] (disjoint) ", a < c, true},
    {"[1,3] > [0,2]   (overlap)  ", b > a, false},
    {"[0,2] < 2       (touching) ", a < interval(2.), false},
    {"[-1,1] > -1     (touching) ", interval(-1., 1.) > interval(-1.), false}};
  bool ok = true;
  for (auto& x : t) { cout << x.what << " got " << x.got << ", must be " << x.must << (x.got == x.must ? "  PASS" : "  FAIL") << endl; ok = ok && x.got == x.must; }
  cout << (ok ? "COMPARISON SEMANTICS OK" : "COMPARISON SEMANTICS WRONG") << endl;
  return ok ? 0 : 1;
}
