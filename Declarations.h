#ifndef VECTORS
#define VECTORS
#include <vector>
#endif
#ifndef MATH
#define MATH
#include <cmath>
#endif
#ifndef FILE_OUTPUT
#define FILE_OUTPUT
#include <fstream>
#endif
#ifndef STRING
#define STRING
#include <string>
#endif
#include <cfloat>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <limits>
#include <functional>
#include <numeric> 
#include <algorithm>
#include <list>
#include <iterator>
#include <queue> 
#include <ctime>
#include <iomanip>

#ifndef SAMPLED_DISTRIBUTION
#define SAMPLED_DISTRIBUTION
#include <random>

namespace myrand
{
  extern std::mt19937 rng;
}

#endif

using namespace std;

struct Transition;
class EventQueue;
struct Propensities;

double DoubleUniformDistribution(double min, double max);
void print_vector(vector<double> *v, ofstream *output, double time, double print_time);
void UpdatePropensities(vector<double> *u, vector<double> *v,Propensities *propensities, EventQueue *events, int voxel, bool erase, double time);
void GrayScottSimulation(vector<ofstream> *outputs);
