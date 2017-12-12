//
//
//
//				Ising model on the three dimensional square lattice
//					class with lattice flips evolution methods and discrete lattice operators to send
//						the boundary effect to infininity and sample a bulk sublattice marginal
//
//
//


#ifndef ISING3D_H
#define ISING3D_H

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <fstream>
#include <bitset>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
using namespace std;


#define J_CRIT 0.221654339

typedef tuple<int, int, int> site; //tuple type holding a site on the lattice
typedef boost::mt19937 RNGType; //our boost random generator


class Ising3d{

public:
int h, v, p; //lattice size heigh, width, depth
int N; //total number of spins
double J; //spin spin coupling
double mag; //average magnetization
double enrg; //average energy
int *** tab; //the array containing the spin configuration
double evoW; // the wolff-SW bonding probability = 1-e^{-2J}
mt19937 gen;
vector<site> VirtualClusterSpins; //stack of spins to be flipped in SW/Wolff evolution steps
RNGType rng;
boost::bernoulli_distribution<> * disB;
boost::variate_generator<RNGType,boost::bernoulli_distribution<> > * linkDice;
boost::uniform_real<double> * dis;
boost::variate_generator<RNGType,boost::uniform_real<double> > * metroDice;

Ising3d(const int, const int, const int); //constructor for critical lattice in hot (disordered) state
Ising3d(const char *); //loads lattice from file
~Ising3d();

void export_configuration(const char *); //exports the configuration into a file
void coolDown(const int ); //cools down the lattice in a vaccum state
void calcMag(); //calculates the average magnetization
void calcEnrg(); //calculates the total energy
int Bnei(const int, const int, const int); //returns the total neighbour magnetic field -- sum of the spin values of the first neighbours
int Edensity(const int, const int, const int); //returns the lattice energy density at specific site location

void flipSWperio(const int); //lattice flips using periodic boundary conditions
void flipSW(const int); //lattice flips using fixed boundary conditions

void dilation(const double); //discrete lattice dilation of the center of the cube sublattice -- uses a random parameter to shift the mesh
void dilationWithPixels(const double); //discrete lattice dilation using duplicates (pixellated version)

double s_distance(const int, const int, const int, const int, const int, const int);
vector<site> * giveMeTheNeighbours( const int);

void cuttingInLattice( const int); //cuts in a central subsection of the lattice
void blowUp();
void dilationHybrid(const double , const double ); //hybrid dilation mixing holes and pixels

void thermalizeLattice( const vector<double> &, const vector<int> &); //thermalizeLattices the system from a number of dilation values and lattice flips numbers

};

#endif
