#include "Ising3d.h"

Ising3d::Ising3d(const int Nx, const int Ny, const int Nz){
	h = Nx;
	v = Ny;
	p = Nz;
	J = J_CRIT; //coupling always set at critical value
	evoW = 1 - exp(-2*J);
	N = h*v*p;
	random_device rd;
	gen.seed(rd());
	rng.seed(rd()); //boost generator seeding
	uniform_int_distribution<int> bitd(0,1);
	tab = new int**[h];
	for(int i = 0; i < h; ++i){
		tab[i] = new int*[v];
		for(int j = 0; j < v; ++j){
			tab[i][j] = new int[p];
			for(int k = 0; k < p; ++k)
				tab[i][j][k] = bitd(gen)*2 - 1; //random spin affected at each lattice site
			}
		}
	calcMag();
	calcEnrg(); //calculating the average energy
	disB = new boost::bernoulli_distribution<> (evoW) ;
	linkDice = new boost::variate_generator<RNGType,boost::bernoulli_distribution<> > (rng,*disB);
	rng.seed(rd()); //reseeding to avoid correlations
	dis = new boost::uniform_real<double>(0,1) ;
	metroDice = new boost::variate_generator<RNGType, boost::uniform_real<double> > (rng,*dis);
}

Ising3d::Ising3d(const char * fname){
	ifstream ifile;
	ifile.open(fname);
	ifile >> h;
	ifile >> v;
	ifile >> p;
	ifile >> J;
	ifile.get();
	N = h*v*p;
	tab = new int**[h];
	char temp;
	for(int i = 0; i < h; i++){
		tab[i] = new int* [v];
		for(int j = 0; j < v; ++j)
			tab[i][j] = new int[p];
		}
	for(int i=0; i < N/8; i++){
		bitset<8> input_byte( ifile.get());
		for(int j=0; j<8; j++){
			int index = i * 8 + j;
			tab[index / (p*v)][ (index / p)%v][index % p] = input_byte[j] == 1 ? +1 : -1;
			}
		}
	ifile.close();
	evoW = 1- exp(-2*J);
	calcMag();
	calcEnrg();
	random_device rd; //initialization of the random number generators
	gen.seed(rd());
	rng.seed(rd()); //boost generator seeding
	disB = new boost::bernoulli_distribution<> (evoW) ;
	linkDice = new boost::variate_generator<RNGType,boost::bernoulli_distribution<> > (rng,*disB);
	rng.seed(rd()); //reseeding to avoid correlations
	dis = new boost::uniform_real<double>(0,1) ;
	metroDice = new boost::variate_generator<RNGType,boost::uniform_real<double> > (rng,*dis);
}

Ising3d::~Ising3d(){
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	delete disB;
	delete linkDice;
	delete dis;
	delete metroDice;
}

void Ising3d::calcMag(){
	mag = 0;
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k)
				mag += tab[i][j][k];
	mag /= N;
}

void Ising3d::calcEnrg(){
	enrg = 0;
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k){
				int B = 0;
				if(i < h-1) //counting 3 out fo the 6 neighbours, each in the "positive direction"
					B += tab[i+1][j][k];
				if(j < v-1)
					B += tab[i][j+1][k];
				if(k < p-1)
					B += tab[i][j][k+1];
				enrg += tab[i][j][k] * B;
				}
	enrg *= -J;
}

void Ising3d::coolDown(int vaccum){ //cools down the lattice to one of its vaccum
	assert(vaccum == +1 or vaccum == -1);
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k)
				tab[i][j][k] = vaccum;
}

int Ising3d::Bnei(const int i, const int j, const int k){
	int B = 0;
	if(i > 0)
		B += tab[i-1][j][k];
	if(i < h-1)
		B += tab[i+1][j][k];
	if(j > 0)
		B += tab[i][j-1][k];
	if(j < v-1)
		B += tab[i][j+1][k];
	if(k > 0)
		B += tab[i][j][k-1];
	if(k < p-1)
		B += tab[i][j][k+1];
	return B;
}


void Ising3d::dilation(const double lambda){
	double r = 1./lambda;
	uniform_real_distribution<> probd(0,1);
	double roundingParameter = probd(gen); //random value for the shifting parameter in the antecedent affectation

	int BORDER_H = (1. - r) * h/2;
	int BORDER_V = (1. - r) * v/2;
	int BORDER_P = (1. - r) * p/2;

	int *** temp = new int ** [h]; //new lattice storing the values of the spins post dilation
	for(int i = 0; i < h; ++i){
		temp[i] = new int * [v];
		for(int j = 0; j < v; ++j){
			temp[i][j] = new int[p];
			for(int k = 0; k < p; ++k)
				temp[i][j][k] = 0; //giving them a null value to not corrupt the assignation prescription on first neighours
			}
		}

	bool *** used = new bool ** [h];
	bool *** assigned = new  bool **[h];
	for(int i = 0; i < h; ++i){
		used[i] = new bool * [v];
		assigned[i] = new bool * [v];
		for(int j = 0; j < v; ++j){
			used[i][j] = new bool[p];
			assigned[i][j] = new bool[p];
			for(int k = 0; k < p; ++k){
				used[i][j][k] = false;
				assigned[i][j][k] = false;
				}
			}
		}


	//first assgining the boundary spins
{
	vector<site> BoundarySites;

	//collecting the sites on the border surfaces
	for(int j = 0; j < v; ++j){
		for(int k = 0; k < p; ++k){
			BoundarySites.push_back( site(0,j,k));
			BoundarySites.push_back( site(h-1,j,k));
			}
		}
	for(int i = 0; i < h; ++i){
		for(int k = 0; k < p; ++k){
			BoundarySites.push_back( site(i,0,k));
			BoundarySites.push_back( site(i,v-1,k));
			}
		}
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j){
			BoundarySites.push_back( site(i,j,0));
			BoundarySites.push_back( site(i,j,p-1));
			}
		}

	shuffle(BoundarySites.begin(),BoundarySites.end(),gen); //RANDOMLY PERMUTTING THEM
	for(auto site : BoundarySites){
		int i = get<0>(site);
		int j = get<1>(site);
		int k = get<2>(site);
		int I = BORDER_H + i * r + roundingParameter;
		int J = BORDER_V + j * r + roundingParameter;
		int K = BORDER_P + k * r + roundingParameter;
		if(not used[I][J][K]){
			temp[i][j][k] = tab[I][J][K];
			used[I][J][K] = true;
			assigned[i][j][k] = true;
			}
		}
}

	//assigning now the spins on the inside
{
	vector<site> InsideSites;

	//collecting them
	for(int i = 1; i < h-1; ++i)
		for(int j = 1; j < v-1; ++j)
			for(int k = 1; k < p-1; ++k)
				InsideSites.push_back( site(i,j,k));
	//shuffling them:
	shuffle(InsideSites.begin(), InsideSites.end(), gen);
	//assigning them:
	for(auto site : InsideSites){
		int i = get<0>(site);
		int j = get<1>(site);
		int k = get<2>(site);
		int I = BORDER_H + i * r + roundingParameter;
		int J = BORDER_V + j * r + roundingParameter;
		int K = BORDER_P + k * r + roundingParameter;
		if(not used[I][J][K]){
			temp[i][j][k] = tab[I][J][K];
			used[I][J][K] = true;
			assigned[i][j][k] = true;
			}
		}
}


//now filling in the holes
{
	vector<site> ToBeAssigned;
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k)
				if(not assigned[i][j][k])
					ToBeAssigned.push_back( site(i,j,k));

	vector< tuple<int, int, int, int>> AssignationData;
	for(auto toBeAssignedSpin : ToBeAssigned){
		int i = get<0>(toBeAssignedSpin);
		int j = get<1>(toBeAssignedSpin);
		int k = get<2>(toBeAssignedSpin);

		//calculating the neighbouring magnetization in the new lattice
		int Bnei = 0;
		if(i > 0)
			Bnei += temp[i-1][j][k];
		if(i < h-1)
			Bnei += temp[i+1][j][k];
		if(j > 0)
			Bnei += temp[i][j-1][k];
		if(j < v-1)
			Bnei += temp[i][j+1][k];
		if(k > 0)
			Bnei += temp[i][j][k-1];
		if(k < p-1)
			Bnei += temp[i][j][k+1];

		if( probd(gen) < 1. / (1 + exp(-2*J*Bnei)) )
			AssignationData.push_back( make_tuple(i,j,k,+1));
		else
			AssignationData.push_back( make_tuple(i,j,k,-1));
		}

	//now assigning the values	
	for(auto Dat : AssignationData){
		int i = get<0>(Dat);
		int j = get<1>(Dat);
		int k = get<2>(Dat);
		int val = get<3>(Dat);
		temp[i][j][k] = val;
		}
}


	//freeing temporary variables
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j){
			delete[] used[i][j];
			delete[] assigned[i][j];
			}
		delete[] used[i];
		delete[] assigned[i];
		}
	delete[] used;
	delete[] assigned;


	//updating spin configuration
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	tab = temp;
	calcMag();
	calcEnrg();
}


double Ising3d::s_distance(const int i1, const int j1, const int k1, const int i2, const int j2, const int k2){
	//euclidean distance between cubic lattice sites
	int deltaI = i1 - i2;
	int deltaJ = j1 - j2;
	int deltaK = k1 - k2;
	return sqrt( deltaI * deltaI + deltaJ * deltaJ + deltaK * deltaK);
}

void Ising3d::flipSWperio(const int Nsteps){

	for(int stepIndex = 0; stepIndex < Nsteps; ++stepIndex){

		bool *** flagged = new bool ** [h];
		for(int i = 0; i < h; ++i){
			flagged[i] = new bool * [v];
			for(int j = 0; j < v; ++j){
				flagged[i][j] = new bool[p];
				for(int k = 0; k < p; ++k)
					flagged[i][j][k] = false;
				}
			}

		{
		int X = 0, Y, Z;
		while(X < h){
			Y = 0;
			while(Y < v){
				Z = 0;
				while(Z < p){
					if( not flagged[X][Y][Z]){
						flagged[X][Y][Z] = true;
						int oldspin = tab[X][Y][Z];
						VirtualClusterSpins.push_back( site(X,Y,Z));
						int stack_n = 1;
						int parsed = 0;
						while(stack_n - parsed){
							int i = get<0>(VirtualClusterSpins[parsed]);
							int j = get<1>(VirtualClusterSpins[parsed]);
							int k = get<2>(VirtualClusterSpins[parsed]);
							//adding the candidate neighbours:
							if( not flagged[(i+1)%h][j][k] and tab[(i+1)%h][j][k] == oldspin and (*linkDice)() ){
								flagged[(i+1)%h][j][k] = true;
								VirtualClusterSpins.push_back( site((i+1)%h,j,k));
								stack_n++;
								}
							if( not flagged[(i-1+h)%h][j][k] and tab[(i-1+h)%h][j][k] == oldspin and (*linkDice)() ){
								flagged[(i-1+h)%h][j][k] = true;
								VirtualClusterSpins.push_back( site((i-1+h)%h,j,k));
								stack_n++;
								}
							if( not flagged[i][(j+1)%v][k] and tab[i][(j+1)%v][k] == oldspin and (*linkDice)() ){
								flagged[i][(j+1)%v][k] = true;
								VirtualClusterSpins.push_back( site(i,(j+1)%v,k));
								stack_n++;
								}
							if( not flagged[i][(j-1+v)%v][k] and tab[i][(j-1+v)%v][k] == oldspin and (*linkDice)() ){
								flagged[i][(j-1+v)%v][k] = true;
								VirtualClusterSpins.push_back( site(i,(j-1+v)%v,k));
								stack_n++;
								}
							if( not flagged[i][j][(k+1)%p] and tab[i][j][(k+1)%p] == oldspin and (*linkDice)() ){
								flagged[i][j][(k+1)%p] = true;
								VirtualClusterSpins.push_back( site(i,j,(k+1)%p));
								stack_n++;
								}
							if( not flagged[i][j][(k-1+p)%p] and tab[i][j][(k-1+p)%p] == oldspin and (*linkDice)() ){
								flagged[i][j][(k-1+p)%p] = true;
								VirtualClusterSpins.push_back( site(i,j,(k-1+p)%p));
								stack_n++;
								}
							parsed++;
							}
						//flipping virtual cluster with half probability now
						if((*metroDice)() < .5){
							for(auto site : VirtualClusterSpins){
								int i = get<0>(site);
								int j = get<1>(site);
								int k = get<2>(site);
								tab[i][j][k] *= -1;
								}
							}
						VirtualClusterSpins.clear();
						}
					Z++;
					}
				Y++;
				}
			X++;
			}
		}

		for(int i = 0; i < h; ++i){
			for(int j = 0; j < v; ++j)
				delete[] flagged[i][j];
			delete[] flagged[i];
			}
		delete[] flagged;
		}
}

void Ising3d::flipSW(const int Nsteps){ //fixed boundary lattice evolutions

	bool *** flagged = new bool**[h];
	for(int i = 0; i < h; ++i){
		flagged[i] = new bool * [v];
		for(int j = 0; j < v; ++j)
			flagged[i][j] = new bool[p];
		}
	for(int stepIndex = 0; stepIndex < Nsteps; ++stepIndex){
		for(int i = 0; i < h; ++i)
			for(int j = 0; j < v; ++j)
				for(int k = 0; k < p; ++k)
					flagged[i][j][k] = i == 0 or i == h-1 or j == 0 or j == v-1 or k == 0 or k == p-1;

		int X = 1, Y, Z;
		while(X < h-1){
			Y = 1;
			while(Y < v-1){
				Z = 1;
				while(Z < p-1){
					if(not flagged[X][Y][Z]){
						flagged[X][Y][Z] = true;
						int oldspin = tab[X][Y][Z];
						VirtualClusterSpins.push_back( site(X,Y,Z));
						int stack_n = 1;
						int parsed = 0;
						while(stack_n - parsed){
							int i = get<0>(VirtualClusterSpins[parsed]);
							int j = get<1>(VirtualClusterSpins[parsed]);
							int k = get<2>(VirtualClusterSpins[parsed]);
							//adding the candidate neighbours:
							if( not flagged[i+1][j][k] and tab[i+1][j][k] == oldspin and (*linkDice)() ){
								flagged[i+1][j][k] = true;
								VirtualClusterSpins.push_back( site(i+1,j,k));
								stack_n++;
								}
							if( not flagged[i-1][j][k] and tab[i-1][j][k] == oldspin and (*linkDice)() ){
								flagged[i-1][j][k] = true;
								VirtualClusterSpins.push_back( site(i-1,j,k));
								stack_n++;
								}
							if( not flagged[i][j+1][k] and tab[i][j+1][k] == oldspin and (*linkDice)() ){
								flagged[i][j+1][k] = true;
								VirtualClusterSpins.push_back( site(i,j+1,k));
								stack_n++;
								}
							if( not flagged[i][j-1][k] and tab[i][j-1][k] == oldspin and (*linkDice)() ){
								flagged[i][j-1][k] = true;
								VirtualClusterSpins.push_back( site(i,j-1,k));
								stack_n++;
								}
							if( not flagged[i][j][k+1] and tab[i][j][k+1] == oldspin and (*linkDice)() ){
								flagged[i][j][k+1] = true;
								VirtualClusterSpins.push_back( site(i,j,k+1));
								stack_n++;
								}
							if( not flagged[i][j][k-1] and tab[i][j][k-1] == oldspin and (*linkDice)() ){
								flagged[i][j][k-1] = true;
								VirtualClusterSpins.push_back( site(i,j,k-1));
								stack_n++;
								}
							parsed++;
							}
					double deltaE = 0;
					for( auto site : VirtualClusterSpins){
						int i = get<0>(site);
						int j = get<1>(site);
						int k = get<2>(site);
						if(i == 1)
							deltaE += oldspin * tab[i-1][j][k];
						if(i == h-2)
							deltaE += oldspin * tab[i+1][j][k];
						if(j == 1)
							deltaE += oldspin * tab[i][j-1][k];
						if(j == v-2)
							deltaE += oldspin * tab[i][j+1][k];
						if(k == 1)
							deltaE += oldspin * tab[i][j][k-1];
						if(k == p-2)
							deltaE += oldspin * tab[i][j][k+1];
						}
					deltaE *= 2*J;
					if( (*metroDice)() < 1./(1+exp(deltaE))){
						for(auto site : VirtualClusterSpins){
							int i = get<0>(site);
							int j = get<1>(site);
							int k = get<2>(site);
							tab[i][j][k] *= -1;
							}
						}
					VirtualClusterSpins.clear();
					}
				Z++;
				}
			Y++;
			}
		X++;
		}

	}

	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] flagged[i][j];
		delete[] flagged[i];
		}
	delete[] flagged;
}

void Ising3d::export_configuration(const char* fname){
	//exporting configuration with serialization
	ofstream ofile;
	ofile.open(fname);
	ofile.precision(9);
	ofile << h << " " << v << " " << p << " " << J << " ";
	for(int i = 0 ; i < N/8; i++){
		bitset<8> the_byte;
		for(int j = 0; j < 8; j++){
			int index = i*8 + j;
			the_byte[j] = (tab[index / (p*v)][(index / p) % v][index % p] == 1);
			}
		ofile << (char)the_byte.to_ulong();
		}
	ofile.close();
}

vector<site> * Ising3d::giveMeTheNeighbours(const int level){
	//returns list of half the neighbours at some specfific level, the half is chosen so that the first coordinate differential is positive
	auto theNeighbours = new vector<site>;
	for(int n1 = 0; n1 <= level; ++n1){
		for(int n2 = -level + n1 + 1; n2 <= level - n1; ++n2){
			theNeighbours->push_back( site(n1, n2, -(level-n1-abs(n2))));
			if( level != n1+abs(n2) )
				theNeighbours->push_back( site(n1, n2, level-n1-abs(n2)));
			}
		}
	return theNeighbours;
}

int Ising3d::Edensity(const int i, const int j, const int k){
	//value of the energy density
	return tab[i][j][k] * Bnei(i,j,k);
}

void Ising3d::cuttingInLattice( const int newSize){
	//cuts in a central subsection of the lattice
	int BORDER = (h - newSize)/2;
	int *** newTab = new int ** [newSize];
	for(int i = 0; i < newSize; ++i){
		newTab[i] =  new int * [newSize];
		for(int j = 0; j < newSize; ++j){
			newTab[i][j] = new int[newSize];
			for(int k = 0; k < newSize; ++k)
				newTab[i][j][k] = tab[BORDER + i][BORDER + j][BORDER +k];
			}
		}
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	tab = newTab;
	h = newSize;
	v = newSize;
	p = newSize;
	N = h*v*p;
	calcMag();
	calcEnrg();
}

void Ising3d::blowUp(){
	//dilation operation, multiplying by 2**3 the size of the cubic lattice
	h *= 2;
	v *= 2;
	p *= 2;
	N = h*v*p;
	int *** newTab = new int ** [h];
	for(int i = 0; i < h; ++i){
		newTab[i] =  new int * [v];
		for(int j = 0; j < v; ++j){
			newTab[i][j] = new int[p];
			for(int k = 0; k < p; ++k)
				newTab[i][j][k] = tab[i/2][j/2][k/2];
			}
		}
	for(int i = 0; i < h/2; ++i){
		for(int j = 0; j < v/2; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	tab = newTab;
	calcMag();
	calcEnrg();
}

void Ising3d::dilationWithPixels(const double lambda){
	//applies a dilation with duplicates/pixels prescription on the border
	double r = 1./lambda;
	int BORDER_H = (1. - r) * h/2;
	int BORDER_V = (1. - r) * v/2;
	int BORDER_P = (1. - r) * p/2;

	//new lattice storing the values of the spins post dilation
	int *** temp = new int ** [h];
	for(int i = 0; i < h; ++i){
		temp[i] = new int * [v];
		for(int j = 0; j < v; ++j){
			temp[i][j] = new int[p];
			for(int k = 0; k < p; ++k)
				temp[i][j][k] = 0;
			}
		}
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k){
				int I = BORDER_H + i * r;
				int J = BORDER_V + j * r;
				int K = BORDER_P + k * r;
				temp[i][j][k] = tab[I][J][K];
				}

	//freeing the temporary variables
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	tab = temp;
	calcMag();
	calcEnrg();
}

void Ising3d::dilationHybrid(const double lambda, const double pHole = 0.5){
	//applies a dilation with hybrid prescription on the border
	assert( pHole >= 0 and pHole <= 1);
	double r = 1./lambda;
	uniform_real_distribution<> probd(0,1);
	double roundingParameter = probd(gen);

	int BORDER_H = (1. - r) * h/2;
	int BORDER_V = (1. - r) * v/2;
	int BORDER_P = (1. - r) * p/2;

	int *** temp = new int ** [h]; //new lattice storing the values of the spins post dilation
	for(int i = 0; i < h; ++i){
		temp[i] = new int * [v];
		for(int j = 0; j < v; ++j){
			temp[i][j] = new int[p];
			for(int k = 0; k < p; ++k)
				temp[i][j][k] = 0; //giving them a null value to not corrupt the assignation prescription on first neighours
			}
		}

	bool *** used = new bool ** [h];
	bool *** assigned = new  bool **[h];
	for(int i = 0; i < h; ++i){
		used[i] = new bool * [v];
		assigned[i] = new bool * [v];
		for(int j = 0; j < v; ++j){
			used[i][j] = new bool[p];
			assigned[i][j] = new bool[p];
			for(int k = 0; k < p; ++k){
				used[i][j][k] = false;
				assigned[i][j][k] = false;
				}
			}
		}


	//first assigning the boundary spins
{
	vector<site> BoundarySites;

	//collecting the sites on the border surfaces
	for(int j = 0; j < v; ++j){
		for(int k = 0; k < p; ++k){
			BoundarySites.push_back( site(0,j,k));
			BoundarySites.push_back( site(h-1,j,k));
			}
		}
	for(int i = 0; i < h; ++i){
		for(int k = 0; k < p; ++k){
			BoundarySites.push_back( site(i,0,k));
			BoundarySites.push_back( site(i,v-1,k));
			}
		}
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j){
			BoundarySites.push_back( site(i,j,0));
			BoundarySites.push_back( site(i,j,p-1));
			}
		}

	shuffle(BoundarySites.begin(),BoundarySites.end(),gen); //RANDOMLY PERMUTTING THEM
	for(auto site : BoundarySites){
		int i = get<0>(site);
		int j = get<1>(site);
		int k = get<2>(site);
		int I = BORDER_H + i * r + roundingParameter;
		int J = BORDER_V + j * r + roundingParameter;
		int K = BORDER_P + k * r + roundingParameter;
		if(not used[I][J][K]){
			temp[i][j][k] = tab[I][J][K];
			used[I][J][K] = true;
			assigned[i][j][k] = true;
			}
		}
}

	//assigning now the spins on the inside
{
	vector<site> InsideSites;

	//collecting them
	for(int i = 1; i < h-1; ++i)
		for(int j = 1; j < v-1; ++j)
			for(int k = 1; k < p-1; ++k)
				InsideSites.push_back( site(i,j,k));
	//shuffling them:
	shuffle(InsideSites.begin(), InsideSites.end(), gen);
	//assigning them:
	for(auto site : InsideSites){
		int i = get<0>(site);
		int j = get<1>(site);
		int k = get<2>(site);
		const int I = BORDER_H + i * r + roundingParameter;
		const int J = BORDER_V + j * r + roundingParameter;
		const int K = BORDER_P + k * r + roundingParameter;
		if(not used[I][J][K]){
			temp[i][j][k] = tab[I][J][K];
			used[I][J][K] = true;
			assigned[i][j][k] = true;
			}
		}
}


	//now filling in the holes
{
	vector<site> ToBeAssigned;
	for(int i = 0; i < h; ++i)
		for(int j = 0; j < v; ++j)
			for(int k = 0; k < p; ++k)
				if(not assigned[i][j][k])
					ToBeAssigned.push_back( site(i,j,k));

	vector< tuple<int, int, int, int>> AssignationData;
	for(auto toBeAssignedSpin : ToBeAssigned){
		int i = get<0>(toBeAssignedSpin);
		int j = get<1>(toBeAssignedSpin);
		int k = get<2>(toBeAssignedSpin);

		if (probd(gen) < pHole){ //assigning by HeatBath
			//calculating the neighbouring magnetization in the new lattice
			int Bnei = 0;
			if(i > 0)
				Bnei += temp[i-1][j][k];
			if(i < h-1)
				Bnei += temp[i+1][j][k];
			if(j > 0)
				Bnei += temp[i][j-1][k];
			if(j < v-1)
				Bnei += temp[i][j+1][k];
			if(k > 0)
				Bnei += temp[i][j][k-1];
			if(k < p-1)
				Bnei += temp[i][j][k+1];

			if( probd(gen) < 1. / (1 + exp(-2*J*Bnei)) )
				AssignationData.push_back( make_tuple(i,j,k,+1));
			else
				AssignationData.push_back( make_tuple(i,j,k,-1));
			}
		else{ //creating a pixel
			const int I = BORDER_H + i * r + roundingParameter;
			const int J = BORDER_V + j * r + roundingParameter;
			const int K = BORDER_P + k * r + roundingParameter;
			const int value = tab[I][J][K];
			AssignationData.push_back( make_tuple( i, j, k, value));
			}
		}

	//now assigning the values	
	for(auto Dat : AssignationData){
		int i = get<0>(Dat);
		int j = get<1>(Dat);
		int k = get<2>(Dat);
		int val = get<3>(Dat);
		temp[i][j][k] = val;
		}
}


	//freeing temporary variables
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j){
			delete[] used[i][j];
			delete[] assigned[i][j];
			}
		delete[] used[i];
		delete[] assigned[i];
		}
	delete[] used;
	delete[] assigned;


	//updating Ising configuration
	for(int i = 0; i < h; ++i){
		for(int j = 0; j < v; ++j)
			delete[] tab[i][j];
		delete[] tab[i];
		}
	delete[] tab;
	tab = temp;
	calcMag();
	calcEnrg();
}

void Ising3d::thermalizeLattice( const vector<double> & dilationVector, const vector<int> & SWflipsVector){
	//MCMC mixing by applying a sequence of dilations and lattice updates
	assert ( dilationVector.size() == SWflipsVector.size()); //sanity checking that both vectors have the same size
	auto dilValue = dilationVector.begin();
	for( auto & swValue : SWflipsVector){
		dilationHybrid( *dilValue);
		flipSW( swValue);
		dilValue++;
		}
}

