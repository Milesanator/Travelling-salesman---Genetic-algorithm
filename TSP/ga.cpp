/************************************************************
*      ga.cpp - GA Program for CSCI964 - Ass2
*      Written by: Koren Ward May 2010
*      Modified by: Miles Tuffs	2016
*      Changes: <Provide details of any changes here>
*************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <vector>
using namespace std;
class EdgeMatrix {
	public:
	int* neighbours;
	int size;
	
	EdgeMatrix() {
		size=0;
		neighbours=new int[4];
	}
	~EdgeMatrix() {
		delete[] neighbours;
	}
	EdgeMatrix(const EdgeMatrix &obj) {
		size = obj.size;
		for(int i=0;i<obj.size;i++){
			neighbours[i] = obj.neighbours[i];
		}
	}
	EdgeMatrix& operator=(EdgeMatrix obj) {
		size = obj.size;
		for(int i=0;i<obj.size;i++){
			neighbours[i] = obj.neighbours[i];
		}
		return *this;
	}
	void add(int x){
		neighbours[size]=x;
		if(size<=4)
			size++;
	}
	void remove(int x) {
		bool removed=false;
		for(int i=0;i<size;i++) {
			if(removed) {
				neighbours[i-1]=neighbours[i];
			}
			if(neighbours[i]==x) {
				neighbours[i]= -1;// dead
				removed=true;
			}
		}
		if(removed) {
			size--;
		}
	}
	bool ifExists(int x){
		for(int i=0;i<size;i++){
			if(neighbours[i]==x){
				return true;
			}
		}
		return false;
	}
};

char	tournamentType[20] = "roulette";
const double cCrossoverRate = 0.75;
const double cMutationRate = 0.001;
const int    cNumGens = 250;
const int    cPopSize = 100; // must be an even number
const int    cTournamentSize = 5;
const int    Seed = 1234;
const double    cTargetFitness=80;
// new vars
int numTowns;
int* townX;
int* townY;
int* townType;
double tollMatrix[4][4] = {
{-1,-1,-1,-1},
{-1,10,7.5,5},
{-1,7.5,5,2.5},
{-1,5,2.5,1},
};
double** costMatrix;

void InitPop(int ***CrntPop,int ***NextPop,double **Fitness,int **BestMenber);
void FreeMem(int **CrntPop,int **NextPop,double *Fitness,int *BestMember);
int Tournament(double *Fitness);
double EvaluateFitness(int *Member);
void Crossover(int *P1,int *P2,int *C1,int *C2);
void Copy(int *P1,int *P2,int *C1,int *C2);
void Mutate(int *Member);
double Rand01();    // 0..1
int RandInt(int n); // 0..n-1
double euclideanDistance(int a, int b); // straight line distance
double calcTolls(int x, int y); // calulate tolls based on road type
void calculateCostMatrix();
void edgeRecombinationCrossover(int* P1, int* P2, int* C);

//===========================================================

int main(int argc,char *argv[]){
  ifstream fin;
  char fileName[25];
  int **CrntPop, **NextPop; // the crnt & next population lives here
  double *Fitness, BestFitness=numeric_limits<double>::max();
  int *BestMember; // fitness vars
  int i;
  bool TargetReached=false;
  // Read data file
  if(argc!=2) {
	cout<<"Enter data filename: \n";
	cin>>fileName; cin.ignore();
  } else {
	strcpy(fileName,argv[1]);
  }
  fin.open(fileName);
  if(!fin.good()){cout<<"File not found\n";exit(1);}
  fin>>numTowns;
  townType=new int[numTowns];townX=new int[numTowns];townY=new int[numTowns];
  for(i=0;i<numTowns;i++){
	fin >> townX[i];
	fin >> townY[i];
	fin >> townType[i];
  }
  fin.close();
	
  InitPop(&CrntPop,&NextPop,&Fitness,&BestMember);
  calculateCostMatrix();
  for(int Gen=0;Gen<cNumGens;Gen++){
    for(i=0;i<cPopSize;i++){

      // Evaluate the fitness of pop members
      Fitness[i]=EvaluateFitness(CrntPop[i]);
      if(BestFitness>Fitness[i]){ // save best member
        BestFitness=Fitness[i];
		// Arr copy crntpop to the BestMember
        for(int j=0;j<numTowns;j++){
			BestMember[j]=CrntPop[i][j];
		}
        if(Fitness[i]<=cTargetFitness){
		  TargetReached=true;
          break;
        }
      }
    }
    if(TargetReached)break;

    // Produce the next population
    for(i=0;i<cPopSize;i+=2){
      int Parent1=Tournament(Fitness);
      int Parent2=Tournament(Fitness);
      if(cCrossoverRate>Rand01()){
		Crossover(CrntPop[Parent1],CrntPop[Parent2],NextPop[i],NextPop[i+1]);
      }else{
        Copy(CrntPop[Parent1],CrntPop[Parent2],NextPop[i],NextPop[i+1]);
	  }
      if(cMutationRate<Rand01()){
	    Mutate(NextPop[i]);
	  }
      if(cMutationRate<Rand01()){
	    Mutate(NextPop[i+1]);
	  }
    }
    int **Tmp=CrntPop; CrntPop=NextPop; NextPop=Tmp;
	//cout << Gen << endl;
    cout<<BestFitness<<endl;
  }
  if(TargetReached) cout<<"Target fitness reached: "<<BestFitness<<"!\n";
  else cout<<"Target fitness not reached: "<<BestFitness<<"!\n";
  FreeMem(CrntPop,NextPop,Fitness,BestMember);
  char s[20];cin.getline(s,20);
  return 0;
}
//===========================================================

void InitPop(int ***CrntPop,int ***NextPop,double **Fitness,int **BestMember){
  int i, j;
  srand(Seed);
  *CrntPop = new int*[cPopSize];
  *NextPop = new int*[cPopSize];
  for(i=0;i<cPopSize;i++){
    (*CrntPop)[i] = new int[numTowns];
    (*NextPop)[i] = new int[numTowns];
  }
  *Fitness    = new double[cPopSize];
  *BestMember = new int[numTowns];
  if(Fitness==NULL||BestMember==NULL)exit(1);
  for(i=0;i<cPopSize;i++){
    for(j=0;j<numTowns;j++){
      (*CrntPop)[i][j] = j;
    }
	random_shuffle(&(*CrntPop)[i][0],&(*CrntPop)[i][numTowns]);
  }
}

void FreeMem(int **CrntPop,int **NextPop,double *Fitness,int *BestMenber){
  for(int i=0;i<cPopSize;i++){
    delete[]CrntPop[i];
    delete[]NextPop[i];
  }
  for(int i=0;i<numTowns;i++){
	delete[]costMatrix[i];
  }
  delete costMatrix;
  delete CrntPop;
  delete NextPop;
  delete Fitness;
  delete BestMenber;
  delete[] townY;
  delete[] townX;
  delete[] townType;
}
//===========================================================

double EvaluateFitness(int *Member){
//Evaluates fitness based on bit pattern
  double dist=0;
  for(int i=1;i<numTowns;i++){
	dist+=costMatrix[Member[i-1]][Member[i]];
  }

  // complete the tour
  dist+=costMatrix[Member[numTowns-1]][Member[0]];
  return dist;
}
//================================================================

int Tournament(double *Fitness){
  if(strcmp(tournamentType,"roulette")==0) {
	int RandNum = RandInt(4000);
	double* cpyFitness = new double[cPopSize];
	for(int i=0;i<cPopSize;i++){
		cpyFitness[i]=Fitness[i];
	}
	double min=numeric_limits<double>::max(),maxNormalised=0;
	for(int i=0;i<cPopSize;i++){
		if(cpyFitness[i]<min){
			min = cpyFitness[i];
		}
	}
	for(int i=0;i<cPopSize;i++) {
		cpyFitness[i] -= min;
	}
	for(int i=0;i<cPopSize;i++){
		if(cpyFitness[i]>maxNormalised) {
			maxNormalised=cpyFitness[i];
		}
	}
	double tempSum=0;
	for(int i=0;i<cPopSize;i++){
		cpyFitness[i]=maxNormalised-cpyFitness[i];
		tempSum+=cpyFitness[i];
	
		if(tempSum>RandNum){
			return i;
		}
	}
	return 0;
  } else {
    int WinFit = -99999, Winner=0;
    for(int i=0;i<cTournamentSize;i++){
      int j = RandInt(cPopSize);
      if(Fitness[j]>WinFit){
		WinFit = Fitness[j];
		Winner = j;
	  }
	}
	return Winner;
  }
}

void Crossover(int *P1,int *P2,int *C1,int *C2){
	edgeRecombinationCrossover(P1, P2, C1);
	edgeRecombinationCrossover(P2, P1, C2);
}

void edgeRecombinationCrossover(int* P1, int* P2, int* C){
//Edge Recombination Crossover - Genetic Algorithms
//This is a crossover techniques for permutation (ordered) chromosomes.
// It strives to introduce the fewest paths possible.
// The idea here is to use as many existing edges, or node-connections, as possible to generate children.
	int i,j, X,Z, left, right, irand;
	int currentNeighbours[4], currentNumNeighbours[4];
	int lowest[4], numLowest;
	lowest[0]=0;lowest[1]=0;lowest[2]=0;lowest[3]=0;

	vector<int> remains;
	for(i=0;i<numTowns;i++) {
		remains.push_back(P1[i]);
	}
//	// Create new EdgeMatrix list
	vector<EdgeMatrix> neighbourList(numTowns);
	for(i=0;i<numTowns;i++) {
		neighbourList[i].size=0;
	}
	// GENERATE NEIGHBOUR LIST
	for(i=0;i<numTowns;i++) {
		left = (i+1)%numTowns;
		right = (i-1)%numTowns;
		if(right==-1 || left == numTowns){
			if(right==-1){
				right=numTowns-1;
			}else {
				left=0;
			}
		}
		neighbourList[P1[i]].add(P1[left]);
		neighbourList[P1[i]].add(P1[right]);
		if(!(neighbourList[P2[i]].ifExists(P2[left]))) {
			neighbourList[P2[i]].add(P2[left]);
		}
		if(!(neighbourList[P2[i]].ifExists(P2[right]))) {
			neighbourList[P2[i]].add(P2[right]);
		}
	}
	// Do edge recombination
	X=P1[0]; // take first inded from Parent 1, as initial start point
	for(i=0;i<numTowns;i++){
		C[i]=X;
		for(j=0;j<numTowns;j++){
			neighbourList[j].remove(X);
			if(remains.size()>j && remains[j]==X){
				remains.erase(remains.begin()+j);
			}
		}
		
		if(neighbourList[X].size==0) {
			irand=RandInt(remains.size());
			Z = remains[irand];
		} else {
			numLowest=0;
			for(int y=0;y<neighbourList[X].size;y++) {
				currentNeighbours[y] = neighbourList[X].neighbours[y];
				currentNumNeighbours[y] = neighbourList[currentNeighbours[y]].size;
			}
			lowest[0] = numTowns+1;
			for(j=0;j<neighbourList[X].size;j++) {
				if(currentNumNeighbours[j]<=lowest[numLowest]) {
					lowest[numLowest]= currentNeighbours[j];
					numLowest++;
					// set this only for comparison
					lowest[numLowest]= currentNeighbours[j];
				}
			}
			if(numLowest!=1) {
				// if more than 1 neighbour have lowest children, choose random
				irand = RandInt(numLowest);
				Z=lowest[irand];
			} else {
				// neighbour with least children
				Z=lowest[0];
			}
		}
		X=Z;
	}
	remains.clear();
	neighbourList.clear();
}

void Mutate(int *Member){
  int Pick1 = RandInt(numTowns);
  int Pick2 = RandInt(numTowns);
  int temp;
  temp=Member[Pick1];
  Member[Pick1]=Member[Pick2];
  Member[Pick2]=temp;
}

void Copy(int *P1,int *P2,int *C1,int *C2){
  for(int i=0;i<numTowns;i++){
    C1[i]=P1[i]; C2[i]=P2[i];
  }
}
//=================================================================

double Rand01(){ // 0..1
  return(rand()/(double)(RAND_MAX));
}

int RandInt(int n){ // 0..n-1
  return int( rand()/(double(RAND_MAX)+1) * n );
}

void calculateCostMatrix(){
	// STEP 2 - PRE-CALCULATE COST MATRIX
	costMatrix= new double*[numTowns];
	for(int i=0;i<numTowns;i++){
		costMatrix[i] = new double[numTowns];
	}
	double dist=0;
	for(int i=0;i<numTowns;i++){
		for(int j=0;j<numTowns;j++){
			if(i!=j){
				dist=euclideanDistance(i, j)*calcTolls(i, j);
				costMatrix[i][j] = dist;
			} else {
				costMatrix[i][j]=-1;
			}
		}
	}
}

double calcTolls(int x, int y) {
	return tollMatrix[townType[x]][townType[y]];
}
// calculate straight line distance between 2 indexes
double euclideanDistance(int a, int b) {
	return sqrt(pow((townX[a]-townX[b]), 2.0) + pow((townY[a]-townY[b]), 2.0));
}


