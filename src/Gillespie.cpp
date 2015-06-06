//============================================================================
// Name        : Gillespie Algorithm for ion channel permeability
// Author      : Haroon Arshad
// Version     : 0.1
// Copyright   : GNU public license
//============================================================================

/*
 * Description:
 * This program uses the Gillespie algorithm for ion channel conformational changes
 * that impact the selectivity of ions
 */

#include <iostream>
#include <cmath>
#include <sys/time.h>
//#include <sys/resource.h>
#include "MersenneTwister.h"
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <numeric>


//#include <>g
using namespace std;

#define FILENAME  "output.txt"
#define SPECIFICATION_ "spec.txt"
FILE *foutput = NULL, *fspec = NULL, *foutput_sep[10];


double getcputime();
//template <size_t n, size_t m>
vector<vector<int > > readInitialReaction(int nr, int ns);
void print(vector<vector<int> > &out);
int minPositivePosition(vector<double> &vec);

/*
 * Initialize random number generator
 */
MTRand _generator;



int main() {

	double _start, _finish;
	//_start = getcputime();
	char tempo, fileend;

	/*
	 * Open file for output
	 */
	char name[200] = FILENAME, name_r[200] = SPECIFICATION_;
	foutput = fopen(name, "w");
	fspec = fopen(name_r, "w");
	//for(int it=0;it<10;++it) {
	//	tempo = (char)(((int)'0')+it);
	//	foutput_sep[it] = fopen(&tempo,"w");
	//}

	/*
	 * Initialize the storage and parameters
	 * nReactions 	- 	number of total reactions occuring in system
	 * nStage 		-	number of stage or state the system can be in (i.e open/closed, calcium binding numbers etc)
	 *
	 * tFinal 		-	time (ms) the system will continue till
	 */
	int nReactions = 12, nStage = 7, reactionPosition, totalIter = 0;//, diagnostic = 0;
	vector<int> output(nStage,0);
	vector<double> prop(nReactions,0), reactionRates(nReactions);
	double tFinal = 12, t=0, total_sum, tau;	// tFinal is in seconds
	double cai = 0.0000001;//;0.0000001;

	/*
	 * Voltage protocol
	 */
	double V = 10;
	typedef enum whatReaction
	{
		chan_cattach1=0,
		chan_cdetach1,
		chan_cattach2,
		chan_cdetach2,
		chan_cattach3,
		chan_cdetach3,
		chan_on1,
		chan_off1,
		chan_on2,
		chan_off2,
		chan_on3,
		chan_off3,
	} whatReaction;

	reactionRates[chan_cattach1] = 20000000*cai; 	// Kon1 (M^-1 s^-1)		:20
	reactionRates[chan_cdetach1] = 50; 	// Koff1 (s^-1)				:50
	reactionRates[chan_cattach2] = 20000000*cai; 	// Kon2 (M^-1 s^-1)		:20
	reactionRates[chan_cdetach2] = 50; 	// Koff2 (s^-1)				:50
	reactionRates[chan_cattach3] = 20000000*cai; 	// Kon3 (M^-1 s^-1)		:20
	reactionRates[chan_cdetach3] = 50;		// Koff2 (s^-1)			:50
	reactionRates[chan_on1] = 75; 			// alpha1 (s^-1)		:10
	reactionRates[chan_on2] = 150; 			// alpha2 (s^-1)		:30
	reactionRates[chan_on3] = 300; 			// alpha3 (s^-1)		:100
	reactionRates[chan_off1] = 10./(1+exp((V-75.0)/-50)); 		// beta1[Voltage]		:?
	reactionRates[chan_off2] = 75./(1+exp((V-120.0)/-50)); 		// beta2[Voltage]		:?
	reactionRates[chan_off3] = 100./(1+exp((V-120.0)/-50)); 		// beta3[Voltage]		:?

	/*
	 * Initial state of ion channel system
	 */
	output[0] = 5000;

	/*
	 * store reaction stage
	 */
	vector<vector<int> > vReaction = readInitialReaction (nStage,nReactions);


	while(t<tFinal)
	{
		fprintf(foutput,"%-5.7f %i %i %i %i %i %i %i ", t, output[0],output[1],output[2],output[3],output[4],output[5],output[6]);
		fprintf(foutput,"%i %i\n", output[0]+output[1]+output[2]+output[3],output[4]+output[5]+output[6]);

		// store each of the 'chemical' reactions that can occur in a propensity vector
		prop[chan_cattach1] = reactionRates[chan_cattach1]*output[0];
		prop[chan_cdetach1] = reactionRates[chan_cdetach1]*output[1];
		prop[chan_cattach2] = reactionRates[chan_cattach2]*output[1];
		prop[chan_cdetach2] = reactionRates[chan_cdetach2]*output[2];
		prop[chan_cattach3] = reactionRates[chan_cattach3]*output[2];
		prop[chan_cdetach3] = reactionRates[chan_cdetach3]*output[3];
		prop[chan_on1] = reactionRates[chan_on1]*output[1];
		prop[chan_off1] = reactionRates[chan_off1]*output[4];
		prop[chan_on2] = reactionRates[chan_on2]*output[2];
		prop[chan_off2] = reactionRates[chan_off2]*output[5];
		prop[chan_on3] = reactionRates[chan_on3]*output[3];
		prop[chan_off3] = reactionRates[chan_off3]*output[6];

		// the total amount from all reaction conditional states
		total_sum = accumulate(prop.begin(),prop.end(),0);


		// produce the cumulative increasing values as new storage format
		for(int i=1; i<nReactions; ++i)	{
			prop[i] = prop[i] + prop[i-1];
		}

		// store the increasing proportional summation values
		for(int i=0; i<nReactions; ++i)	{
			prop[i] = (prop[i]/total_sum);//-randStore;
		}

		// generate time interval of possible reaction taking place
		tau = (log(1/_generator.rand()))/total_sum;
		t = t + tau;

		reactionPosition = minPositivePosition(prop);

		for (int i=0; i<nStage; ++i)	{
			/**///if(vReaction[i][reactionPosition]!=0){
			output[i] = output[i] + vReaction[i][reactionPosition];
			/**///fprintf(foutput_sep[i],"%-5.7f %i\n",t,output[i]);}
			/**///else{}
		}
		//fprintf(foutput_sep[7],"%-5.7f %i\n",t,output[0]+output[1]+output[2]+output[3]);
		//fprintf(foutput_sep[8],"%-5.7f %i\n",t,output[4]+output[5]+output[6]);

		//printf("time - %g\n",t);
		totalIter++;

	}
	// Here is the original consecutive "large-file" print
	fprintf(foutput,"%-5.7f %i %i %i %i %i %i %i ", t, output[0],output[1],output[2],output[3],output[4],output[5],output[6]);
	fprintf(foutput,"%i %i\n", output[0]+output[1]+output[2]+output[3],output[4]+output[5]+output[6]);
	// The seperated file output is in the above loop


	fclose(foutput); fclose(fspec);
	/*for(int it=0;it<10;++it) {
		fclose(foutput_sep[it]);
	}
	*/
	//print(vReaction);

	//_finish = getcputime();
	//cout << "Time took for calculation is: " << _finish - _start << " seconds computing " << totalIter << " Gillespie simulation iterations.";

	return 0;
}




/*
 * Computation times for code analysis purposes
 */
/*
double getcputime(void)
{
	struct timeval tim;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	tim=ru.ru_utime;
	double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
	tim=ru.ru_stime;
	t+=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
	return t;
}
*/

/*
 * Read reaction vectors for state changes written in another file
 */
vector<vector<int> > readInitialReaction(int ns, int nr)
{
	vector<vector<int> > output(ns,vector<int>(nr));
	ifstream filein("reactionstates.dat");
	if(!filein)	{
		printf("ERROR: file not found. Make sure file is in the correct location");
		system("pause");
	}
	for (int row=0; row<ns; ++row)	{
		for (int col=0; col<nr; ++col)	{
			filein >> output[row][col];
		}
	}
	return output;
}


/*
 * Print 2D vector to screen (visual checking purpose)
 */
void print(vector<vector<int> > &out)
{
	for(size_t i=0; i<out.size(); ++i)	{
		for(size_t j=0; j<out[0].size(); ++j)	{
			printf(" %i",out[i][j]);
		}
		printf("\n");
	}
}


int minPositivePosition(vector<double> &vec)
{
	double randNumber = _generator.rand();
    size_t i = 0;
    while(i<vec.size()-1)
    {
    	if((i>0)&&(vec[i-1]!=vec[i])&&(randNumber < vec[i]))
    	{
    		return i;
    		break;
    	}
    	else if(vec[i]==1)
    	{
    		return i;
    		break;
    	}
    	else
    	{
    		if(randNumber < vec[0])	{
    			return i;
    			break;
    		}
    	}
      ++i;
    }
    return vec.size()-1;
}
