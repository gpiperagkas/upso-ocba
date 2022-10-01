/*
UPSO implementation with OCBA for simulation of stochasticsystems. 
*/

/*	Copyright (C) 2022  Grigorios Piperagkas

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/



#include<iostream>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<ctime>
#include<iomanip>
#include<cstdio>

using namespace std;



//Algorithm parameter structs 


struct problem{
	int dim ;
	int bench; 

};

struct algoparam{
	int swarmSize;
	int iterations;
	int k;
	double c1 ; 
	double c2 ;
	double xi ;
	double unifactor ;
	int radius;

};

struct simulator{
	int T;
	int n0;
	int k;
	int delta;


};

//********************************

#include "subroutines.cpp"
#include "writefiles.cpp"
//#include "testfunctions.cpp"

int main (int argc, char *argv[])
{
	int i,tests; 
	int j,it,l,m,maxfevals;
	int *fevals; 

	
	double xmini[]={-100,-30,-5.12,-600,-20,-8,-10,-10,-4,-10,-10,-100,-65.536,-500};
	double xmaxi[]={100,30,5.12,600,30,8,10,10,5,10,10,100,65.536,500};

	double finalvalue[100];
	double fe[100];
	double mean;
	double std;
	double goal;	
	double unifactor;	
	double suc,fev;
	
	FILE * fout;
	char foutname[30];
	char fouttmp[20];
	
	algoparam algo;
	problem prob; 
	simulator ocba;

	prob.bench=atoi(argv[1]); //for simulation set value 100
	prob.dim=atoi(argv[2]); // set dimensions of simulator input vector
	unifactor=atoi(argv[3]);  //1; // set global pso variant as default
	//ocba.T=30;


	tests= 100; //set tests equal to one for one optimization procedure	
	algo.swarmSize= 5* prob.dim;
	maxfevals= 1000* prob.dim;//does not play any role in simulation	


	double temp3[prob.dim];
	double temp4;

	double *Xmin=NULL;
	Xmin = new double [prob.dim];

	double *Xmax=NULL;
	Xmax = new double [prob.dim];

	for (i=0;i<prob.dim;i++){ //set bounds correctly from simulator data
			Xmin[i]=xmini[prob.bench];
			Xmax[i]=xmaxi[prob.bench];
	}
	
	algo.c1=2.05;
	algo.c2=2.05;
	algo.xi=0.729;
	algo.k=2;
	algo.radius=5;
	
	double **particles=NULL;
	particles = new double *[algo.swarmSize];
	for (i=0;i<algo.swarmSize;i++){
		particles[i]= new double[prob.dim];
	}

	double **pbest=NULL;
	pbest = new double *[algo.swarmSize];
	for (i=0;i<algo.swarmSize;i++){
		pbest[i] = new double[prob.dim];
	}
	
	double **lbest=NULL;
	lbest = new double *[algo.swarmSize];
	for (i=0;i<algo.swarmSize;i++){
		lbest[i] = new double[prob.dim];
	}

	double *gbest=NULL ;
	gbest = new double[prob.dim];

	double **velocity=NULL;
	velocity = new double *[algo.swarmSize];
	for (i=0; i<algo.swarmSize;i++){
		velocity[i]=new double[prob.dim];
	}
	
    	double **lvel=NULL;
    	lvel = new double *[algo.swarmSize];
	for (i=0; i<algo.swarmSize;i++){
	    lvel[i]=new double[prob.dim];
	}

	double **gvel=NULL;
	gvel = new double *[algo.swarmSize];
	for (i=0; i<algo.swarmSize;i++){
		gvel[i]=new double[prob.dim];
	}

	double *obj_particles=NULL;
	obj_particles = new double[algo.swarmSize];

	double *obj_pbest=NULL;
	obj_pbest = new double[algo.swarmSize];

	double *obj_lbest=NULL;
	obj_lbest = new double[algo.swarmSize];

	double *obj_gbest=NULL;
	obj_gbest = new double[1];

	double *uf=NULL;
	uf = new double[algo.swarmSize];

	srand48(time(NULL));
	fevals=(int *)malloc(sizeof(int));

	mean = 0 ;

	if (unifactor<=1.0)
		for (i=0;i<algo.swarmSize;i++)
			uf[i]=unifactor;


	for (l=0;l<tests;l++){ 				// run PSO $tests times in order to find mean and std 
		cout << "********   Start of " << l << " test of PSO \n";
		*fevals = 0;
		initialization(particles, pbest, lbest, velocity, Xmin, Xmax, prob, algo);
		evaluateBests(pbest, lbest, gbest, obj_pbest, obj_lbest, obj_gbest, fevals, prob, algo);
		while (*fevals < maxfevals ){
			flyGlobal(particles, pbest, gbest, gvel, velocity, Xmin, Xmax, prob, algo);//global variant
			flyLocal(particles, pbest, lbest, lvel, velocity, prob, algo);//local variant
			evalposition(particles, lvel, gvel , velocity, uf, Xmin, Xmax, prob, algo);//position update
			positions_feval(particles, obj_particles, fevals, prob, algo, ocba); // Objective functions evaluations
			evaluatePbest(particles, pbest, obj_particles, obj_pbest, prob, algo);//personal best evaluation
			evaluateLocalbest(particles, lbest, obj_particles, obj_lbest, prob, algo);//neighbourhood best evaluation
			evaluateGlobalbest(pbest, gbest, obj_pbest, obj_gbest, prob,algo);//global guide evaluation
			//cout << " obj_gbest :  "<< obj_gbest[0] << "\n";

		}
		fe[l]=*fevals;
		finalvalue[l]=obj_gbest[0];
		cout << "Final value of "<< l <<"th test is " << finalvalue[l] <<"\n";
		mean = mean + finalvalue[l];
	cout << "TEST finalvalue tou test "<< l<<": " << finalvalue[l] <<"\n";
	}
//	m=0;
//	suc=0;
//	for (i=0;i<100;i++){
//		if (fe[i]<1000000){
//			m++;
//			fev=fev+fe[i];	
//		}
//		else{
//			suc++;
//		}
//	}
//	fev=fev/m; //mean value of function evaluations
//	suc=(tests - suc)/tests; //success ratio

	mean = mean/tests;
	cout << "\n" ;
	cout << "Mean value of tests is "<< mean ;
	temp4 =0; 
	for  (m=0;m<tests;m++){
			temp4=temp4+ ((finalvalue[m]-mean)*(finalvalue[m]-mean));
	}
	temp4=temp4/tests;
	std=sqrt(temp4);

	cout << "\n";
	cout << "Std value of tests is "<< std <<"\n";

	for (i=0;i<tests;i++){
		cout << "final value of "<<i<<"th test is :"<< finalvalue[i]<<"in functions evaluations : "<< fe[i]<<"\n";
	}
	cout << "Function evaluations : "<< *fevals << "\n";
	
	strcpy(foutname,"Results__TestProb_");
	sprintf(fouttmp,"%d",prob.bench);
	strcat(foutname,fouttmp);

	fout = fopen(foutname,"a");
	outfile(fout,mean,std,uf,unifactor,prob,algo);

	fclose(fout);


}


