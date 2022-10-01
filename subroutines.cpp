/*************************/
/* MAIN UPSO SUBROUTINES */
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

#include "testfunctions.cpp"

//****************************************************************************
// OPTIMAL COMPUTING BUDGET ALLOCATION PROCEDURE 
//---------------------------------------------------------------------
// BASED ON: C. CHEN, L. SHI, LH. LEE, ´Stochastic systems simulation 
// optimization´, Front. Electr. Electron. Eng. China 2011, 6(3): 468-480 
//****************************************************************************

void ocbasim(double **particles, double *obj_particles, problem prob, algoparam algo, simulator ocba){
	int i,j,t,l,ind;
	double temp[prob.dim];
	int rep[algo.swarmSize];//current replications
	int repb[algo.swarmSize];//previous allocation replications
	double tobj[algo.swarmSize][ocba.T];//replication values for each particle	
	double mean[algo.swarmSize];//mean value of replications
	double std[algo.swarmSize];//standard deviation of replications
	double tmp1,tmp2;
	int sysb; //index of the best individual
	double min,sumrep;
	double rat[algo.swarmSize]; //ratio value for each particle	
	

	//initialize temporary archive of replications
	for (i=0;i<algo.swarmSize;i++)
		for (j=0;j<ocba.T;j++)
			tobj[i][j]=0;
	
	//PERFORM INITIAL REPLICATIONS
	l=0;
	for (i=0;i<algo.swarmSize;i++){
		tmp1=0;
		tmp2=0;
		rep[i]=ocba.n0;
		repb[i]=rep[i]; 
		for (j=0;j<prob.dim;j++)
			temp[j]=particles[i][j];
		for (t=0;t<ocba.n0;t++){
			tobj[i][t]= Objective(temp,prob);//replication for each particle
			tmp1=tmp1+tobj[i][t];		
		}
		mean[i]=tmp1/ocba.n0;// compute mean and standard deviation of initial replications
		for (t=0;t<ocba.n0;t++){
			tmp2=tmp2+pow(tobj[i][t]-mean[i],2);
		}
		std[i]=sqrt(tmp2/ocba.n0);
	}
	
	//MAIN LOOP FOR EXTRA REPLICATIONS AND BUDGET ALLOCATION
	sumrep = ocba.n0*algo.swarmSize;
	while ( sumrep < ocba.T ){


		for (i=0;i<algo.swarmSize;i++)
			repb[i]=rep[i];	

		//find minimum solution untill now
		min=mean[0];
		sysb=0;
		for (i=0;i<algo.swarmSize;i++){
			if (mean[i]<min){
				min=mean[i];
				sysb=i;
			}
		}
		
		//INCREASE COMPUTING BUDGET BY DELTA AND ALLOCATE ACCORDING TO RATIOS
		if (sysb!=0)
			ind=0;
		else
			ind=1;
		for (i=0;i<algo.swarmSize;i++){
			if (i!=sysb){
				rat[i]=((std[i]*pow(mean[sysb]-mean[ind],2))/(std[ind]*pow(mean[sysb]-mean[i],2)));	
			}else{
				rat[i]=0;
			}
		}
		tmp1=0;
		
		for (i=0;i<algo.swarmSize;i++){
			tmp1=tmp1+pow(rat[i]/std[i],2);	
		}
		
		rat[sysb]=std[sysb]*sqrt(tmp1);
		tmp2=0;
		for (i=0;i<algo.swarmSize;i++)
				tmp2=tmp2+rat[i];
		rep[ind]=(sumrep+ocba.delta)/tmp2;//integer division

		for (i=1;i<algo.swarmSize;i++){
			rep[i]=rat[i]*rep[ind];
			if (rep[i]-repb[i]>0){
				for (j=0;j<prob.dim;j++)
					temp[j]=particles[i][j];
				for (t=0;t<(rep[i]-repb[i]);t++)
					tobj[i][repb[i]+t]=Objective(temp,prob);
				
			}	
		}
		for (i=0;i<algo.swarmSize;i++){
			tmp1=0;
			tmp2=0;
			for (t=0;t<rep[i];t++)
				tmp1=tmp1+tobj[i][t];
			mean[i]=tmp1/rep[i];
			for (t=0;t<rep[i];t++)
                        	tmp2=tmp2+pow(tobj[i][t]-mean[i],2);
                	std[i]=sqrt(tmp2/rep[i]);
		}
		sumrep=sumrep+ocba.delta;
	}	
}






//*******************************
// INITIALIZATION   *************
//*******************************

void initialization (double **particles,double **pbest, double **lbest, double **velocity,double *Xmin, double *Xmax, problem prob,algoparam algo ){
	int i,j;
	for (i=0;i<algo.swarmSize; i++){
		for (j=0;j<prob.dim;j++){
			particles[i][j] = Xmin[j] + ( Xmax[j] - Xmin[j] )*drand48();
			pbest[i][j]=particles[i][j];
			lbest[i][j]=pbest[i][j];
			velocity[i][j] = drand48();
		}
	}
}

//*******************************
// EVALUATE INITIAL BESTs  ******
//*******************************

void evaluateBests (double **pbest, double **lbest, double *gbest, double *obj_pbest, double *obj_lbest,double *obj_gbest, int *fevals, problem prob, algoparam algo){
	int i,j,k;
	double temp1[prob.dim];
    double temp2[prob.dim];	
	
	for (i=0;i<prob.dim;i++){
		gbest[i]=pbest[0][i];
		temp1[i]=gbest[i];
	}
	obj_pbest[0]=Objective(temp1,prob);
	obj_gbest[0]=obj_pbest[0];
	obj_lbest[0]=obj_gbest[0];
	(*fevals)++;
	for (i=1;i<algo.swarmSize;i++){
		for (j=0;j<prob.dim;j++){
			temp1[j] = pbest[i][j];
			temp2[j] = gbest[j];
		}
		obj_pbest[i]= Objective (temp1, prob);
		(*fevals)++;
		obj_lbest[i]= obj_pbest[i];
		if (obj_pbest[i]< obj_gbest[0]){
			for (k=0;k<prob.dim;k++)
				gbest[k]=pbest[i][k];
			obj_gbest[0] = obj_pbest[i];
		}
	}

}


//*******************************
// FLY GLOBAL *******************
//*******************************

void flyGlobal(double **particles, double **pbest, double *gbest, double **gvel, double **velocity, double *Xmin, double *Xmax, problem prob, algoparam algo){

	int i, j;
	double r1,r2;
	double vmax[prob.dim];

	for (i=0;i<prob.dim;i++){
		vmax[i]=(Xmax[i]-Xmin[i])/algo.k ;
	}
	for (i=0;i<algo.swarmSize; i++){
		for (j=0;j<prob.dim;j++){
			r1=drand48();
			r2=drand48();
			gvel[i][j] = algo.xi* (velocity[i][j] + algo.c1*r1*(pbest[i][j]-particles[i][j]) + algo.c2*r2*(gbest[j]-particles[i][j]));
			if (gvel[i][j]> vmax[j] ){
				gvel[i][j]=vmax[j];
			}
			if (gvel[i][j]< -vmax[j]){
				gvel[i][j]= -vmax[j];
			}
		}
	}



}


//*********************************
// FLY LOCAL       ****************
//*********************************

void flyLocal(double **particles, double **pbest, double **lbest, double **lvel, double **velocity, problem prob, algoparam algo){
	int i,j,k;
	double r1,r2; 
	double temp1[prob.dim];
	double temp2[prob.dim];
	
	for (i=0;i<algo.swarmSize;i++){
		for (j=0;j<prob.dim;j++){
			r1=drand48();
			r2=drand48();
			lvel[i][j] = algo.xi* (velocity[i][j] + algo.c1*r1*(pbest[i][j] - particles[i][j]) + algo.c2*r2*(lbest[i][j])-particles[i][j]);
		
		}
	}

}


//**********************************************************
// EVALUATE UNIFICATION FACTOR - CHAOTIC MAP (LOGISTIC) ****
//**********************************************************

void evalUnifactor(double *uf, algoparam algo){
	int m=4;
	int i;

	for (i=0;i<algo.swarmSize;i++)
		uf[i]=m*uf[i]*(1-uf[i]);//logistic map

}


//****************************************************
//EVALUATE VELOCITY AND POSITION *********************
//****************************************************

void evalposition(double **particles, double **lvel, double **gvel, double **velocity, double *uf, double *Xmin, double *Xmax, problem prob, algoparam algo){
	int i,j; 
	double temp[prob.dim];
	double obj1;
	for (i=0;i<algo.swarmSize;i++){
		for (j=0;j<prob.dim;j++){

			velocity[i][j]= uf[i] * gvel[i][j] + (1-uf[i])*lvel[i][j];
			particles[i][j] = particles[i][j] + velocity[i][j];	

			if (particles[i][j]>Xmax[j]){
				particles[i][j]=Xmax[j];
			}
			if (particles[i][j]<Xmin[j]){
				particles[i][j]=Xmin[j];
			}
		}
	}


}

//******************************************************
// FUNCTION EVALUATIONS  *******************************
//******************************************************

void positions_feval(double **particles, double *obj_particles,int *fevals, problem prob, algoparam algo, simulator ocba){
	int i,j;
	double temp[prob.dim];

	if (prob.bench==100){
		ocbasim(particles, obj_particles, prob, algo, ocba);
		
	}else{
		for (i=0;i<algo.swarmSize;i++){
			for (j=0;j<prob.dim;j++){
				temp[j]=particles[i][j];
			}
			obj_particles[i]= Objective(temp,prob);//Position objective evaluation
			(*fevals)++;//function evaluations counter
		}
	}

}


//******************************************************
// EVALUATE PERSONAL BESTS *****************************
//******************************************************

void evaluatePbest(double **particles, double **pbest, double *obj_particles, double *obj_pbest, problem prob, algoparam algo){
    int i,k;
	double temp1[prob.dim];

	for (i=0;i<algo.swarmSize;i++){
		if (obj_particles[i]<obj_pbest[i]){ // if new pbest found, update objective_pbest
			for (k=0;k<prob.dim;k++){
				pbest[i][k]=particles[i][k];
			}
			obj_pbest[i]=obj_particles[i];
		}
	}
}


//***************************************************************
// EVALUATE GLOBAL BEST *****************************************
//***************************************************************

void evaluateGlobalbest(double **pbest, double *gbest, double *obj_pbest, double *obj_gbest, problem prob, algoparam algo) {
	int i,j;


	for (i=0;i<algo.swarmSize;i++){//compare gbest with the personal bests
		if (obj_pbest[i]< obj_gbest[0]){
			for (j=0;j<prob.dim;j++)
				gbest[j]=pbest[i][j];
			obj_gbest[0]=obj_pbest[i];
		}
	}

}


//***********************************************************************
// EVALUATE LOCAL BEST **************************************************
//***********************************************************************

void evaluateLocalbest (double **particles, double **lbest, double *obj_particles, double *obj_lbest, problem prob, algoparam algo){
	int i,j; 
	int k,l;
	double temp1[prob.dim];


	for (i=0;i<algo.swarmSize;i++){
		for (j=(i-algo.radius);j<=(i+ algo.radius);j++ ){
			k=j;
			if (j<0)//start of the ring topology
				k=algo.swarmSize +j;
			else if (j>algo.swarmSize - 1)//end of ring topology
				k=j-algo.swarmSize;
			if (obj_particles[k]<obj_lbest[k]){//if new lbest found, evaluate new lbest objective
				for (l=0;l<prob.dim;l++)
					lbest[k][l]=particles[k][l];
				obj_lbest[k]=obj_particles[k];
				//cout << "obj_lbest = " << obj_lbest[k] << "\n";
			}
		}
	}

}


