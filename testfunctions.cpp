// BENCHMARK FUNCTIONS //
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

double Objective (double X[], problem prob){
	
	int  i,j ;
	double tmp1, tmp2[prob.dim];
	double pi=3.14159265358979;
	double sol[prob.dim];
	double F; 

//FOR SIMULATION, GIVE prob.bench=100, and put inside if-loop the simulator caller, with input the 
// decision vector X[] and output the Objective value F.

	if ( prob.bench==100 ){

		// call simulation

	}



    // SPHERE, X=[0,...,0], F=0, [-100,100]
	//--------------------------------------
	if (prob.bench == 0) {
	  F = 0.0;
	  for (i=0; i<prob.dim; i++)
	 	 F = F + pow(X[i],2);
	  for (i=0; i<prob.dim; i++)
		 sol[i] = 0;
    }

	// ROSENBROCK, X=[1,...,1], F=0, [-30,30]
	//----------------------------------------
	else if (prob.bench == 1) {
			F = 0.0;
			for (i=0; i<prob.dim-1; i++)
				F = F + 100.0*pow((X[i+1]-X[i]*X[i]),2) + pow((X[i]-1.0),2);
			for (i=0; i<prob.dim; i++)
				sol[i] = 1;
	}

    // RASTRIGIN, X=[0,...,0], F=0, [-5.12,5.12]
	//-------------------------------------------
    else if (prob.bench == 2) {
		F = 0.0;
		for (i=0; i<prob.dim; i++)
			F = F + pow(X[i],2)  + 10 - 10.0*cos(2*pi*X[i]);
		for (i=0; i<prob.dim; i++)
			sol[i] = 0;
	}

	// GRIEWANK, X=[0,...,0], F=0, [-600,600]
	//---------------------------------------
	else if (prob.bench == 3) {
		tmp1 = 0.0;
	    for (i=0; i<prob.dim; i++)
			tmp1 = tmp1 + pow(X[i],2);
		F = tmp1/4000.0;
		tmp1 = 1.0;
		for (i=0; i<prob.dim; i++)
			tmp1 = tmp1*cos(X[i]/sqrt((double) i + 1.0));
		F = F - tmp1 + 1.0;
		for (i=0; i<prob.dim; i++)
			sol[i] = 0;
    }


    // ACKLEY, X=[0,...,0], F=0, [-20,30]
	//------------------------------------
	else if (prob.bench == 4) {
		F   = 20.0 + exp(1.0);
		tmp1 = 0.0;
		for (i=0; i<prob.dim; i++)
			tmp1 = tmp1 + X[i]*X[i];
		tmp1 = sqrt(tmp1/prob.dim);
		F   = F - 20.0*exp(-0.2*tmp1);
		tmp1 = 0.0;
		for (i=0; i<prob.dim; i++)
			tmp1 = tmp1 + cos(2.0*pi*X[i]);
		tmp1 = tmp1/prob.dim;
		F   = F - exp(tmp1);
		for (i=0; i<prob.dim; i++)
			sol[i] = 0;
	}

    // FREUDENSTEIN-ROTH, X=[5,4], F=0
	//---------------------------------
	else if (prob.bench == 5) {
		tmp1 = ( -13.0 + X[0] + ( (5.0-X[1])*X[1] - 2.0 ) * X[1] );
		F   = tmp1 * tmp1;
		tmp1 = ( -29.0 + X[0] + ( (X[1]+1.0)*X[1] - 14.0) * X[1]);
		F   = F + tmp1 * tmp1;
	}

    // LEVY No.3, 18 global minimizers with F=-176.542
	//-------------------------------------------------
	else if (prob.bench == 6) {
		tmp1 = 0.0;
		for (i=0; i<5; i++)
			tmp1 = tmp1 + (i+1.0) * cos(i*X[0]+i+1.0);
		F   = tmp1;
		tmp1 = 0.0;
		for (i=0; i<5; i++)
			tmp1 = tmp1 + (i+1.0) * cos((i+2.0)*X[1]+i+1.0);
		F   = F * tmp1;
	}

    // LEVY, X=[1,...,1], F=0, DIM=3,4,5,8,10
	//----------------------------------------
 	else if (prob.bench == 7) {
		for (i=0; i<prob.dim; i++)
			tmp2[i] = 1.0 + (X[i]-1.0)/4.0;
		tmp1 = 0.0;
		F   = sin(pi*tmp2[0])*sin(pi*tmp2[0]) + (tmp2[prob.dim-1]-1.0)*(tmp2[prob.dim-1]-1.0);
		for (i=0; i<(prob.dim-1); i++) {
			tmp1 = tmp2[i]-1.0;
			tmp1 = tmp1 * tmp1 * ( 1.0 + 10.0*sin(pi*tmp2[i+1])*sin(pi*tmp2[i+1]));
			F   = F + tmp1;
		}
	}

    // POWELL BADLY SCALED, X=[1.098*10^(-5),9.106], F=0
	//---------------------------------------------------
	else if (prob.bench == 8) {
		tmp1 = 10000.0*X[0]*X[1] - 1.0;
		F   = tmp1 * tmp1;
		tmp1 = exp(-X[0]) + exp(-X[1]) - 1.0001;
		F   = F + tmp1 * tmp1;
	}


    // WOOD, X=[1,1,1,1], F=0
	//------------------------
	else if (prob.bench == 9) {
		tmp1 = 10.0*(X[1]-X[0]*X[0]);
		F   = tmp1 * tmp1;
		tmp1 = (1-X[0]);
		F   = F + tmp1 * tmp1;
		tmp1 = (X[3]-X[2]*X[2]);
		F   = F + 90.0*tmp1*tmp1;
		tmp1 = (1-X[2]);
		F   = F + tmp1 * tmp1;
		tmp1 = (X[1]+X[3]-2.0);
		F   = F + 10.0*tmp1*tmp1;
		tmp1 = (X[1]-X[3]);
		F   = F + 0.1*tmp1*tmp1;
	}

    // TRIGONOMETRIC, F=0
	//--------------------
	else if (prob.bench == 10) {
		tmp1 = 0.0;
		for (i=0; i<prob.dim; i++)
			tmp1 = tmp1 + cos(X[i]);
		F = 0.0;
		for (i=0; i<prob.dim; i++)
			tmp2[i] = prob.dim - tmp1 + (i+1.0)*(1.0-cos(X[i])) - sin(X[i]);
		for (i=0; i<prob.dim; i++)
			F = F + tmp2[i]*tmp2[i];
	}

    // SCHAFFER F6, X=[0,...,0], F=0, [-100,100]
	//-------------------------------------------
	else if (prob.bench == 11) {
		tmp1 = sin( sqrt(X[0]*X[0]+X[1]*X[1]) );
		tmp1 = tmp1 * tmp1 - 0.5;
		F   = tmp1;
		tmp1 = 1.0 + 0.001*(X[0]*X[0]+X[1]*X[1]);
		tmp1 = tmp1 * tmp1;
		F   = 0.5 + F/tmp1;
	}

    // SCHWEFEL'S DOUBLE SUM, X=[0,...,0], F=0, [-65.536, 65.536]
	//------------------------------------------------------------
	else if (prob.bench == 12) {
		F = 0.0;
		for (i=0; i<prob.dim; i++) {
			tmp1 = 0.0;
			for (j=0; j<=i; j++)
				tmp1 = tmp1 + X[j];
			F = F + tmp1*tmp1;
		}
	}

    // SCHWEFEL FUNCTION
    //-------------------
	else if (prob.bench == 13) {
		F = 418.9829*prob.dim;
		tmp1=0.0;
		for (j=0; j<prob.dim; j++)
			tmp1 += X[j]*sin(sqrt(abs(X[j])));
		F = F + tmp1;
	}

	//for (i=0; i<prob.dim; i++)
	//		tmp2[i] = X[i]-sol[i];
	//*DIST = norms(prob.dim, tmp2, 1);

	return F;

}	
