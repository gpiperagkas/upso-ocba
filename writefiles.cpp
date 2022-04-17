/********************************/
// OUTPUT FILES FUNCTIONS //
/********************************/
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

void outfile(FILE *fout, double mean, double std, double *uf, double unifactor, problem prob, algoparam algo){



if ((prob.dim==2)&&(uf[0]==0.0)&&(uf[1]==0.0)&&(uf[2]==0.0)){
	fprintf(fout,"******************************************\n");
	fprintf(fout,"PROBLEM PARAMETERS:\n");
    fprintf(fout,"******************************************\n");
	fprintf(fout,"TEST PROBLEM: %d\n",prob.bench);
	fprintf(fout,"SWARM SIZE: %d\n", algo.swarmSize);
	fprintf(fout,"LOCAL PSO RADIUS: %d\n",algo.radius);
	fprintf(fout,"PSO PARAMETERS x, c1, c2 : %5.3f %4.2f %4.2f\n", algo.xi, algo.c1, algo.c2);
	fprintf(fout,"EXPERIMENTS 100\n" );
	fprintf(fout,"_____________RESULTS_______________\n");
	fprintf(fout,"DIMs - UniFact - Mean Value - Stand. Dev. \n");
}

	fprintf(fout, "%d  - ", prob.dim);
	fprintf(fout, "%2.1f - ",unifactor);
	fprintf(fout, "%le - ",mean);
	fprintf(fout, "%le \n",std);

//	fprintf(fout,"MEAN FUNCTION EVALUATIONS: %6.1f\n",fev);
//	fprintf(fout,"SUCCESS RATE: %3.2f\n",suc);
//	fprintf(fout,">>>>>>>>>>>>>>>>>>>>>MEAN OBJECTIVE VALUE: %8.5f\n",mean);
//	fprintf(fout,">>>>>>>>>>>>>>>>>STANDARD DEVIATION VALUE: %6.4f\n",std);
//	fprintf(fout,"******************************************\n");

}
