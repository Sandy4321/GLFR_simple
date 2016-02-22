/*    
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/  
/* 
 * File:   main.c
 * Author: Ba-Dung
 *
 * Created on October 6, 2015, 3:48 PM
 * for network sizes<5000 nodes
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N_max 20000 // maximum number of nodes
#define n_max 500 // maximum number of communities 
#define m_max 250000 // maximum number of edges

typedef struct edgepair{
	unsigned int end1;
	unsigned int end2;
	bool flag;
};
struct edgepair lsgraph[m_max]; // edge list for the graph
unsigned long m, ml; // the expected number of edges in the graph, the number of edges created

unsigned int N=1000; // number of nodes
unsigned char k_avg=20, kmin = 1, kmax=50, smin=20, smax=100;
double gam_ma = 2.0, beta = 2.0; // exponents of the power law distributions for node degree and community size
double mu=0.5, delta=0.5; // mean and half-range of the uniform distribution for the fractions of external links of communities
float interval= 0.025; // minimum distance between any two different fractions of external links of communities

unsigned int n; // number of communities
unsigned char s[n_max]; // community sizes
double mu_c[n_max]; // fractions of external links assigned for communities

unsigned int Kc[n_max], Kc_int[n_max]; // total internal degree of communities, total internal degree of communities
double Kc_int_assigned[n_max]; // total expected internal degree of communities

unsigned int membership[N_max]; // community memberships of nodes

unsigned char k[N_max], k_int[N_max]; // node degrees, internal degrees of nodes
double k_int_assigned[N_max]; // expected degrees of nodes

// check if two edges have the same end points
bool duplicated_edge(struct edgepair edge1, struct edgepair edge2){
	if ((edge1.end1 == edge2.end1) && (edge1.end2 == edge2.end2)) return true;
	if ((edge1.end1 == edge2.end2) && (edge1.end2 == edge2.end1)) return true;
	return false;
}

// check if one edge is a multi-edge in a list in which each two consequence elements represent two end points of an edge
bool multi_edge(struct edgepair edge, unsigned int stubs[], unsigned int pos1, unsigned int pos2){	
	unsigned int i = pos1;
	while (i<pos2){
		struct edgepair edge2;
		edge2.end1 = stubs[i];
		edge2.end2 = stubs[i + 1];
		if (duplicated_edge(edge, edge2)) return true;		
		i = i + 2;
	}
	return false;
}

void generate_internal_links(){

        unsigned int stubs[2*m_max]; // list of stubs to generate links for each community
        unsigned int nstubs; // number of stubs in the stub list
        unsigned int n_multi_edges = 0, n_self_edges = 0;
        
	printf("generating links within communities...\n");

	// creating a list of stubs for each community
        for (unsigned int i = 0; i < n; i++){		
		nstubs = 0;
		for (unsigned int j = 0; j < N; j++){
			if (membership[j] == i){
				//add a number of stubs to the stub list, each stub is the id of the node, the number of stubs is the degree of the node
				for (unsigned int i2 = 0; i2 < k_int[j]; i2++){
					stubs[nstubs] = j;
					nstubs++;
				}
			}
		}

                // randomly reorder the list of the stubs in which each stub is reordered at least once
                unsigned int pos;
                unsigned int temp;
                for (unsigned int i = 0; i < nstubs; i++){
			pos = rand() % nstubs;
			temp = stubs[i];
			stubs[i] = stubs[pos];
			stubs[pos] = temp;
                }
                
		// rewire links to avoid self-edges
		for (unsigned int i = 0; i < nstubs-1; i=i+2){
			if (stubs[i] == stubs[i+1]){
				pos = 0;
				while (((stubs[i] == stubs[pos + 1]) || (stubs[pos] == stubs[i + 1])) && (pos<nstubs-1)){
					pos = pos+2;
				}
				// if we can find a stub to change an end point of the self-edge
				if (pos < nstubs-1){
					temp = stubs[i+1];
					stubs[i+1] = stubs[pos+1];
					stubs[pos+1] = temp;
				}
			}
		}

		// rewire links to avoid multi-edges
		for (unsigned int i = 0; i < nstubs-1; i=i+2){
			struct edgepair edge;
			edge.end1 = stubs[i];
			edge.end2 = stubs[i+1];

			// check if this edge is multi-edge
			if (multi_edge(edge, stubs, 0, i - 1) || multi_edge(edge, stubs, i + 2, nstubs - 1)){
				// search for an edge to rewire that can decrease the number of multi-edges
				pos = 0;
				struct edgepair edge_temp1, edge_temp2;
				edge_temp1.end1 = stubs[i];
				edge_temp1.end2 = stubs[pos + 1];
				edge_temp2.end1 = stubs[pos];
				edge_temp2.end2 = stubs[i + 1];
				while ((multi_edge(edge_temp1, stubs, 0, nstubs - 1) || multi_edge(edge_temp2, stubs, 0, nstubs - 1) || (stubs[i]==stubs[pos+1]) || (stubs[pos]==stubs[i+1])) && (pos < nstubs - 1)){
					pos = pos + 2;
					edge_temp1.end1 = stubs[i];
					edge_temp1.end2 = stubs[pos + 1];
					edge_temp2.end1 = stubs[pos];
					edge_temp2.end2 = stubs[i + 1];
				}

				// if we can find an edge to rewrite
				if (pos < nstubs - 1){
					temp = stubs[i + 1];
					stubs[i + 1] = stubs[pos + 1];
					stubs[pos+1] = temp;
				}
			}		
		}
                
		// connecting these stubs to create links within the community
		unsigned int mltemp;
		mltemp = ml;

		unsigned int i = 0;

		while (i< (nstubs-1)){		
			
			lsgraph[ml].end1 = stubs[i];			
			lsgraph[ml].end2 = stubs[i+1];
			lsgraph[ml].flag = true;
			i=i+2;
			ml++;
		}                 	
                
		//mark and count the number of self-edges and multi-edges in each community                
		for (unsigned int i = mltemp; i < ml; i++){
                    
			if (lsgraph[i].end1 == lsgraph[i].end2){
				lsgraph[i].flag = false; // this edge will not be counted
				n_self_edges++;
			}

			if (lsgraph[i].flag == true){
				for (unsigned int i2 = i + 1; i2 < ml; i2++){
					if (duplicated_edge(lsgraph[i], lsgraph[i2])){
						lsgraph[i2].flag = false;
						n_multi_edges++;
					}
				}
			}
		}	
	}

	//printf("total number of edges of the internal graphs= %d\n", ml);
	//printf("number of self-edges and number of multi-edges within communities= (%d, %d)\n", n_self_edges, n_multi_edges);
}
 
void generate_external_links(){		

	printf("generating links between communities...\n");

	// creating stubs for the external graph
	unsigned int stubs[2*m_max]; // list of stubs to generate links between communities
        int nstubs = 0;
        
	for (unsigned int i = 0; i < N; i++){		
		if (k[i]-k_int[i]>0) {
			for (unsigned int i2 = 0; i2 < k[i]-k_int[i]; i2++){
				stubs[nstubs] = i;
				nstubs++;					
			}
		}
	}
	
	if (nstubs == 0) return;

	// randomly reorder stubs in the list in which each stub is reordered at least once
	unsigned int pos;
	unsigned int temp;
	
	for (unsigned int i = 0; i < nstubs; i++){
		pos = rand() % nstubs;
		temp = stubs[i];
		stubs[i] = stubs[pos];
		stubs[pos] = temp;
	}	

	// remove edges that connect nodes in the same community
	for (unsigned int i = 0; i < nstubs - 1; i = i + 2){
		if (membership[stubs[i]] == membership[stubs[i + 1]]){
			pos = 0;
			while (((membership[stubs[i]] == membership[stubs[pos + 1]]) || (membership[stubs[i + 1]] == membership[stubs[pos]])) && (pos < nstubs - 1)){
				pos = pos + 2;
			}
			// if we can find an edge to rewire
			if (pos < nstubs - 1){
				temp = stubs[i + 1];
				stubs[i + 1] = stubs[pos + 1];
				stubs[pos + 1] = temp;
			}
		}
	}

	// remove multi-edges	
		for (unsigned int i = 0; i < nstubs - 1; i = i + 2){
			struct edgepair edge;
			edge.end1 = stubs[i];
			edge.end2 = stubs[i + 1];

			// check if this edge is a multi-edge	
			if (multi_edge(edge, stubs, 0, i - 1) || multi_edge(edge, stubs, i + 2, nstubs - 1)){
				// search for an edge to rewire that can decrease the number of multi-edges
				pos = 0;
				struct edgepair edge_temp1, edge_temp2;
				edge_temp1.end1 = stubs[i];
				edge_temp1.end2 = stubs[pos + 1];
				edge_temp2.end1 = stubs[pos];
				edge_temp2.end2 = stubs[i + 1];

				while ((multi_edge(edge_temp1, stubs, 0, nstubs - 1) || multi_edge(edge_temp2, stubs, 0, nstubs - 1) || (membership[stubs[i]] == membership[stubs[pos + 1]]) || (membership[stubs[pos]] == membership[stubs[i + 1]])) && (pos < nstubs - 1)){
					pos = pos + 2;
					edge_temp1.end1 = stubs[i];
					edge_temp1.end2 = stubs[pos + 1];
					edge_temp2.end1 = stubs[pos];
					edge_temp2.end2 = stubs[i + 1];
				}

				// if we can find an edge to rewrite
				if (pos < nstubs - 1){
					temp = stubs[i + 1];
					stubs[i + 1] = stubs[pos + 1];
					stubs[pos + 1] = temp;
				}
			}
		} // end for
	
	unsigned int mltemp = ml;
	unsigned int i = 0;
        
	// connecting these stubs to create the external graph
	while (i< (nstubs-1)){		
		lsgraph[ml].end1 = stubs[i];
		lsgraph[ml].end2 = stubs[i+1];
		lsgraph[ml].flag = true;		
		i = i+2;
		ml++;
	}
	
	// calculate the number of self-edges and multi-edges in the external graph
	unsigned int n_self = 0, n_multi = 0;
        
	for (unsigned int i = mltemp; i < ml; i++){

		if (membership[lsgraph[i].end1] == membership[lsgraph[i].end2]){
			lsgraph[i].flag = false;
			n_self++;
		}

		if (lsgraph[i].flag == true){			
			for (unsigned int i2 = i + 1; i2 < ml; i2++){
				if (lsgraph[i2].flag == true){
					if (duplicated_edge(lsgraph[i], lsgraph[i2])){
						n_multi++;
						lsgraph[i2].flag = false;
					}
				}
			}
		}
	}
        //printf("total number of edges of the external graph= %d\n", ml-mltemp);
	//printf("number of self-edges and number of multi-edges between communities= (%d,%d)\n", n_self, n_multi);
}

// this part of code is based on the source code of the LFR-benchmark at https://sites.google.com/site/andrealancichinetti/files
double integral(double a, double b) {

	if (fabs(a + 1)>1e-10)
		return (1 / (a + 1)*pow(b, a + 1));
	else
		return (log(b));
}

double average_degree(double dmax, double dmin, double gamma) {

	return (1 / (integral(gamma, dmax) - integral(gamma, dmin)))*(integral(gamma + 1, dmax) - integral(gamma + 1, dmin));

}

//bisection method to find the inferior limit, in order to have the expected average degree
double solve_dmin(double dmax, double dmed, double gamma) {

	double dmin_l = 1;
	double dmin_r = dmax;
	double average_k1 = average_degree(dmin_r, dmin_l, gamma);
	double average_k2 = dmin_r;

	while (fabs(average_k1 - dmed)>1e-7) {

		double temp = average_degree(dmax, ((dmin_r + dmin_l) / 2.), gamma);
		if ((temp - dmed)*(average_k2 - dmed)>0) {

			average_k2 = temp;
			dmin_r = ((dmin_r + dmin_l) / 2.);

		}
		else {

			average_k1 = temp;
			dmin_l = ((dmin_r + dmin_l) / 2.);

		}
	}

	return dmin_l;
}
// end of the part of the source code


double round1e7(double val){
	return (double)round(val * 1e7) / 1e7;
}

void generate_network(unsigned int N, unsigned char k_avg, unsigned char kmax, unsigned char smin, unsigned char smax, double mu, double delta, double gamma, double beta){
        
	// display the network information        
	printf("network information:\n");        
	printf("number of nodes= %d\n", N); // =1000 by default
	printf("average node degree= %d\n", k_avg); // =20 by default
	printf("maximum node degree= %d\n", kmax); // =50 by default
	printf("exponent for the node degree distribution= %.2f\n", -gamma); // =2.0 by default
	printf("exponent for the community size distribution= %.2f\n", -beta);  // =2.0 by default
	printf("minimum community size= %d\n", smin); // =20 by default
	printf("maximum community size= %d\n", smax); // =100 by default
        printf("average of the fractions of external links of communities= %.2f\n", mu); // =0.5 by default
        printf("maximum up and down from average for the fractions= %.2f\n", delta); // =0.5 by default
        printf("----------------------------------------\n");    
        
	kmin = solve_dmin(kmax, k_avg, -gamma) + 1; // adjust kmin to keep the average degree as expected by k_avg
                
        // generate node degrees following a power law distribution		
	m = 0;
	for (unsigned int i = 0; i < N; i++){
		double y = (double)(rand() % N + 1) / N;
                double exp;
                if (gamma==1){
                    // node degrees follow a uniform distribution
                    exp = rand() % (kmax-kmin+1) + kmin;
                }
                else{
                    // node degrees follow a power law distribution
                    exp = pow((pow(kmax, -gamma + 1) - pow(kmin, -gamma + 1))*y + pow(kmin, -gamma + 1), (double)1 / (-gamma + 1));
                }    
		k[i] = round(exp);
                //printf("k[%d]= %d\n", i, k[i]);
		m = m + k[i];
	}

	// the total number of edges in the network is the half of the total number of node degrees
	m = m / 2;	

	// generate community sizes randomly that the community sizes follow a power law distribution
	unsigned int Ntemp = N;
	unsigned int i = 0;

	while (Ntemp > 0){
			if (Ntemp <= smax){
				s[i] = Ntemp;
			}
			else{
				int max_c = N / smin;
				double y = (double)(rand() % max_c + 1) / max_c;  
                                double exp;
                                if (beta==1){
                                    // community sizes follow a uniform distribution
                                    exp = rand() % (smax-smin+1) + smin;
                                }
                                else{
                                    // community size follow a power law distribution
                                    exp = pow((pow(smax, -beta + 1) - pow(smin, -beta + 1))*y + pow(smin, -beta + 1), (double)1 / (-beta + 1));
                                }
				s[i] = round(exp);
                                //printf("s[%d]= %d\n", i, s[i]);
			}
			Ntemp = Ntemp - s[i];
			i++;
	}        
	n = i;	

	unsigned int index = 0;

	// assign nodes into communities
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < s[i]; j++){
			membership[j + index] = i;
		}
		index = index + s[i];
	}

	// if a node have the degree that is larger than its community size, change community of the node if possible
	for (i = 0; i < N; i++){
		if (k[i]>s[membership[i]]){
			unsigned int i2 = 0;
			bool found = false;
			while ((i2 < N) && (found == false)){
				if ((k[i2] <= s[membership[i]]) && (k[i] <= s[membership[i2]])){
					unsigned int temp = membership[i2];
					membership[i2] = membership[i];
					membership[i] = temp;
					found = true;
				}
				i2++;
			}
		}
	}

	// calculate total degree Kc for each community
	for (i = 0; i < n; i++){
		Kc[i] = 0;
	}

	for (i = 0; i < N; i++){
		Kc[membership[i]] = Kc[membership[i]] + k[i];
	}

	unsigned int Kc_max = Kc[0];
	for (i = 1; i < n; i++){
		if (Kc_max < Kc[i])  Kc_max = Kc[i];
	}	

        // set the minimum valid value and the maximum valid value for the fractions of external links of communities
	double left = 0.025, right;
	right = floor((double)(2 * m - Kc_max) / (2 * m) / 0.025)*0.025;

	//printf("minimum valid value= %.3f, maximum valid value= %.3f\n", left, right);	

        // adjusted mu if it is invalid
        if (mu <left) mu = left;
	if (mu > right) mu = right;
        
        // generate the fractions of external links assigned for communities		        
	double mixing_min, mixing_max, mu_avg, mu_actual= 0;		
        
        // calculate the left bound and the right bound of the range of the uniform distribution
	mixing_min = mu - delta;
	mixing_max = mu + delta;

	// adjust the range if it contains invalid values for the fractions
	if (mixing_min < left){
		mixing_min = left;
		mixing_max = mu + mu - mixing_min;
	}	
	
	if (mixing_max >right){
		mixing_max = right;
		mixing_min = mu - (mixing_max - mu);
	}
        
        // assign the fractions of external links for communities
        // the fractions are uniformly distributed in the range from mixing_min to mixing_max
	for (i = 0; i < n; i++){
		int i_intterval = rand() % ((int)round((mixing_max - mixing_min) /interval)+1);
		mu_c[i] = i_intterval*interval + mixing_min;
		mu_avg = mu_avg + mu_c[i];
	}		
		
        // if the average of the fractions assigned for communities is different from mu
        // reduce the gap between the generated value and the expected value to control the ambiguity of the community structure				
		bool exit = false;

		while ((fabs(mu - mu_avg / n) >= 0.0045) && (!exit)){					
			int pos = rand() % n;
			if ((mu - mu_avg / n) >= 0.0045){
				if (round1e7(mu_c[pos]) < round1e7(mixing_max)){	
							mu_c[pos] = mu_c[pos] + 0.025;	
							mu_avg = mu_avg + 0.025;	
							if ((mu - mu_avg / n) <0){
								exit = true;
							}
				}						
			}
			else{
				if (round1e7(mu_c[pos]) > round1e7(mixing_min)){
							mu_c[pos] = mu_c[pos] - 0.025;
							mu_avg = mu_avg - 0.025;
							if ((mu - mu_avg / n) >0){
								exit = true;
							}
				}						
			}                    
		}	
			
	// at this stage, we already have individual fractions of external links assigned for communities
			
	mu_avg = mu_avg / n;				

	// calculating the internal degree and external degree for nodes
	for (i = 0; i < N; i++){

            	if (k[i] == 0) {
			printf("error: the degree of node %d is zero!\n");
			return;
		}
                        
		// calculated k_int_assigned from mu_c of communities			
		k_int_assigned[i] = k[i] - (mu_c[membership[i]] * k[i]);
                k_int[i] = round(k_int_assigned[i]);

		// ensure the condition for the existence of the community that p_in_i > p_out_i
		if (Kc[membership[i]] == 0) {
			printf("error: the community %d of the node %d have no link!\n", membership[i], i);
			return;
                }

		if (2 * m - Kc[membership[i]] !=0){
                        while ((double)k_int[i] / Kc[membership[i]] < (double)(k[i] - k_int[i]) / (2 * m - Kc[membership[i]])){
                            	k_int[i]++;
                        }
                }

                mu_actual = mu_actual + (double)(k[i] - k_int[i]) / k[i];
	}

	mu_actual = mu_actual / N;

	// create links for the graph	
	ml = 0;

	generate_internal_links();
	generate_external_links();	

	// sort the list so that end1 < end2 to export graph
	for (i = 0; i < ml; i++){
		if (lsgraph[i].end1 > lsgraph[i].end2){
			unsigned int temp = lsgraph[i].end1;
			lsgraph[i].end1 = lsgraph[i].end2;
			lsgraph[i].end2 = temp;
		}
	}

	// sort lsgraph according to end1	
	for (unsigned int i = 0; i < ml; i++){
		for (unsigned int j = i + 1; j < ml; j++){
			if ((lsgraph[j].end1 < lsgraph[i].end1) || ((lsgraph[j].end1 == lsgraph[i].end1) && (lsgraph[j].end2 < lsgraph[i].end2))){
				struct edgepair temp;
				temp = lsgraph[i];
				lsgraph[i] = lsgraph[j];
				lsgraph[j] = temp;
			}
		}
	}			
	
        // export the list of links and community memberships of nodes into files
        // set file names
	char fname[250];

        FILE *fcommunity, *fnetwork, *fstatistic;
        char fname_community[250], fname_network[250], fname_statistic[250];

	sprintf(fname, "N%dk%dkmax%dsmin%dsmax%dmu%.2fdelta%.2f", N, k_avg, kmax, smin, smax, mu, delta);

	sprintf(fname_community, "%s%s", fname, "_community.dat");
	sprintf(fname_network, "%s%s", fname, "_network.dat");
	sprintf(fname_statistic, "%s%s", fname, "_statistic.dat");

	// export only the edges that have flag=true to create the network
	fcommunity = fopen(fname_community, "w");
	for (i = 0; i < N; i++){
		fprintf(fcommunity, "%d %d\n", i + 1, membership[i] + 1);
	}
	fclose(fcommunity);
        
        int mltemp=0;
	fnetwork = fopen(fname_network, "w");
	for (i = 0; i < ml; i++){
		if (lsgraph[i].flag == true) {
                    fprintf(fnetwork, "%d %d\n", lsgraph[i].end1 + 1, lsgraph[i].end2 + 1);
                    mltemp++;
                }
	}
	fclose(fnetwork);
        
        printf("number of nodes in the network: %d\n", N);        
        printf("number of edges in the network: %d\n", mltemp);
        printf("average node degree= %.2f\n", (float)2*m/N);
        printf("average of the fractions of external links of nodes= %.2f\n", mu_actual);
}

int parameter(char *st){
    if (strcmp(st, "-N")==0) return 1;
    if (strcmp(st, "-k")==0) return 2;
    if (strcmp(st, "-kmax")==0) return 3;
    if (strcmp(st, "-smin")==0) return 4;
    if (strcmp(st, "-smax")==0) return 5;
    if (strcmp(st, "-mu")==0) return 6;
    if (strcmp(st, "-delta")==0) return 7;
    if (strcmp(st, "-gamma")==0) return 8;
    if (strcmp(st, "-beta")==0) return 9;
    return 0;
}

int main(int argc, char** argv) {    
    
    time_t t;
    srand((unsigned)time(&t));
    
    int i_arg=1;
    while (i_arg< argc){
        //printf("argv[%d]= %s\n", i, argv[i]);        
        switch(parameter(argv[i_arg])){
            case 1:
                N = atoi(argv[i_arg+1]);
                break;
            case 2:
                k_avg = atoi(argv[i_arg+1]);
                break;
            case 3:
                kmax = atoi(argv[i_arg+1]);
                break;
            case 4:
                smin = atoi(argv[i_arg+1]);
                break;
            case 5:
                smax = atoi(argv[i_arg+1]);
                break;
            case 6:
                mu = atof(argv[i_arg+1]);
                break;
            case 7:
                delta = atof(argv[i_arg+1]);
                break;
            case 8:
                gam_ma = atof(argv[i_arg+1]);
                break;
            case 9:
                beta = atof(argv[i_arg+1]);
                break;
            default:
                printf("%s is a wrong input parameter, see the Readme.txt file for more information!\n", argv[i_arg]);
                return (EXIT_SUCCESS);
        }
        i_arg = i_arg+2;
    }
    
    generate_network(N, k_avg, kmax, smin, smax, mu, delta, gam_ma, beta);
    
    return (EXIT_SUCCESS);
}

