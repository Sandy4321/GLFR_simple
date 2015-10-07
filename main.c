/* 
 * File:   main.c
 * Author: Ba-Dung
 *
 * Created on October 6, 2015, 3:48 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define N_max 5000 // maximum number of nodes
#define n_max 500 // maximum number of communities 
#define m_max 200000 // maximum number of edges

typedef struct edgepair{
	unsigned int end1;
	unsigned int end2;
	bool flag;
};
struct edgepair lsgraph[N_max]; // edge list for the graph

unsigned int N=1000;
unsigned char k_avg=20, kmin = 1, kmax=50, smin=20, smax=100;
double gamma = 2.0, beta = 2.0; // exponents of the power law distributions for node degree and community size
double mu=0.5, delta=0.5; // mean and half-range of the uniform distribution for the fractions of external links of communities
float interval= 0.025; // minimum distance between any two different fractions of external links of communities

unsigned int n; // number of communities
unsigned char s[n_max]; // community sizes
double mixing[n_max]; // fractions of external links assigned for communities

unsigned int Kc[n_max], Kc_in[n_max]; // total internal degree of communities, total internal degree of communities
double Kc_in_expected[n_max]; // total expected internal degree of communities

unsigned int membership[N_max]; // community memberships of nodes

unsigned char k[N_max], k_in[N_max]; // node degrees, internal degrees of nodes
double k_in_expected[N_max]; // expected degrees of nodes

unsigned int stubs[2*m_max]; // list of stubs to generate edges
unsigned int nstubs = 0; // number of stubs in the stub list
unsigned long m, ml; // the number of edges in the graph, the number of edges in each subgraph

bool duplicated_edge(struct edgepair edge1, struct edgepair edge2){
	if ((edge1.end1 == edge2.end1) && (edge1.end2 == edge2.end2)) return true;
	if ((edge1.end1 == edge2.end2) && (edge1.end2 == edge2.end1)) return true;
	return false;
}

bool duplicated_edge2(struct edgepair edge, unsigned int pos1, unsigned int pos2){	
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

	//printf("generating internal links...\n");

	// generating edges within each community
	unsigned int n_dup = 0, n_selfloop = 0;

	for (unsigned int i = 0; i < n; i++){
		// creating a list of stubs for each community
		nstubs = 0;
		for (unsigned int j = 0; j < N; j++){
			if (membership[j] == i){
				//add the number of internal degree of vertex as stubs to the stub list								
				for (unsigned int i2 = 0; i2 < k_in[j]; i2++){
					stubs[nstubs] = j;
					nstubs++;
				}
			}
		}
	
		unsigned int pos;
		unsigned int temp;

		// now randomly reorder the list in which each stub is reordered at least once
		for (unsigned int i = 0; i < nstubs; i++){
			pos = rand() % nstubs;
			temp = stubs[i];
			stubs[i] = stubs[pos];
			stubs[pos] = temp;
		}
                
		// rewire edges if selfloop
		for (unsigned int i = 0; i < nstubs-1; i=i+2){
			if (stubs[i] == stubs[i+1]){
				pos = 0;
				while (((stubs[i] == stubs[pos + 1]) || (stubs[pos] == stubs[i + 1])) && (pos<nstubs-1)){
					pos = pos+2;
				}
				// if we can find an edge to rewire
				if (pos < nstubs-1){
					temp = stubs[i+1];
					stubs[i+1] = stubs[pos+1];
					stubs[pos+1] = temp;
				}
			}
		}

		// rewire edges if duplicated
		for (unsigned int i = 0; i < nstubs-1; i=i+2){
			struct edgepair edge;
			edge.end1 = stubs[i];
			edge.end2 = stubs[i+1];

			// check if this edge is duplicated	
			// rewrire if this edge is a duplicated edge

			if (duplicated_edge2(edge, 0, i - 1) || duplicated_edge2(edge, i + 2, nstubs - 1)){
				// search for an edge to rewire that can decrease the number of duplicated edges
				pos = 0;
				struct edgepair edge_temp1, edge_temp2;
				edge_temp1.end1 = stubs[i];
				edge_temp1.end2 = stubs[pos + 1];
				edge_temp2.end1 = stubs[pos];
				edge_temp2.end2 = stubs[i + 1];
				while ((duplicated_edge2(edge_temp1, 0, nstubs - 1) || duplicated_edge2(edge_temp2, 0, nstubs - 1) || (stubs[i]==stubs[pos+1]) || (stubs[pos]==stubs[i+1])) && (pos < nstubs - 1)){
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
                
		// connecting these stubs to create edges
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
                
		//mark and count the selfloop and duplicated edges in each community again after rewiring

		for (unsigned int i = mltemp; i < ml; i++){
			if (lsgraph[i].end1 == lsgraph[i].end2){
				lsgraph[i].flag = false;
				n_selfloop++;
			}

			if (lsgraph[i].flag == true){
				for (unsigned int i2 = i + 1; i2 < ml; i2++){
					if (duplicated_edge(lsgraph[i], lsgraph[i2])){
						lsgraph[i2].flag = false;
						n_dup++;
					}
				}
			}
		}	
	}

	printf("total number of edges of the internal graphs= %d\n", ml);
	printf("total number of self-edges and multi-edges in the internal graphs= (%d, %d)\n", n_selfloop, n_dup);	
}
 
void generate_external_links(){	

	//generating edges between different communities by attaching stubs and rewire links

	//printf("generating external links...\n");

	// creating external stubs for communities
	nstubs = 0;
	for (unsigned int i = 0; i < N; i++){		
		//nstubs = nstubs + k_out[i];
		if (k[i]-k_in[i]>0) {
			for (unsigned int i2 = 0; i2 < k[i]-k_in[i]; i2++){
				stubs[nstubs] = i;
				nstubs++;					
			}
		}
	}
	
	if (nstubs == 0) return;

	// randomly reorder the list of external stubs
	unsigned int pos;
	unsigned int temp;

	// now randomly reorder the list in which each stub is reordered at least once
	for (unsigned int i = 0; i < nstubs; i++){
		pos = rand() % nstubs;
		temp = stubs[i];
		stubs[i] = stubs[pos];
		stubs[pos] = temp;
	}	

	// remove self-edges and edges connecting nodes in the same community
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

			// check if this edge is duplicated	
			// rewrire if this edge is a duplicated edge

			if (duplicated_edge2(edge, 0, i - 1) || duplicated_edge2(edge, i + 2, nstubs - 1)){
				// search for an edge to rewire that can decrease the number of duplicated edges
				pos = 0;
				struct edgepair edge_temp1, edge_temp2;
				edge_temp1.end1 = stubs[i];
				edge_temp1.end2 = stubs[pos + 1];
				edge_temp2.end1 = stubs[pos];
				edge_temp2.end2 = stubs[i + 1];

				while ((duplicated_edge2(edge_temp1, 0, nstubs - 1) || duplicated_edge2(edge_temp2, 0, nstubs - 1) || (membership[stubs[i]] == membership[stubs[pos + 1]]) || (membership[stubs[pos]] == membership[stubs[i + 1]])) && (pos < nstubs - 1)){
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
	unsigned int n_dup = 0, n_selfloop = 0;
	for (unsigned int i = mltemp; i < ml; i++){

		if (membership[lsgraph[i].end1] == membership[lsgraph[i].end2]){
			lsgraph[i].flag = false;
			n_selfloop++;
		}

		if (lsgraph[i].flag == true){			
			for (unsigned int i2 = i + 1; i2 < ml; i2++){
				if (lsgraph[i2].flag == true){
					if (duplicated_edge(lsgraph[i], lsgraph[i2])){
						n_dup++;
						lsgraph[i2].flag = false;
					}
				}
			}
		}
	}
        printf("total number of edges of the external graph= %d\n", ml-mltemp);
	printf("number of self-edges and multi-edges in the external graph= (%d,%d)\n", n_selfloop, n_dup);
}

// this part of code is based on the source code of the LFR-benchmark at https://sites.google.com/site/andrealancichinetti/files
double integral(double a, double b) {

	if (fabs(a + 1.)>1e-10)
		return (1. / (a + 1.)*pow(b, a + 1.));
	else
		return (log(b));
}

double average_degree(double dmax, double dmin, double gamma) {

	return (1. / (integral(gamma, dmax) - integral(gamma, dmin)))*(integral(gamma + 1, dmax) - integral(gamma + 1, dmin));

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

double round1000(double val){
	return (double)round(val * 1000) / 1000;
}

void generate_network(unsigned int N, unsigned char k_avg, unsigned char kmax, unsigned char smin, unsigned char smax, double mu, double delta, double gamma, double beta){
        
	// display the network information
	printf("\nnetwork information:\n\n");
        printf("----------------------------------------\n");
	printf("number of nodes= %d\n", N); // =1000 by default
	printf("average node degree= %d\n", k_avg); // =20 by default
	printf("maximum node degree= %d\n", kmax); // =50 by default
	printf("exponent for the degree distribution= %.2f\n", -gamma); // =2.0 by default
	printf("exponent for the community size distribution= %.2f\n", -beta);  // =2.0 by default
	printf("minimum community size= %d\n", smin); // =20 by default
	printf("maximum community size= %d\n", smax); // =100 by default
        printf("average of the fractions of external links of communities= %.2f\n", mu); // =0.5 by default
        printf("maximum up and down from the average for the fractions= %.2f\n", delta); // =0.5 by default
	printf("----------------------------------------\n");
        
	kmin = solve_dmin(kmax, k_avg, -gamma) + 1; // adjust kmin to keep the average degree as expected by k_avg
                
        // generate node degrees following a power law distribution		
	m = 0;
	for (unsigned int i = 0; i < N; i++){
		double y = (double)(rand() % N + 1) / N;
		double exp = pow((pow(kmax, -gamma + 1) - pow(kmin, -gamma + 1))*y + pow(kmin, -gamma + 1), (double)1 / (-gamma + 1));
		k[i] = round(exp);
		m = m + k[i];
	}

	// the total number of edges in the network is the half of the total number of node degrees
	m = m / 2;	

	// generate community sizes following a power law distribution
	unsigned int Ntemp = N;
	unsigned int i = 0;

	while (Ntemp > 0){
			if (Ntemp <= smax){
				s[i] = Ntemp;
			}
			else{
				int max_c = N / smin;
				double y = (double)(rand() % max_c + 1) / max_c;
				double exp = pow((pow(smax, -beta + 1) - pow(smin, -beta + 1))*y + pow(smin, -beta + 1), (double)1 / (-beta + 1));
				s[i] = round(exp);
			}
			Ntemp = Ntemp - s[i];
			i++;
	}

	n = i;	

	unsigned int index = 0;

	// randomly generate community membership
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < s[i]; j++){
			membership[j + index] = i;
		}
		index = index + s[i];
	}

	// if the vertex degree is larger than its community size, it is exchanged with another vertex
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

        // set the minimum valid value and the maximum valid value for the fractions
	double left = 0.025, right;
	right = floor((double)(2 * m - Kc_max) / (2 * m) / 0.025)*0.025;

	//printf("left= %.3f, right= %.3f\n", left, right);

	// generate the fractions of external links assigned for communities		
        
	double mixing_min, mixing_max, mu_expected=0, mu_actual = 0;		

	// generate mixing for each community randomly from mu +/- delta			

        if (mu <left) mu = left;
	if (mu > right) mu = right;
        
        // calculate the left bound and the right bound of the range of the uniform distribution
	mixing_min = mu - delta;
	mixing_max = mu + delta;

	// adjust the range if it contains the values that less then the minimum valid value
	if (mixing_min < left){
		mixing_min = left;
		mixing_max = mu + mu - mixing_min;
	}	

	// adjust the range if it contains values that greater than the maximum valid value
	if (mixing_max >right){
		mixing_max = right;
		mixing_min = mu - (mixing_max - mu);
	}		

        // assign the fractions of external links for communities
        // the fractions are uniformly distributed in the range from mixing_min to mixing_max
	for (i = 0; i < n; i++){
		int i_interval = rand() % ((int)round((mixing_max - mixing_min) /0.025)+1);
		mixing[i] = i_interval*0.025 + mixing_min;
		mu_expected = mu_expected + mixing[i];
	}		
		
        // if the average of the fractions generated is different from the expected average defined by mu, reduce the gap	
				
		bool exit = false;

		while ((fabs(mu - mu_expected / n) >= 0.0045) && (!exit)){					
			int pos = rand() % n;
			if ((mu - mu_expected / n) >= 0.0045){
				if (round1000(mixing[pos]) < round1000(mixing_max)){	
							mixing[pos] = mixing[pos] + 0.025;	
							mu_expected = mu_expected + 0.025;	
							if ((mu - mu_expected / n) <0){
								exit = true;
							}
				}						
			}
			else{
				if (round1000(mixing[pos]) > round1000(mixing_min)){
							mixing[pos] = mixing[pos] - 0.025;
							mu_expected = mu_expected - 0.025;
							if ((mu - mu_expected / n) >0){
								exit = true;
							}
				}						
			}                    
		}	
			
	// now we already have individual mixing for communities
			
	mu_expected = mu_expected / n;			

	//printf("mu= %.3f, mixing_min= %.3f, mixing_max= %.3f, mu_expected= %.3f\n",mu, mixing_min, mixing_max, mu_expected);			
                
        for (i = 0; i < n; i++){
		Kc_in[i] = 0;
		Kc_in_expected[i] = 0;
        }

	// calculating the internal degree and external degree for nodes
	for (i = 0; i < N; i++){

		// calculated k_in_expected from mu_expected of communities			
		k_in_expected[i] = k[i] - (mixing[membership[i]] * k[i]);
		Kc_in_expected[membership[i]] = Kc_in_expected[membership[i]] + k_in_expected[i];
                k_in[i] = round(k_in_expected[i]);

		// ensure the condition of community
		// conflict with the condition of existing community -> increase k_in	

		if (Kc[membership[i]] == 0) {
			printf("the number of internal links of the community %d is zero!\n", i);
			return;
                }

		if (2 * m - Kc[membership[i]] == 0) {
			printf("the number of external links of the community %d is zero!\n", i);
			return;
		}

		while ((double)k_in[i] / Kc[membership[i]] < (double)(k[i] - k_in[i]) / (2 * m - Kc[membership[i]])){
			k_in[i]++;
		}			

		Kc_in[membership[i]] = Kc_in[membership[i]] + k_in[i];
		if (k[i] == 0) {
			printf("the degree of node %d is zero!\n");
			return;
		}
		
                mu_actual = mu_actual + (double)(k[i] - k_in[i]) / k[i];				
	}

	mu_actual = mu_actual / N;

	// create edges by connecting internal and external stubs of each vertices
	// generating internal edges
	ml = 0;

	generate_internal_links();
	generate_external_links();
	
	unsigned int n_removed = 0;
	for (unsigned int i = 0; i < ml; i++){
		if (lsgraph[i].flag == false) n_removed++;
	}

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
			
        // set file names
	char fname[250];

        FILE *fcommunity, *fnetwork, *fstatistic;
        char fname_community[250], fname_network[250], fname_statistic[250];

	sprintf(fname, "N%d_k%d_kmax%d_smin%d_smax%d_mu%.2f_fhalf-range%.2f", N, k_avg, kmax, smin, smax, mu, delta);

	sprintf(fname_community, "%s%s", fname, "_community.dat");
	sprintf(fname_network, "%s%s", fname, "_network.dat");
	sprintf(fname_statistic, "%s%s", fname, "_statistic.dat");

	// exporting the network
	fcommunity = fopen(fname_community, "w");
	for (i = 0; i < N; i++){
		fprintf(fcommunity, "%d %d\n", i + 1, membership[i] + 1);
	}
	fclose(fcommunity);

	fnetwork = fopen(fname_network, "w");
	for (i = 0; i < ml; i++){
		if (lsgraph[i].flag == true) fprintf(fnetwork, "%d %d\n", lsgraph[i].end1 + 1, lsgraph[i].end2 + 1);
	}
	fclose(fnetwork);
        
        printf("-------------------------------------------\n");
        printf("total number of nodes in the network: %d\n", N);
        printf("total number of edges in the network: %d\n", ml);
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
    
    int i=1;
    while (i< argc){
        //printf("argv[%d]= %s\n", i, argv[i]);        
        switch(parameter(argv[i])){
            case 1:
                N = atoi(argv[i+1]);
                break;
            case 2:
                k_avg = atoi(argv[i+1]);
                break;
            case 3:
                kmax = atoi(argv[i+1]);
                break;
            case 4:
                smin = atoi(argv[i+1]);
                break;
            case 5:
                smax = atoi(argv[i+1]);
                break;
            case 6:
                mu = atof(argv[i+1]);
                break;
            case 7:
                delta = atof(argv[i+1]);
                break;
            case 8:
                gamma = atof(argv[i+1]);
                break;
            case 9:
                beta = atof(argv[i+1]);
                break;
            default:
                printf("%s is a wrong input parameter, a default value may be used!\n", argv[i]);
        }
        i = i+2;
    }
    
    generate_network(N, k_avg, kmax, smin, smax, mu, delta, gamma, beta);
    
    return (EXIT_SUCCESS);
}

