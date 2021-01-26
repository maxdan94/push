/*
Maximilien Danisch
Septembre 2017
http://bit.ly/danisch
maximilien.danisch@gmail.com

to compile:
gcc RootedPageRank.c -o RootedPageRank -O9

to execute:
./RoutedPageRank net.txt source eps res.txt
- net.txt should contain on each line two unsigned separated by a space: "source target\n" that is the input directed graph.
- source: the id of the root node: the random walk restarts from that node with probability 0.15
- eps is the precision: infinite norm of pr_{new}-pr_{old}<eps to stop
- res.txt will contain an approximation of the pagerank (30 iterations using the power iteration method). "nodeID PageRankValue\n" on each line.

to sort the outputvby value:
LC_NUMERIC=C sort -gr -k2,2 res.txt >resSORT.txt
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <time.h>

#define NLINKS 10000000 //maximum number of nonzero entries, will increase if needed
#define ALPHA 0.15 //restart probability of pagerank

// graph datastructure:
typedef struct {
	unsigned long s;//source node
	unsigned long t;//target node
} edge;

//sparse graphe structure
typedef struct {
	unsigned long n;//number of nodes
	unsigned long long e;//number of edges
	edge *el;//edge list
	unsigned long *dout;//outdegree
} sparse;

//compute the maximum of three unsigned long long
unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the weighted sparsematrix
sparse* readedgelist(char* edgelist){
	unsigned long long i, e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	g->el=malloc(e1*sizeof(edge));
	FILE *file;
	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	while (fscanf(file,"%lu %lu\n", &(g->el[g->e].s), &(g->el[g->e].t))==2) {
		g->n=max3(g->n,g->el[g->e].s,g->el[g->e].t);
		if (++(g->e)==e1) {
			e1+=NLINKS;
			g->el=realloc(g->el,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->el=realloc(g->el,g->e*sizeof(edge));

	//computing out-degrees
	g->dout=calloc(g->n,sizeof(unsigned long));
	for (i=0;i<g->e;i++){
		g->dout[g->el[i].s]++;
	}

	return g;
}

//free the graph stucture
void freegraph(sparse *g){
	free(g->el);
	free(g->dout);
	free(g);
}

//one step of random walk on the graph: initial probalities in v1 and result stored in v2
void onestep(sparse* g, double* v1, double* v2){
	unsigned long long i;
	bzero(v2,sizeof(double)*g->n);
	for (i=0;i<g->e;i++){
		v2[g->el[i].t]+=v1[g->el[i].s]/((double)(g->dout[g->el[i].s]));//
	}
}

//returns an approximation of the pagerank (restar probability=a, number of iterations=k)
double* pagerank(sparse* g, unsigned long source, double alpha, double eps){
	double *v1,*v2,*v3;
	double s;
	unsigned i;
	unsigned long j,n=g->n;
	double prec=1.;
	unsigned it=0;

	v1=calloc(n,sizeof(double));

	v1[source]=1.;
	v2=malloc(n*sizeof(double));

	while(prec>eps){
		it++;
		prec=0;
		onestep(g,v1,v2);//one step of random walk stored in v2
		for (j=0;j<g->n;j++){
			v2[j]*=(1.-alpha);
		}
		v2[source]+=alpha;
		for (j=0;j<g->n;j++){
			prec += (v1[j]>v2[j]) ? v1[j]-v2[j] : v2[j]-v1[j];
		}
		v3=v2,v2=v1,v1=v3;
	}

	printf("Number of iterations: %u\n",it);

	free(v2);
	return v1;
}

void printres(FILE* file,unsigned long n, double* vect){
	unsigned long i,n0=0;
	double sum=0;
	for (i=0;i<n;i++){
		sum+=vect[i];
		if (vect[i]>0){
			fprintf(file,"%lu %le\n",i,vect[i]);
			n0++;
		}
	}
	printf("Number of nonzero entries: %lu\n",n0);
	printf("Sum of values: %le\n",sum);
}

int main(int argc,char** argv){
	unsigned long source;
	double eps;
	sparse* g;
	double* pr;
	FILE* file;

	time_t t0,t1,t2;

	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[1]);

	g=readedgelist(argv[1]);

	printf("Number of nodes = %lu\n",g->n);
	printf("Number of edges = %llu\n",g->e);

	source=atof(argv[2]);
	printf("source node = %lu\n",source);

	eps=atof(argv[3]);
	printf("precision = %le\n",eps);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing approximation of pagerank\n");

	pr=pagerank(g,source,ALPHA,eps);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Printing results to file %s\n",argv[4]);
	file=fopen(argv[4],"w");
	printres(file,g->n,pr);
	fclose(file);

	freegraph(g);

	t2=time(NULL);
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}

