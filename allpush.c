/*
Maximilien Danisch
Novembre 2020
http://bit.ly/danisch
maximilien.danisch@gmail.com

Info:
Implementation of the push method to approximate the rooted pagerank.
Algorithm described page 6 here: http://www.leonidzhukov.net/hse/2015/networks/papers/andersen06localgraph.pdf

to compile:
gcc allpush.c -o allpush -O9

to execute:
./allpush net.txt eps pagerank.txt
- net.txt should contain on each line twos unsigned separated by a space: "source target\n" that is the input directed graph.
- find an approximation of the rooted pagerank for each node
- eps precision
- res.txt will contain an approximation of the pagerank with a restart probability of 0.15. each line coresponds to a node: "nodeID1 PageRankValue1 nodeID2 PageRankValue2 nodeID3 PageRankValue3..." (contains only nonzero values).

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <time.h>

#define NLINKS 10000000 //maximum number of edges, will increase if needed
#define ALPHA 0.15 //restart probability of pagerank

typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

//directed graph datastructure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long long e;//number of edges
	edge *edges;//list of edges
	unsigned long *d;//d[i]=out-degree of node i
	unsigned long long *cd;//cumulative out-degree cd[0]=0 length=n+1
	unsigned long *adj;//concatenated lists of out-neighbors of all nodes
} adjlist;

//compute the maximum of three unsigned long long
unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
adjlist* readedgelist(char* input){
	unsigned long long e1=NLINKS;
	FILE *file=fopen(input,"r");

	adjlist *g=malloc(sizeof(adjlist));
	g->n=0;
	g->e=0;
	g->edges=malloc(e1*sizeof(edge));//allocate some RAM to store edges

	while (fscanf(file,"%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (++(g->e)==e1) {//increase allocated RAM if needed
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);

	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//building the adjacency matrix
void mkadjlist(adjlist* g){
	unsigned long i,s,t;
	unsigned long long j;

	g->d=calloc(g->n,sizeof(unsigned long));

	for (j=0;j<g->e;j++) {
		g->d[g->edges[j].s]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned long long));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+g->d[i-1];
		g->d[i-1]=0;
	}

	g->adj=malloc(g->e*sizeof(unsigned long));

	for (j=0;j<g->e;j++) {
		s=g->edges[j].s;
		t=g->edges[j].t;
		g->adj[ g->cd[s] + g->d[s]++ ]=t;
	}

	free(g->edges);
}

//freeing memory
void free_adjlist(adjlist *g){
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g);
}

//dictionary-like datastruture
typedef struct {
	unsigned long nmax;//maximum number of element in dict
	unsigned long n;//number of elements in dict
	unsigned long* list;//elements in dict
	double* val;//elements in dict
} Dict ;

Dict* allocdict(unsigned long n){
	Dict* dict=malloc(sizeof(Dict));
	dict->nmax=n;
	dict->n=0;
	dict->list=malloc(n*sizeof(unsigned long));
	dict->val=calloc(n,sizeof(double));
	return dict;
}

void cleandict(Dict* dict){
	unsigned long i;
	for (i=0;i<dict->n;i++){
		dict->val[dict->list[i]]=0;
	}
	dict->n=0;
}

void free_dict(Dict *dict){
	free(dict->list);
	free(dict->val);
	free(dict);
}

//The heart of the code
void push(adjlist* g, double alpha, unsigned long source, double eps, Dict* r, Dict* p){
	unsigned long n=0;
	unsigned long long i;
	unsigned long u,v;
	double val;
	static unsigned long *list=NULL;//stores the nodes with r->val[v]>eps*g->d[v]
	if (list==NULL){
		list=malloc(g->n*sizeof(unsigned long));
	}

	r->val[source]=1;
	r->list[(r->n)++]=source;
	if (1.>eps*g->d[source]){
		list[n++]=source;
	}

	while (n>0){
		u=list[--n];
		val=r->val[u];
		r->val[u]=0;

		if (p->val[u]==0){
			p->list[(p->n)++]=u;
		}
		p->val[u]+=alpha*val;
		for (i=g->cd[u];i<g->cd[u+1];i++){
			v=g->adj[i];
			if (r->val[v]==0){//add v to set
				r->list[(r->n)++]=v;
			}
			if (r->val[v]<eps*g->d[v]){
				r->val[v]+=(1-alpha)*val/g->d[u];
				if (r->val[v]>eps*g->d[v]){//add v to list of nodes to consider for push
					list[n++]=v;
				}
			}
			else{
				r->val[v]+=(1-alpha)*val/g->d[u];
			}
		}
	}

}


void print_dict(FILE* file,Dict* dict){
	unsigned long i,u;
	fprintf(file,"%lu",dict->n);
	for (i=0;i<dict->n;i++){
		u=dict->list[i];
		fprintf(file," %lu %le",u,dict->val[u]);
	}
	fprintf(file,"\n");
}

int main(int argc,char** argv){
	unsigned long i;
	double eps;
	adjlist *g;
	Dict *r, *p;
	FILE *file;

	time_t t0,t1,t2;

	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[1]);

	g=readedgelist(argv[1]);

	printf("Number of nodes = %lu\n",g->n);
	printf("Number of edges = %llu\n",g->e);

	mkadjlist(g);

	eps=atof(argv[2]);
	printf("epsilon = %le\n",eps);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing approximation of pagerank\n");
	printf("Printing results to file %s\n",argv[3]);

	r=allocdict(g->n);
	p=allocdict(g->n);

	file=fopen(argv[3],"w");
	for (i=0;i<g->n;i++){
		cleandict(r);
		cleandict(p);
		push(g,ALPHA,i,eps,r,p);
		print_dict(file,p);
	}
	fclose(file);
	free_dict(r);
	free_dict(p);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	free_adjlist(g);


	t2=time(NULL);
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}

