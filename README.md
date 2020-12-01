# push
Push method for the computation of rooted page rank

## Info:
push.c is an implementation of the push method to approximate the rooted pagerank.  
The algorithm is described page 6 here: https://papers-gamma.link/paper/174  
Contrarily to the paper, we do not consider the lazy random walk, but plain random walk.  
We approximately solve $pr(s)=\alpha s + T pr(s)$, where $T_ij=A_ij/dout_j$ is the transition matrix and vector s should be such that $s[u]=1$ if u is the root node and 0 otherwise.  
allpush.c runs the push method for each node in the graph. The overall running time is in $O(m+n/(\epsilon \alpha))$.  
RootedPageRank.c is an implementation of the power iteration method.

## To compile:
- gcc push.c -o push -O9
- gcc allpush.c -o allpush -O9
- gcc RootedPageRank.c -o RootedPageRank -O9


## To execute:

./push net.txt source eps pagerank.txt
- net.txt should contain on each line two unsigned separated by a space: "source target\n" that is the input directed graph.
- find the rooted pagerank of node source
- eps precision
- res.txt will contain an approximation of the pagerank with a restart probability of 0.15. "nodeID PageRankValue\n" on each line (contains only nonzero values).

./allpush net.txt eps pagerank.txt
- net.txt should contain on each line two unsigned separated by a space: "source target\n" that is the input directed graph.
- eps precision
- res.txt will contain an approximation of the pagerank with a restart probability of 0.15: each line coresponds to a node: "k nodeID1 PageRankValue1 nodeID2 PageRankValue2... nodeIDk PageRankValuek" (contains only nonzero values and k is the number of nonzero values).


./RootedPageRank net.txt source eps res.txt
- net.txt should contain on each line two unsigned separated by a space: "source target\n" that is the input directed graph.
- source: the id of the root node: the random walk restarts from that node with probability 0.15
- will print only entries of the pagrank vector that are larger than eps
- res.txt will contain an approximation of the pagerank (30 iterations using the power iteration method). "nodeID PageRankValue\n" on each line.
