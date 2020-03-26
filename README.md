# k Shortest Simple Paths

Implementation in C++/Cython of algorithms for computing the k shortest simple paths from a source to a destination in a weighted directed graph.

## Authors
- [Ali Al Zoobi](http://www-sop.inria.fr/members/Ali.Al-Zoobi)
- [David Coudert](http://www-sop.inria.fr/members/David.Coudert)
- [Nicolas Nisse](http://www-sop.inria.fr/members/Nicolas.Nisse)


## Installation
Todo


## Usage
Todo

## License
Todo

## Existing implementations:

* In [NetworkX](https://networkx.github.io), method [shortest_simple_paths](https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/simple_paths.html) implements in Python the Yen's algorithm [8]
* In [JGraphT](https://jgrapht.org), method [YenKShortestPath](https://jgrapht.org/javadoc/org/jgrapht/alg/shortestpath/YenKShortestPath.html).
  In [JGraphT](https://jgrapht.org) also offers an implementation of [Eppstein k shortest paths algorithm](https://jgrapht.org/javadoc/org/jgrapht/alg/shortestpath/EppsteinKShortestPath.html) [2]
* Implementation in Python for [NetworkX](https://networkx.github.io) of the algorithm described in [6]: [k-shortest-path](https://github.com/dsaidgovsg/k-shortest-path)
* Implementations in C++, Jaca, Scala, etc. of the algorithm proposed in [6] can be found [here](http://thinkingscale.com/k-shortest-paths-cpp-version/) (warning: [possible memory leak](https://stackoverflow.com/questions/6709066/c-k-shortest-paths-algorithm))


## References
1. D. Ajwani, E. Duriakova, N. Hurley, U. Meyer and A. Schickedanz. An Empirical Comparison of k-Shortest Simple Path Algorithms on Multicores. In Proceedings of the 47th International Conference on Parallel Processing (ICPP 2018), pages 78:1--78:12, ACM, 2018. DOI:10.1145/3225058.3225075
2. D. Eppstein. "Finding the k Shortest Paths" (PDF). SIAM J. Comput. 28 (2): 652–673, 1998. doi:10.1137/S0097539795290477.
3. G. Feng. Finding k shortest simple paths in directed graphs: A node classification algorithm. Networks, 64(1):6–17, 2014. doi:10.1002/net.21552
4. J. Hershberger, M. Maxel, S. Suri. Finding the k shortest simple paths: A new algorithm and its implementation. ACM Transactions on Algorithms (TALG), 3:4.45, November 2007. DOI:10.1145/1290672.1290682
5. D. Kurz and P. Mutzel. A sidetrack-based algorithm for finding the k shortest simple paths in a directed graph. In International Symposium on Algorithms and Computation (ISAAC), volume 64 of LIPIcs, pages 49:1–49:13. Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik, December 2016.doi:10.4230/LIPIcs.ISAAC.2016.49.
6. Q. V. Martins, Ernesto and M. B. Pascoal, Marta. (2003). A new implementation of Yen's ranking loopless paths algorithm. Quarterly Journal of the Belgian, French and Italian Operations Research Societies. 1. 121-133. DOI:10.1007/s10288-002-0010-2
7. J. Y. Yen. Finding the K shortest loopless paths in a network. Management Science, 17:712-716, 1971.
8. J. Y. Yen. Another algorithm for finding the K shortest loopless network paths. In Proc. of 41st Mtg. Operations Research Society of America, volume 20, 1972.
