
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <cstdlib>
#include "../include/easy_digraph.h"

using namespace std::chrono;
using namespace directed_graph;

EasyDirectedGraph<size_t, uint32_t, uint32_t> *load_digraph(std::string filename, bool verbose=true)
{
  // Read graph
  if (verbose)
    {
      std::cout << "# graph (" << filename << ") : ";
      auto start = high_resolution_clock::now();
      EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = new EasyDirectedGraph<size_t, uint32_t, uint32_t>(filename);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start); // can also be cast to microseconds
      std::cout << G->n << " vertices and " << G->m << " edges : " << duration.count() << " ms\n" << std::endl;
      return G;
    }
  else
    return new EasyDirectedGraph<size_t, uint32_t, uint32_t>(filename);
}

/*
 * Build a list a x random pairs of vertices in V
 */
std::vector<std::pair<size_t, size_t> > random_pairs(std::vector<size_t> V, uint32_t x)
{
  size_t n = V.size();
  std::vector<std::pair<size_t, size_t> > pairs;
  pairs.clear();
  std::srand(std::time(nullptr));
  while (pairs.size() < x)
    {
      size_t u = n * (rand() / (RAND_MAX + 1.0));
      size_t v = u;
      while (u == v)
	v = n * (rand() / (RAND_MAX + 1.0));
      pairs.push_back(std::make_pair(V[u], V[v]));
      if (pairs.size() == x)
	{ // remove possible duplicates
	  std::sort(pairs.begin(), pairs.end());
	  auto new_end = std::unique(pairs.begin(), pairs.end());
	  pairs.erase(new_end, pairs.end());
	}
    }
  return pairs;
}


/*
 * Run specified algorithm from source to target.
 * Display timing for given list of k values
 */
void run_algo1(EasyDirectedGraph<size_t, uint32_t, uint32_t> *G,
	      size_t source, size_t target,
	      std::string algorithm, size_t version,
	      std::vector<size_t> k_values)
{
  G->reset_kssp();
  auto start_init = high_resolution_clock::now();
  G->init_kssp(algorithm, source, target, version, false);
  auto stop_init = high_resolution_clock::now();
  auto duration_init = duration_cast<milliseconds>(stop_init - start_init);
  size_t total_duration = duration_init.count();
  
  size_t cpt = 0;
  for (auto const& k: k_values)
    {
      auto start = high_resolution_clock::now();
      while (cpt < k)
	{
	  auto sol = G->next_path();
	  if (sol.second == 0)
	    break;
	  cpt++;
	}
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      total_duration += duration.count();
      
      std::cout << algorithm << (version == 1 ? "" : "*");
      std::cout << "\t" << source << "\t" << target << "\t" << k;
      std::cout << "\t" << cpt;
      std::cout << "\t" << G->used_trees();
      std::cout << "\t" << total_duration;
      std::cout << std::endl;
	    
      if (cpt != k)
	break;
    }
}



/*
 * For the input graph, solves the problem for
 * - SB, SB* and PSB
 * - k = 100, 200, 250, ..., 1000
 * - 1000 randomly selected pairs
 *
 */
void experience1(std::string filename)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"SB", 1}, {"SB", 3}, {"PSB", 3}};
  std::vector<size_t> ks = {100, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000};

  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename);
  
  std::vector<std::pair<size_t, size_t> > pairs = random_pairs(G->int_to_vertex, 1000);
  
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& p: pairs)
    for (auto const& algo: algorithms)
      run_algo1(G, p.first, p.second, algo.first, algo.second, ks);

  delete G;
}


/*
 * For the input graph, solves the problem for
 * - SN, SN* and AC*
 * - k = 100, 200, 250, ..., 1000
 * - 1000 randomly selected pairs
 *
 */
void experience6(std::string filename)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"NC", 1}, {"SB", 1}, {"SB", 3}, {"PSB", 3}};
  std::vector<size_t> ks = {100, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000};

  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename);
  
  std::vector<std::pair<size_t, size_t> > pairs = random_pairs(G->int_to_vertex, 1000);
  
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& p: pairs)
    for (auto const& algo: algorithms)
      run_algo1(G, p.first, p.second, algo.first, algo.second, ks);

  delete G;
}



/*
 * Run specified algorithm from source to target.
 * k = 1000
 */
void run_algo2(EasyDirectedGraph<size_t, uint32_t, uint32_t> *G,
	      size_t source, size_t target,
	       std::string algorithm, size_t version)
{
  G->reset_kssp();
  auto start = high_resolution_clock::now();
  G->init_kssp(algorithm, source, target, version, false);
  size_t cpt = 0;
  while (cpt < 1000)
    {
      auto sol = G->next_path();
      if (sol.second == 0)
	break;
      cpt++;
    }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  size_t total_duration = duration.count();

  std::cout << algorithm << (version == 1 ? "" : "*");
  std::cout << "\t" << source << "\t" << target << "\t" << 1000;
  std::cout << "\t" << cpt;
  std::cout << "\t" << G->used_trees();
  std::cout << "\t" << total_duration;
  std::cout << std::endl;
}

/*
 * For the input graph, solves the problem for
 * - SB, SB* and PSB
 * - all pairs
 * - k = 1000
 *
 */
void experience2(std::string filename)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"SB", 1}, {"SB", 3}, {"PSB", 3}};

  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename);
    
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& source: G->int_to_vertex)
    for (auto const& target: G->int_to_vertex)
      if (source != target)
	for (auto const& algo: algorithms)
	  run_algo2(G, source, target, algo.first, algo.second);

  delete G;
}




/*
 * For the input graph, solves the problem for
 * - SB, SB* and PSB
 * - k = 100, 200, 250, ..., 1000, 5000, 10000
 * - 1000 randomly selected pairs
 *
 */
void experience3(std::string filename)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"SB", 1}, {"SB", 3}, {"PSB", 3}};
  std::vector<size_t> ks = {100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 10000};

  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename);
  
  std::vector<std::pair<size_t, size_t> > pairs = random_pairs(G->int_to_vertex, 10000);
  
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& p: pairs)
    for (auto const& algo: algorithms)
      run_algo1(G, p.first, p.second, algo.first, algo.second, ks);

  delete G;
}








/*
 * For the input graph, solves the problem for
 * - SB, SB* and PSB
 * - k = 100, 200, 250, ..., 1000, 5000, 10000
 * - 1000 randomly selected pairs
 *
 */
void experience4(std::string filename, size_t source, size_t target)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"SB", 1}, {"SB", 3}, {"PSB", 3}};

  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename);
  
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& algo: algorithms)
    {
      G->reset_kssp();
      auto start_init = high_resolution_clock::now();
      G->init_kssp(algo.first, source, target, algo.second, false);
      auto stop_init = high_resolution_clock::now();
      auto duration_init = duration_cast<milliseconds>(stop_init - start_init);
      size_t total_duration = duration_init.count();
  
      size_t cpt = 0;
      bool stop = false;
      while (not stop)
	{
	  auto start = high_resolution_clock::now();
	  for (int i = 0; i < 100; i++)
	    {
	      auto sol = G->next_path();
	      if (sol.second == 0)
		{
		  stop = true;
		  break;
		}
	      cpt++;
	    }
	  auto stop = high_resolution_clock::now();
	  auto duration = duration_cast<milliseconds>(stop - start);
	  total_duration += duration.count();
      
	  std::cout << algo.first << (algo.second == 1 ? "" : "*");
	  std::cout << "\t" << source << "\t" << target << "\t" << cpt;
	  std::cout << "\t" << cpt;
	  std::cout << "\t" << G->used_trees();
	  std::cout << "\t" << total_duration;
	  std::cout << std::endl;
	}
      G->report();
    }

  delete G;
}


/*
 */
void experience5(std::string path, size_t source, size_t target)
{
  std::vector<std::pair<std::string, size_t> > algorithms = {{"SB", 1}, {"SB", 3}, {"PSB", 3}};

  std::cout << "i\tSB\tSB*\tPSB" << std::endl;
  for (int i = 0; i < 2880; i++)
    {
      std::cout << i;
      std::string filename = path + "g" + std::to_string(i) + ".gr";
      EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = load_digraph(filename, false);
      std::vector<int> res;
      res.clear();
      for (auto const& algo: algorithms)
	{
	  G->reset_kssp();
	  G->init_kssp(algo.first, source, target, algo.second, false);
	  size_t cpt = 0;
	  while (true)
	    {
	      auto sol = G->next_path();
	      if (sol.second == 0)
		break;
	      cpt++;
	    }

	  std::cout << "\t" << G->used_trees();
	  res.push_back(G->used_trees());
	}
      delete G;
      if (res[0] < res[1])
	std::cout << "\t<----";
      std::cout << std::endl;
    }
}





int main(int argc, char* argv[])
{
  // if (argc < 2)
  //   print_usage();

  switch (std::stoi(argv[1])) {
  case 1:
    experience1(argv[2]);
    break;
  case 6:
    experience6(argv[2]);
    break;
  case 2:
    experience2(argv[2]);
    break;
  case 3:
    experience3(argv[2]);
    break;
  case 4:
    experience4(argv[2], atoi(argv[3]), atoi(argv[4]));
    break;
  case 5:
    experience5(argv[2], atoi(argv[3]), atoi(argv[4]));
    break;
  default:
    break;
  }
}
