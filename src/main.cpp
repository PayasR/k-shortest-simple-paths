
#include <iostream>
#include <chrono>
#include "../include/easy_digraph.h"

using namespace std::chrono;
using namespace directed_graph;


class ParserData {
public:
  std::string filename;
  size_t source;
  size_t target;
  size_t k;
  std::vector<std::pair<std::string, size_t> > algorithms;
  size_t repetitions;
  bool show_weights;
  bool show_paths;

  explicit ParserData(int argc, char* argv[]):
    k(10), repetitions(1), show_weights(false), show_paths(false)
  {
    if (argc < 6)
      {
	usage(argc, argv);
	ERROR("Wrong number of parameters");
      }

    filename = argv[1];
    source = std::stoi(argv[2]);
    target = std::stoi(argv[3]);

    for (int i = 4; i < argc; )
      {
	if (strcmp(argv[i], "-k") == 0)
	  {
	    k = std::stoi(argv[++i]);
	  }
	else if (strcmp(argv[i], "-r") == 0)
	  {
	    repetitions = std::stoi(argv[++i]);
	  }
	else if (strcmp(argv[i], "-a") == 0)
	  {
	    i++;
	    if (strcmp(argv[i], "PSB") == 0)      // Parsimonious Sidetrack Based
	      algorithms.push_back(std::make_pair("PSB", 3));
	    else if (strcmp(argv[i], "SB") == 0)  // Sidetrack Based
	      algorithms.push_back(std::make_pair("SB", 1));
	    else if (strcmp(argv[i], "SB*") == 0) // Sidetrack Based with updates
	      algorithms.push_back(std::make_pair("SB", 3));
	    else if (strcmp(argv[i], "NC") == 0)  // Feng
	      algorithms.push_back(std::make_pair("NC", 1));
	    else if (strcmp(argv[i], "Y") == 0)   //"Yen";
	      algorithms.push_back(std::make_pair("Y", 1));
	    else
	      ERROR("Unkown algorithm " << argv[i]);
	  }
	else if (strcmp(argv[i], "-w") == 0)
	  {
	    show_weights = true;
	  }
	else if (strcmp(argv[i], "-p") == 0)
	  {
	    show_paths = true;
	  }
	else if (strcmp(argv[i], "-h") == 0)
	  {
	    usage(argc, argv);
	    std::exit(0);
	  }
	else
	  ERROR("Unknown option " << argv[i]);
	i++;
      }
  }

  void usage(int argc, char* argv[])
  {
    std::cout << "----\n";
    std::cout << "USAGE: " << argv[0] << " <graph_filename> <source> <target> -k <k> -r <r> -a <algo> -a <algo> ...\n";
    std::cout << "----\noptions:\n";
    std::cout << " -h: print this help and exit\n";
    std::cout << " -k: number of paths to compute (default: 10)\n";
    std::cout << " -r: number of repetitions to average running time (default: 1)\n";
    std::cout << " -a: algorithm to run among\n";
    std::cout << "     Y: basic algorithm as proposed by Yen\n";
    std::cout << "     NC: node classification algorithm proposed by Feng\n";
    std::cout << "     SB: sidetrack based algorithm proposed by Kurz and Mutzel\n";
    std::cout << "     SB*: improvement of SB using shortest path tree updates\n";
    std::cout << "     PSB: parsimonious sidetrack based proposed by Al Zoobi, Coudert and Nisse\n";
    std::cout << " -w: whether to display the weight of computed paths (default: false)\n";
    std::cout << " -p: whether to display computed paths (default: false)\n";
    std::cout << std::endl;
  }
};


void run_algorithm(EasyDirectedGraph<size_t, uint32_t, uint32_t> *G,
		   size_t source, size_t target, size_t k,
		   std::string algorithm, size_t version, ParserData P)
{
  size_t total_duration = 0;
  size_t cpt;

  for (size_t j = 0; j < P.repetitions; j++)
    {
      G->reset_kssp();
      std::vector<uint32_t> res;
      res.clear();
      auto start = high_resolution_clock::now();
      G->init_kssp(algorithm, source, target, version, false);
      for (cpt = 0; cpt < k; cpt++)
	{
	  auto sol = G->next_path();
	  if (sol.second == 0)
	    break;
	  res.push_back(sol.second);
	  if (P.show_paths)
	    {
	      std::cout << sol.second << "\t[" << sol.first[0];
	      for (size_t i = 1; i < sol.first.size(); i++)
		std::cout << ", " << sol.first[i];
	      std::cout << "]" << std::endl;
	    }
	}
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stop - start);
      total_duration += duration.count();

      if (P.show_weights)
	{
	  std::cout << algorithm << ": path lengths = [" << res[0];
	  for (size_t i = 1; i < res.size(); i++)
	    std::cout << ", " << res[i];
	  std::cout << "]" << std::endl;
	}
    }

  std::cout << algorithm << (version == 1 ? "" : "*");
  std::cout << "\t" << source << "\t" << target << "\t" << k;
  std::cout << "\t" << cpt;
  std::cout << "\t" << G->used_trees();
  std::cout << "\t" << total_duration / P.repetitions;
  std::cout << std::endl;
  G->reset_kssp();
}



int main(int argc, char* argv[])
{
  ParserData P = ParserData(argc, argv);

  // Read graph
  std::cout << "# graph (" << P.filename << ") : ";
  auto start = high_resolution_clock::now();
  EasyDirectedGraph<size_t, uint32_t, uint32_t> *G = new EasyDirectedGraph<size_t, uint32_t, uint32_t>(P.filename);
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start); // can also be cast to microseconds
  std::cout << G->n << " vertices and " << G->m << " edges : " << duration.count() << " ms\n" << std::endl;

  // Run algorithms
  std::cout << "Algo\tsource\ttarget\tk\tk*\ttrees\ttime (ms)" <<std::endl;
  for (auto const& algo: P.algorithms)
    run_algorithm(G, P.source, P.target, P.k, algo.first, algo.second, P);

  delete G;
}
