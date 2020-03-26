/*
 * Tools
 *
 * useful methods
 */

#ifndef TOOLS_H
#define TOOLS_H

#include <cstdint>
#include <vector>
#include <iostream>

namespace tools_for_dijkstra {

  // This method helps using priority_queue as a min queue !
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  struct CompareQueueNode {
    bool operator()(std::pair<TI,TV> const &lhs, std::pair<TI,TV> const &rhs) {
      return lhs.second > rhs.second;
    }
  };

}


namespace tools_for_kssp {

  // Storage of candidate path
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class CandidatePath
  {
  public:
    std::vector<TI> path;
    TV weight;
    size_t dev_idx;

    CandidatePath() = default;

    // Constructor
    explicit CandidatePath(std::vector<TI> path, const TV weight, const size_t deviation_index);

    // Destructor
    virtual ~CandidatePath();
  };

  // Constructor
  template<typename TI, typename TV>
  CandidatePath<TI, TV>::CandidatePath(std::vector<TI> ppath, const TV weight, const size_t deviation_index):
    path(ppath), weight(weight), dev_idx(deviation_index)
  {
    /*
    path.clear();
    for (TI u: ppath)
      path.push_back(u);
    */
  }

  // Destructor
  template<typename TI, typename TV>
  CandidatePath<TI, TV>::~CandidatePath()
  {
  }

  // if used with priority_queue
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  struct CompareCandidatePath {
    bool operator()(CandidatePath<TI, TV> const &lhs, CandidatePath<TI, TV> const &rhs) {
      return lhs.weight > rhs.weight;
    }
  };


  // ---------------------------------------------------------------------------
  // Data structure to store candidate paths
  // Used by SidetrackBased and ParsimoniousSidetrackBased
  // ---------------------------------------------------------------------------
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class CandidateSimple
  {
  public:
    TV weight;
    std::size_t dev_idx;
    std::vector<std::pair<TI,TI> > sidetrack_sequence;
    std::vector<std::size_t> tree_index_sequence;

    CandidateSimple() = default;

    // Constructor
    explicit CandidateSimple(const TV weight,
			     const std::size_t deviation_index,
			     std::vector<std::pair<TI,TI> > sidetrack_seq,
			     std::vector<std::size_t> tree_index_seq);

    // Copy constructor
    CandidateSimple(CandidateSimple const *other);

    // Destructor
    virtual ~CandidateSimple();
  };

  // Constructor
  template<typename TI, typename TV>
  CandidateSimple<TI, TV>::CandidateSimple(const TV weight,
					   const std::size_t deviation_index,
					   std::vector<std::pair<TI,TI> > sidetrack_seq,
					   std::vector<std::size_t> tree_index_seq):
    weight(weight), dev_idx(deviation_index),
    sidetrack_sequence(sidetrack_seq), tree_index_sequence(tree_index_seq)
  {
  }

  
  // Copy constructor
  template<typename TI, typename TV>
  CandidateSimple<TI, TV>::CandidateSimple(CandidateSimple<TI, TV> const *other):
    weight(other->weight), dev_idx(other->dev_idx),
    sidetrack_sequence(other->sidetrack_sequence),
    tree_index_sequence(other->tree_index_sequence)
  {
  }


  // Destructor
  template<typename TI, typename TV>
  CandidateSimple<TI, TV>::~CandidateSimple()
  {
  }

  
  // if used with priority_queue
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  struct CompareCandidateSimple {
    bool operator()(std::pair<TV, CandidateSimple<TI,TV> *> &lhs, std::pair<TV, CandidateSimple<TI,TV> *> &rhs)
    {
      return lhs.first > rhs.first;
    }
  };


  // ---------------------------------------------------------------------------
  // Data structure to store non-simple candidate paths
  // Used by ParsimoniousSidetrackBased
  // ---------------------------------------------------------------------------
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class CandidateHard
  {
  public:
    TV weight;
    std::size_t idx_min;
    CandidateSimple<TI, TV> candidate;
    std::pair<TI,TI> sidetrack;  // used for Kurz and Mutzel
    std::vector<std::pair<std::size_t,TI> > sidetracks; // used for Al Zoobi and Coudert
    std::size_t tree_index;
    bool tree_used;

    CandidateHard() = default;

    // Constructor for SidetrackBased
    explicit CandidateHard(const TV weight,
			   const std::size_t idx,
			   CandidateSimple<TI, TV> candidate,
			   std::pair<std::size_t,TI> bad_sidetrack,
			   std::size_t tree_idx);

    // Constructor for ParsimoniousSidetrackBased
    explicit CandidateHard(const TV weight,
			   const std::size_t idx,
			   CandidateSimple<TI, TV> candidate,
			   std::vector<std::pair<std::size_t,TI> > bad_sidetracks,
			   std::size_t tree_idx);

    // Destructor
    virtual ~CandidateHard();
  };

  // Constructor for SidetrackBased
  template<typename TI, typename TV>
  CandidateHard<TI, TV>::CandidateHard(const TV weight,
				       const std::size_t idx,
				       CandidateSimple<TI, TV> candidate_path,
				       std::pair<std::size_t,TI> bad_sidetrack,
				       std::size_t tree_idx):
    weight(weight), idx_min(idx), candidate(candidate_path), sidetrack(bad_sidetrack),
    tree_index(tree_idx), tree_used(false)
  {
  }

  // Constructor for ParsimoniousSidetrackBased
  template<typename TI, typename TV>
  CandidateHard<TI, TV>::CandidateHard(const TV weight,
				       const std::size_t idx,
				       CandidateSimple<TI, TV> candidate_path,
				       std::vector<std::pair<std::size_t,TI> > bad_sidetracks,
				       std::size_t tree_idx):
    weight(weight), idx_min(idx), candidate(candidate_path), sidetracks(bad_sidetracks),
    tree_index(tree_idx), tree_used(false)
  {
  }

  // Destructor
  template<typename TI, typename TV>
  CandidateHard<TI, TV>::~CandidateHard()
  {
  }


  // if used with priority_queue
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  struct CompareCandidateHard {
    bool operator()(std::pair<TV, CandidateHard<TI,TV> *> &lhs, std::pair<TV, CandidateHard<TI,TV> *> &rhs)
    {
      return lhs.first > rhs.first;
    }
  };

}
#endif
