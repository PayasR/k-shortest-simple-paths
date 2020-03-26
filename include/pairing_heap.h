/*
 * Pairing heap
 *
 * Implements a pairing heap data structure as described in [1]. See also [2]
 * for more details.
 * 
 * This implementation is templated by the type TI of items and the type TV of
 * the value associated with an item. The top of the heap is the item with
 * smallest value, i.e., this is a min heap data structure. The number of items
 * in the heap is not fixed. Furthermore, this implementation supports both
 * decrease key operations. More precisely, it supports the following
 * operations:
 *
 * - empty(): return true if the heap is empty, and false otherwise.
 *
 * - push(item, value): push an item to the heap with specified value.
 *
 * - top(): access the pair (item, value) at the top of the heap, i.e., with
 *   smallest value in time O(1).
 *   This operation assumes that the heap is not empty.
 *
 * - top_item(): access the item at the top of the heap in time O(1)
 *   This operation assumes that the heap is not empty.
 *
 * - top_value(): access the value of the item at the top of the heap in O(1)
 *   This operation assumes that the heap is not empty.
 *
 * - pop(): remove top item from the heap in amortize time O(log(n))
 *
 * - decrease(item, value): change the value associated with the item to the new
 *   specified value in time o(log(n)). The new value must be smaller than the
 *   previous one. Otherwise the structure of the heap is no longer guaranteed.
 *   This operation assumes that the item is already in the heap.
 *
 * - contains(item): check whether specified item is in the heap in time O(1).
 *
 * - value(item): return the value associated with the item in the heap.
 *   This operation assumes that the item is already in the heap.
 *
 * References:
 *
 * [1] M. L. Fredman, R. Sedgewick, D. D. Sleator, and R. E. Tarjan.
 *     "The pairing heap: a new form of self-adjusting heap".
 *     Algorithmica. 1 (1): 111-129, 1986. doi:10.1007/BF01840439.
 *
 * [2] https://en.wikipedia.org/wiki/Pairing_heap
 *
 * Author:
 * - David Coudert <david.coudert@inria.fr>
 *
 */


#ifndef PAIRING_HEAP_H
#define PAIRING_HEAP_H

#include <iostream>
#include <map>
#include <iterator>

namespace heap {

  template<
    typename TI,  // type of items stored in the node
    typename TV   // value associated with the stored item
    >
  struct PairingHeapNode {
    //template<typename, typename> friend class PairingHeap;
    //protected:
    TI item;   // item contained in the node
    TV value;  // value associated with the item

    PairingHeapNode<TI, TV> *prev;  // Previous sibling of the node or parent
    PairingHeapNode<TI, TV> *next;  // Next sibling of the node
    PairingHeapNode<TI, TV> *child; // First child of the node

    explicit PairingHeapNode(const TI &some_item, const TV &some_value):
      item(some_item), value(some_value), prev(nullptr), next(nullptr), child(nullptr) {
    }

    bool operator<(PairingHeapNode const& other) const {
      return value <= other.value;
    }
  };

  
  template<
    typename TI,  // type of items stored in the node
    typename TV   // value associated with the stored item
    >
  class PairingHeap
  {
  public:
    // Constructor
    explicit PairingHeap();

    // Copy constructor
    PairingHeap(PairingHeap<TI, TV> const *other);

    // Destructor
    virtual ~PairingHeap();
    
    // Return true if the heap is empty, else false
    bool empty();

    // Insert an item into the heap with specified value (priority)
    void push(const TI &some_item, const TV &some_value);

    // Return the top pair (item, value) of the heap
    std::pair<TI, TV> top();

    // Return the top item of the heap
    TI top_item();

    // Return the top value of the heap
    TV top_value();

    // Remove the top element from the heap
    void pop();

    // Decrease the value of specified item
    void decrease(const TI &some_item, const TV &new_value);

    // Check if specified item is in the heap
    bool contains(TI const& some_item);
    
    // Return the value associated with the item
    TV value(const TI &some_item);

    // Return the number of items in the heap
    size_t size();

  private:
    // Pointer to the top of the heap
    PairingHeapNode<TI, TV> *root;

    // Map used to access stored items
    std::map<TI, PairingHeapNode<TI, TV> *> nodes;

    // Pair list of heaps and return pointer to the top of resulting heap
    PairingHeapNode<TI, TV> *_pair(PairingHeapNode<TI, TV> *p);

    // Merge 2 heaps and return pointer to the top of resulting heap
    PairingHeapNode<TI, TV> *_merge(PairingHeapNode<TI, TV> *a, PairingHeapNode<TI, TV> *b);

    // Make b a child of a
    static void _link(PairingHeapNode<TI, TV> *a, PairingHeapNode<TI, TV> *b);

    // Remove p from its parent children list
    static void _unlink(PairingHeapNode<TI, TV> *p);
    
  };

  // Constructor
  template<typename TI, typename TV>
  PairingHeap<TI, TV>::PairingHeap():
    root(nullptr)
  {
    nodes.clear();
  }

  
  // Copy constructor
  template<typename TI, typename TV>
  PairingHeap<TI, TV>::PairingHeap(PairingHeap<TI, TV> const *other):
    root(nullptr)
  {
    nodes.clear();
    for (auto const& it:other->nodes)
      push(it.first, it.second->value);
  }

  
  // Destructor
  template<typename TI, typename TV>
  PairingHeap<TI, TV>::~PairingHeap()
  {
    for (auto const& it: nodes)
      delete it.second;
    root = nullptr;
    nodes.clear();
  }

  
  // Return true if the heap is empty, else false
  template<typename TI, typename TV>
  inline bool PairingHeap<TI, TV>::empty()
  {
    return root == nullptr;
  }


  // Insert an item into the heap with specified value (priority)
  template<typename TI, typename TV>
  void PairingHeap<TI, TV>::push(const TI &some_item, const TV &some_value)
  {
    PairingHeapNode<TI, TV> *p = new PairingHeapNode<TI, TV>(some_item, some_value);
    nodes[some_item] = p;
    root = root == nullptr ? p : _merge(root, p);
  }

  
  // Return the top pair (value, item) of the heap
  template<typename TI, typename TV>
  inline std::pair<TI, TV> PairingHeap<TI, TV>::top()
  {
    return std::make_pair(root->item, root->value);
  }
    
  // Return the top item of the heap
  template<typename TI, typename TV>
  inline TI PairingHeap<TI, TV>::top_item()
  {
    return root->item;
  }

  // Return the top value of the heap
  template<typename TI, typename TV>
  inline TV PairingHeap<TI, TV>::top_value()
  {
    return root->value;
  }


  // Remove the top element from the heap
  template<typename TI, typename TV>
  void PairingHeap<TI, TV>::pop()
  {
    PairingHeapNode<TI, TV> *p = root->child;
    nodes.erase(root->item);
    delete root;
    root = _pair(p);
  }


  // Decrease the value of specified item
  // If the item is not in the heap, push it
  template<typename TI, typename TV>
  void PairingHeap<TI, TV>::decrease(const TI &some_item, const TV &new_value)
  {
    if (contains(some_item))
      {
	PairingHeapNode<TI, TV> *p = nodes[some_item];
	p->value = new_value;
	if(p->prev != nullptr)
	  {
	    _unlink(p);
	    root = _merge(root, p);
	  }
      }
    else
      push(some_item, new_value);
  }

  
  // Check if specified item is in the heap
  template<typename TI, typename TV>
  inline bool PairingHeap<TI, TV>::contains(TI const& some_item)
  {
    return nodes.count(some_item) > 0;
  }
    
  // Return the value associated with the item
  template<typename TI, typename TV>
  inline TV PairingHeap<TI, TV>::value(const TI &some_item)
  {
    return nodes[some_item]->value;
  }

  
  // Return the number of items in the heap
  template<typename TI, typename TV>
  inline size_t PairingHeap<TI, TV>::size()
  {
    return nodes.size();
  }


  // Pair list of heaps and return pointer to the top of resulting heap
  template<typename TI, typename TV>
  inline PairingHeapNode<TI, TV> *PairingHeap<TI, TV>::_pair(PairingHeapNode<TI, TV> *p)
  {
    if(p == nullptr)
      {
	return nullptr;
      }

    /*
     * Move toward the end of the list, counting elements along the way.
     * This is done in order to:
     * - know whether the list has odd or even number of nodes
     * - speed up going-back through the list
     */
    size_t children = 1;
    PairingHeapNode<TI, TV> *it = p;
    while(it->next != nullptr)
      {
	it = it->next;
	children++;
      }

    PairingHeapNode<TI, TV> *result;

    if(children % 2 == 1)
      {
	PairingHeapNode<TI, TV> *a = it;
	it = it->prev;
	a->prev = a->next = nullptr;
	result = a;
      }
    else
      {
	PairingHeapNode<TI, TV> *a = it;
	PairingHeapNode<TI, TV> *b = it->prev;
	it = it->prev->prev;
	a->prev = a->next = b->prev = b->next = nullptr;
	result = _merge(a, b);
      }

    for(size_t i = 0; i < (children - 1) / 2; i++)
      {
	PairingHeapNode<TI, TV> *a = it;
	PairingHeapNode<TI, TV> *b = it->prev;
	it = it->prev->prev;
	a->prev = a->next = b->prev = b->next = nullptr;
	result = _merge(_merge(a, b), result);
      }

    return result;
  }

  
  // Merge 2 heaps and return pointer to the top of resulting heap
  template<typename TI, typename TV>
  inline PairingHeapNode<TI, TV> *PairingHeap<TI, TV>::_merge(PairingHeapNode<TI, TV> *a, PairingHeapNode<TI, TV> *b)
  {
    if(*a < *b)  // Use comparison method of PairingHeapNode
      {
	_link(a, b);
	return a;
      }
    else
      {
	_link(b, a);
	return b;
      }
  }

  // Make b a child of a
  template<typename TI, typename TV>
  inline void PairingHeap<TI, TV>::_link(PairingHeapNode<TI, TV> *a, PairingHeapNode<TI, TV> *b)
  {
    if(a->child != nullptr)
      {
	b->next = a->child;
	a->child->prev = b;
      }
    b->prev = a;
    a->child = b;
  }

  // Remove p from its parent children list
  template<typename TI, typename TV>
  inline void PairingHeap<TI, TV>::_unlink(PairingHeapNode<TI, TV> *p)
  {
    if(p->prev->child == p)
      {
	p->prev->child = p->next;
      }
    else
      {
	p->prev->next = p->next;
      }
    if(p->next != nullptr)
      {
	p->next->prev = p->prev;
      }
    p->prev = nullptr;
    p->next = nullptr;
  }

}

#endif
