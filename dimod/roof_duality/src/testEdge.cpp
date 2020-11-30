
template <typename capacity_t> class ImplicationEdge {
public:
  typedef capacity_t capacity_type;

  ImplicationEdge(int from_vertex, int to_vertex, capacity_t capacity,
                  capacity_t reverse_capacity)
      : from_vertex(from_vertex), to_vertex(to_vertex), residual(capacity) {
    assert((!capacity || !reverse_capacity) &&
           "Either capacity or reverse edge capacity must be zero.");
    _encoded_capacity = (!capacity) ? -reverse_capacity : capacity;
  }

  // This empty constructor must be inlined for having a low overhead for std::vector resize that will be used in the parallel version of implicit network creation. We do not want anything to be done during the call to resize for the vectors of ImplicatioEdge, which will call the constructor of ImplicationEdge.
  inline ImplicationEdge() __attribute__((always_inline)) { };

  inline fillData(int from_vertex_, int to_vertex_, capacity_t capacity_,
                  capacity_t reverse_capacity_, int reverse_edge_index_, int symmetric_edge_index_)  __attribute__((always_inline))
  {
    assert((!capacity || !reverse_capacity) &&
           "Either capacity or reverse edge capacity must be zero.");
    from_vertex = from_vertex_;
    to_vertex = to_vertex_;
    residual = capacity_;
    _encoded_capacity = (!capacity_) ? -reverse_capacity_ : capacity_;
    reverse_edge_index = reverse_edge_index_;
    symmetric_edge_index = symmetric_edge_index_;
  }


  void print() {
    std::cout << std::endl;
    std::cout << from_vertex << " --> " << to_vertex << std::endl;
    std::cout << "Capacity : " << getCapacity() << std::endl;
    std::cout << "Residual : " << residual << std::endl;
    std::cout << "Reverse Edge Capaciy : " << getReverseEdgeCapacity()
              << std::endl;
    std::cout << "Reverse Edge Residual : " << getReverseEdgeResidual()
              << std::endl;
  }

  // The from_vertex is redundant, since we use an adjacency list, but we are
  // not wasting any space as of now, since the compiler uses the same amount of
  // storage for padding when residual is of a type that takes 8 bytes and in
  // our implmementation we use long long int (8 byte type).
  int from_vertex;
  int to_vertex;
  int reverse_edge_index;
  int symmetric_edge_index;
  capacity_t residual;

private:
  // The value of capacity is not needed per se to compute a max-flow but is
  // needed for the purpose of verifying it. We use the encoded capacity to
  // store the capacity of an edge, but when it is a residual/reverse edge,
  // instead of saving 0 we save the negative of the original edge's capacity,
  // this way we can both verify max-flow and also return the residual capacity
  // of the edge itself and also its reverse/residual without hopping through
  // memory. This is not the best software engineering practice, but is needed
  // to save memory. Current size is 32 byte when long long int is used for
  // capacity and 2 such edges fit in a typical cache line, adding 8 bytes will
  // not allow that to be possible.
  capacity_t _encoded_capacity;

public:
  inline capacity_t getFlow() {
    return ((_encoded_capacity > 0) ? (_encoded_capacity - residual)
                                    : -residual);
  }
  inline capacity_t getCapacity() {
    return ((_encoded_capacity > 0) ? _encoded_capacity : 0);
  }

  inline capacity_t getReverseEdgeCapacity() {
    return ((_encoded_capacity > 0) ? 0 : -_encoded_capacity);
  }

  inline capacity_t getReverseEdgeResidual() {
    return ((_encoded_capacity > 0) ? (_encoded_capacity - residual)
                                    : (-_encoded_capacity - residual));
  }

  // Needed for the purpose of making residual network symmetric.
  void scaleCapacity(int scale) { _encoded_capacity *= scale; }
};
