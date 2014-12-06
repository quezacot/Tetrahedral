#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <vector>
using std::vector;
#include <cassert>

#include <iostream>
using std::cout;
using std::endl;

#include "CS207/Util.hpp"
#include "CS207/Point.hpp"

//constant length in problem 1
#define LENGTH

//Debug codes turn off
//#define Dbg_graph
//#define Dbg_edge_graph
#define Dbg
//#define Dbg_edge_it


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:
  struct internal_node;
  struct internal_edge;


 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;
  // Type of user defined node.
  typedef V node_value_type;
  // Type of user defined edge.
  typedef E edge_value_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph() {
    num_edge = 0;
  }
  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    num_edge = 0;
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position modifiable.*/
    Point& position(){
      return graph_->nodes[idx_].p;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[idx_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ( graph_ == n.graph_ && idx_ == n.idx_ );
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if( *this == n )
        return false;
      if( graph_ == n.graph_ ){// nodes in the same Graph
        if( idx_ < n.idx_ ){
          return true;
        }
      }else{ // nodes not in the same Graph
        if( graph_ < n.graph_ ){
          return true;
        }
      }
      return false;
    }

#ifdef Dbg_graph
    void printnode(){
      const Point& point = graph_->nodes[idx_].p;
      cout << "x:" << point.x << " y:" << point.y << " z:" << point.z<< endl;//" value: "<<graph_->nodes[idx_].node_value;
    }
#endif
#ifdef Dbg
    friend std::ostream& operator<< (std::ostream& stream, const Node& n){
      const Point& point = n.position();
      const internal_node& ina = n.graph_->nodes[n.index()];
      stream << "["<<n.idx_<<"] x:" << point.x << " y:" << point.y << " z:" << point.z<< endl;//" value: "<<n.value() << " number of link edge:" << ina.num_l_edge << endl;
      return stream;
    }
#endif
    /** Obtain the user defined type V stored in this node.
     * @return the @a node_value as a reference
     *
     * Complexity: O(1).
     */
    node_value_type& value(){
      return graph_->nodes[idx_].node_value;
    }
    /** Obtain the user defined type V stored in this node.
     * @return the @a node_value as a const reference
     *
     * Complexity: O(1).
     */
    const node_value_type& value() const {
      return graph_->nodes[idx_].node_value;
    }

    /** Obtain the number of incidents edges connected to this node
     * @return the number of edges connected to this node
     * s.t. 0 <= degree() < num_nodes()
     *
     * Complexity: O(1).
     */
    size_type degree() const {
      return graph_->nodes[idx_].link_edge.size();
    }

    /** Obtain an incident iterator pointing to the first incident edge of this node.
     * @return an incident_iterator pointing to the first edge connecting to this node.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const {
      return IncidentIterator( graph_, idx_, 0 );
    }
    /** Obtain an incident iterator represents the end of incident iterator.
     * @return an incident_iterator with idx_ equal to the number of incident edges.
     *  If there is no edge connecting to this node, this iterator will equal to the begin iterator.
     *  i.e. @a graph_->nodes[idx_].link_edge.size() == 0
     *
     *  Complexity: O(1).
     */
    incident_iterator edge_end() const {
      return  IncidentIterator( graph_, idx_, graph_->nodes[idx_].link_edge.size() );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type idx_;
    Node( const Graph* g, size_type index )
        : graph_( const_cast<Graph*>(g) ), idx_(index) {
    }
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& v = node_value_type ()) {
    size_type added_idx = nodes.size();
    nodes.push_back(internal_node(position,v));
    return Node(this, added_idx);
  }


  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if( this == n.graph_ && nodes.size() > n.index() )
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }

  /** Remove a node in the graph.
   * @param[in] Node The node we want to remove
   * @return the index of removed node if @a n was in this graph and removed.
   *         return the size() if @a n is not in this graph;
   * @post new size() == old size() - 1 if @a n is removed.
   *       new size() == old size() if @a n is not in this graph.
   *       all former created Node objects and node iterators may be invalidated after a node is removed.
   *
   * Complexity: O(num_nodes()^2).
   */
  size_type remove_node(const Node& n){
    if( !has_node(n) )
      return size();
    num_edge -= nodes[n.index()].link_edge.size();
    // Remove the edges in other nodes' adjacent list
    for( auto eit = n.edge_begin(); eit != n.edge_end(); ++eit ){
      Edge e((*eit).graph_, (*eit).idx2_, (*eit).idx1_ );
      remove_oneedge(e);
    }

    // Preserve the index to modify
    // the node with index pre_index is changed to post_index
    size_type pre_index = nodes.size() - 1;
    size_type post_index = n.index();

    // Remove this node
    if( vector_swap_erase(nodes, post_index) ){
      // Relabel the node index that changed due to remove (pre_index => post_index)
      for( IncidentIterator eit = Node(this, post_index).edge_begin(); eit != Node(this, post_index).edge_end(); ++eit ){
        Edge e((*eit).graph_, (*eit).idx2_, pre_index );
        // idx1 > idx2 => idx1 < idx2
        if( post_index < (*eit).idx2_ ){
          // move edge_value
          nodes[post_index].link_edge[(*eit).pos2_].edge_value = nodes[e.idx1_].link_edge[e.pos2_].edge_value;
          // change num_l_edge
          ++nodes[post_index].num_l_edge;
          --nodes[e.idx1_].num_l_edge;
        }
        // label the new index (pre_index => post_index)
        nodes[e.idx1_].link_edge[e.pos2_].node2 = post_index;
      }
    }
    /* node erase shift down version.
    nodes.erase(nodes.begin() + n.index());
    for( unsigned index1 = 0; index1 < n.graph_->size(); ++index1){
      for( unsigned index2 = 0; index2 < n.graph_->nodes[index1].link_edge.size(); ++index2){
        internal_edge& edg = n.graph_->nodes[index1].link_edge[index2];
        if( edg.node2 > n.index() )
          edg.node2 -= 1;
      }
    }*/
    return post_index;
  }

  /** Remove a node in the graph.
   * @param[in] NodeIterator The iterator pointing to the node we want to remove.
   * @pre @a n_it can be dereferenced.
   * @return the iterator pointing to the next node of removed node if @a *n_it was in graph and removed.
   *         return end() if @a n_it is not pointing to a node in this graph;
   * @post new size() == old size() - 1 if @a *n_it is removed.
   *       new size() == old size() if @a *n_it is not in this graph.
   *       all former created Node objects and node iterators may be invalidated after a node is removed.
   *
   * Complexity: O(num_nodes()^2).
   */
  node_iterator remove_node(node_iterator n_it){
    size_type index = remove_node(*n_it);
    return NodeIterator( (*n_it).graph_, index );
  }

    /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node( graph_, idx1_ );
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node( graph_, idx2_ );
    }
#ifdef LENGTH
    /** Initial lenght*/
    double length() const {
      return graph_->nodes[idx1_].link_edge[pos2_].length;
    }
#endif
    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& e) const {
      if( graph_ == e.graph_ &&
          node1() == e.node1() && node2() == e.node2() )
        return true;
      else
        return false;
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& e) const {
      if( *this == e )
        return false;
      if( graph_ == e.graph_ ){// edges are in the same Graph
        if( idx1_< e.idx1_ ){
          return true;
        }else if ( idx1_ == e.idx1_ && idx2_ < e.idx2_ ){
          return true;
        }
      }else{ // edges are not in the same Graph
        if( graph_ < e.graph_ ){
          return true;
        }
      }
      return false;
    }

    /** Obtain the user defined type E stored in this edge.
     * @return the @a edge_value as a reference
     *
     * Complexity: O(num_nodes()).
     */
    edge_value_type& value(){
      if( idx2_ < idx1_ )
        return Edge(graph_, idx2_, idx1_ ).value();
      return graph_->nodes[idx1_].link_edge[pos2_].edge_value;
    }
    /** Obtain the user defined type E stored in this edge.
     * @return the @a edge_value as a const reference
     *
     * Complexity: O(num_nodes()).
     */
    const edge_value_type& value() const {
      if( idx2_ < idx1_ )
        return Edge(graph_, idx2_, idx1_ ).value();
      return graph_->nodes[idx1_].link_edge[pos2_].edge_value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type idx1_, idx2_, pos2_;
    Edge( const Graph* g, size_type index1, size_type index2 )
        : graph_( const_cast<Graph*>(g) ), idx1_(index1), idx2_(index2) {
      assert(idx1_<graph_->nodes.size());
      for(unsigned i = 0; i < graph_->nodes[idx1_].link_edge.size(); ++i)
        if(graph_->nodes[idx1_].link_edge[i].node2 == idx2_)
          pos2_ = i;
    }
    Edge( const Graph* g, size_type index1, size_type pos2, bool indexing )
        : graph_( const_cast<Graph*>(g) ), idx1_(index1), pos2_(pos2) {
      if(indexing)
        idx2_ = graph_->nodes[idx1_].link_edge[pos2].node2;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
#ifdef Dbg
    size_type total_edges = 0;
    for(unsigned i = 0; i < nodes.size(); ++i)
      total_edges += nodes[i].num_l_edge;
    assert(total_edges == num_edge);
#endif
    return num_edge;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    size_type idx1 = a.index();
    size_type idx2 = b.index();

    assert( idx1 < nodes.size() && idx2 < nodes.size() );

    if( has_edge( a, b ) )
      return Edge( this, a.index(), b.index() );
#ifdef Dbg_edge_graph
    if( has_edge( a, b ) )
      cout<< "Edge repeated, not added! " << "node a: " << Node( this, a.index() ) << endl << " node b: " << Node( this, b.index() ) << endl;
#endif
    internal_edge ab, ba;
    ab.node2 = b.index();
    ba.node2 = a.index();
#ifdef LENGTH
    ab.length = ba.length = norm(a.position() - b.position());
#endif
    nodes[a.index()].link_edge.push_back(ab);
    nodes[b.index()].link_edge.push_back(ba);

    if( a.index() < b.index() )
      ++nodes[a.index()].num_l_edge;
    else
      ++nodes[b.index()].num_l_edge;

    ++num_edge;
#ifdef Dbg_edge_graph
    cout<< "Edge added! " << "node a: " << Node( this, a.index() ) << endl << " node b: " << Node( this, b.index() ) << endl;
#endif
    return Edge( this, a.index(), b.index() );
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // assert( a != b );
    if( a == b )
      return false;
    // Search the edge in the smaller @a list_edge
    //
    const vector<internal_edge>& edgea = nodes[a.index()].link_edge;
    const vector<internal_edge>& edgeb = nodes[b.index()].link_edge;
    if( edgea.size() < edgeb.size() ){
      for( unsigned i = 0; i < edgea.size(); ++i ){
        if( edgea[i].node2 == b.index() )
          return true;
      }
    }else{
      for( unsigned i = 0; i < edgeb.size(); ++i ){
        if( edgeb[i].node2 == a.index() )
          return true;
      }
    }
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type index) const {
    assert( index < num_edge );
    unsigned a,b;
#ifdef Dbg_edge_graph
    cout<<"look for edge "<<index<<endl;
#endif
    for( unsigned i = 0, cidx = 0; cidx <= index && i < nodes.size(); ++i){
      const internal_node& ina = nodes[i];
      if( cidx + ina.num_l_edge < index ){ // the edge is not in this link_edge, skip to next
        cidx += ina.num_l_edge;
      }else{
        const vector<internal_edge>& edgea = nodes[i].link_edge;
        for( unsigned j = 0; cidx <= index && j < edgea.size(); ++j ){
          if( i < edgea[j].node2 ){ // skip the edge(i,j) if j>i
            ++cidx;
            a = i;
            b = edgea[j].node2;
          }
        }
      }
    }
#ifdef Dbg_edge_graph
    cout<<"found edge( "<<a<<", "<<b<<" )"<<endl;
#endif
    return Edge( this, a, b );
  }

  /** Remove an edge in the graph.
   * @param[in] Edge The Edge we want to remove
   * @return 0 if @a e is removed.
   *         return num_edge() if the edge is not in this graph.
   * @post new num_edge() == old num_edge() - 1 if edge was in this graph and removed.
   *       new num_edge() == old num_edge() if @a e is not in this graph.
   *       all former created Edge objects and edge iterators may be invalidated after an edge is removed.
   *
   * Complexity: O(num_nodes()).
   */
  size_type remove_edge(const Edge& e){
    if( !e.graph_->has_edge(Node(e.graph_, e.node1().index()), Node(e.graph_, e.node2().index() )) )
      return num_edge;
    --num_edge;
    remove_oneedge(e);
    remove_oneedge( Edge(e.graph_, e.node2().index(), e.node1().index()) );
    return 0;
  }

  /** Remove an edge in the graph.
   * @param[in] Node The nodes @a a and @a b connecting the edge we want to remove
   * @pre @a a and @a b are in the same graph.
   * @return 0 if the edge(@a a,@a b) is removed.
   *         return num_edge() if the edge is not in the graph.
   * @post new num_edge() == old num_edge() - 1 if edge(@a a,@a b) was in this graph and removed.
   *       new num_edge() == old num_edge() if the edge(@a a,@a b) is not in this graph.
   *       all former created Edge objects and edge iterators may be invalidated after an edge is removed.
   *
   * Complexity: O(num_nodes()).
   */
  size_type remove_edge(const Node& a, const Node& b){
    assert(a.graph_ == b.graph_);
    return remove_edge( Edge(a.graph_, a.index(), b.index()) );
  }

  /** Remove an edge in the graph.
   * @param[in] EdgeIterator The iterator @a e_it pointing to the edge we want to remove
   * @pre e_it can be dereferenced.
   * @return the next edge of @a e_it if *e_it is removed.
   *         return @a e_it if the edge is not in the graph.
   * @post new num_edge() == old num_edge() - 1 if *e_it was in this graph and removed.
   *       new num_edge() == old num_edge() if *e_it is not in this graph.
   *       all former created Edge objects and edge iterators may be invalidated after an edge is removed.
   *
   * Complexity: O(num_nodes()).
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    if( e_it.node2_pos_ == e_it.graph_->nodes[e_it.node1_idx_].link_edge.size() )
      ++e_it;
    else if( e_it.node2_pos_ < e_it.graph_->nodes[e_it.node1_idx_].link_edge.size() &&
             e_it.node1_idx_ > e_it.graph_->nodes[e_it.node1_idx_].link_edge[e_it.node2_pos_].node2 )
      ++e_it;
    return e_it;
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator  : private equality_comparable<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Obtain the abstract node this iterator pointing.
     * @pre @a node_idx_ < num_nodes()
     * @return Node with this graph's pointer and index of this node
     *
     * Complexity: O(1).
     */
    Node operator*() const {
      return Node( graph_, node_idx_ );
    }
    /**Increment NodeIterator and return the next position.
     * @post the @node_idx_ increase by 1, may point to an invalid position.
     * @return the modified NodeIterator.
     *
     * Complexity: O(1).
     */
    NodeIterator& operator++() {
      ++node_idx_;
      return *this;
    }
    /** Test the equality of NodeIterator:
     * @return true if the two NodeIterators belong to the same graph, same node and pointing to the same index.
     *
     */
    bool operator==(const NodeIterator& target) const {
      return ( graph_ == target.graph_&& node_idx_ == target.node_idx_ );
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_idx_;
    NodeIterator( const Graph* g, size_type index )
        : graph_( const_cast<Graph*>(g) ), node_idx_(index) {
    }
  };

  /** Obtain a node_iterator pointing to the start of the graph's nodes.
   * @return a node_iterator at the beginning position of the graph's nodes, it could be invalid if there is no node in the graph.
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator( this, 0 );
  }
  /** Obtain a node_iterator representing the end of the graph's nodes.
   * @return a node_iterator with index = @a nodes.size()
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator( this, nodes.size() );
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereference the edge iterator
     * @pre @a node1_idx_ < @a nodes.size()
     * @pre @a node2_pos_ < @a nodes[node1_idx_].link_edge.size()
     * @return Edge connecting the nodes index @a node1_idx_ and @a nodes[node1_idx_].link_edge[node2_pos_]
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return Edge( graph_, node1_idx_, node2_pos_, true );
    }
    /** Increase the edge iterator
     * Increase the iterator to the next edge in @a link_edge of node(@a node1_idx).
     * If it is the last one, then point to the first edge of next node.
     * To deal with the duplicated edges stored in the @a link_edges. (Both edge(i,j) and edge(j,i) are stored)
     * Only count the edge(i,j) if i < j, skip the ones that j < i.
     * @return the EdgeIterator advanced to next position, or be end.
     *
     * Complexity: O(1).
     */
    EdgeIterator& operator++() {
#ifdef Dbg_edge_it
      cout<<"++0 node1: "<<node1_idx_<<"   node2: "<<graph_->nodes[node1_idx_].link_edge[node2_pos_].node2<<endl;
#endif
      const vector<internal_edge>& edgev = graph_->nodes[node1_idx_].link_edge;
      // First, traverse the edges in the @a link_edge of current node.
      //
      ++node2_pos_;
      for( ; node2_pos_ < edgev.size(); ++node2_pos_ ){
        if( node1_idx_ < edgev[node2_pos_].node2 ){ // skip the edge(i,j) if j>i
          return *this;// found next one
        }
      }
      // Then traverse to the @a link_edge of next node if the edge did not found in the @a link_edge of previous node.
      //
      ++node1_idx_;
      for( ; node1_idx_ < graph_->nodes.size(); ++node1_idx_){
        const vector<internal_edge>& edgea = graph_->nodes[node1_idx_].link_edge;
        for( node2_pos_ = 0; node2_pos_ < edgea.size(); ++node2_pos_ ){
#ifdef Dbg_edge_it
          cout<<"++1 node1: "<<node1_idx_<<"   node2: "<<graph_->nodes[node1_idx_].link_edge[node2_pos_].node2<<endl;
#endif
          if( node1_idx_ < edgea[node2_pos_].node2 ){ // skip the edge(i,j) if j>i
            return *this;// found next one
          }
        }
      }
#ifdef Dbg_edge_it
      cout<<"++2 node1: "<<node1_idx_<<endl;
#endif
      // When @a node1_idx increase to @a graph_->nodes.size(), it is essentially at the end no matter what value @a node2_pos_ is.
      // Set @a node2_pos_ to 0 for consistency, so it will equal to the end iterator.
      //
      if( node1_idx_ == graph_->nodes.size() )
        node2_pos_ = 0; // set to end
      return *this; // node1_idx_ = nodes.size()
    }
    /** Test the equality of EdgeIterator.
     * @param[in] target EdgeIterator
     * @return True if both EdgeIterators are in the same graph and connecting the same two nodes.
     */
    bool operator==(const EdgeIterator& target ) const {
      return ( graph_ == target.graph_ &&  node1_idx_ == target.node1_idx_ && node2_pos_ == target.node2_pos_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_idx_;
    size_type node2_pos_;
    EdgeIterator( const Graph* g, size_type index1, size_type pos2 )
        : graph_( const_cast<Graph*>(g) ), node1_idx_(index1), node2_pos_(pos2) {
    }
  };

  /** Obtain the begin iterator of edge iterator
   * @return the first edge in the first non-empty link_edge.
   * The begin iterator will equal to end iterator if there is no any edge in this graph.
   *
   * Complexity: O(num_nodes()).
   */
  edge_iterator edge_begin() const {
#ifdef Dbg_edge_it
    cout<<"edge_begin"<<endl;
#endif
    // no nodes at all
    if( nodes.empty() ){
      return EdgeIterator( this, 0, 0 );
    }
    // find the first edge
    for( size_type idx = 0; idx < nodes.size(); ++idx ){
      const vector<internal_edge>& edgea = nodes[idx].link_edge;
      for( size_type pos = 0; pos < edgea.size(); ++pos ){
        if( idx < edgea[pos].node2 ){ // skip the edge(i,j) if j>i
          return EdgeIterator( this, idx, pos );
        }
      }
    }
    // No edge was found.
    // Set the begin iterator equals to end iterator
    //
    return EdgeIterator( this, nodes.size(), 0 );
  }
  /** Obtain the end of edge iterator
   * @return the end iterator, which is defined as @a node1_idx_ = @a nodes.size() and @a node2_pos_ = 0.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
#ifdef Dbg_edge_it
    cout<<"edge_end"<<endl;
#endif
    return EdgeIterator( this, nodes.size(), 0 );
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Dereference the incident iterator
     * @pre NodeIterator != node.edge_end().
     * @return the Edge connecting nodes @a node1_idx_ and @a graph_->nodes[node1_idx_].link_edge[node2_pos_]
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return Edge( graph_, node1_idx_, node2_pos_, true );
    }
    /** Increase the incident iterator.
     * @post Increase @a node2_pos_ to the next index in the link_edge of current node.
     * @post @a node2_pos_ will not increase if incident iterator equals to the end iterator.
     * @return the advanced IncidentIterator, may be valid or invalid position.
     *
     * Complexity: O(1).
     */
    IncidentIterator& operator++() {
      if( node2_pos_ < graph_->nodes[node1_idx_].link_edge.size() )
        ++node2_pos_;
      return *this;
    }
    /** Compare the equality of IncidentIterator.
     * @return True if the two IncidentIterators are in the same graph, have the same index of two sides of nodes.
     */
    bool operator==(const IncidentIterator& target) const {
      return ( graph_ == target.graph_ && node1_idx_ == target.node1_idx_ && node2_pos_ == target.node2_pos_ );
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_idx_;
    size_type node2_pos_;
    IncidentIterator( const Graph* g, size_type node1, size_type node2 )
        : graph_( const_cast<Graph*>(g) ), node1_idx_(node1),  node2_pos_(node2){
    }
  };

 private:
  // Internal use. Only remove one edge.
  void remove_oneedge(const Edge& e){
    vector<internal_edge>& edges = e.graph_->nodes[e.idx1_].link_edge;
    if( e.idx1_ < e.idx2_ )
      nodes[e.idx1_].num_l_edge -= 1;
    vector_swap_erase(edges, e.pos2_);
  }
  // Swap the element we want to remove and then pop_back
  template <typename T>
  bool vector_swap_erase(vector<T>& v, const size_type index){
    bool swapped = true;
    if(index == v.size()-1)
      swapped = false;
    else
      v[index] = v.back();
    v.pop_back();
    return swapped;
  }
  // The internal structure for real nodes and edges.
  // @a link_edge indicates the edges connecting to this node.
  // store the node index of the other side of an edge.
  // @a num_l_edge recording the number of edges with larger node index.
  // boost the speed of find an edge by index
  //
  struct internal_node {
    Point p;
    vector<internal_edge> link_edge;
    size_type num_l_edge;
    node_value_type node_value;
    internal_node():num_l_edge(0) {}
    internal_node(Point p, node_value_type node_value) :
      p(p), num_l_edge(0), node_value(node_value) {}
  };

  struct internal_edge {
      size_type node2;
#ifdef LENGTH
      double length;
#endif
      edge_value_type edge_value;
  };

  vector<internal_node> nodes;
  size_type num_edge;


};

#endif
