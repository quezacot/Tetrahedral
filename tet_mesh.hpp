/*
 * CS 207: Final Project Tetrahedral Mesh
 * Yung-jen Cheng HUID #20947802
 * Jeffrey (Shih-kai) Shen HUID #70949288
 */

#include<iostream>
#include<vector>
using std::vector;
#include <cmath>
using std::sqrt;
#include<algorithm>
using std::swap;
using std::sort;
using std::unique;

#include "Graph.hpp"
#include "CS207/Util.hpp"
#include "CS207/Point.hpp"

const unsigned NUM_TET_ADJ_TET = 4;

#pragma once
/** @file tet_mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and tetrahedrals such that:
 *  -- All tetrahedrals have four nodes and six edges.
 *  -- All edges belong to at least one tetrahedral.
 */

/** @class Mesh
 * @brief A template for 3D tetrahedral meshes.
 *
 * Users can add tetrahedrals and retrieve nodes, edges, and tetrahedrals.
 */
template <typename N, typename E, typename T>
class Mesh {

  private:
    struct internal_edge;
    struct internal_tetrahedral;

  public:
    /** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
    typedef N node_value_type;
    typedef E edge_value_type;
    typedef T tet_value_type;
    typedef unsigned size_type;
    class Node;
    class Edge;
    class Tetrahedral;
    typedef Graph<node_value_type, internal_edge> g_real_type;
    typedef Graph<internal_tetrahedral, bool> g_tet_type;

    typedef Node node_type;
    typedef Edge edge_type;

    /** Type of node iterators, which iterate over all mesh nodes. */
    class NodeIterator;
    /** Synonym for NodeIterator */
    typedef NodeIterator node_iterator;

    /** Type of edge iterators, which iterate over all mesh edges. */
    class EdgeIterator;
    /** Synonym for EdgeIterator */
    typedef EdgeIterator edge_iterator;
    /** Type of incident iterators, which iterate incident edges to a node. */
    class IncidentIterator;
    /** Synonym for IncidentIterator */
    typedef IncidentIterator incident_iterator;

    /** Type of edge iterators, which iterate over all mesh edges. */
    class TetrahedralIterator;
    /** Synonym for EdgeIterator */
    typedef TetrahedralIterator tet_iterator;

    /** Return the number of nodes in the mesh. */
    size_type num_nodes() const {
      return g_real_.num_nodes();
    }

    /** Return the number of edges in the mesh. */
    size_type num_edges() const {
      return g_real_.num_edges();
    }

    /** Return the number of tetrahedrals in the mesh. */
    size_type num_tetrahedral() const {
      return g_tet_.num_nodes();
    }

    class Node : private totally_ordered<Node> {
      public:

        /** Return this node's position.
         * @return The node's Point object
         * Complexity O(1)
         */
        Point& position() {
          return m_->g_real_.node(node_uid_).position();
        }

       /** Return this node's position as a constant.
    	 * @return The node's Point object
    	 * Complexity O(1)
    	 */
        const Point& position() const {
          return m_->g_real_.node(node_uid_).position();
        }

        /** Return this node's index, a number in the range [0, g_real.graph_size).
         * @return The node's index i, s.t. 0 <= i < g_real.num_nodes()
         * Complexity O(1)
         * */
        size_type index() const {
          return node_uid_;
        }

        /** Get this node's value (modifiable).
         * @return This node's node_value_type value as a reference.
         * Complexity O(1)
         */
        node_value_type& value(){
          return m_->g_real_.node(node_uid_).value();
        }

        /** Get this node's value (non-modifiable).
         * @return This node's node_value_type value as a constant.
         * Complexity O(1)
         */
        const node_value_type& value() const {
          return m_->g_real_.node(node_uid_).value();
        }

        /**Return a vector of tetrahedrals adjacent to the Node
         * @pre Valid Node.
         * @post return 0 <= result_vector.size() <= num_tetrahedrals()
         * @return vector containing Tetrahedrals
         *
         * Complexity: O(g_real_.node(node_uid_).degree())
         */
        vector<Tetrahedral> nodeAdjTetrahedral() const {
          vector<Tetrahedral> result;
          // Traverse all the incident edges
          for( typename g_real_type::IncidentIterator icit = m_->g_real_.node(node_uid_).edge_begin();
              icit != m_->g_real_.node(node_uid_).edge_end(); ++icit ){
            //  Collect all the tetrahedrals from the adjacent tetrahedrals of the edges
            const vector<Tetrahedral>& tets = Edge( m_, (*icit).node1().index(), (*icit).node2().index() ).edgeAdjTetrahedral();
            for( auto it = tets.begin(); it != tets.end(); ++it ){
              result.push_back((*it));
            }
          }

          //Delete duplicate tetrahedrals
          sort( result.begin(), result.end() );
          auto last = unique( result.begin(), result.end() );
          result.erase(last, result.end());
          return result;
        }

        /** Test whether this node and @a x are equal.
         * @param[in] @a x is a node
         * @return True if this node has the same mesh pointer and uid;
         * otherwise False.
         *
         * Complexity O(1)
         */
        bool operator==(const Node& x) const {
        	return ( (m_ == x.m_) && (node_uid_ == x.node_uid_) );
        }

        /** Test whether this node is less than @a x in the global order.
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         * Complexity O(1)
         */
        bool operator<(const Node& n) const {
          if( *this == n )
            return false;
          if( m_ == n.m_ ){// nodes in the same Mesh
            if( node_uid_ < n.node_uid_ ){
              return true;
            }
          }else{ // nodes not in the same Mesh
            if( m_ < n.m_ ){
              return true;
            }
          }
          return false;
        }

        /** Obtain an incident iterator pointing to the first incident edge of this node.
         * @return an incident_iterator pointing to the first edge connecting to this node.
         * Complexity: O(1).
         */
        incident_iterator edge_begin() const {
          return IncidentIterator( m_, m_->g_real_.node(node_uid_).edge_begin() );
        }
        /** Obtain an incident iterator represents the end of incident iterator.
         * @return an incident_iterator with idx_ equal to the number of incident edges.
         *  Complexity: O(1).
         */
        incident_iterator edge_end() const {
          return  IncidentIterator( m_, m_->g_real_.node(node_uid_).edge_end() );
        }

      private:
        /** Private constructor for Mesh to construct nodes
         * @param[in] mesh This new node's parent mesh.
         * @param[in] uid The uID that @a mesh uses to uniquely identify this node. @a uid >= 0.
         */
        Node( const Mesh* m, size_type node_uid )
            : m_( const_cast<Mesh*>(m) ), node_uid_(node_uid) {
        }

        //Private member variables
        /* Representative Invariants
         * m_ != nullptr
         * 0 <= uid_ < g.real_.num_nodes()
         */
        Mesh* m_; //Pointer to the parent mesh
        size_type node_uid_; //uid for the node in the parent mesh

        // Allow Mesh to access Node's private member data and functions.
        friend class Mesh;
    };

    class Edge : private totally_ordered<Edge> {
      public:

    	/** Return one of the two edge's nodes with uid with @a i.
    	 * @pre 0 <= @a i < 2
    	 * @post result_node.index() == node_uid1_ if node_uid1_ < node_uid2_
    	 *       else, result_node.index() == node_uid2_
    	 * @return Node such that if node_uid1_ < node_uid2_ returns node.index() == node_uid1_
    	 *                        else, returns node.index() == node_uid2_
    	 * Complexity: O(1).
    	 */
        Node node(size_type i) const {
        	assert( 0 <= i && i < 2 );
        	size_type node1, node2;
        	if( node_uid1_ < node_uid2_ ) {
        		node1 = node_uid1_;
        		node2 = node_uid2_;
        	}
        	else {
        		node1 = node_uid2_;
        		node2 = node_uid1_;
        	}
        	if( i == 0 ) {
        		return Node( m_, m_->g_real_.node(node1).index() );
        	}
        	else {
        		return Node( m_, m_->g_real_.node(node2).index() );
        	}
        }

        /** Return the node with the smaller uid of the edges 2 nodes.
         * @pre Valid Edge of the Mesh
         * @post result_node.index() == node_uid1_ if node_uid1_ < node_uid2_
         *       else, result_node.index() == node_uid2_
         * @return Node such that if node_uid1_ < node_uid2_ returns node.index() == node_uid1_
         *                        else, returns node.index() == node_uid2_
         * Complexity: O(1).
         */
        Node node1() const {
          return node(0);
        }

        /** Return the node with  the greater uid of the edges 2 nodes.
         * @pre Valid Edge of the Mesh
         * @post result_node.index() == node_uid2_ if node_uid1_ < node_uid2_
         *       else, result_node.index() == node_uid1_
         * @return Node such that if node_uid1_ < node_uid2_ returns node.index() == node_uid2_
         *                        else, returns node.index() == node_uid1_
         * Complexity: O(1).
         */
        Node node2() const {
          return node(1);
        }

        /** Retrieve the Edge's value (Modifiable)
         * @pre Valid Edge.
         * @return reference to this Edge's value.
         *
         * Complexity: same as g_real_.edge.value()
         */
        edge_value_type& value(){
          return m_->getEdgefrom2Nodes( node_uid1_, node_uid2_ ).value().edge_value;
        }

        /** Retrieve the Edge's value (Cannot be modified)
         * @pre Valid Edge.
         * @return reference to this Edge's value.
         *
         * Complexity: same as g_real_.edge.value()
         */
        const edge_value_type& value() const {
          return m_->getEdgefrom2Nodes( node_uid1_, node_uid2_ ).value().edge_value;
        }

        /** Return the length of the Edge
         * @pre Both nodes have valid positions.
         * @return Double length between the two nodes by Euclidean distance formula.
         * Complexity: O(1).
         */
        double length() const {
        	return norm( node1().position() - node2().position() );
        }

        /**Return a vector of tetrahedrals adjacent to the Edge
         * @pre Valid Edge.
         * @post return 1<= vector.size()
         * @return vector containing Tetrahedrals
         *
         * Complexity: O(d) //From getEdgefrom2Nodes which uses the underlying graph's
         * incident iterator
         */
        vector<Tetrahedral> edgeAdjTetrahedral() const {
          vector<Tetrahedral> result;
          const typename g_real_type::Edge& e = m_->getEdgefrom2Nodes( node_uid1_, node_uid2_ );

          //Loop through tet_uid and it to result
          for( unsigned int i = 0; i < e.value().tet_uid.size(); ++i) {
            result.push_back( m_->tetrahedral( e.value().tet_uid[i] ) );
          }
          return result;
        }

        /** Test whether this edge and @a x are equal.
         * @param[in] x Edge in a mesh
         * @return True if this Edge's mesh pointer is the same as @a x's mesh pointer &&
         * both nodes' uids match.
         *
         * Equal edges are from the same mesh and have the same nodes.
         *
         * Complexity: O(1).
         */
        bool operator==(const Edge& x) const {
          return ( (m_ == x.m_) &&
                 ( ( ( node_uid1_ == x.node_uid1_ ) && ( node_uid2_ == x.node_uid2_ ) ) ||
                   ( ( node_uid2_ == x.node_uid1_ ) && ( node_uid1_ == x.node_uid2_ ) ) ) );
        }

        /** Test whether this edge is less than @a x in the global order.
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         */
        bool operator<(const Edge& e) const {
          if( *this == e )
            return false;
          if( m_ == e.m_ ){// edges are in the same Mesh
            if( node_uid1_< e.node_uid1_ ){
              return true;
            }else if ( node_uid1_ == e.node_uid1_ && node_uid2_ < e.node_uid2_ ){
              return true;
            }
          }else{ // edges are not in the same Mesh
            if( m_ < e.m_ ){
              return true;
            }
          }
          return false;
        }

      private:

        /** Private constructor for Mesh to construct edges.
         * @param[in] m This new edge's parent mesh.
         * @param[in] node_uid1 The uID that @a mesh uses to access this Edge's first node in the adjacency list. @a node_uid1_ >= 0.
         * @param[in] node_uid2 The uID that @a mesh uses to access this Edge's second node in the adjacency list. @a node_uid2_ >= 0.
         */
        Edge(const Mesh* m, size_type node_uid1, size_type node_uid2 )
        		: m_(const_cast<Mesh*>(m)), node_uid1_(node_uid1), node_uid2_(node_uid2) {}

        // Private Member Variables
        /* Representative Invariants
         * m_ != nullptr
         * 0 <= node_uid1_, node_uid2 < g.real_.num_nodes()
         * No self edges: edge s.t. node_uid1_ != node_uid2_
         * Every edge belongs to at least one triangle:
         * 	∀i,j,Edge(i,j) in g_real_,∃ a tetrahedral tri from g_tri_ s.t.
         * 	g_tri.hasedge(tet.node(0).index(), tet.node(1).index()) ,
         * 	|| g_tet_.hasedge(tet.node(0).index(), tet.node(2).index()),
         * 	|| g_tet_.hasedge(tet.node(1).index(), tet.node(2).index())
         */
        Mesh* m_; //Pointer to the parent mesh for the edge
        size_type node_uid1_, node_uid2_; //uid for the node in the parent mesh

        // Allow Mesh to access Node's private member data and functions.
        friend class Mesh;
    };

    class Tetrahedral : private totally_ordered<Edge> {
      public:

    	/** Return one of the four tetrahedral's nodes with uid with @a i.
    	 * @pre 0 <= @a i < 4
    	 * @return Node with the @ith smallest index among the four nodes of this tetrahedral
    	 * Complexity: O(1).
    	 */
        Node node(size_type i) const {
          assert( 0 <= i && i <= 3 );
          const size_type& id = m_->g_tet_.node(tet_uid_).value().node_uid_[i];
          return Node( m_, id);
        }

        /** Return one of the six tetrahedral's edges with uid with @a i.
         * @pre 0 <= @a i < 3, 0 <= @a j < 3
         * @return Edge (node(i), node(j))
         * Complexity: O(1).
         */
        Edge edge(size_type i, size_type j) const {
          assert( i != j );
          assert( 0 <= i && i <= 3 );
          assert( 0 <= j && j <= 3 );
          return Edge(m_, node(i).index(), node(j).index());
        }

        /** Return one of the four tetrahedral's nodes with uid with @a i.
         * @pre 0 <= @a i < 3
         * @post edge(0) is (node0, node1),
         *       edge(1) is (node0, node2),
         *       edge(2) is (node0, node3),
         *       edge(3) is (node1, node2),
         *       edge(4) is (node1, node3)
         *       edge(5) is (node2, node3)
         *       where node0.index() < node1.index() < node2.index() < node3.index()
         * @return Edge such that @a i is the ordering of listed in @post
         *
         * Complexity: O(1).
         */
        Edge edge(size_type i) const {
        	assert( 0 <= i && i < 6 );
        	const size_type& id0 = m_->g_tet_.node(tet_uid_).value().node_uid_[0];
        	const size_type& id1 = m_->g_tet_.node(tet_uid_).value().node_uid_[1];
        	const size_type& id2 = m_->g_tet_.node(tet_uid_).value().node_uid_[2];
        	const size_type& id3 = m_->g_tet_.node(tet_uid_).value().node_uid_[3];

        	switch(i){
        	case 0:
        		return Edge(m_, id0, id1);
        	case 1:
        		return Edge(m_, id0, id2);
        	case 2:
        		return Edge(m_, id0, id3);
        	case 3:
        		return Edge(m_, id1, id2);
        	case 4:
        		return Edge(m_, id1, id3);
        	case 5:
        		return Edge(m_, id2, id3);
        	default:
        		assert(false);
        		return Edge(m_, id0, id1);
        	}
        }

        /** Return this tetrahedral's uid, a number in the range [0, g_tet_.num_nodes()).
         * @return The tetrahedral's uid i, s.t. 0 <= i < g_tet_.num_nodes()
         * Complexity O(1)
         * */
        size_type index() const {
        	return tet_uid_;
        }

        /** Get this tetrahedral's value (modifiable).
         * @return This tetrahedral's node_value_type value as a reference.
         * Complexity O(1)
         */
        tet_value_type& value(){
          return m_->g_tet_.node(tet_uid_).value().tet_value;
        }

        /** Get this tetrahedral's value (Non-modifiable).
         * @return This tetrahedral's node_value_type value as a reference.
         * Complexity O(1)
         */
        const tet_value_type& value() const {
          return m_->g_tet_.node(tet_uid_).value().tet_value;
        }

        /** Return a vector of Tetrahedrals adjacent to the Tetrahedral
         * @pre Valid Tetrahedral.
         * @post return 0 <= result_vector.size() <= 4
         * @return vector containing Tetrahedral
         *
         * Complexity: O(g_tet_.node(tet_uid_).degree())
         */
        vector<Tetrahedral> tetAdjTetrahedral() const {
          vector<Tetrahedral> adj;
          //Use incident iterator to iterate through all adjacent tetrahedral in g_tet_
          for( typename g_tet_type::IncidentIterator it = m_->g_tet_.node(tet_uid_).edge_begin();
               it != m_->g_tet_.node(tet_uid_).edge_end(); ++it ){
            adj.push_back( Tetrahedral(m_, (*it).node2().index() ) );
          }
          return adj;
        }

        /** Return the volume of this tetrahedral
         * @pre Valid Tetrahedral.
         * @return Double volume of this tetrahedral.
         *         The volume will be positive if its sign is the same as the original volume.
         *         The volume will be negative if its sign is different as the original volume.
         * Complexity: O(1)
         */
        double volume() const {
          double v( dot( (node(1).position() - node(0).position()),
                cross ( node(2).position() - node(0).position(), node(3).position() - node(0).position() )  )/6.0 );
          if( m_->g_tet_.node(tet_uid_).value().volume_sign ){
            return v;
          } else {
            return -1.0 * v;
          }

//        	return dot( (node(1).position() - node(0).position()),
//        			cross ( node(2).position() - node(0).position(), node(3).position() - node(0).position() )  )/6.0;

        }

        /** Return true if tetrahedral is on the surface
         * @pre Valid Tetrahedral.
         * @return bool if this tetrahedral is on the surface
         * Complexity: O(1)
         */
        bool isSurface() const {
        	return m_->g_tet_.node(tet_uid_).degree() < NUM_TET_ADJ_TET;
        }


        /** Test whether this Tetrahedral and @a x are equal.
         * @param[in] x Tetrahedral in a mesh
         * @return Equal edges are from the same mesh and have the same tetrahedral uids.
         * Complexity: O(1).
         */
        bool operator==(const Tetrahedral& t) const {
        	return ( ( m_ == t.m_ ) && ( tet_uid_ == t.tet_uid_ ) );
        }

        /** Test whether this Tetrahedral is less than @a x in the global order.
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         */
        bool operator<(const Tetrahedral& t) const {
          if( *this == t )
            return false;
          if( m_ == t.m_ ){// tetrahedral in the same Mesh
            if( tet_uid_ < t.tet_uid_ ){
              return true;
            }
          }else{ // tetrahedral not in the same Mesh
            if( m_ < t.m_ ){
              return true;
            }
          }
          return false;
        }

        /** Calculate the triangle's specific edge's outward normal vector
         * @param[in] Valid edge of this triangle
         * @return Point the edge's outward normal vector
         *
         * Complexity: O(1).

        Point outward(const Edge& e) const {
          // initialize to invalid index
        	size_type uid3 = m_->g_real_.num_nodes();
        	for( unsigned int i = 0; i < NUM_TET_ADJ_TET; ++i ) {
        		const size_type& node_uid = m_->g_tet_.node(tet_uid_).value().node_uid_[i];
        		if( node_uid != e.node(0).index() &&
        		    node_uid != e.node(1).index() ) {
        			uid3 = node_uid;
        			break;
        		}
        	}
        	assert(uid3 < m_->g_real_.num_nodes());

        	// Vector math
        	const Point& A = e.node(0).position();
        	const Point& B = e.node(1).position();
        	const Point& C = m_->g_real_.node(uid3).position();
        	Point BA = A-B;
        	Point normBA = Point( -1.0*BA.y, BA.x, BA.z );
        	double D = dot( normBA, C-B );

        	if(D > 0) {
        		return -1.0 * normBA;
        	}
        	else {
        		return normBA;
        	}
        }*/

      private:
        /** Private constructor for Mesh to construct tetrahedrals
         * @param[in] mesh This new tetrahedral's parent mesh.
         * @param[in] uid The uID that @a mesh uses to uniquely identify this tetrahedral. @a uid >= 0.
         */
        Tetrahedral( const Mesh* m, size_type tet_uid )
            : m_( const_cast<Mesh*>(m) ), tet_uid_(tet_uid) {
        }

        //Private member variables
        /* Representative Invariants
         * m_ != nullptr
         * 0 <= uid_ < g.tet_.num_nodes()
         * No repeated node_uids
         * There exists a valid edge between all of the node_uids
         */
        Mesh* m_; //Pointer to the parent mesh
        size_type tet_uid_; //uid for the tetrahedral in the parent mesh

        // Allow Mesh to access Node's private member data and functions.
        friend class Mesh;

    };// End of Tetrahedral class

    ///////////////
    // Iterators //
    ///////////////

    /** @class Mesh::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : private equality_comparable<NodeIterator> {
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
       * @return Node corresponding to the node in @a g_real_.
       *
       * Complexity: same as g_real_type::NodeIterator operator*(), probably O(1).
       */
      Node operator*() const {
        return Node( m_, (*nit_).index() );
      }
      /**Increment NodeIterator and return the next position.
       * @post the @a nit_ increase by 1, may point to an invalid position.
       * @return the modified NodeIterator.
       *
       * Complexity: same as g_real_type::NodeIterator operator++(), probably O(1).
       */
      NodeIterator& operator++() {
        ++nit_;
        return *this;
      }
      /** Test the equality of NodeIterator:
       * @return true if the two NodeIterators belong to the same mesh, have the same @a nit_.
       *
       */
      bool operator==(const NodeIterator& target) const {
        return ( m_ == target.m_ && nit_ == target.nit_ );
      }

     private:
      friend class Mesh;
      Mesh* m_;
      typename g_real_type::NodeIterator nit_;
      NodeIterator( const Mesh* m, typename g_real_type::NodeIterator nit )
          : m_( const_cast<Mesh*>(m) ), nit_(nit) {
      }
    };

    /** Obtain a node_iterator pointing to the start of the mesh's nodes.
     * @return a node_iterator with @a nit_ = @a g_real_.node_begin(),
     *           it could be invalid if there is no node in the mesh.
     *
     * Complexity: same as g_real.node_begin(), probably O(1).
     */
    node_iterator node_begin() const {
      return NodeIterator( this, g_real_.node_begin() );
    }
    /** Obtain a node_iterator representing the end of the mesh's nodes.
     * @return a node_iterator with @a nit_ = @a g_real_.node_end()
     *
     * Complexity: same as g_real_.node_end(), probably O(1).
     */
    node_iterator node_end() const {
      return NodeIterator( this, g_real_.node_end() );
    }
    //=====End of NodeIterator related functions=====//

    /** @class Mesh::EdgeIterator
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
       * @return the Edge corresponding to the edge in g_real_.
       *
       * Complexity: same as g_real_type::EdgeIterator operator*(), probably O(num_nodes()).
       */
      Edge operator*() const {
        return Edge( m_, (*eit_).node1().index(), (*eit_).node2().index() );
      }
      /** Increase the edge iterator
       * @post the @a eit_ increase by 1, may point to an invalid position.
       * @return the modified EdgeIterator.
       *
       * Complexity: same as g_real_type::EdgeIterator operator++(), probably O(num_nodes()).
       */
      EdgeIterator& operator++() {
        ++eit_;
        return *this;
      }
      /** Test the equality of EdgeIterator.
       * @param[in] target EdgeIterator
       * @return True if both EdgeIterators are in the same mesh and have the same eit_.
       */
      bool operator==(const EdgeIterator& target ) const {
        return ( m_ == target.m_ && eit_ == target.eit_ );
      }

     private:
      friend class Mesh;
      Mesh* m_;
      typename g_real_type::EdgeIterator eit_;
      EdgeIterator( const Mesh* m, typename g_real_type::EdgeIterator eit )
          : m_( const_cast<Mesh*>(m) ), eit_(eit) {
      }
    };

    /** Obtain the begin iterator of edge iterator
     * @return the first edge iterator with @a eit_ = g_real_.edge_begin(),
     *         it could be invalid if there is no edge in the mesh.
     *
     * Complexity: same as g_real_.edge_begin(), probably O(num_nodes()).
     */
    edge_iterator edge_begin() const {
      return EdgeIterator( this, g_real_.edge_begin() );
    }
    /** Obtain the end of edge iterator
     * @return the end edge iterator, set eit_ = g_real_.edge_end()
     *
     * Complexity: same as g_real_.edge_end(), probably O(1).
     */
    edge_iterator edge_end() const {
      return EdgeIterator( this, g_real_.edge_end() );
    }
    //=====End of EdgeIterator related functions=====//

    /** @class Mesh::TetrahedralIterator   TetrahedralIterator tet_iterator
     * @brief Iterator class for Tetrahedrals. A forward iterator. */
    class TetrahedralIterator : private equality_comparable<TetrahedralIterator> {
     public:
      // These type definitions help us use STL's iterator_traits.
      /** Element type. */
      typedef Tetrahedral value_type;
      /** Type of pointers to elements. */
      typedef Tetrahedral* pointer;
      /** Type of references to elements. */
      typedef Tetrahedral& reference;
      /** Iterator category. */
      typedef std::input_iterator_tag iterator_category;
      /** Difference between iterators */
      typedef std::ptrdiff_t difference_type;

      /** Construct an invalid TetrahedralIterator. */
      TetrahedralIterator() {
      }

      /** Obtain the abstract tetrahedral this iterator pointing.
       * @pre @a tit_ < num_nodes()
       * @return Node with this mesh's pointer and index of this node
       *
       * Complexity: g_tet_type::NodeIterator operator*(), probably O(1).
       */
      Tetrahedral operator*() const {
        return Tetrahedral( m_, (*tit_).index() );
      }
      /**Increment TetrahedralIterator and return the next position.
       * @post the @a tet_ increase by 1, may point to an invalid position.
       * @return the modified TetrahedralIterator.
       *
       * Complexity: g_tet_type::NodeIterator operator++(), probably O(1).
       */
      TetrahedralIterator& operator++() {
        ++tit_;
        return *this;
      }
      /** Test the equality of TetrahedralIterator:
       * @return true if the two TetrahedralIterator belong to the same mesh, and have the same tit_.
       */
      bool operator==(const TetrahedralIterator& target) const {
        return ( m_ == target.m_ && tit_ == target.tit_ );
      }

     private:
      friend class Mesh;
      Mesh* m_;
      typename g_tet_type::NodeIterator tit_;
      TetrahedralIterator( const Mesh* m, typename g_tet_type::NodeIterator tit )
          : m_( const_cast<Mesh*>(m) ), tit_(tit) {
      }
    };

    /** Obtain a tet_iterator pointing to the start of the mesh's tetrahedral.
     * @return a tet_iterator at the beginning position of the mesh's tetrahedrals,
     *         it could be invalid if there is no tetrahedral in this mesh.
     *
     * Complexity: same as g_tet_.node_begin(), probably O(1).
     */
    tet_iterator tetrahedral_begin() const {
      return TetrahedralIterator( this, g_tet_.node_begin() );
    }
    /** Obtain a tet_iterator representing the end of the mesh's tetrahedral.
     * @return a tet_iterator with index = @a num_tetrahedrals()
     *
     * Complexity: same as g_tet_.node_end(), probably O(1).
     */
    tet_iterator tetrahedral_end() const {
      return TetrahedralIterator( this, g_tet_.node_end() );
    }
    //=====End of TetrahedralIterator related functions=====//

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
        return Edge( m_, (*icit_).node1().index(), (*icit_).node2().index() );
      }
      /** Increase the incident iterator.
       * @post Increase @a node2_pos_ to the next index in the link_edge of current node.
       * @post @a node2_pos_ will not increase if incident iterator equals to the end iterator.
       * @return the advanced IncidentIterator, may be valid or invalid position.
       *
       * Complexity: O(1).
       */
      IncidentIterator& operator++() {
        ++icit_;
        return *this;
      }
      /** Compare the equality of IncidentIterator.
       * @return True if the two IncidentIterators are in the same graph, have the same index of two sides of nodes.
       */
      bool operator==(const IncidentIterator& target) const {
        return ( icit_ == target.icit_ );
      }

     private:
      friend class Mesh;
      Mesh* m_;
      typename g_real_type::IncidentIterator icit_;
      IncidentIterator( const Mesh* m, typename g_real_type::IncidentIterator icit )
          : m_( const_cast<Mesh*>(m) ), icit_(icit) {
      }
    };
    //=======================End of IncidentIterator=============================//

    /** Return the total number of nodes in the mesh.
     * Complexity: same as g_real_.num_nodes(), probably O(1)
     */
    Node node(size_type i) const {
      return Node( this, g_real_.node(i).index() );
    }

    /** Return the total number of tetrahedrals in the mesh.
     * Complexity: same as g_tet_.size(), probably O(1)
     */
    Tetrahedral tetrahedral(size_type i) const {
      return Tetrahedral( this, g_tet_.node(i).index() );
    }

    /** Add a node to the mesh, returning the added node.
     * @param[in] p The new node's position
     *            node_value user defined node_value
     * @post new g_real_.size() == old g_real_.size() + 1
     * @post result_node.index() == old g_real_.size()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& p, const node_value_type& node_value){
    	typename g_real_type::Node n = g_real_.add_node(p, node_value);
    	return Node( this, n.index() );
    }

    /**Add a tetrahedral to the mesh, return the added tetrahedral.
     * @param[in] n0, n1, n2, n3: The new tetrahedral's four nodes
     *            tet_value: User defined tet_value
     * @return a Tetrahedral object tet with tet.node(0) is min(@a n0, @a n1, @a n2, @n3)
     *           tet.node(1) is the 2nd smallest among (@a n0, @a n1, @a n2),
     *           tet.node(2) is the 3rd smallest among (@a n0, @a n1, @a n2),
     *           and tet.node(3) is max(@a n0, @a n1, @a n2, @n3)
     * @pre @a n0, @a n1, @a n2, @a n3 are valid Mesh::Node
     * @pre Tetrahedral compsed of @a n0, @a n1, @a n2, @a n3
     * @post new g_tet_.size() == old g_tet_.size() + 1
     *       result_tet.index() == old g_tet_.size()
     * Complexity: same as g_real_.add_node()
     */
    Tetrahedral add_tetrahedral(const Node& n0, const Node& n1, const Node& n2, const Node& n3, const tet_value_type& tet_value = tet_value_type()) {
      // all three nodes have to be existed
      assert( n0.index() < g_real_.size() && n1.index() < g_real_.size() &&
              n2.index() < g_real_.size() && n3.index() < g_real_.size() );

      size_type tet_uid = g_tet_.size();

      // sort n0, n1, n2, n3
      size_type ni[4] = { n0.index(), n1.index(), n2.index(), n3.index() };
      if( ni[0] > ni[1] ) swap( ni[0], ni[1] );
      if( ni[2] > ni[3] ) swap( ni[2], ni[3] );
      if( ni[0] > ni[2] ) swap( ni[0], ni[2] );
      if( ni[1] > ni[3] ) swap( ni[1], ni[3] );
      if( ni[1] > ni[2] ) swap( ni[1], ni[2] );

      // add edges in g_real_
      vector<typename g_real_type::Edge> edges(6);
      edges[0] = g_real_.add_edge( g_real_.node(ni[0]), g_real_.node(ni[1]) );
      edges[1] = g_real_.add_edge( g_real_.node(ni[0]), g_real_.node(ni[2]) );
      edges[2] = g_real_.add_edge( g_real_.node(ni[0]), g_real_.node(ni[3]) );
      edges[3] = g_real_.add_edge( g_real_.node(ni[1]), g_real_.node(ni[2]) );
      edges[4] = g_real_.add_edge( g_real_.node(ni[1]), g_real_.node(ni[3]) );
      edges[5] = g_real_.add_edge( g_real_.node(ni[2]), g_real_.node(ni[3]) );

      // add tetrahedral uid to edges
      for( unsigned i = 0; i < edges.size(); ++i ){
        edges[i].value().tet_uid.push_back(tet_uid);
      }

      // add node in g_tet_
      internal_tetrahedral itri( ni[0], ni[1], ni[2], ni[3], tet_value );
      g_tet_.add_node( Point(), itri );

      // check the sign of volume and store it
      Tetrahedral tet( this, tet_uid );
      if ( tet.volume() < 0 ){
        g_tet_.node(tet_uid).value().volume_sign = false;
      }

      // add edges in g_tet_ to record the adjacent tetrahedrals of a tetrahedral.
      add_tet_edges(ni[0], ni[1], ni[2], ni[3], tet_uid);
      add_tet_edges(ni[2], ni[3], ni[0], ni[1], tet_uid);

      return tet;
    }


  private:

    /* Get the Graph::Edge from two nodes
     * @pre: the edge in g_real_ need to be existed.
     * Complexity: O(g_real_.node(@ id1).degree() or g_real_.node(@ id2).degree())
     */
    typename g_real_type::Edge getEdgefrom2Nodes ( const size_type& id1, const size_type& id2 ) const{
      // traverse the node with smaller degree
      size_type i1, i2;
      if( g_real_.node(id1).degree() < g_real_.node(id2).degree() ){
        i1 = id1;
        i2 = id2;
      } else {
        i1 = id2;
        i2 = id1;
      }

      const typename g_real_type::Node& n1 = g_real_.node(i1);
      const typename g_real_type::Node& n2 = g_real_.node(i2);
      for( typename g_real_type::incident_iterator it = n1.edge_begin(); it != n1.edge_end(); ++it ){
        if( (*it).node2() == n2 ){
          return (*it);
        }
      }
      // check the edge exists before call this function!
      assert(false);
      return g_real_.edge(0);
    }

    /* Add edges in the g_tet_
     * @param[in] en0, en1: two nodes' indices of tetrahedral @ti forming an edge
     *            nx, ny: The other two nodes of tetrahedral @ti
     *            ti: index of the tetrahedral
     * @post new g_tet_.num_edges() - old g_tet_.num_edges() 0 or 1 or 2
     * Complexity: same as g_real_.add_node()
     *
     */
    void add_tet_edges( size_type en0, size_type en1, size_type nx, size_type ny, size_type ti ){
      const vector<Tetrahedral> tets = Edge(this, en0, en1).edgeAdjTetrahedral();
      // loop through the adjacent tetrahedrals of edge(en0, en1)
      for( auto itet = tets.begin(); itet != tets.end(); ++itet ){
        // skip self tetrahetral
        if( (*itet).index() != ti ){
          // check the nodes of this tetrahedral.
          // if it has node nx or ny, then it is an adjacent tetrahedral of tetrahedral ti
          for( unsigned j = 0; j < 4; ++j ){
            if ( (*itet).node(j).index() == nx || (*itet).node(j).index() == ny ){
              g_tet_.add_edge( g_tet_.node((*itet).index()), g_tet_.node(ti) );
            }
          }
        }
      }
    }

    //internal_edge data structure to store tet_uid and user defined edge_values
    struct internal_edge {
      vector<size_type> tet_uid;
      edge_value_type edge_value;
    };

    /* internal_tetrahedral data structure to store tetrahedral related data,
     * tetrahedral's node uid, and user defined tet_value
     */
    struct internal_tetrahedral {
      size_type node_uid_[NUM_TET_ADJ_TET];
      tet_value_type tet_value;
      bool volume_sign;
      internal_tetrahedral(size_type id0, size_type id1, size_type id2, size_type id3, tet_value_type tet_value ) :
        tet_value(tet_value) {
        node_uid_[0] = id0;
        node_uid_[1] = id1;
        node_uid_[2] = id2;
        node_uid_[3] = id3;
        volume_sign = true;
      }
    };


    /* Representative Invariants
     * for all i, node_i in g_tet_, node_i.degree() <= 3
     */

    //Graph that stores phyiscal nodes and edges
    Graph<node_value_type, internal_edge> g_real_;

    //Graph to store tetrahedral as nodes
    // bool is not used
    Graph<internal_tetrahedral, bool> g_tet_;

};
