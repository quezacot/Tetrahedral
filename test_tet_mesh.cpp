/*
 * CS 207: Final Project Tetrahedral Mesh
 * Yung-jen Cheng HUID #20947802
 * Jeffrey (Shih-kai) Shen HUID #70949288
 */

/**
 * @file test_tet_mesh.cpp
 * Tetrahedral Mesh testing file
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Tetrahedrals (one per line) defined by 4 indices into the point list
 */

#include <fstream>
#include <cmath>
using std::pow;
using std::exp;
#include <iostream>
using std::cout;
using std::endl;
#include <map>
using std::map;
#include <string>
using std::string;

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "CS207/Point.hpp"

#include "tet_mesh.hpp"


//DEBUG
//#define TETRAHEDRAL_OUTPUT
//#define DEBUG_MESH

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;
struct NodeData;
struct EdgeData;
struct TetrahedralData;
typedef Mesh<NodeData, EdgeData, TetrahedralData> MeshType;
typedef typename MeshType::Node Node;
typedef typename MeshType::Edge Edge;

/** Custom structure of data to store with Nodes */
struct NodeData {
	Point velocity;  //< Node velocity
	double mass;     //< Node mass
};

/** Custom structure of data to store with Edges */
struct EdgeData {
	double length;   //< Edge length
};

/** Custom structure of data to store with Tetrahedral */
struct TetrahedralData {
	double initialVolume;   //< Initial Tetrahedral Volume
};

struct NullConstraint;

/** Change a mesh's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] m      Mesh
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]	  constrain Constrain object defining the constraints for the simulation
 *                Defaulted to FixedConstraint where the Points (0,0,0) and (1,0,0) are fixed
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the mesh and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename M, typename F, typename C = NullConstraint>
double symp_euler_step(M& m, double t, double dt, F force, C constrain = C() ) {

	//Set constraints
	constrain(m, t);

	// Compute the {n+1} node positions
	for (auto it = m.node_begin(); it != m.node_end(); ++it) {
		auto n = *it;

		// Update the position of the node according to its velocity
		// x^{n+1} = x^{n} + v^{n} * dt
		n.position() += n.value().velocity * dt;
	}

	// Compute the {n+1} node velocities
	for (auto it = m.node_begin(); it != m.node_end(); ++it) {
		auto n = *it;
		// v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
		n.value().velocity += force(n, t) * (dt / n.value().mass);
	}



	return t + dt;
}


/** Dashpot Force Functor that returns the Dashpot Force
 */
struct DashpotForce {
	double K_;
	double C_;
	/** DashpotForce Constructor.
	 * @param[in] K Spring constant in N/m
	 * @param[in] K Damping constant in in N*s/m
	 */
	DashpotForce(double K = 100, double C=100) : K_(K), C_(C) {}

	/** Calculates Dashpot Force
	 * @param[in] n Valid node.
	 * @param[in] t Valid time.
	 * @return Point object that represents the dashpot force.
	 */
	Point operator()(Node n, double t) {
		(void) t;     // silence compiler warnings
		Point dashpotForce(0,0,0);
		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
			Node incidentNode = (*it).node2();
			if( (*it).node1().index() == n.index() ) {
				incidentNode = (*it).node2();
			} else {
				incidentNode = (*it).node1();
			}

			double distance = norm(n.position() - incidentNode.position());
			//Spring component
			double spring_comp = K_ * (distance - (*it).value().length);
			//Damping component
			Point c_comp = C_* (n.value().velocity - incidentNode.value().velocity) *
					       (n.position() - incidentNode.position()) / distance;
			Point unitVector = (n.position() - incidentNode.position()) / distance;

			dashpotForce += (-1.0) * (spring_comp + c_comp) * unitVector;
		}
		return dashpotForce;
	}
};

/** Volume Penalty Force Functor that returns the Volume Penalty Force
 */
struct VolumePenaltyForce {
	MeshType* m_;
	double K_;
	/** VolumePenaltyForce Constructor.
	 * @param[in] g Gravity in m/s^2.
	 */
	VolumePenaltyForce(MeshType* m, double K) : m_(m), K_(K) {}

	/** Calculates VolmePenalty Force
	 * @param[in] n Valid node.
	 * @param[in] t Valid time.
	 * @return Point object that represents the volume penalty force.
	 */
	Point operator()(Node n, double t) {
		(void) t;     // silence compiler warnings
		Point volumePenaltyForce(0,0,0);
		auto AdjTetrahedral = n.nodeAdjTetrahedral();

		// Loop through all of the node's adj tetrahedrals
		for(unsigned k = 0; k < AdjTetrahedral.size(); ++k) {
			MeshType::Tetrahedral tet = AdjTetrahedral[k];
			// Mass weighted centroid
			Point baryCenter(0,0,0);
			double totalMass = 0;

			// Loop through each of the tet's node to calculate baryCenter
			for(unsigned i = 0; i < NUM_TET_ADJ_TET; ++i) {
				baryCenter += tet.node(i).position() * tet.node(i).value().mass;
				totalMass += tet.node(i).value().mass;
			}
			baryCenter = baryCenter/totalMass;

			Point unitVector = (n.position() - baryCenter)/norm(n.position() - baryCenter);
			volumePenaltyForce += -K_ * (tet.volume() - tet.value().initialVolume) * unitVector;

		}
		return volumePenaltyForce;
	}
};

/** Gravity Force Functor that returns the Gravity Force
 */
struct GravityForce {
	double gravity_;
	/** GravityForce Constructor.
	 * @param[in] g Gravity in m/s^2.
	 */
	GravityForce(double g = grav) : gravity_(g) {}

	/** Calculates Gravity Force
	 * @param[in] n Valid node.
	 * @param[in] t Valid time.
	 * @return Point object that represents the gravity force calculated as m*g.
	 */
	Point operator()(Node n, double t) {
		(void) t;     // silence compiler warnings
		return n.value().mass * Point(0,0,-1.0 * gravity_);
	}
};

/** Combine Force Functor that returns a combination of forces
 * @param[in] Two valid forces in f1 and f2.
 */
template <typename Force1, typename Force2>
struct CombineForces {
	Force1 f1_;
	Force2 f2_;

	/** CombineForces Constructor.
	 * @param[in] f1 First valid force.
	 * @param[in] f2 Second valid force.
	 */
	CombineForces(Force1 f1, Force2 f2)
			: f1_(f1), f2_(f2) {}

	/** Calculates Combine Forces
	 * @param[in] n Valid node.
	 * @param[in] t Valid time.
	 * @return Point object that represents the combination of forces of @a f1_ and @a f2_.
	 */

	Point operator() (Node n, double t){
		return f1_(n,t) + f2_(n,t);
	}
};

/** Combine Force Function that returns a combination of forces
 * @param[in] Two valid forces in f1, and f2.
 * @pre Valid forces that input a node and time with a function call of f(n,t).
 * @return A CombineForces object that adds up the forces.
 */
template <typename Force1, typename Force2>
CombineForces<Force1, Force2> make_combined_force(Force1 f1, Force2 f2){
	return CombineForces<Force1, Force2> (f1, f2);
}

/** Combine Force Function that returns a combination of forces
 * @param[in] Three valid forces in f1, f2, and f3.
 * @pre Valid forces that input a node and time with a function call of f(n,t).
 * @return A CombineForces object that adds up the forces.
 */
template <typename Force1, typename Force2, typename Force3>
CombineForces<CombineForces<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3){
	return CombineForces<CombineForces<Force1, Force2>, Force3> (CombineForces<Force1, Force2>(f1,f2), f3);
}

/** Null Constraint
 */
struct NullConstraint {

	/** Null Constraint Setter
	 * @param[in] g Valid mesh.
	 * @param[in] t Valid time.
	 * @return Point object that represents the combination of forces of @a f1_ and @a f2_.
	 */
	void operator()(MeshType& m, double t) {
		(void) t;     // silence compiler warnings
		(void) m;
	}

};

/** Fixed Constraint where you can specify some points to be static
 */
struct FixedConstraint {
  vector<Point> cpoints;
  // the vector of Points that you want to constrain
  FixedConstraint(const vector<Point>& v) : cpoints(v) {}

	/** Fixed Constraint Setter
	 * @param[in] g Valid mesh.
	 * @param[in] t Valid time.
	 * @post The velocity of Point(0,0,0) and Point(1,0,0) are 0.
	 */
	void operator()(MeshType& m, double t) {
			(void) t;     // silence compiler warnings
			for(auto it = m.node_begin(); it != m.node_end(); ++it) {
			  for( auto p = cpoints.begin(); p != cpoints.end(); ++p ){
			    if( (*it).position() == (*p) )
			      (*it).value().velocity =  Point(0,0,0);
			  }
			}
		}
};

/** Horizontal Plane Constraint that models an impassable plane.
 */
struct HPlaneConstraint {
double z_constraint_; //coordinate for the horizontal plane

	/** HPlaneConstraint Constructor.
	 * @param[in] z_constraint Sets the z-coordinate to define the horizontal plane.
	 */
	HPlaneConstraint(double z_constraint) : z_constraint_(z_constraint) {}

	/** Horizontal Constraint Setter
	 * @param[in] g Valid mesh.
	 * @param[in] t Valid time.
	 * @post The velocity of @a z_constraint_ is 0 and is set to the closest point to
	 *       the horizontal plane as defined by @a z_constraint_.
	 */
	void operator()(MeshType& m, double t) {
		(void) t;     // silence compiler warnings
		for(auto it = m.node_begin(); it != m.node_end(); ++it) {
			Node n = (*it);
			if(n.position().z < z_constraint_) {
				n.position().z = z_constraint_;
				n.value().velocity.z = 0;
				// To move the ball
//				n.value().velocity.x *= -0.03;
//				n.value().velocity.y *= -0.03;
			}
		}
	}
};

/** Combine Constraints Functor that returns a combination of constraints
 * @param[in] Two valid constraints in c1 and c2.
 */
template <typename Constraint1, typename Constraint2>
struct CombineConstraints {
	Constraint1 c1_;
	Constraint2 c2_;

	CombineConstraints(Constraint1 c1, Constraint2 c2)
			: c1_(c1), c2_(c2) {}

	void operator() (MeshType& m, double t){
		c1_(m,t);
		c2_(m,t);
	}
};

/** Combine Constraints Function that returns a combination of constraints
 * @param[in] Two valid constraints in @a c1, and @a c2.
 * @pre Valid constraints that input a node and time with a function call of c(g,t).
 * @return A CombineConstraints object that executes the constraints.
 */
template <typename Constraint1, typename Constraint2>
CombineConstraints<Constraint1, Constraint2> make_combined_constraints(Constraint1 c1, Constraint2 c2){
	return CombineConstraints<Constraint1, Constraint2> ({c1, c2});
}

/** Combine Constraints Function that returns a combination of constraints
 * @param[in] Valid constraints in @a c1, @a c2, and @a c3.
 * @pre Valid constraints that input a node and time with a function call of c(g,t).
 * @return A CombineConstraints object that executes the constraints.
 */
template <typename Constraint1, typename Constraint2, typename Constraint3>
CombineConstraints<CombineConstraints<Constraint1, Constraint2>, Constraint3> make_combined_constraints(Constraint1 c1, Constraint2 c2, Constraint3 c3){
	return CombineConstraints<CombineConstraints<Constraint1, Constraint2>, Constraint3> ({{c1, c2}, c3});
}

#ifdef TETRAHEDRAL_OUTPUT
void printAllTetrahedralInformation(const MeshType::Mesh& m, double time) {
	cout << endl;
	for(unsigned int j = 0; j < m.num_tetrahedral(); ++j) {
		cout << endl;
		MeshType::Tetrahedral tet = m.tetrahedral(j);
		cout << "Tetrahedral: " << tet.index() << "@" << time << endl;
		cout << "Surface: " << tet.isSurface() << endl;
		cout << "Volume: " << tet.volume() << endl;
		cout << "Value initialVolume: " << tet.value().initialVolume << endl;

		cout << "Adj Tetrahedral from this tetrahedral: ";
		auto AdjTetrahedral = tet.tetAdjTetrahedral();
		for(unsigned int l = 0; l < AdjTetrahedral.size(); ++l) {
			cout << AdjTetrahedral.at(l).index() << ", ";
		}
		cout << endl;

		for(unsigned int i = 0; i < 6; ++i) {
			cout << "Edge " << i << " : " <<
					tet.edge(i).node(0).position() << ", " <<
					tet.edge(i).node(1).position() << endl;
			cout << tet.edge(i).node(0).index() << ", " <<
					tet.edge(i).node(1).index() << endl;

			// Adj Tetrahedral from edge
			cout << "Adj Tetrahedral from edge: ";
			for(unsigned int k = 0; k < tet.edge(i).edgeAdjTetrahedral().size(); ++k) {
				cout << tet.edge(i).edgeAdjTetrahedral().at(k).index() << ", ";
			}
			cout << endl;
			cout << endl;
		}
		cout << "Node 0: " << tet.node(0).position() << endl;
		cout << "Node 1: " << tet.node(1).position() << endl;
		cout << "Node 2: " << tet.node(2).position() << endl;
		cout << "Node 3: " << tet.node(3).position() << endl;
	}
	cout << endl;
}
#endif

int main(int argc, char* argv[])
{
	// Check arguments
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
		exit(1);
	}

	cout << "Start" << endl;

  //Declare mesh
  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;

  while (CS207::getline_parsed(nodes_file, p)) {
	  //Add nodes
      mesh_node.push_back(mesh.add_node(p, NodeData()));
//      cout << p << endl;
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,4> t;
  while (CS207::getline_parsed(tris_file, t)) {
	//Initialize each tetrahedral's initial value to be the average of its nodes
    mesh.add_tetrahedral(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]], mesh_node[t[3]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_tetrahedral() << std::endl;

  // Initialization of mass and length
  // Set initial conditions for your nodes, if necessary.
  // Construct Forces/Constraints

  //Zero initial velocity
  //Set mass
  for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
	  auto n = *it;
	  n.value().velocity = Point(0,0,0);
	  n.value().mass = 1.0/mesh.num_nodes();
	}

  //To set rest length for all of the Edges to their initial length
  for (auto ei = mesh.edge_begin(); ei != mesh.edge_end(); ++ei ) {
	  (*ei).value().length = (*ei).length();
  }

  //To set initial Volume for all of the Tetrahedral to their initial Volume
  for (auto ti = mesh.tetrahedral_begin(); ti != mesh.tetrahedral_end(); ++ti ) {
	  (*ti).value().initialVolume = (*ti).volume();
  }


  // Launch the SDLViewer
  CS207::SDLViewer viewer;

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.launch();
  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 10.0;

  for (double t = t_start; t < t_end; t += dt) {
	  //std::cout << "t = " << t << std::endl;

	  //Setting Forces
	  GravityForce gravity(grav);
	  DashpotForce dashpot(100,0.01); // 1000 0.001
	  VolumePenaltyForce volumePenalty(&mesh, 1);

	  //Final Force
	  auto f = make_combined_force( gravity, dashpot, volumePenalty);
	  //Setting Constraints
	  HPlaneConstraint hplane(-2);
	  FixedConstraint fc( vector<Point>{Point(0,0,0), Point(1,0,0)} );

	  //Final Constraints
	  auto c = make_combined_constraints( hplane, fc, NullConstraint());

	  //Final Version - Add constraints
	  symp_euler_step(mesh, t, dt, f, c);

	  //printAllTetrahedralInformation(mesh, 0);

	  // Update viewer with nodes' new positions
	  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
	  viewer.set_label(t);

	  // These lines slow down the animation for small meshs, like grid0_*.
	  // Feel free to remove them or tweak the constants.
	  if (mesh.num_nodes() < 100)
		  CS207::sleep(0.001);
  }



#ifdef DEBUG_MESH
	cout << endl;
	cout << "Testing Tetrahedrals" << endl;
	printAllTetrahedralInformation(mesh, 0);

	cout << endl;
	cout << "Testing Edges" << endl;
	int i = 0;
	for(auto ei = mesh.edge_begin(); ei != mesh.edge_end(); ++ei) {
		cout << "Edge " << i << " : " <<
				(*ei).node(0).position() << ", " <<
				(*ei).node(1).position() << endl;
		++i;
	}

	cout << endl;
	cout << "Testing Nodes" << endl;
	int j = 0;
	for(auto ni = mesh.node_begin(); ni != mesh.node_end(); ++ni) {
		cout << "Node " << j << " : ";
		const MeshType::Node & n = (*ni);
		for(unsigned k = 0; k < n.nodeAdjTetrahedral().size(); ++k) {
			cout << n.nodeAdjTetrahedral().at(k).index() << ", ";
		}
		cout << endl;
		++j;
	}

	cout << "End" << endl;
#endif




  return 0;
}

