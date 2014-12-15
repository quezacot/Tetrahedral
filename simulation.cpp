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
#include "Graph.hpp"
#include "CollisionDetector.hpp"


//DEBUG
//#define TETRAHEDRAL_OUTPUT
//#define DEBUG_MESH
//#define SurfaceColor
//#define MESH2

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;
struct NodeData;
struct EdgeData;
struct TetrahedralData;
typedef Mesh<NodeData, EdgeData, TetrahedralData> MeshType;
typedef typename MeshType::Node Node;
typedef typename MeshType::Edge Edge;
typedef CS207::SDLViewer<CS207::ViewerCallback> ViewerType;
typedef Graph<double, double> GraphType;

/** Custom structure of data to store with Nodes */
struct NodeData {
	Point velocity;  //< Node velocity
	double mass;     //< Node mass
	double color;    //< Node color
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
	MeshType* m_; // mesh pointer
	double K_; // mass spring constant
	/** VolumePenaltyForce Constructor.
	 * @param[in] m mesh pointer.
	 * @param[in] K mass spring constant.
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

/** Drag Force Functor that returns the force generated by mouse motion event
 */
struct DragForce {
  // The drag force
  Point dforce;
  // Multiply the force. Need tuning for different scenario.
  double coeff;

  /** Default DragForce Constructor.
   */
  DragForce(double coeff): dforce(Point(0,0,0)), coeff(coeff){}

  /** Calculates Gravity Force
   * @param[in] n Valid node.
   * @param[in] t Valid time.
   * @return Point object that represents the drag force generated by mouse motion.
   */
  Point operator()(Node n, double t) {
    (void) t;     // silence compiler warnings
    //cout<<"DragForce: x:"<<dforce.x<<" y:"<<dforce.y<<" z:"<<dforce.z<<endl;
    return n.value().mass * dforce * coeff;
  }
};

/** Wind Force
 *  Code written by Tian Lan and Xide Xia
 */
struct WindForce {
	Point w;
	WindForce(Point wind): w(wind) {}

	template <typename NODE>
	Point operator()(NODE n, double t) {
		double c = 0.00004;
		auto normal = Point(0,0,0);

		auto AdjTetrahedral = n.nodeAdjTetrahedral();

		for(unsigned k = 0; k < AdjTetrahedral.size(); ++k) {
			MeshType::Tetrahedral tet = AdjTetrahedral[k];
			if(tet.isSurface() == false) {
				return Point(0,0,0);
			}
			Point tnorm;
			tnorm = cross((tet).node(0).position()-(tet).node(1).position(),
					(tet).node(0).position()-(tet).node(2).position());
			tnorm = tnorm/norm(tnorm);
			normal = normal + tnorm;
		}
		(void) t;
		return c*dot((w-n.value().velocity),normal)*normal;
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

/** Null Constraint. No real constraint in this functor
 */
struct NullConstraint {

	/** Null Constraint Setter
	 * @param[in] g Valid mesh.
	 * @param[in] t Valid time.
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
	/* Contructor of Fixed Constraint
	 * @param[in] v The vector of Points that you want to constrain
	 */
	FixedConstraint(const vector<Point>& v) : cpoints(v) {}

	/** Fixed Constraint Setter
	 * @param[in] g Valid mesh.
	 * @param[in] t Valid time.
	 * @post The velocity of Points in @a cpoints are 0.
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

/** Hold Constraint to stop objects moving when mouse button pressing on them
 */
struct HoldConstraint {
  // @a hold means the object is prevent from moving.
  // @a pressed means the mouse button event of this object is pressed but not yet released.
  bool hold, pressed;
  /** Default constructor of Hold Constraint.
   *  default setting of @a hold and @a pressed are false.
   */
  HoldConstraint() : hold(false), pressed(false) {}

  /** Fixed Constraint Setter
   * @param[in] g Valid mesh.
   * @param[in] t Valid time.
   * @post The velocity of all points of this object are 0 when @a hold == true.
   */
  void operator()(MeshType& m, double t) {
    (void) t;     // silence compiler warnings
    if( hold ){
      for(auto it = m.node_begin(); it != m.node_end(); ++it) {
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
				// Simulate friction to move the ball
				n.value().velocity.x *= 0;
				n.value().velocity.y *= 0;
			}
		}
	}
};

/** Box Constraint that constructs an impassable box
 */
struct BoxConstraint {
	double constraint_; // coordinate for the box
	double frictionCoef_; // friction coef

	/** BoxConstraint Constructor.
	 * @param[in] constraint Sets the coordinate to define the the box.
	 * @param[in] frictionCoef friction coefficient for the friction force exerted by the box
	 */
	BoxConstraint(double constraint, double frictionCoef) : constraint_(constraint), frictionCoef_(frictionCoef) {}

	/** Horizontal Constraint Setter
	 * @param[in] m Valid mesh.
	 * @param[in] t Valid time.
	 */
	void operator()(MeshType& m, double t) {
		(void) t;     // silence compiler warnings
		for(auto it = m.node_begin(); it != m.node_end(); ++it) {
			Node n = (*it);

			// Constraint on z axis
			if(n.position().z > constraint_) {
				n.position().z = constraint_;
				n.value().velocity.z = 0;
				// Simulate friction to move the ball
				n.value().velocity.x *= frictionCoef_;
				n.value().velocity.y *= frictionCoef_;
			}
			else if (n.position().z < -constraint_) {
				n.position().z = -constraint_;
				n.value().velocity.z = 0;
				// Simulate friction to move the ball
				n.value().velocity.x *= frictionCoef_;
				n.value().velocity.y *= frictionCoef_;
			}
			// Constraint on x axis
			if(n.position().x > constraint_) {
				n.position().x = constraint_;
				n.value().velocity.x = 0;
				// Simulate friction to move the ball
				n.value().velocity.y *= frictionCoef_;
				n.value().velocity.z *= frictionCoef_;
			}
			else if (n.position().x < -constraint_) {
				n.position().x = -constraint_;
				n.value().velocity.x = 0;
				// Simulate friction to move the ball
				n.value().velocity.y *= frictionCoef_;
				n.value().velocity.z *= frictionCoef_;
			}
			// Constraint on y axis
			else if(n.position().y > constraint_) {
				n.position().y = constraint_;
				n.value().velocity.y = 0;
				// Simulate friction to move the ball
				n.value().velocity.x *= frictionCoef_;
				n.value().velocity.z *= frictionCoef_;
			}
			else if(n.position().y < -constraint_) {
				n.position().y = -constraint_;
				n.value().velocity.y = 0;
				// Simulate friction to move the ball
				n.value().velocity.x *= frictionCoef_;
				n.value().velocity.z *= frictionCoef_;
			}
		}
	}
};

/** Models collision between two balls
 * @param[in] Two valid constraints in @a c1, and @a c2.
 * @pre Valid g1, g2
 *
 * Code written by Tian Lan and Xide Xia
 */
template<typename M>
void CollisionConstraint(M& g1,M& g2, std::vector<unsigned> list1,std::vector<unsigned> list2){
  Point center1 = Point(0, 0, 0);
  for (auto it=g1.node_begin(); it != g1.node_end(); ++it){
    center1 += (*it).position()/g1.num_nodes();
  }
  Point center2 = Point(0, 0, 0);
  for (auto it=g2.node_begin(); it != g2.node_end(); ++it){
    center2 += (*it).position()/g2.num_nodes();
  }
  Point center0 = (center1+center2)/2;
  Point n1 = (center0-center1)/norm(center0-center1);
  Point n2 = (center0-center2)/norm(center0-center2);

  for (auto it1 = list1.begin(); it1 != list1.end(); ++it1){
    Node node = g1.node(*it1);
    Point p = node.position();
    Point p1 = center0-p;
    Point v = node.value().velocity;
    Point v1 = center0-v;
    node.position() = dot(n1,p1)*n1+p;
    node.value().velocity = v -dot(n1,v1)*(-n1);
  }

  for (auto it2 = list2.begin(); it2 != list2.end(); ++it2){
    Node node = g2.node(*it2);
    Point p = node.position();
    Point p1 = center0-p;
    Point v = node.value().velocity;
    Point v1 = center0-v;
    node.position() = dot(n2,p1)*n2+p;
    node.value().velocity = v -dot(n2,v1)*(-n2);
  }
}

// Color Functor
struct makePatterns {
  float longestPath_ = 1.0;
  makePatterns(const float& longestPath) : longestPath_(longestPath) {};

  CS207::Color operator()(const MeshType::Node& n) const {
    if(int(n.position().z*100) % 4 == 0) {
      return CS207::Color::make_hsv(0.2,1,1);
    }
    else if (int(n.position().z*100) % 4 == 1) {
      return CS207::Color::make_hsv(0.4,1,1);
    }
    else if (int(n.position().z*100) % 4 == 2) {
      return CS207::Color::make_hsv(0.6,1,1);
    }
    else {
      return CS207::Color::make_hsv(0.9,1,1);
    }
  }
};

/** Creates a functor to color the ball with 2 colors.
 */
struct twoColor {
  CS207::Color operator()(const MeshType::Node& n) const {
    return CS207::Color::make_hsv(n.value().color,1,1);
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

/* Mouse left key pressed and released listener
 * Catch mouse button events to interact with meshes.
 */
class MouseLeftClickCallback : public CS207::ViewerCallback {
  private:
    MeshType* m_;
    ViewerType* viewer_;
    DragForce* dforce_;
    HoldConstraint* hc_;

  public:
    /** Constructor of MouseLeftClickCallback
     * @param[in] m the mesh that will be listened
     * @param[in] v the viewer that will interact
     * @param[in] df the drag force impact on the mesh @a m
     * @param[in] hc the hold constraint constrains the mesh @a m
     */
    MouseLeftClickCallback(MeshType& m, ViewerType& v, DragForce& df, HoldConstraint& hc):
                            m_(&m), viewer_(&v), dforce_(&df), hc_(&hc) {}
  /**Functor executed when mouse button event happens.
   * Monitor the pressing and releasing of mouse left key.
   * @param[in] SDL_EVENT The event that occurred.
   */
  void operator()(const SDL_Event& event) {
    // Mouse left key pressed.
    if (event.type == SDL_MOUSEBUTTONDOWN && event.button.button == SDL_BUTTON_LEFT) {
      if ( event.button.state == SDL_PRESSED ) {
        // Check if the position clicked is on the mesh @a m_. It is only approximated.
        Point newcenter(viewer_->refpoint(event.button.x, event.button.y));
        Point maxp((*m_->node_begin()).position()), minp((*m_->node_begin()).position());
        // find out the box contained the mesh @a m_
        for( auto nit = m_->node_begin(); nit != m_->node_end(); ++nit ){
          const Point& np = (*nit).position();
          if (np.x > maxp.x ) maxp.x = np.x;
          if (np.y > maxp.y ) maxp.y = np.y;
          if (np.z > maxp.z ) maxp.z = np.z;
          if (np.x < minp.x ) minp.x = np.x;
          if (np.y < minp.y ) minp.y = np.y;
          if (np.z < minp.z ) minp.z = np.z;
        }
        // Check if the clicked position is in this box.
        if( maxp.x > newcenter.x && maxp.y > newcenter.y && maxp.z > newcenter.z &&
            minp.x < newcenter.x && minp.y < newcenter.y && minp.z < newcenter.z){
          cout<<"catched!"<<endl;
          hc_->hold = true;
          hc_->pressed = true;
        }
      }
    } else if (event.type == SDL_MOUSEBUTTONUP && event.button.button == SDL_BUTTON_LEFT) {
      // Mouse left key released.
      // Stop the drag force and let the mesh @a m_ free to move.
      if ( event.button.state == SDL_RELEASED ) {
        //cout<<"release hold"<<endl;
        dforce_->dforce = Point(0,0,0);
        hc_->pressed = false;
        hc_->hold = false;
      }
    }
  }
};

/* Mouse left key motion listener
 * Catch mouse motion event to interact the meshes.
 */
class MouseLeftDragCallback : public CS207::ViewerCallback {
  private:
    ViewerType* viewer_;
    DragForce* dforce_;
    HoldConstraint* hc_;

  public:
    /** Constructor of MouseLeftDragCallback
     * @param[in] v the viewer that will interact
     * @param[in] df the drag force impact on the mesh @a m
     * @param[in] hc the hold constraint constrains the mesh @a m
     */

    MouseLeftDragCallback(ViewerType& v, DragForce& df, HoldConstraint& hc):
                          viewer_(&v), dforce_(&df), hc_(&hc) {}

    /**Functor executed when mouse motion event happens.
     * Monitor the mouse motion when mouse left key is pressed.
     * @param[in] SDL_EVENT The event that occurred.
     */
  void operator()(const SDL_Event& event) {
    if (event.type == SDL_MOUSEMOTION ) {
      if ( event.motion.state == SDL_BUTTON(1) ) {
        // Apply the force generated by mouse motion to the mesh. It is only approximated.
        if( hc_->pressed ){
          hc_->hold = false;
          dforce_->dforce = viewer_->moveinline(event.motion.xrel, event.motion.yrel);
        }
      }
    }
  }
};

/* Wind Listener
 * Code adapted from Ruitao Du and Yingzhuo Zhang
 */
struct Listener_Wind : public CS207::ViewerCallback {
  WindForce& wind; // Wind Force
  Point pre_level; // Wind force level
  double increment_; // Wind force increment

  // Wind listener constructor
  Listener_Wind(WindForce& w, double increment): wind(w), pre_level(wind.w), increment_(increment){}

  void operator()(const SDL_Event& event) {
    if (event.type == SDL_KEYDOWN) {
      // Stop Wind
      if (event.key.keysym.sym == SDLK_s) {
        std::cout<<"Stop Wind"<<std::endl;
        wind.w=Point(0,0,0);
        std::cout << "Wind F:" << wind.w << std::endl;
      }
      // Resume Wind
      if (event.key.keysym.sym == SDLK_w) {
        std::cout<<"Resume Wind"<<std::endl;
        wind.w=pre_level;
        std::cout << "Wind F:" << wind.w << std::endl;
      }
      // Keyboard d to increase wind
      if (event.key.keysym.sym == SDLK_d){
        std::cout<<"Increase Wind"<<std::endl;
        wind.w+=increment_;
        pre_level=wind.w;
        std::cout << "Wind F:" << wind.w << std::endl;
      }
      // Keyboard a to decrease wind
      if (event.key.keysym.sym == SDLK_a){
        std::cout<<"Decrease Wind"<<std::endl;
        wind.w-=increment_;
        pre_level=wind.w;
        std::cout << "Wind F:" << wind.w << std::endl;
      }
    }
  }
};

/** Creates a functor to color the skull to show where the shortest paths are.
 * @param[in] graph Valid graph
 * @param[in] length Length of each box's edge
 * @post Add 8 nodes to the graph and 12 edges to form a box to be displayed in the viewer.
 */
void addBox(GraphType& graph, double length) {
	// Add nodes
	GraphType::Node n0 = graph.add_node(Point(length,length,-length));
	GraphType::Node n1 = graph.add_node(Point(-length,length,-length));
	GraphType::Node n2 = graph.add_node(Point(-length,-length,-length));
	GraphType::Node n3 = graph.add_node(Point(length,-length,-length));
	GraphType::Node n4 = graph.add_node(Point(length,length,length));
	GraphType::Node n5 = graph.add_node(Point(-length,length,length));
	GraphType::Node n6 = graph.add_node(Point(-length,-length,length));
	GraphType::Node n7 = graph.add_node(Point(length,-length,length));

	// Add edges
	graph.add_edge(n0, n1);
	graph.add_edge(n1, n2);
	graph.add_edge(n2, n3);
	graph.add_edge(n3, n0);
	graph.add_edge(n4, n5);
	graph.add_edge(n5, n6);
	graph.add_edge(n6, n7);
	graph.add_edge(n7, n4);
	graph.add_edge(n0, n4);
	graph.add_edge(n1, n5);
	graph.add_edge(n2, n6);
	graph.add_edge(n3, n7);
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
	// Box edge length
	// Used for adding Box and Box constraint
	double boxLength = 4;

	// Declare graph
	GraphType graph;
	std::vector<typename GraphType::node_type> graph_node;
	addBox(graph, boxLength);

	// Declare mesh
	MeshType mesh;
	std::vector<typename MeshType::node_type> mesh_node;
#ifdef MESH2
	MeshType mesh2;
	std::vector<typename MeshType::node_type> mesh_node2;
#endif
	// Read all Points and add them to the Mesh
	std::ifstream nodes_file(argv[1]);
	Point p;

	// Add in first mesh
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
#ifdef MESH2
	// Add in second mesh
	std::ifstream nodes_file2(argv[1]);
	while (CS207::getline_parsed(nodes_file2, p)) {
		p = p - Point(0,0,3);
		//Add nodes
		mesh_node2.push_back(mesh2.add_node(p, NodeData()));
	}

	// Read all mesh triangles and add them to the Mesh
	std::ifstream tris_file2(argv[2]);
	std::array<int,4> t2;
	while (CS207::getline_parsed(tris_file2, t2)) {
		//Initialize each tetrahedral's initial value to be the average of its nodes
		mesh2.add_tetrahedral(mesh_node2[t2[0]], mesh_node2[t2[1]], mesh_node2[t2[2]], mesh_node2[t2[3]]);
	}
#endif
	// Print out the stats
	std::cout << mesh.num_nodes() << " "
			<< mesh.num_edges() << " "
			<< mesh.num_tetrahedral() << std::endl;

	// Initialization of mass and length
	// Set initial conditions for your nodes, if necessary.

	// Initialize Colors
	for (auto it = mesh.node_begin(); it != mesh.node_end(); ++it) {
		auto n = *it;
		if(n.position().z > 0) {
			n.value().color = .8; // .8
		}
		else{
			n.value().color = .4; //.4
		}
	}

#ifdef SurfaceColor
	for (auto tIt = mesh.tetrahedral_begin(); tIt != mesh.tetrahedral_end(); ++tIt) {
		if((*tIt).isSurface() == true) {
			for (unsigned i = 0; i < 4; ++i) {
				(*tIt).node(i).value().color = 1;
			}
		}
	}
#endif

	// Initialize velocity and mass
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

#ifdef MESH2
	// Initialize mesh 2
	// Initialize velocity and mass
	for (auto it = mesh2.node_begin(); it != mesh2.node_end(); ++it) {
		auto n = *it;
		n.value().velocity = Point(0,0,0);
		n.value().mass = 1.0/mesh.num_nodes();
	}

	//To set rest length for all of the Edges to their initial length
	for (auto ei = mesh2.edge_begin(); ei != mesh2.edge_end(); ++ei ) {
		(*ei).value().length = (*ei).length();
	}

	//To set initial Volume for all of the Tetrahedral to their initial Volume
	for (auto ti = mesh2.tetrahedral_begin(); ti != mesh2.tetrahedral_end(); ++ti ) {
		(*ti).value().initialVolume = (*ti).volume();
	}
#endif


  // Launch the SDLViewer
  ViewerType viewer;

	//Initiate color functors
	makePatterns colorFunctor = makePatterns(1);

	// Add in first mesh
	auto node_map = viewer.empty_node_map(mesh);
	viewer.add_nodes(mesh.node_begin(), mesh.node_end(), colorFunctor, node_map);
	viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

#ifdef MESH2
	// Add in second mesh
	auto node_map2 = viewer.empty_node_map(mesh2);
	viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), node_map2);
	viewer.add_edges(mesh2.edge_begin(), mesh2.edge_end(), node_map2);
#endif

	// Add in graph nodes
	auto node_map_graph = viewer.empty_node_map(graph);
	viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map_graph);
	viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map_graph);

  // Setting Forces
  GravityForce gravity(grav);
#ifndef MESH2
  DashpotForce dashpot(1000,0.1); // 1000 0.001 // one ball
  VolumePenaltyForce volumePenalty(&mesh, 1);// one ball
#endif
  WindForce wind(Point(0,-1000,0));
  DragForce dforce(boxLength*2);
#ifdef MESH2
  DragForce dforce2(boxLength*2);
  DashpotForce dashpot(3000,0.1); // 1000 0.001 //two balls
  VolumePenaltyForce volumePenalty(&mesh, 400);// two balls
#endif

  // Setting Constraints
  HPlaneConstraint hplane(-3);
  //FixedConstraint fc( vector<Point>{Point(0,0,0), Point(2,0,0), Point(1,1.732,0)} );
  //HPlaneConstraint hplane(-3);
  BoxConstraint boxconstraint(boxLength,0.1);
  HoldConstraint hc1;
#ifdef MESH2
  HoldConstraint hc2;
#endif

  // Add Listeners
  MouseLeftClickCallback mlp(mesh, viewer, dforce, hc1);
  viewer.register_callback(&mlp);
  MouseLeftDragCallback mld(viewer, dforce, hc1);
  viewer.register_callback(&mld);
#ifdef MESH2
  MouseLeftDragCallback mld2(viewer, dforce2, hc2);
  viewer.register_callback(&mld2);
  MouseLeftClickCallback mlp2(mesh2, viewer, dforce2, hc2);
  viewer.register_callback(&mlp2);
#endif
  // Initialize Wind Listener
  Listener_Wind lwind(wind, 1000);
  viewer.register_callback(&lwind);

  // Launch viewer
	viewer.launch();
	viewer.center_view();

	// Begin the mass-spring simulation
	double dt = 0.001;
	double t_start = 0.0;
	double t_end   = 100.0;

	for (double t = t_start; t < t_end; t += dt) {
		//std::cout << "t = " << t << std::endl;

	  // Combine forces
	  auto sub_f = make_combined_force( gravity, dashpot, wind );
	  auto f1 = make_combined_force(sub_f, volumePenalty, dforce );
    // Combine constraints
	  auto c1 = make_combined_constraints( hc1, boxconstraint, NullConstraint());

		symp_euler_step(mesh, t, dt, f1, c1);

#ifdef MESH2
    auto f2 = make_combined_force(sub_f, volumePenalty, dforce2 );
    auto c2 = make_combined_constraints( hc2, boxconstraint, NullConstraint());
		symp_euler_step(mesh2, t, dt, f2, c2);
#endif

		//printAllTetrahedralInformation(mesh, 0);

#ifdef MESH2
		// Collision detection
		CollisionDetector<MeshType> col;
		col.add_object(mesh);
		col.add_object(mesh2);
		col.check_collisions();
		std::vector<unsigned> collision;
		std::vector<unsigned> collision2;

		for (auto it=col.begin(); it!= col.end(); ++it){
			auto boom = *it;
			Node n = boom.n1;
			if (boom.mesh1 == &mesh)
				collision.push_back(n.index());
			if (boom.mesh1 == &mesh2)
				collision2.push_back(n.index());
		}

		// Call Collision Detection
		CollisionConstraint(mesh, mesh2, collision, collision2);
#endif

		// Update viewer with nodes' new positions
		viewer.add_nodes(mesh.node_begin(), mesh.node_end(), twoColor(), node_map);

#ifdef MESH2
		viewer.add_nodes(mesh2.node_begin(), mesh2.node_end(), node_map2);
#endif
		viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map_graph);

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

