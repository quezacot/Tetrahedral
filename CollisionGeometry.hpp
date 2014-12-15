#include<iostream>
#include<string>
#include<cmath>
#include<cassert>
#include<vector>
#include<unordered_map>

#include "CS207/Point.hpp"
#include "BoundingBox.hpp"

// functions
bool equal(Point, Point);
bool equal(double, double, double epsilon=0.005);
Point plane_normal(Point, Point, Point);
bool on_same_side(Point, Point, Point, Point);
Point plane_line_intersect(Point, Point, Point, Point, Point);
bool is_inside_triangle(Point, Point, Point, Point);


/** Create a bounding box around set of objects positioned in space
 * @pre The type of *IT must implement the "position" concept; that
 * 	is the function (*IT).position() must be defined and return a
 *	Point
 * @pre @a first and @a last must define a valid range
 * @param[in] @a first is an iterator to the first object in the
 * 	range of objects
 * @param[in] @a last is an iterator to the last object in the
 * 	range
 * @returns a bounding box around the set of objects in the range
 * 	@a first to @a last
 */
template<typename IT>
BoundingBox build_bb(IT first, IT last) {
	BoundingBox b = BoundingBox();
	for(auto it = first; it != last; ++it) {
		b |= BoundingBox((*it).position());
	}
	return b;
}

/** Determines whether a line segment passes through triangle
 * @param[in] Triangle is defined by t1, t2, t3
 * @param[in] Line segment is defined by a and b
 * @returns true if the line segement intersects plane of triangle
 */
bool is_colliding(Point t1, Point t2, Point t3, Point p0,
					 Point p1) {
	// Declare variables
	Point n;
	double t;

	// r(t) = p0 + t*(p1 - p0) is the equation of the line; find
	// t at the point of intersection; t between 0 and 1 indicates
	// that the intersection point lies on the line segment
	n = plane_normal(t1, t2, t3);
	t = (dot(n, t1) - dot(p0, n)) / (dot(n, p1 - p0));
	if( 0 <= t <= 1 ) {
		Point intersect = plane_line_intersect(t1, t2, t3, p0, p1);
		return is_inside_triangle(t1, t2, t3, intersect);
	}
	else {
		return false;
	}
}

/** Returns true if a given point is inside the triangle
 * @param[in] Triangle defined by t1, t2, t3
 * @param[in] Point to check is p
 * @returns true is point is inside triangle, false otherwise
 */
bool is_inside_triangle(Point t1, Point t2, Point t3, Point p) {
	return on_same_side(t1, t2, t3, p) &&
		   on_same_side(t1, t3, t2, p) &&
		   on_same_side(t2, t3, t1, p);
}

/** Returns true if two points are on the same side of a line
 * @param[in] Line is defined by points a and b
 * @param[in] The two points are p0, and p1
 * @returns true if p0 and p1 are on the same side of line defined
 * 	by b - a
 */
bool on_same_side(Point a, Point b, Point p0, Point p1) {
	Point cross1 = cross(b - a, p0 - a);
	Point cross2 = cross(b - a, p1 - a);
	return dot(cross1, cross2) >= 0;
}


/* Return the area of the triangle formed by t1, t2, t3 */
double triangle_area(Point t1, Point t2, Point t3) {
	return 0.5 * norm(cross(t3 - t1, t3 - t2));
}

/** Compute the intersection of a plane and a line
 * @pre p0 and p1 must form a valid line, p0 != p1 &&
 * 	dot(p0, p1) != 0
 * @pre line formed by p0 and p1 is not parallel to the plane
 * @pre p0 and p1 must form a line not on the plane
 * @param[in] Line is defined by the two points p0 and p1
 * @param[in] Plane is defined by three points t1, t2, t3
 * @return a Point of intersection between the point and plane
 */
Point plane_line_intersect(Point t1, Point t2, Point t3,
						   Point p0, Point p1) {
	// Declare variables
	Point n;
	double t;

	n = plane_normal(t1, t2, t3);
	// Check preconditions
	assert( !equal(p0, p1) );
	assert( !equal(dot(p0, p1), 0) );
	assert( !equal(dot(n, p1 - p0), 0) );

	// r(t) = p0 + t*(p1 - p0) is the equation of the line; find
	// t at the point of intersection; then compute the intersect
	t = (dot(n, t1) - dot(p0, n)) / (dot(n, p1 - p0));
	return p0 + t * (p1 - p0);
}

/** Determines whether plane and line will intersect
 * @returns a bool if there is an intersection between line and plane
 */
bool is_plane_line_intersect(Point v1, Point v2, Point v3, Point a, Point b) {
	// Make sure that a and b form a valid line
	if( equal(a, b) || equal(dot(a, b), 0) )
		return false;
	// Make sure that v1, v2, v3 form a valid plane
	else if( equal(triangle_area(v1, v2, v3), 0) )
		return false;
	// Make sure that line is not parallel to the plane
	else if( equal(dot(plane_normal(v1, v2, v3), a-b), 0) )
		return false;
	else 
		return true;
}

/** Compute the normal to the plane defined by the three
 *	points p1, p2, and p3.
 */
Point plane_normal(Point p1, Point p2, Point p3) {
	Point v1 = p3 - p1;
	Point v2 = p2 - p1;
	return cross(v1, v2);
}

/* Return true if the point is on the plane defined by v1, v2, v3*/
bool is_on_plane(Point v1, Point v2, Point v3, Point p) {
	// Check if the vectors joining the point p to every vertex
	// on the plane dotted with the normal vector is 0
	Point n = plane_normal(v1, v2, v3);
	return equal(dot(n, v1 - p), 0) &&
		   equal(dot(n, v2 - p), 0) &&
		   equal(dot(n, v3 - p), 0);
}

/* Return true if the point is on the line defined by a and b */
bool is_on_line(Point a, Point b, Point p) {
	// Declare variables
	Point v1, v2, crpr;
	// If b and p jut out from a in the same or opposite
	// directions then p lies on the line formed by a and b
	v1 = b - a;
	v2 = p - a;
	crpr = cross(v1, v2);
	return equal(crpr, Point(0, 0, 0));
}

/* Returns true if the two doubles are within epsilon */
bool equal(double a, double b, double epsilon) {
	return std::abs(a - b) < epsilon;
}

/* Returns true if the two points x, y, z coords satisfy equal */
bool equal(Point a, Point b) {
	return equal(a.x, b.x) && equal(a.y, b.y) && equal(a.z, b.z);
}

/* returns area of tet defined by 4 nodes */
double tet_area(Point a, Point b, Point c, Point d) {
	return abs(dot(a-d, cross(b-d, c-d)))/6.0;
}
/* returns area of tet defined by 4 nodes */
template <typename N>
double tet_area(N a, N b, N c, N d) {
	return tet_area(a.position(), b.position(), c.position(), d.position());
}




