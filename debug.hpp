#include "CS207/Point.hpp"
#include<iostream>
#include<string>

// basic debugs
template <typename S1, typename S2 = std::string>
void db(const S1 s1,const S2 s2 = "") {
	std::cout << s1 << " " << s2 << std::endl;
}
template <typename S1>
void db(const S1 s1, const Point p) {
	std::cout << s1 << " " << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}


/* Debug color
 * c: 30 - 37
 * black, red, green, yellow, blue, magenta, cyan, white
 */
template <typename S1>
void dbc(const S1 s1, int c) {
	std::cout << "\x1b[" << c << "m" << s1 << "\x1b[0m" << std::endl;
}
template <typename S1> // debug red
void dbr(const S1 s1) {
	std::cout << "\x1b[31m" << s1 << "\x1b[0m" << std::endl;
}
template <typename S1> // debug green
void dbg(const S1 s1) {
	std::cout << "\x1b[32m" << s1 << "\x1b[0m" << std::endl;
}

// collection debug
template <typename S, typename M>
void db_map(const S s, const M m) {
	db(s);
	for ( auto& x: m )
	    std::cout << x.first << ": " << x.second.index() << std::endl;
}
template <typename S, typename V>
void db_vec(const S s, const V v) {
	db(s);
	for ( auto& x: v)
	    std::cout << x << ", ";
	std::cout << std::endl;
}
