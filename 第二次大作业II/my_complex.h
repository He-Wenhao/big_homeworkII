#pragma once
#include<math.h>
#include<iostream>
#include<iomanip>

struct my_complex {
	double x;
	double y;
	my_complex() :x{ 0 }, y{ 0 } {};
	my_complex(double _x, double _y=0) :x{ _x }, y{ _y } {};
};

my_complex conj(const my_complex & a);
double abs(const my_complex &a);
my_complex sqrt(const my_complex& a);

//operator
  my_complex operator+(const my_complex &a, const my_complex &b) {
	return my_complex(a.x + b.x, a.y + b.y);
}
  my_complex operator+(const my_complex &a, double b) {
	return my_complex(a.x + b, a.y);
}
  my_complex operator-(const my_complex &a, const my_complex &b) {
	return my_complex(a.x - b.x, a.y - b.y);
}
  my_complex operator-(const my_complex &a, double b) {
	  return my_complex(a.x - b, a.y);
  }
  my_complex operator*(const my_complex &a, const my_complex &b) {
	return my_complex(a.x*b.x - a.y*b.y, a.x * b.y + a.y*b.x);
}
  my_complex operator*(const double &a, const my_complex &b) {
	return my_complex(a*b.x,a*b.y);
}
my_complex operator/(const my_complex&a, const my_complex&b) {
	return (1 / (b*conj(b)).x)*a*conj(b);
}

std::ostream& operator<<(std::ostream&os, my_complex a) {
	os << std::setprecision(5) << a.x << std::showpos << a.y << "i";
	return os;
}

//conjugate
  my_complex conj(const my_complex & a) {
	return my_complex(a.x, -a.y);
}

//abslute value
double abs(const my_complex &a) {
	return sqrt((a*conj(a)).x);
}

//sqrt
  my_complex sqrt(const my_complex& a) {
	return my_complex(sqrt((a.x + sqrt(a.x*a.x + a.y*a.y)) / 2), sqrt((-a.x + sqrt(a.x*a.x + a.y*a.y)) / 2));
}
