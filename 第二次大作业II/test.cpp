#include"my_complex.h"
#include<iostream>

using namespace std;
typedef my_complex c;

c equation(c t) {
	return t*t+2;
}

c Muller_sol_eq(double precision,int N0,c x0, c x1, c x2, c (*f)(c)) {
	c h0 = x0 - x2;
	c h1 = x1 - x2;
	c delta0 = (f(x0) - f(x2)) / h0;
	c delta1 = (f(x1) - f(x2)) / h1;
	c d = (delta0 - delta1) / (h0 - h1);
	c E,b,D,h,p;
	for (int i = 0; i < N0; i++) {
		b = delta0 - h0 * d;
		D = sqrt(b*b - 3 * f(x2)*d);
		if (abs(b - D) < abs(b + D)) {
			E = b + D;
		}
		else {
			E = b - D;
		}
		h = -2 * f(x2) / E;
		p = x2 + h;
		if (abs(h) < precision) {
			return p;
		}
		cout << p << endl;//!!!!!!!!!
		x0 = x1;
		x1 = x2;
		x2 = p;
		h0 = x0 - x2;
		h1 = x1 - x2;
		delta0 = (f(x0) - f(x2)) / h0;
		delta1 = (f(x1) - f(x2)) / h1;
		d = (delta0 - delta1) / (h0 - h1);
	}
	cout << "failed after" << N0 << "iterations"<<endl;
	return p;
}


//derivative on real axes
//利用 5 点中心差分 公式
c derive(double h,c (*f)(c),double x) {
	return (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / (12 * h);
}



int main() {
	c x0(0, 0), x1(0,-2), x2(2, 0);
	cout << Muller_sol_eq(0.0000001, 70, x0,x1,x2, equation)<<endl;
	system("pause");
}