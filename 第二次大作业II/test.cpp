#include<complex>
#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>
#include<fstream>

using namespace std;
typedef complex<double> c;
constexpr double PI = 3.1415926535897932384;
constexpr double tf = 882.5598585646;
constexpr double A0 = 2.65094122459;
constexpr double omega = 0.0142385476661;
constexpr double Ip = 0.5;
//S
c S(c t, double p) {
	return 3 * pow(A0, 2)*t / 32. + pow(p, 2)*t / 2. + Ip * t +
		pow(A0, 2) / (960 * omega)*(-240.*sin(omega*t / 2.) +
			40.*sin(3 * omega*t / 2.) + 24.*sin(5 * omega*t / 2.) +
			15.*sin(omega*t) - 45.*sin(2 * omega*t) -
			5.*sin(3 * omega*t))
		+ (160 * A0*p / (960 * omega))*(3.*cos(omega*t / 2.) + cos(3 * omega*t / 2. ) - 3. * cos(omega*t)) - A0 * p / 6. / omega;
}

//S的一阶导数(鞍点方程)
c S_1(c t,double p) {
	return 0.5*pow((p + 2.65094122459*sin(0.0142385476661*t)*pow(sin(0.0142385476661*t / 4.), 2)), 2) + 0.5;
}

//S的二阶导数
c S_2(c t, double p) {
	return A0 * omega*sin(omega*t / 4.)*(p + A0 * sin(omega*t)*pow(sin(omega*t / 4.), 2))*(cos(omega*t)*sin(omega*t / 4.) + 0.5*sin(omega*t)*cos(omega*t / 4.));
}


//对鞍点方程的muller算法
//precision=精度,N0=最大迭代次数,x0,x1,x2为初始位置,p是鞍点方程参数p
c Muller(double precision,int N0,c x0, c x1, c x2, double p) {
	c h0 = x0 - x2;
	c h1 = x1 - x2;
	//割线斜率delta0 delta1
	c delta0 = (S_1(x0,p) - S_1(x2,p)) / h0;
	c delta1 = (S_1(x1,p) - S_1(x2,p)) / h1;
	//拟合抛物线的参数d
	c d = (delta0 - delta1) / (h0 - h1);
	c E,b,D,h,p0;
	for (int i = 0; i < N0; i++) {
		//计算近似解p0
		b = delta0 - h0 * d;
		D = sqrt(b*b - 3. * S_1(x2,p)*d);
		if (abs(b - D) < abs(b + D)) {
			E = b + D;
		}
		else {
			E = b - D;
		}
		h = -2. * S_1(x2,p) / E;
		p0 = x2 + h;
		//判断精度是否达到
		if (abs(S_1(p0,p)) < precision) {
			return p0;
		}
		//更新拟合的参数,进行下一次循环
		x0 = x1;
		x1 = x2;
		x2 = p0;
		h0 = x0 - x2;
		h1 = x1 - x2;
		delta0 = (S_1(x0,p) - S_1(x2,p)) / h0;
		delta1 = (S_1(x1,p) - S_1(x2,p)) / h1;
		d = (delta0 - delta1) / (h0 - h1);
	}
	//cout << "failed after" << N0 << "iterations"<<endl;
}


//解第(1)问方程
vector<c> solve_S(double precision, int N0, double p ,double deltat) {
	vector<c> List;
	double ty = 0;
	//在规定范围内搜索根
	while (List.end()-List.begin()<6)
	{
		for (double tx = 0; tx < tf&&List.end() - List.begin() < 6; tx += deltat) {
			c x0(tx, ty);
			c x1(tx + deltat, ty);
			c x2(tx, ty + deltat);
			//求出局部的根
			c root = Muller(precision, N0, x0, x2, x1, p);
			//判断这个根是否存在
			bool in_list = false;
			for (int i = 0; i < List.end() - List.begin(); i++) {
				if (abs(List[i] - root) < 100 * precision) {
					in_list = true;
				}
			}
			//判断这个根是否符合要求
			if (in_list == false&&root.imag()>0&&root.real()>=0&&root.real()<=tf) {
				List.push_back(root);
			}
		}
		ty += deltat;
	}
	return List;
}

//鞍点法电离几率幅(第(2)问)
c Mp_Saddle(double p) {
	vector<c> list = solve_S(1e-13, 20, p, 1);
	c result = 0;
	for (int i = 0; i < 6; i++) {
		c t = list.at(i);
		result += (cos(S(t, p)) + c(0, 1)*sin(S(t, p))) / S_2(t, p);
	}
	return -pow(2 * Ip, 5. / 4.) / sqrt(2)*result;
}

//电场E
c E(c t,double p) {
	return -A0 * omega*sin(omega*t / 4.)*(cos(omega*t)*sin(omega*t / 4.) + 0.5*sin(omega*t)*cos(omega*t / 4.));
}
//机械动量q
c q(c t,double p) {
	return p + A0 * sin(omega*t)*pow(sin(omega*t / 4.), 2);
}
//第(3)问被积函数I
c I(c t, double p) {
	return pow(2, 7. / 2.)*(2 * Ip)*q(t, p)*E(t, p) / PI / pow(2 * Ip + pow(q(t, p), 2), 3)*(cos(S(t, p)) + c(0, 1)*sin(S(t, p)));
}


//直接积分法电离几率幅(第(3)问)
//对I进行述职积分 用变步长积分,im是最大次数
c Mp_integral(double p, int im) {
	double max = tf;//max 是积分上限
	double delta = max;//delta是步长
	c T = (I(max,p) + I(0,p)) / 2.;//T是积分的数值近似值
	int n = 1;//n是区间个数
	c T_before;//记录前一次迭代值
	for (int i = 0; i < im; i++) {
		//利用公式T_2n=0.5T_n+h/2*(sigma x_k+0.5 from 0 to n-1)
		//先求sigma的那一项
		T_before = T;
		c sigma = 0;
		for (int j = 0; j < n; j++) {
			sigma += I(0.5*delta + j * delta,p);
		}
		T = T / 2. + delta * sigma / 2.;
		delta = delta / 2;
		n = n * 2;
		//判断精度是否符合要求
		if (abs(T - T_before) / abs(T) < 1e-8) {
			return T;
		}
	}
	return T ;
}


void test1() {
	vector<c> result = solve_S(1e-13, 70, 1, 1);
	for (int i = 0; i < 6; i++) {
		cout <<setprecision(10)<< result[i] << endl;
	}
}

void test2() {
	fstream os;
	os.open("temp2.txt");
	for (double p = 0.01; p < 2 + 0.01; p += 0.01) {
		cout << p<<endl;//!!!!!!
		os << p << "\t" << norm(Mp_Saddle(p)) << endl;
	}
}

void test3() {
	fstream os;
	os.open("temp3.txt");
	for (double p = 0.01; p < 2 + 0.01; p += 0.01) {
		cout << p << endl;//!!!!!!
		os << p << "\t" << norm(Mp_integral(p,20)) << endl;
	}
}

int main() {
	
	test2();
	system("pause");
}