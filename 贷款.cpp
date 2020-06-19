#include<iostream>
#include<cmath>
using namespace std;
class sd//商贷
{
public:
	void time(int a);
	void calculate(double m, double n, double p,int b,double w);
	void calculate(double m, int n,double r);
	void calculate1(double m, double n, double p, int b, double w);
	void calculate1(double m, int n, double r);
	double zj=0;//总价
	double dj=0;//单价
	double mj=0;//面积
	double jc=0;//揭成
	int jn=0;//揭年
	double nll=0;//年利率
	double yll=0;//月利率
	double sf=0;//首付
	double dk=0;//贷款总额
	double lx=0;//利息
	double hk=0;//还款总额
	int mon=0;//还款月数
	double yhk;//每月还款额
	double ylx[3600] = { 0 };//每月还款利息
	double ybj[3600] = { 0 };//每月还款本金
	double sybj[3600] = { 0 };//剩余本金
};
void sd::time(int n)
{
	mon = 12 * n;
	cout << "还款月数:" << mon << "月" << endl;
}
void sd::calculate(double m, double n, double p,int b,double w)
{
	nll = w*0.01;
	yll = nll*1.0 / 12;
	zj = m * n;
	sf = zj * 0.1 * (10 - p);
	dk = zj - sf;
	int i = 0;
	sybj[0] = dk;
	for (i = 0; i < 12 * b; i++)
	{
		yhk = dk*((yll * pow((1 + yll),b))*1.0/ (pow((1 + yll) , b) - 1));
		ylx[i] = sybj[i] * yll;
		ybj[i] = yhk - ylx[i];
		sybj[i + 1] = sybj[i] - ybj[i];
		lx = lx+ ylx[i];
	}
	hk = lx + dk;
	cout << "首付:" << sf << "元" << endl;
	cout << "贷款总额:" << dk << "元" << endl;
	cout << "支付利息:" << lx << "元" << endl;
	cout << "还款总额:" << hk << "元" << endl;
}
void sd::calculate(double m, int n,double r)
{
	mon = 12 * n;
	yll = r * 1.0 / 12*0.01;
	dk = m;
	sybj[0] = dk;
	int i = 0;
	for (i = 0; i < 12 * mon; i++)
	{
		yhk = dk * ((yll * pow((1 + yll), mon)) * 1.0 / (pow((1 + yll), mon) - 1));
		ylx[i] = sybj[i] * yll;
		ybj[i] = yhk - ylx[i];
		sybj[i + 1] = sybj[i] - ybj[i];
		lx = lx + ylx[i];
	}
	hk = lx + dk;

	cout << "贷款总额:" << dk << "元" << endl;
	cout << "还款月数:" << mon<< "月" << endl;
	cout << "支付利息:" << lx << "元" << endl;
	cout << "还款总额:" << hk << "元" << endl;
}

void sd::calculate1(double m, double n, double p, int b, double w)
{
	nll = w*0.01;
	yll = nll * 1.0 / 12;
	zj = m * n;
	sf = zj * 0.1 * (10 - p);
	dk = zj - sf;
	mon = 12 * b;
	lx = (mon + 1) * dk * yll * 1.0 / 2;
	hk = lx + dk;
	cout << "首付:" << sf << "元" << endl;
	cout << "贷款总额:" << dk << "元" << endl;
	cout << "支付利息:" << lx << "元" << endl;
	cout << "还款总额:" << hk << "元" << endl;
	cout << "还款月数:" << mon << "月" << endl;
}
void sd::calculate1(double m, int n, double r)
{
	mon = 12 * n;
	yll = r * 1.0 / 12 * 0.01;
	dk = m;
	lx = (mon + 1) * dk * yll * 1.0 / 2;
	hk = lx + dk;
	cout << "贷款总额:" << dk << "元" << endl;
	cout << "还款月数:" << mon << "月" << endl;
	cout << "支付利息:" << lx << "元" << endl;
	cout << "还款总额:" << hk << "元" << endl;
}
class gd:public sd//公贷
{
public:
	void d(double m , double n , double p, int b , double w );
	void b(double a, double b, double c, int d, double e);
	void s(double m, int n, double p);
	void p(double a, int b, double c);
};
void gd::d(double m , double n , double p , int b , double w )
{
	time(m);
	calculate(m, n, p, b, w);
}
void gd::b(double a, double b, double c, int d, double e)
{
	calculate1(a, b, c, d, e);
}
void gd::s(double m, int n, double p)
{
	calculate(m, n, p);
}
void gd::p(double a, int b, double c)
{
	calculate1(a, b, c);
}
class zd//组贷
{
public:
	double yll2;
	double dk2;
	double lx2;
	double hk2;
	double mon2;
	double yhk2;//每月还款额
	double yll3;
	double dk3;
	double lx3;
	double hk3;
	double yhk3;//每月还款额
	double ylx2[3600] = { 0 };//每月还款利息
	double ybj2[3600] = { 0 };//每月还款本金
	double sybj2[3600] = { 0 };//剩余本金
	double ylx3[3600] = { 0 };//每月还款利息
	double ybj3[3600] = { 0 };//每月还款本金
	double sybj3[3600] = { 0 };//剩余本金
	void calculate2(double a, double b, double c, double d, int y);
	void calculate3(double a, double b, double c, double d, int y);
};
void zd::calculate2(double a, double b, double c, double d, int y)
{
	mon2 = 12 * y;
	yll2 = c * 0.01 / 12;
	yll3 = d * 0.01 / 12;
	dk3 = b;
	dk2 = a;
	sybj2[0] = dk2;
	int i = 0;
	for (i = 0; i < 12 * mon2; i++)
	{
		yhk2 = dk2 * ((yll2 * pow((1 + yll2), mon2)) * 1.0 / (pow((1 + yll2), mon2) - 1));
		ylx2[i] = sybj2[i] * yll2;
		ybj2[i] = yhk2 - ylx2[i];
		sybj2[i + 1] = sybj2[i] - ybj2[i];
		lx2 = lx2 + ylx2[i];
		hk2 = dk2 + lx2;
	}
	sybj3[0] = dk3;
	for (i = 0; i < 12 * mon2; i++)
	{
		yhk3 = dk3 * ((yll3 * pow((1 + yll3), mon2)) * 1.0 / (pow((1 + yll3), mon2) - 1));
		ylx3[i] = sybj3[i] * yll3;
		ybj3[i] = yhk3 - ylx3[i];
		sybj3[i + 1] = sybj2[i] - ybj3[i];
		lx3 = lx3 + ylx3[i];
		hk3 =hk2+ dk3 + lx3;
	}
	lx3 = lx3 + lx2;
	cout << "贷款总额:" << a + b << "元" << endl;
	cout << "支付利息:" << lx3 << "元" << endl;
	cout << "还款金额:" << hk3 << "元" << endl;
	cout << "还款月数:" << mon2 << "元" << endl;
}
void zd::calculate3(double a, double b, double c, double d, int y)
{
	mon2 = 12 * y;
	lx2 = (mon2 + 1) * (a*c*0.01/12+b*d*0.01/12) * 1.0 / 2;
	hk2 = lx2 + a + b;
	cout << "贷款总额:" << a + b << "元" << endl;
	cout << "支付利息:" << lx2 << "元" << endl;
	cout << "还款金额:" << hk2 << "元" << endl;
	cout << "还款月数:" << mon2 << "元" << endl;
}
int main()
{
	sd s;
	gd m;
	zd n;
	int choice;
	cout << "请选择贷款方式（1：商业贷款，2：公积金贷款，3：组合型贷款）："; 
	cin >> choice;
	if (choice == 1)
	{
		int choice1;
		cout << "请选择还款方式（1：等额本息，2：等额本金）";
		cin >> choice1;
		if (choice1 == 1)
		{
			int choice2;
			cout << "请选择贷款方式（1：根据面积、单价计算，2：根据贷款总额计算）";
			cin >> choice2;
			if (choice2 == 1)
			{
				double a;//单价
				double b;//面积
				double c;//揭成
				int d;//揭年
				double e;//年利率
				cin >> a >> b >> c >> d >> e;
				s.calculate(a, b, c, d, e);
			}
			else if (choice2 == 2)
			{
				double f;//贷款总额
				int year;//揭年数
				double yell;//年利率
				cin >> f >> year >> yell;
				s.calculate(f, year, yell);
			}
			else
				return 0;
		}
		else if (choice1 == 2)
		{
			int choice3;
			cout << "请选择贷款方式（1：根据面积、单价计算，2：根据贷款总额计算）";
			cin >> choice3;
			if (choice3 == 1)
			{
				double a;//单价
				double b;//面积
				double c;//揭成
				int d;//揭年
				double e;//年利率
				cin >> a >> b >> c >> d >> e;
				s.calculate1(a, b, c, d, e);
				s.time(d);
			}
			else if (choice3 == 2)
			{
				double f;//贷款总额
				int year;//揭年数
				double yell;//年利率
				cin >> f >> year >> yell;
				s.calculate1(f, year, yell);
			}
			else
				return 0;
		}
		else
			return 0;
	}
	else if (choice == 2)
	{
		int choice4;
		cout << "请选择还款方式（1：等额本息，2：等额本金）";
		cin >> choice4;
		if (choice4 == 1)
		{
			int choice5;
			cout << "请选择贷款方式（1：根据面积、单价计算，2：根据贷款总额计算）";
			cin >> choice5;
			if (choice5 == 1)
			{
				double a1;//单价
				double b1;//面积
				double c1;//揭成
				int d1;//揭年
				double e1;//年利率
				cin >> a1 >> b1 >> c1 >> d1 >> e1;
				m.d(a1, b1, c1, d1, e1);
			}
			else if (choice5 == 2)
			{
				double f1;//贷款总额
				int year1;//揭年数
				double yell1;//年利率
				cin >> f1 >> year1 >> yell1;
				m.s(f1, year1, yell1);
			}
			else
				return 0;
		}
		else if (choice4 == 2)
		{
			int choice6;
			cout << "请选择贷款方式（1：根据面积、单价计算，2：根据贷款总额计算）";
			cin >> choice6;
			if (choice6== 1)
			{
				double a2;//单价
				double b2;//面积
				double c2;//揭成
				int d2;//揭年
				double e2;//年利率
				cin >> a2 >> b2 >> c2 >> d2 >> e2;
				m.b(a2, b2, c2, d2, e2);
			}
			else if (choice6 == 2)
			{
				double f2;//贷款总额
				int year2;//揭年数
				double yell2;//年利率
				cin >> f2 >> year2 >> yell2;
				m.p(f2, year2, yell2);
			}
			else
				return 0;
		}
		else
			return 0;
	}
	
	else if (choice == 3)
	{
	double x,y,h,k;
	int year3;
	cin >> x >> y >> h >> k >> year3;
	int choice7;
	cout << "请选择贷款方式（1：等额本息，2：等额本金）";
	cin >> choice7;
	if (choice7 == 1)
		n.calculate2(x, y, h, k, year3);
	else
		n.calculate3(x, y, h, k, year3);
	}
	else
		return 0;
}