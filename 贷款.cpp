#include<iostream>
#include<cmath>
using namespace std;
class sd//�̴�
{
public:
	void time(int a);
	void calculate(double m, double n, double p,int b,double w);
	void calculate(double m, int n,double r);
	void calculate1(double m, double n, double p, int b, double w);
	void calculate1(double m, int n, double r);
	double zj=0;//�ܼ�
	double dj=0;//����
	double mj=0;//���
	double jc=0;//�ҳ�
	int jn=0;//����
	double nll=0;//������
	double yll=0;//������
	double sf=0;//�׸�
	double dk=0;//�����ܶ�
	double lx=0;//��Ϣ
	double hk=0;//�����ܶ�
	int mon=0;//��������
	double yhk;//ÿ�»����
	double ylx[3600] = { 0 };//ÿ�»�����Ϣ
	double ybj[3600] = { 0 };//ÿ�»����
	double sybj[3600] = { 0 };//ʣ�౾��
};
void sd::time(int n)
{
	mon = 12 * n;
	cout << "��������:" << mon << "��" << endl;
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
	cout << "�׸�:" << sf << "Ԫ" << endl;
	cout << "�����ܶ�:" << dk << "Ԫ" << endl;
	cout << "֧����Ϣ:" << lx << "Ԫ" << endl;
	cout << "�����ܶ�:" << hk << "Ԫ" << endl;
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

	cout << "�����ܶ�:" << dk << "Ԫ" << endl;
	cout << "��������:" << mon<< "��" << endl;
	cout << "֧����Ϣ:" << lx << "Ԫ" << endl;
	cout << "�����ܶ�:" << hk << "Ԫ" << endl;
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
	cout << "�׸�:" << sf << "Ԫ" << endl;
	cout << "�����ܶ�:" << dk << "Ԫ" << endl;
	cout << "֧����Ϣ:" << lx << "Ԫ" << endl;
	cout << "�����ܶ�:" << hk << "Ԫ" << endl;
	cout << "��������:" << mon << "��" << endl;
}
void sd::calculate1(double m, int n, double r)
{
	mon = 12 * n;
	yll = r * 1.0 / 12 * 0.01;
	dk = m;
	lx = (mon + 1) * dk * yll * 1.0 / 2;
	hk = lx + dk;
	cout << "�����ܶ�:" << dk << "Ԫ" << endl;
	cout << "��������:" << mon << "��" << endl;
	cout << "֧����Ϣ:" << lx << "Ԫ" << endl;
	cout << "�����ܶ�:" << hk << "Ԫ" << endl;
}
class gd:public sd//����
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
class zd//���
{
public:
	double yll2;
	double dk2;
	double lx2;
	double hk2;
	double mon2;
	double yhk2;//ÿ�»����
	double yll3;
	double dk3;
	double lx3;
	double hk3;
	double yhk3;//ÿ�»����
	double ylx2[3600] = { 0 };//ÿ�»�����Ϣ
	double ybj2[3600] = { 0 };//ÿ�»����
	double sybj2[3600] = { 0 };//ʣ�౾��
	double ylx3[3600] = { 0 };//ÿ�»�����Ϣ
	double ybj3[3600] = { 0 };//ÿ�»����
	double sybj3[3600] = { 0 };//ʣ�౾��
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
	cout << "�����ܶ�:" << a + b << "Ԫ" << endl;
	cout << "֧����Ϣ:" << lx3 << "Ԫ" << endl;
	cout << "������:" << hk3 << "Ԫ" << endl;
	cout << "��������:" << mon2 << "Ԫ" << endl;
}
void zd::calculate3(double a, double b, double c, double d, int y)
{
	mon2 = 12 * y;
	lx2 = (mon2 + 1) * (a*c*0.01/12+b*d*0.01/12) * 1.0 / 2;
	hk2 = lx2 + a + b;
	cout << "�����ܶ�:" << a + b << "Ԫ" << endl;
	cout << "֧����Ϣ:" << lx2 << "Ԫ" << endl;
	cout << "������:" << hk2 << "Ԫ" << endl;
	cout << "��������:" << mon2 << "Ԫ" << endl;
}
int main()
{
	sd s;
	gd m;
	zd n;
	int choice;
	cout << "��ѡ����ʽ��1����ҵ���2����������3������ʹ����"; 
	cin >> choice;
	if (choice == 1)
	{
		int choice1;
		cout << "��ѡ�񻹿ʽ��1���ȶϢ��2���ȶ��";
		cin >> choice1;
		if (choice1 == 1)
		{
			int choice2;
			cout << "��ѡ����ʽ��1��������������ۼ��㣬2�����ݴ����ܶ���㣩";
			cin >> choice2;
			if (choice2 == 1)
			{
				double a;//����
				double b;//���
				double c;//�ҳ�
				int d;//����
				double e;//������
				cin >> a >> b >> c >> d >> e;
				s.calculate(a, b, c, d, e);
			}
			else if (choice2 == 2)
			{
				double f;//�����ܶ�
				int year;//������
				double yell;//������
				cin >> f >> year >> yell;
				s.calculate(f, year, yell);
			}
			else
				return 0;
		}
		else if (choice1 == 2)
		{
			int choice3;
			cout << "��ѡ����ʽ��1��������������ۼ��㣬2�����ݴ����ܶ���㣩";
			cin >> choice3;
			if (choice3 == 1)
			{
				double a;//����
				double b;//���
				double c;//�ҳ�
				int d;//����
				double e;//������
				cin >> a >> b >> c >> d >> e;
				s.calculate1(a, b, c, d, e);
				s.time(d);
			}
			else if (choice3 == 2)
			{
				double f;//�����ܶ�
				int year;//������
				double yell;//������
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
		cout << "��ѡ�񻹿ʽ��1���ȶϢ��2���ȶ��";
		cin >> choice4;
		if (choice4 == 1)
		{
			int choice5;
			cout << "��ѡ����ʽ��1��������������ۼ��㣬2�����ݴ����ܶ���㣩";
			cin >> choice5;
			if (choice5 == 1)
			{
				double a1;//����
				double b1;//���
				double c1;//�ҳ�
				int d1;//����
				double e1;//������
				cin >> a1 >> b1 >> c1 >> d1 >> e1;
				m.d(a1, b1, c1, d1, e1);
			}
			else if (choice5 == 2)
			{
				double f1;//�����ܶ�
				int year1;//������
				double yell1;//������
				cin >> f1 >> year1 >> yell1;
				m.s(f1, year1, yell1);
			}
			else
				return 0;
		}
		else if (choice4 == 2)
		{
			int choice6;
			cout << "��ѡ����ʽ��1��������������ۼ��㣬2�����ݴ����ܶ���㣩";
			cin >> choice6;
			if (choice6== 1)
			{
				double a2;//����
				double b2;//���
				double c2;//�ҳ�
				int d2;//����
				double e2;//������
				cin >> a2 >> b2 >> c2 >> d2 >> e2;
				m.b(a2, b2, c2, d2, e2);
			}
			else if (choice6 == 2)
			{
				double f2;//�����ܶ�
				int year2;//������
				double yell2;//������
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
	cout << "��ѡ����ʽ��1���ȶϢ��2���ȶ��";
	cin >> choice7;
	if (choice7 == 1)
		n.calculate2(x, y, h, k, year3);
	else
		n.calculate3(x, y, h, k, year3);
	}
	else
		return 0;
}