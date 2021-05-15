#include "pch.h"
#include <iostream>
#include<stdio.h>
#include<cmath>
using namespace std;
bool isComposite[1000005];
int primeNumber;
int primeSet[1000005];
int stimulateNumber; //模拟光子的次数
double xrayNumber;
double XrayEnergy[10][5][5]; //存储被激发物质的特征X射线能量信息,p[i][0]为能量值,p[i][1]为比例
double R = 1382, H = 500; //探测物质底面半径和高
double standE[20] = {1, 1.5, 1.839, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60}; //单位:keV
double lnE[20];
double pAbsorption[20] = { 1.567e3, 5.333e2, 3.071e2, 2.775e3, 9.767e2, 4.514e2, 2.438e2, 1.458e2, 63.79, 33.15, 9.848, 4.089, 1.161, 0.4687, 0.2307, 0.1288 };
double cAbsorption[20] = { 1.317e-2, 2.393e-2,3.081e-2,3.388e-2, 4.962e-2, 6.135e-2, 7.110e-2, 7.983e-2, 9.507e-2, 0.1076, 0.1289, 0.1402, 0.1501, 0.1534, 0.1538, 0.1526 };
double plnAbsorption[20], clnAbsorption[20]; //两种效应的截面取对数后的值
long long randNum;
int channel[1030];
const int Cu = 0, Rb = 1, Mo = 2, Ag = 3, Ba = 4, Tb = 5;
const int ka1 = 0, ka2 = 1, kb1 = 2, kb2=3;
const double pi = 3.1415926; //圆周率
const double h = 6.626e-34; //普朗克常数
const double e = 1.602e-19; //电子电荷
const double c = 3e8; //光速(m/s)
const double density = 2.33; //Si密度(g/cm^3)
const double resolution = 139.0 / 5900; //Si(PIN)探测器分辨率
const double CmToUm = 1e4; //cm到um的转换系数
const int Photoelectric = 1;
const int Compton = 0;
const double electronComptonWavelength = 2.426e-11; //电子的康普顿波长(m)
const long long BaseNumber = 1ll << 42;
int scatterCount[100] = { 0 }; //测试用，统计散射次数
void getPrimeSet()
{
	for (int i = 2; i <= 1000000; i++)
		for (long long j = 1ll * i * i; j <= 1000000; j += i)
			isComposite[j] = 1;
	for (int i = 2; i <= 1000000; i++)
		if (!isComposite[i])
			primeSet[primeNumber++] = i;
}
void getlnAbsorption()
{
	for (int i = 0; i < 16; i++)
	{
		lnE[i] = log(standE[i]);
		plnAbsorption[i] = log(pAbsorption[i]);
		clnAbsorption[i] = log(cAbsorption[i]);
	}
}
/****生成一个新的随机数****/
double getRandomNumber()
{
	for (int i = 0; i < 17; i++)
		randNum = randNum * 5 % BaseNumber;
	return 1.0 * randNum / BaseNumber;
}
/****插值法计算截面，输入（能量、反应类型），输出（反应截面，单位:cm^2/g）****/
double getSection(double E, int reactionType)
{
	int i = 0;
	for (; E > standE[i]; i++);
	i--;
	double A1, A2;
	if (reactionType == Photoelectric)
		A1 = plnAbsorption[i], A2 = plnAbsorption[i + 1];
	else if (reactionType == Compton)
		A1 = clnAbsorption[i], A2 = clnAbsorption[i + 1];
	double E2 = lnE[i + 1], E1 = lnE[i], E0 = log(E);
	double logAns = (A1 * (E2 - E0) + A2 * (E0 - E1)) / (E2 - E1);
	return exp(logAns);
}
/****判断光子是否在探测物质内****/
bool isInMaterial(double *r)
{
	return r[0] * r[0] + r[1] * r[1] <= R && r[2] >= 0 && r[2] <= H;
}
/********返回向量相加r+e*len的结果********/
void addVector(double *r, double *e, double len)
{
	for (int i = 0; i < 3; i++)
		r[i] += e[i] * len;
}
/********叉乘********/
void CrossProduct(double *z, double *x, double *y)
{
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}
/********光子在后表面的反射********/
void reflection(double *r)
{
	r[2] -= 2 * (r[2] - H);
}
/********给能量(keV)算波长(m)********/
double getwavelength(double E)
{
	return h * c / (E * e * 1000);
}
/********给波长(m)算能量(keV)********/
double getenergy(double W)
{
	return h * c / (W*e * 1000);
}
double calculusForCompton(double a, double b)
{
	return (b*b - a * a) / (2 * a) + a * log(b / a) - (b - a)*(b - a) + (b - a)*(b - a)*(b - a) / 3;
}
double differentialForCompton(double x0, double x)
{
	return x / x0 + x0 / x + 2 * (x0 - x) + (x - x0)*(x - x0);
}
/********计算发生康普顿散射后的 能量 & 方向********/
double GaussExtend(double u0)
{
	double x1 = getRandomNumber();
	double x2 = getRandomNumber();
	double y = sqrt(-2 * log(x1))*cos(2 * pi*x2); 
	return u0 + u0 * y*resolution / 3;
}
double ComptonEffect(double E, double *direction)
{
	/********get wavelength(Wn1) after compton scatter********/
	double e0 = getRandomNumber(); //cout << e0 << endl;
	double Wn = getwavelength(E) / electronComptonWavelength;
	double x, x1, y1, y2;
	x = Wn;
	do
	{
		x1 = x;
		y1 = calculusForCompton(Wn, x) - e0 * calculusForCompton(Wn, Wn + 2);
		y2 = differentialForCompton(Wn, x1);
		x = x1 - y1 / y2; //cout << x << ' ' << x1 << endl;
	} while (fabs(x - x1) >= 1e-5);
	double Wn1 = x1;
	/********get new direction of photon after compton scatter********/
	double newAngle = acos(Wn + 1 - Wn1);
	double xnew = getRandomNumber();
	double ynew = getRandomNumber();
	double znew = getRandomNumber();
	double u[3] = { xnew,ynew,znew };
	double uv[3];
	CrossProduct(uv, u, direction);
	double length = 0;
	for (int i = 0; i < 3; i++)
		length += uv[i] * uv[i];
	for (int i = 0; i < 3; i++)
		uv[i] = uv[i] / length;
	double newDirection[3];
	double alpha = Wn + 1 - Wn1;
	double belta = sqrt(1 - alpha * alpha);
	for (int i = 0; i < 3; i++)
		newDirection[i] = direction[i] * alpha + uv[i] * belta;
	//cout << alpha << ' ' << direction[0] * newDirection[0] + direction[1] * newDirection[1] + direction[2] * newDirection[2] << endl;
	for (int i = 0; i < 3; i++)
			direction[i] = newDirection[i];
	/*double angel1 = pi * getRandomNumber();
	double angel2 = 2 * pi * getRandomNumber();
	double u[3] = { sin(angel1) * cos(angel2) , sin(angel1)*sin(angel2), cos(angel1) };
	double newDirection[3];
	CrossProduct(newDirection, direction, u);
	double length = 0;
	for (int i = 0; i < 3; i++)
		length += newDirection[i] * newDirection[i];
	for (int i = 0; i < 3; i++)
		newDirection[i] = newDirection[i] / sqrt(length) * tan(newAngle);
	for (int i = 0; i < 3; i++)
		direction[i] += newDirection[i];*/
	/********返回散射后的能量值********/
	return getenergy(Wn1 * electronComptonWavelength);
}
/********单个光子的运输模拟********/
void singlePhotonStimulate(double energy, double *r, double *direction)
{
	double DepositEnergy = 0; //沉积的能量
	double e0 = 0; //随机数
	int scatter = 0;
	while(true)
	{
		double pAbsorb = density*getSection(energy, Photoelectric); //光电效应吸收系数(cm^-1)
		double cAbsorb = density*getSection(energy, Compton); //康普顿效应吸收系数(cm^-1)
		double totalAbsorb = pAbsorb + cAbsorb; //发生反应的总吸收系数(cm^-1)
		/********Determine the location of reaction********/
		e0 = getRandomNumber();
		double len = -log(1 - e0) / totalAbsorb;
		len *= CmToUm; //单位转换
		addVector(r, direction, len);
		if (r[2] > H) reflection(r);
		if (!isInMaterial(r))
			break;
		/********Determine the type of reaction********/
		e0 = getRandomNumber();
		bool reactionType = (e0 < pAbsorb / totalAbsorb); 
		/********If photoelectric effect happens********/
		if (reactionType == Photoelectric)
		{
			DepositEnergy += energy;
			break;
		}
		/********If compton effect happens********/
		else /****reactionType == Compton****/
		{
			double newEnergy = ComptonEffect(energy, direction);
			DepositEnergy += energy - newEnergy;
			energy = newEnergy;
		}
	} //cout << DepositEnergy << ' ';
	DepositEnergy = GaussExtend(DepositEnergy); //cout << DepositEnergy << endl;
	int ch = DepositEnergy / 55 * 1024; 
	channel[ch]++;
}
/********得到分辨率(输入成型时间)*********/
double getResolution(double formTime)
{
	double FWHM = 14.7 + 17.5*exp(-0.48*formTime);
	return FWHM / 685.96;
}
/********得到不同元素的特征X射线信息(预处理)********/
void getXrayEnergy()
{
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 2; k++) 
				cin >> XrayEnergy[i][j][k];
}
void stimulate(int Z)
{
	double s = 0;
	for (int i = 0; i < 4; i++)
		s +=XrayEnergy[Z][i][1]; //cout << s << endl;
	double e0;
	for (int i = 0; i < stimulateNumber; i++)
	{
		/********定义光子入射方向********/
		double r0[] = { 0,0,0 };
		double direction0[] = { 0,0,1 };
		/********决定这次入射的光子能量********/
		e0 = getRandomNumber(); //cout << e0 << endl;
		double s0 = 0;
		int j = 0;
		for (j = 0; j < 4; j++)
		{
			s0 += XrayEnergy[Z][j][1];
			if (e0 < s0 / s)
				break;
		} //cout << j << endl;
		singlePhotonStimulate(XrayEnergy[Z][j][0], r0, direction0);
	}
}
/****测试用，最终版会删去****/
void test()
{
	randNum = 1;
	stimulateNumber = 100000;
	stimulate(Cu);
	for (int i = 1; i <= 1024; i++)
		printf("%d\n", channel[i]);
}
/********输入********/
void init()
{
	int x;
	do
	{
		cout << "input a number between 1 and 78498 to determine initial random number:\n";
		cin >> x;
	} while (!(x >= 1 && x <= 78498));
	randNum = primeSet[x - 1];
	cout << "input the The atomic number of the excited substance\n";
	cout << "input the Number of simulations\n";
	cin >> stimulateNumber;
}
int main()
{
	FILE *fs = freopen("xrayinfo.txt", "r", stdin);
	FILE *fp = freopen("xray.txt", "w", stdout);
	getPrimeSet();
	getlnAbsorption();
	getXrayEnergy();
	//init();
	test();
    cout << "Hello World!\n"; 
	fclose(fs);
	fclose(fp);
	return 0;
}
