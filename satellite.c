/*
输入输出均为角度制，目前alpha还需重新计算
*/

#include<stdio.h>
#include<math.h>
#define pi 3.1415926535897932384626
#define G 6.67259e-11  //G为引力常数，m为地球质量
#define m 5.965e24
double alpha = 40; //当前时刻的格林尼治恒星时角alpha，即春分点与大地坐标系的x轴的夹角（格林尼治恒星时角需要根据当前时间另外计算）
double L, B; //星下点的经度L，纬度B

void location(double a, double e, double i, double omega, double w, double M0, double t){
	i = i / 90 * (pi / 2);
	w = w / 90 * (pi / 2);
	omega = omega / 90 * (pi / 2);
	M0 = M0 / 90 * (pi / 2);
	alpha = alpha / 90 * (pi / 2); //将角度都转换为弧度制
	double M, n; //t时刻的平近点角M，卫星平均角速度n
	n = sqrt(G * m / pow(a, 3));
	M = M0 + n * t;
	double E = 1, E0 = 0; //初始化偏近点角E
	while(fabs(E0 - E) > 1e-6){
		E = E0 -(E0 - M - e * sin(E0))/(1 - e * cos(E0));
		E0 = E;
	} //迭代求解开普勒方程，计算偏近点角E
	double xx, yy, zz; //卫星在轨道坐标系中的坐标
	double b = a * sqrt(1 - e * e);
	xx = a * cos(E) - a * e;
	yy = b * sin(E);
	zz = 0;
	double x, y, z; //大地坐标系，从轨道坐标系转换到大地坐标系
	x = xx * (cos(w) * cos(omega - alpha) - sin(w) * cos(i) * sin(omega - alpha)) + yy * (- sin(w) * cos(omega - alpha) - cos(w) * cos(i) * sin(omega - alpha)) + zz * (sin(i) * sin(omega - alpha));
	y = xx * (cos(w) * sin(omega - alpha) + sin(w) * cos(i) * cos(omega - alpha)) + yy * (- sin(w) * sin(omega - alpha) + cos(w) * cos(i) * cos(omega - alpha)) + zz * (- sin(i) * cos(omega - alpha));
	z = xx * (sin(w) * sin(i)) + yy * (cos(w) * sin(i)) + zz * cos(i);
	double theta = atan((z * a)/(b * sqrt(x * x + y * y)));
	double e_square = (a * a - b * b)/(a * a), ee_square = (a * a - b * b)/(b * b);
	L = atan(y / x) / (pi / 2) * 90; //得到经纬度
	B = atan((z + ee_square * b * pow(sin(theta), 3)) / (sqrt(x * x + y * y) - e_square * a * pow(cos(theta), 3))) / (pi / 2) * 90;
}

int main(){
	double a, e, i, omega, w, M0, t; //轨道半长轴a，轨道偏心率e，轨道倾角i，升交点赤经omega，近地点幅角w，t0时的平近点角M0，时刻t
	scanf("%lf %lf %lf %lf %lf %lf %lf", &a, &e, &i, &omega, &w, &M0, &t); //输入轨道六根数和时刻t
	location(a, e, i, omega, w, M0, t);
	printf("%lf %lf\n", L, B);
	return 0;
}