#include <iostream>
#include <math.h>
#include <fstream>
#include "class_vector.h"
#include "move.h"
#include "Const.h"
#include "Bunch.h"
#include "diffr.h"
#include "Any_dist.h"


long double dt = 1.0l;
long double k1 = 0.1l;
long double k2 = k1 / 1000.0l;
long double t = 0.0l;
long double mass_of_electron = 9.1e-31;
vector g(0.0l, 10.0l, 0.0l);
vector r1(0.0l, 0.0l, 0.0l);
vector v01(1.0l, 0.0l, 0.0l); // m/s
/*particle a1(mass_of_electron, r1, v01 * mass_of_electron, q_sgs);*/
vector F;
std::ofstream outF("data.txt");

//Функция для аппроксимации 
long double new_func(long double x) // задаем функцию
{
	return pow(cosh(x), 2);
}

int main()
{
	//1. Попрыгунчик
	/*while (t < 16) {
		outF << t << " " << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;

		vector e_v = a1.p.norm() / mass_of_electron;
		F = e_v * (-(k1 * v01.abs() + k2 * v01.abs() * v01.abs()));

		a1.r = a1.r + (a1.p / a1.m) * dt;
		a1.p = a1.p + (g * a1.m + F) * dt;

		t += dt;
	}
	outF << t <<" " << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;
	*/
	//2. Если просто 20 частиц летят параллельно друг другу с некоторым фиксированным расстоянием
	/*std::ofstream outF_r("20.txt");
	long double N = 20.0l;
	double particles[20] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 };
	vector r_i(0.0l, 0.0l, particles[0]);
	particle a1(mass_of_electron, r_i, v01 * mass_of_electron, q_sgs);
	for(int i = 0.0l;i<=N; i++){
		vector r_i(0.0l,0.0l ,particles[i]);
		particle a1(mass_of_electron, r_i, v01 * mass_of_electron, q_sgs);
		a1.r = a1.r + (a1.p / a1.m) * dt;
		a1.p = a1.p + g  * dt;
		outF_r << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;
	}
	outF_r << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;*/
	//3. Если интересует движение вдоль плоскости под углом
	/*::ofstream outF_angle("angle.txt");
	vector r_a(1.0l*cos(45), 0.0l, 0.0l);
	particle a1(mass_of_electron, r_a, v01 * mass_of_electron, q_sgs);
	while (t < 10) {
		outF_angle << t << " " << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;

		vector e_v = a1.p.norm() / mass_of_electron;
		F = e_v * (-(k1 * v01.abs() + k2 * v01.abs() * v01.abs()));

		a1.r = a1.r + (a1.p / a1.m) * dt;
		a1.p = a1.p + (g * a1.m + F) * dt;

		t += dt;
	}
	outF_angle << t << " " << a1.r.x << " " << a1.r.y << " " << a1.r.z << std::endl;*/
	//4. Введем кулоновскую силу и промоделируем движение частиц с ней 
	/*long double charge = -1.6e-19;
	vector F_k_1 = {};
	vector F_k_2 = {};
	long double r_q = 0.0000001l;
	vector v_a(1.0l, 0.0l, 0.0l); // m/s
	vector v_b(2.0l, 0.0l, 0.0l); // m/s
	vector r_a(0.0l, -0.000001l, 0.0l); // m
	vector r_b(0.0l, 0.000001l, 0.0l);
	particle a(mass_of_electron, r_a, v_a * mass_of_electron, q_sgs);
	particle b(mass_of_electron, r_b, v_b * mass_of_electron, q_sgs*(+1));
	std::ofstream outF_Kulon_a("Kulon_a.txt");
	std::ofstream outF_Kulon_b("Kulon_b.txt");
	while (t < 10) {
		outF_Kulon_a << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << std::endl;
		outF_Kulon_b << t << " " << b.r.x << " " << b.r.y << " " << b.r.z << std::endl;
		F_k_1 = ((a.r - b.r) / (a.r - b.r).abs()) * ((q * q) / ((a.r - b.r).abs() * (a.r - b.r).abs())) * k;
		F_k_2 = ((b.r - a.r) / (b.r - a.r).abs()) * ((q * q) / ((b.r - a.r).abs() * (b.r - a.r).abs())) * k;
		a.r = a.r + (a.p / a.m) * dt;
		a.p = a.p + F_k_1 * dt;

		b.r = b.r + (b.p / b.m) * dt;
		b.p = b.p - F_k_1 * dt;

		t += dt;

	}
	outF_Kulon_a << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << std::endl;
	outF_Kulon_b << t << " " << b.r.x << " " << b.r.y << " " << b.r.z << std::endl;*/
	
	//5. Неподвижный ион и подвижный электрон
/*	vector F_k_1 = {};
	vector v_a(1.0l, 0.0l, 0.0l); // m/s
	vector v_b(2.0l, 0.0l, 0.0l); // m/s
	vector r_a(0.0l, 0.0l, 0.0l); // m
	vector r_b(0.0l, -10.0l, 0.0l);
	particle a(mass_of_electron, r_a, v_a * mass_of_electron, q_sgs);
	particle b(mass_of_electron, r_b, v_b * mass_of_electron, q_sgs*(5));
	std::ofstream outF_electron_a("electron.txt");
	std::ofstream outF_Ion_b("Ion.txt");
	while (t < 105)
	{
		outF_electron_a << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << std::endl;
		outF_Ion_b << t << " " << b.r.x << " " << b.r.y << " " << b.r.z << std::endl;
		F_k_1 = ((a.r - b.r) / (a.r - b.r).abs()) * (( q * q) / ((a.r - b.r).abs() * (a.r - b.r).abs())) * k;
		a.r = a.r + (a.p / a.m) * dt;
		a.p = a.p + F_k_1 * dt;

		b.r = b.r;
		b.p = b.p - F_k_1 * dt;
		
		if (a.r.y == b.r.y)
		{
			break;
		}
		t += dt;
	}
	outF_electron_a << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << std::endl;
	outF_Ion_b << t << " " << b.r.x << " " << b.r.y << " " << b.r.z << std::endl;
*/
//6. Метод Эйлера с пересчетом для вычисления силы Лоренца
/*vector v_a(1.0l, 0.0l, 0.0l); // m/s
vector v_b(2.0l, 0.0l, 0.0l); // m/s
vector r_a(0.0l, 0.0l, 0.0l); // m
vector r_b(0.0l, -10.0l, 0.0l);
vector r_Q (0.0l, -10.0l, 0.0l);
vector r_q (0.0l, 10.0l, 0.0l);
vector R = r_Q - r_q;
double Q = v_a.ugol(R);
particle a(mass_of_electron, r_a, v_a * mass_of_electron, q_sgs);
particle b(mass_of_electron, r_b, v_b * mass_of_electron, q_sgs*(5));
particle a_pr = {};
a_pr.m = a.m;
a_pr.q = a.q;

long double gamma= sqrt(1 / (1 - (v_a.pscal(v_a) / c_sgs * c_sgs)));
vector p = v_a*m_e_sgs * gamma;
vector E_field = (R * q_sgs / (R.abs()*R.abs()*R.abs()))*(1 - (v_a.pscal(v_a) / c_sgs * c_sgs)) / (1 - (v_a.pscal(v_a) / c_sgs * c_sgs)*sin(Q)*sin(Q));
vector B_field = v_a.pvec(E_field)*(1 / c_sgs);
vector sila = (E_field + (p.pvec(B_field))/(gamma*m_e_sgs*c_sgs))*q_sgs;
while (t < 10)
{
	a_pr.r = (a.p / (a.m * a.gamma())) * dt + a.r;
	a_pr.p = sila * dt + a.p;
	a.r = a.r + ((a_pr.p / a_pr.gamma() + a.p / a.gamma()) / (2 * a.m) * dt);
	a.p = a.p + ((sila + sila) / 2) * dt;
	t += dt;

	// outF0 << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << " " << a.p.x << " " << a.p.y << " " << a.p.z << std::endl; //смотрим на траекторию 1 частицы
}*/
//7. Бокс-Мюллера
/*long double sigm = 1; // отклонение
long double mu = 0; // мат ожидание

long double eps0 = {}; // распределение 1
long double z0 = {}; // случайная величина 1

std::ofstream outF_BM("Box-M.txt");
std::mt19937 mersenne{ static_cast <std::mt19937::result_type> (std::time(nullptr)) };

std::uniform_real_distribution <> die{ -1, 1 };

for (int count{ 1 }; count <= 1000; ++count)
{
	long double x = die(mersenne);
	long double y = die(mersenne);

	long double s = x * x + y * y;

	if (s > 0 and s <= 1)
	{
		z0 = x * sqrt((-2 * log(s)) / s);

		eps0 = mu + sigm * z0;

		outF_BM << count << " "<< eps0 << std::endl;
	}

}*/
//9. Бокс-Мюллера
/*long double sigm = 1; // отклонение
long double mu = 0; // мат ожидание

long double eps0 = {}; // распределение 1
long double z0 = {}; // случайная величина 1

std::ofstream outF_BM("Box-M.txt");
std::mt19937 mersenne{ static_cast <std::mt19937::result_type> (std::time(nullptr)) };

std::uniform_real_distribution <> die{ -1, 1 };

for (int count{ 1 }; count <= 1000; ++count)
{
	long double x = die(mersenne);
	long double y = die(mersenne);

	long double s = x * x + y * y;

	if (s > 0 and s <= 1)
	{
		z0 = x * sqrt((-2 * log(s)) / s);

		eps0 = mu + sigm * z0;

		outF_BM << count << " "<< eps0 << std::endl;
	}

}*/
//10. Распределение 
/*bunch b(1000, 1, 0.5, 1, 1, vector(1, 2, 3), 0.2, 1);
b.generate();

std::ofstream ofile("bunch.txt");
for (int i = 0; i < 1000; i++) {
	ofile << b.mass_ch[i].get_W() << " " << b.mass_ch[i].r.x << " " << b.mass_ch[i].r.y << " " << b.mass_ch[i].r.z << std::endl;
}*/
//11. Задача с магнитами
/*const long double magn_h = 0.5; // ширина магнита
const long double magn_a = 1.5; // расстояние между 1 и 2, 3 и 4
const long double magn_b = 3.5; // расстояние между 2 и 3
const long double B = 30; // модуль магнитной индукции
vector magn_pole(vector& r) // определяем положение и магнитную индукцию 4 магнитов
{
	if (0 <= r.z and r.z <= magn_h)
	{
		return vector(-B, 0, 0);
	};
	if ((magn_h + magn_a) <= r.z and r.z <= (magn_h * 2 + magn_a))
	{
		return vector(B, 0, 0);
	};
	if ((magn_h * 2 + magn_a + magn_b) <= r.z and r.z <= (magn_h * 3 + magn_a + magn_b))
	{
		return vector(B, 0, 0);
	};
	if ((magn_h * 3 + magn_a * 2 + magn_b) <= r.z and r.z <= (magn_h * 4 + magn_a * 2 + magn_b))
	{
		return vector(-B, 0, 0);
	};
	return vector(0, 0, 0);
};
vector sila(particle& part) // predict
{
	return part.p.vec(magn_pole(part.r)) * (part.q / (c_sgs * part.gamma() * part.m));
};

long double W_zover_xyz(long double x, long double y, long double z) // зависимость энергии от координаты z (прямая, убывающая)
{
	return -z;
}
bunch b_magn(1000, 1.5 * m_e_sgs * c_sgs * c_sgs, 0.008 * m_e_sgs * c_sgs * c_sgs, m_e_sgs, -q_sgs, vector(0, 0, -1), 0.001, lz, 0); // создается пучок
b_magn.generate(W_zover_xyz);

long double min_z = {};
long double max_z = {};

for (int i = 0; i < 1000; i++) // ищем минимум и максимум конечноых z координат пучка
{
	particle a = b_magn.mass_ch[i]; // particle a = b_magn.mass_ch[0] для выведения траектории
	particle a_pr = {};
	a_pr.m = a.m;
	a_pr.q = a.q;

	std::ofstream outF0("data_magn.txt");

	// outF0 << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << " " << a.p.x << " " << a.p.y << " " << a.p.z << std::endl; //смотрим на траекторию 1 частицы

	t = 0;
	// метод Эйлера с пересчетом
	while (t < t_end)
	{
		a_pr.r = (a.p / (a.m * a.gamma())) * dt + a.r;
		a_pr.p = sila(a) * dt + a.p;
		a.r = a.r + ((a_pr.p / a_pr.gamma() + a.p / a.gamma()) / (2 * a.m) * dt);
		a.p = a.p + ((sila(a_pr) + sila(a)) / 2) * dt;

		t += dt;

		// outF0 << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << " " << a.p.x << " " << a.p.y << " " << a.p.z << std::endl; //смотрим на траекторию 1 частицы
	}

	// написать сравнение минимумов и максимумов z
	if (i == 0 || a.r.z < min_z) // || - это or
	{
		min_z = a.r.z;
	}

	if (i == 0 || a.r.z > max_z) // || - это or
	{
		max_z = a.r.z;
	}
}
std::cout << "min_z = " << min_z << std::endl;
std::cout << "max_z = " << max_z << std::endl;
std::cout << "delta_z = " << (max_z - min_z) << std::endl;
std::cout << "l_z = " << lz << std::endl;*/
//12. Апроксимация
long double a = -2.0l;
long double b = 2.0l; // границы, в которых необходима аппроксимация 
//посчитаем константу с для нахождения плотности распределения 
const int N = 500; // число точек 
long double dx = (b - a) / N; //шаг по сетке 
long double C = 1 / (integ(new_func, dx, a, b)); // константа из условия непрерывности, понадобится для вычисления значения функции
int r = N + 1;
int l = 0; //необходимо для сравнения с случайным числом и корректной аппроксимации
long double x_i[N];
int F_x_i[N];
int k_i[N];
int b_i[N];
F_x_i[0] = 0;
F_x_i[N - 1] = 1;
//массивы значений точки, в которой определена функция, самой функции, а также коэффициентов для прямой 
for (int i = 0; i <= N; i++) // массив х
{
	x_i[i] = a + dx * i;
}
for (int i = 1; i <= N - 1; ++i)
{
	F_x_i[i] = F_x_i[i - 1] + integ_sr(new_func, (x_i[i] - x_i[i - 1]) / N, x_i[i - 1], x_i[i]) * C;
}//это то, как для точек будут вычисляться значения функции
//теперь коэффициенты для уравнения прямой
for (int i = 0; i <= N - 1; ++i)
{
	k_i[i] = (x_i[i + 1] - x_i[i]) / (F_x_i[i + 1] - F_x_i[i]);
	b_i[i] = x_i[i] - k_i[i] * F_x_i[i];
}
// напишем функцию, которая позволит выводить случайные числа в заданном отрезке 
std::mt19937 mersenne{ static_cast <std::mt19937::result_type> (std::time(nullptr)) };
std::uniform_real_distribution <> die{ 0, 1 };
long double random_point = die(mersenne);
while (r - l > 1) {
	int m = (r + l) / 2;
	if (random_point < F_x_i[m]) {
		r = m;
		k_i[l] * random_point + b_i[l];
	}
	else {
		l = m;
		k_i[l] * random_point + b_i[l];
	}
}
std::ofstream outF("appr.txt");

any_dist<long double> dist(new_func, a, b);
for (int i = 0; i <= dist.N; ++i)
{
	outF << dist.x_i[i] << " " << dist.dist_func[i] << " " << new_func(dist.x_i[i]) << std::endl;
}

// выводим случайные числа и строим гистограмму
std::ofstream outF1("appr1.txt");
//13. аппроксимация с эммитансом (мега-функция)
/*
bunch b_1(100, 1.5*m_e_sgs * c_sgs * c_sgs, 0.15 * m_e_sgs * c_sgs * c_sgs, m_e_sgs, 100e-9, vector(1, 2, 3), 0.2, 1); // создается пучок
b_1.generate_dist(poln_dist);
std::ofstream outF06("poln_dist.txt");
outF06.precision(10);
for (int i = 0; i < b_1.N; i++)
{
	outF06 << b_1.mass_ch[i].r.x << " " << b_1.mass_ch[i].r.y << " " << b_1.mass_ch[i].r.z << " " << b_1.mass_ch[i].p.x << " " << b_1.mass_ch[i].p.y << " " << b_1.mass_ch[i].p.z << std::endl;
}
*/
//14. Мультиполь +мультипольное взаимодействие 
// мультипольный момент (электрический мультиполь)
	// здесь моделировали один мультиполь и строили его поле
	/*
	multy o(3, q_sgs, vector(2, 2.5, 0), vector(-2, 1, 0), vector(1, 1, 0)); // для одного мультиполя поле

	std::ofstream ofile1("multy.txt");

	for (int i = 0; i < o.N; i++) // выводили заряд и расположение частиц мультиполя
	{
		ofile1 << o.mass_toch[i].q << " " << o.mass_toch[i].r.x << " " << o.mass_toch[i].r.y << " " << o.mass_toch[i].r.z << std::endl;
	}

	long double min_x = 0.01;
	long double max_x = 5.01;
	long double min_y = 0.01;
	long double max_y = 6.01;
	long double shag = 0.25;
	int n_x = floor((max_x - min_x) / shag + 1);
	int n_y = floor((max_y - min_y) / shag + 1);

	ofile1 << n_x << " " << n_y << " 0 0 0" << std::endl;
	for (int i = 0; i < n_x; i++)
	{
		vector tn(min_x + shag * i, 0, 0); // определяем положение точки по x

		for (int j = 0; j < n_y; j++)
		{
			tn.y = min_y + shag * j; // определяем положение точки по y

			vector E = o.get_E_static(tn); // определяем поле в каждой полученной точке

			ofile1 << tn.x << " " << tn.y << " " << E.x << " " << E.y <<  " " << E.abs() << std::endl;
		}
	}
	*/

	// здесь моделировали два мультиполя и строили их силу взаимодействия
/*std::ofstream ofile2("multy_x_2.txt");
vector r01(2, 2.5, 0);
vector r02(5, 2, 0);

multy o1(3, q_sgs, r01, vector(-2, 1, 0), vector(1, 1, 0)); // силу между двумя мультиполями
multy o2(2, q_sgs, r01, vector(2, 1, 0), vector(1, 1, 0));

int n_oo = 100;
vector dl = (r02 - r01) / (n_oo - 1);

for (int i = 0; i < n_oo; i++) // смещаем один мультиполь относительно второго, первый стоит на месте
{
	vector r = r01 + dl * i; // очередное положение второго мультиполя
	o2.set_r0(r);

	vector F = o2.get_F_static(&o1);
	ofile2 << dl.abs() * i << " " << F.abs() << std::endl;
}

// магнитный квадруполь

bunch b(1, 50e6, 0.15 * m_e_sgs * c_sgs * c_sgs, m_e_sgs, 100e-9, vector(1, 2, 3), 1e-6, 1e-6); // создается пучок
b.generate();


particle a = b_magn.mass_ch[0]; // particle a = b_magn.mass_ch[0] для выведения траектории
particle a_pr = {};
a_pr.m = a.m;
a_pr.q = a.q;

std::ofstream outF05("data_magn.txt");

outF05 << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << " " << a.p.x << " " << a.p.y << " " << a.p.z << std::endl; //смотрим на траекторию 1 частицы

t = 0;
t_end = 6e-9l;
// метод Эйлера с пересчетом
while (t < t_end)
{
	a_pr.r = (a.p / (a.m * a.gamma())) * dt + a.r;
	a_pr.p = sila2(&a, 0) * dt + a.p;
	a.r = a.r + ((a_pr.p / a_pr.gamma() + a.p / a.gamma()) / (2 * a.m) * dt);
	a.p = a.p + ((sila2(&a_pr, 0) + sila2(&a, 0)) / 2) * dt;

	t += dt;

	outF05 << t << " " << a.r.x << " " << a.r.y << " " << a.r.z << " " << a.p.x << " " << a.p.y << " " << a.p.z << std::endl; //смотрим на траекторию 1 частицы
}
*/
//15. БАК(мега-функция)
/* коллайдер
//1. задаем функцию опис. магнитное поле
//НАДО вывести Максимальное отклонение пучка от Z в те моменты времени, когда первая частица влетает в магнит
int N_1 = 100;//количество частиц
long double W = 5e7*ev_sgs;
long double delta_W = W / 1000;
lz = 1e-4;
long double R = 1e-4;
t_max = 1e-8;
dt = t_max / 1e4;
bunch b_collad(N_1, W, delta_W, m_e_sgs, 1000 * q_sgs, vector(0, 0, 0), R, lz, 0);
b_collad.generate();
mega_func(dt, N_1, b_collad.mass_ch, 0, nullptr, kulon_force_relat, sila2, false, t_max, &n_iter, &mass_time, &mass_r, nullptr, true, stop_cond_for_collad);
std::cout << "N_iter = " << n_iter << std::endl;
std::ofstream outF08("collad.txt");
for (int i = 0; i < n_iter; i++)
{
	outF08 << mass_time[i];
	for (int j = 0; j < 10; j++)
	{
		outF08 << " " << mass_r[j][i].x << " " << mass_r[j][i].y << " " << mass_r[j][i].z;
	}
	outF08 << std::endl;
}*/


	return 0;
}

