#pragma once
#include "class_vector.h"
#include "move.h"
#include <random> // Для вихря
#include <ctime> // Для вихря
#include "const.h"
#include "C:\Users\User\Desktop\lesiv2.0\lesiv\function\function\Any_dist.h"

class bunch
{
public:

	int N = {}; // Число частиц
	long double W = {}; // средняя энергия частиц
	long double WW = {}; //разброс энергии

	vector r0 = {}; // средний вектор
	long double R = {}; //радиус пучка
	long double lz = {}; //длина сгустка по оси z
	long double pp_xy = {}; // разброс импульса по осям x y

	long double p = {}; //импульс по z
	long double pp = {}; //отклонение импульса (добавка к импульсу)

	long double m0 = {}; // Масса 1 частицы
	long double q0 = {}; // Заряд 1 частицы

	particle* mass_ch; // указатель на массив

	std::mt19937 mersenne{ static_cast <std::mt19937::result_type> (std::time(nullptr)) };

	std::uniform_real_distribution <> die{ -1, 1 };

	long double norm_rasp(long double mu, long double sigm) // функция, возвращает нормальное распределение
	{
		
		while (true)
		{
			long double x = die(mersenne);
			long double y = die(mersenne);

			long double s = x * x + y * y;

			if (s > 0 and s <= 1)
			{
				long double z0 = x * sqrt((-2 * log(s)) / s);

				return (mu + sigm * z0);
			}
		}
	}

	bunch(int N, long double W, long double WW, long double m0, long double q0, vector r0, long double R, long double lz, long double pp_xy = 0.1)
	{
		this->N = N;
		this->W = W;
		this->WW = WW;
		this->m0 = m0;
		this->q0 = q0;
		this->r0 = r0;
		this->R = R;
		this->lz = lz;
		this->pp_xy = pp_xy;
		p = sqrt((W * W) / (c_sgs * c_sgs) - (m0 * m0) * (c_sgs * c_sgs));
		pp = sqrt(((W + WW) * (W + WW)) / (c_sgs * c_sgs) - (m0 * m0) * (c_sgs * c_sgs)) - p;
	
		mass_ch = new particle[N];
	}

	void generate(void) // генерирует пучок с нормальными распределениями по осям
	{
		for (int i = 0; i < N; ++i)
		{
			mass_ch[i].m = m0;
			mass_ch[i].q = q0;
			mass_ch[i].r.x = norm_rasp(r0.x, R / 3.5);
			mass_ch[i].r.y = norm_rasp(r0.y, R / 3.5);
			mass_ch[i].r.z = norm_rasp(r0.z, lz / 7);
			mass_ch[i].p.x = norm_rasp(0, pp_xy / 7);
			mass_ch[i].p.y = norm_rasp(0, pp_xy / 7);
			mass_ch[i].p.z = norm_rasp(p, pp / 7);	
		}
	}

	void generate(long double (* W_zover_xyz)(long double x, long double y, long double z)) // генерирует пучок с нормальными распределениями по координатам, px, py, а pz зависит от z
	{
		for (int i = 0; i < N; ++i)
		{
			mass_ch[i].m = m0;
			mass_ch[i].q = q0;
			mass_ch[i].r.x = norm_rasp(r0.x, R / 3.5);
			mass_ch[i].r.y = norm_rasp(r0.y, R / 3.5);
			mass_ch[i].r.z = norm_rasp(r0.z, lz / 7);
			mass_ch[i].p.x = norm_rasp(0, pp_xy / 7);
			mass_ch[i].p.y = norm_rasp(0, pp_xy / 7);
			long double W_zz = (WW / 2) * W_zover_xyz((mass_ch[i].r.x - r0.x) / R, (mass_ch[i].r.y - r0.y) / R, (mass_ch[i].r.z - r0.z) / (0.5 * lz)) + W; // энергия конкретной частицы со случайными координатами
			mass_ch[i].p.z = sqrt((W_zz * W_zz) / (c_sgs * c_sgs) - (m0 * m0) * (c_sgs * c_sgs));
		}
	}

	void generate_dist(long double(*func)(long double)) // генерирует пучок с распределение с плотностью вероятности
	{
		any_dist<long double> dist(func, -1, 1);
		for (int i = 0; i < N; ++i)
		{
			mass_ch[i].m = m0;
			mass_ch[i].q = q0;
			mass_ch[i].r.x = dist.get()*R + r0.x;
			mass_ch[i].r.y = dist.get()*R + r0.y;
			mass_ch[i].r.z = dist.get()*lz/2 + r0.z;
			mass_ch[i].p.x = dist.get()*pp_xy/2;
			mass_ch[i].p.y = dist.get()*pp_xy/2;
			mass_ch[i].p.z = dist.get()*pp / 2 + p;
		}
	}

	~bunch() 
	{
		delete[] mass_ch;
	}
};