#pragma once
#include <math.h>
#include <random> // Для вихря
#include <ctime> // Для вихря
#include "diffr.h"



// класс для аппроксимации(генерация случайных чисел с заданным распределением)	
template <typename T>
class any_dist
{
public:
	any_dist(T(*dens_func)(T), T a, T b, int N = 500) {
		this->dens_func = dens_func;
		this->a = a;
		this->b = b;
		this->N = N;

		x_i = new T[N + 1]; // массив х
		dist_func = new T[N + 1]; // массив значений функции
		appr_k = new T[N]; // массив для k
		appr_c = new T[N]; // массив для c

		delta_x = (b - a) / N; // шаг для значений от а до в

		for (int i = 0; i <= N; ++i) // массив х
		{
			x_i[i] = a + delta_x * i;
		}

		make_dist_func();
	};
	~any_dist() {
		delete[] x_i;
		delete[] dist_func;
		delete[] appr_c;
		delete[] appr_k;
	};

	void make_dist_func(void) { // массив значений функции
		T c = 1 / integ(dens_func, (b - a) / 1e5l, a, b);
		dist_func[0] = 0;
		dist_func[N] = 1;
		for (int i = 1; i <= N - 1; ++i)
		{
			dist_func[i] = dist_func[i - 1] + integ_sr(dens_func, (x_i[i] - x_i[i - 1]) / 500, x_i[i - 1], x_i[i]) * c;
		}

		// задаем коэффициенты для участка апроксимации
		for (int i = 0; i <= N - 1; ++i)
		{
			appr_k[i] = (x_i[i + 1] - x_i[i]) / (dist_func[i + 1] - dist_func[i]);
			appr_c[i] = x_i[i] - appr_k[i] * dist_func[i];
		}
	}
	//генерация случайных чисел на определенном отрезке
	T get(void) {
		T rnd = die(mersenne);
		int l = 0, r = N;
		while (r - l > 1) {
			int m = (r + l) / 2;
			if (rnd < dist_func[m])
				r = m;
			else
				l = m;
		}
		return appr_k[l] * rnd + appr_c[l];
	}

	T(*dens_func)(T);
	T a, b, delta_x;
	int N;
	T* x_i; // обращаемся к массиву x
	T* dist_func; // обращаемся к массиву функции
	std::mt19937 mersenne{ static_cast <std::mt19937::result_type> (std::time(nullptr)) }; //задаем вихрь Мерсена
	std::uniform_real_distribution <> die{ 0, 1 }; //задаем диапазон случайных чисел
	T* appr_k; // задаем k для апроксимации, здесь обращаемся к массиву
	T* appr_c; // задаем c для апроксимации, здесь обращаемся к массиву
};