#pragma once
#include <math.h>

/* Суммирование дифференциалов*/
long double diff_sym(long double f(long double), long double x0)
{
	long double dx = 1.503e-7l;
	long double tmp = (f(x0 + dx) - f(x0 - dx)) / (2 * dx);
	return tmp;
}

/*Диффернцирование правой и левой части*/
long double diff_pr(long double f(long double), long double x0, long double dx)
{
	long double tmp = (f(x0 + dx) - f(x0)) / (dx);
	return tmp;
}

long double diff_sl(long double f(long double), long double x0, long double dx)
{
	long double tmp = (f(x0) - f(x0 - dx)) / (dx);
	return tmp;
}

/* Правильная производная с учетом коэффициентов точности (их находили в википедии)*/
long double diff_koef(long double f(long double), long double x0)
{
	long double dx = 1.503e-7l;
	long double tmp = (1.0l / 280.0l) * (f(x0 - 4.0l * dx));
	tmp += (-4.0l / 105.0l) * (f(x0 - 3.0l * dx));
	tmp += (1.0l / 5.0l) * (f(x0 - 2.0l * dx));
	tmp += (-4.0l / 5.0l) * (f(x0 - 1.0l * dx));
	tmp += (-1.0l / 280.0l) * (f(x0 + 4.0l * dx));
	tmp += (4.0l / 105.0l) * (f(x0 + 3.0l * dx));
	tmp += (-1.0l / 5.0l) * (f(x0 + 2.0l * dx));
	tmp += (4.0l / 5.0l) * (f(x0 + 1.0l * dx));
	long double dif = tmp / dx;
	return dif;
}

/* Интегрирование*/
template <typename T>
T integ (T f(T), T dx, T a, T b) {
	T sum = 0;
	T N = (b - a) / dx;
	for (int i = 0; i <= N; ++i) {
		sum += f(a + dx * i) + f(a + dx * (i + 1));
	}

	return sum * dx / 2;
}

/* Интегрирование второй вариант для проверки ошибки*/
template <typename T>
T integ_mod (T f(T), T dx, T a, T b) {
	T sum = 0;
	T N = (b - a) / dx;
	for (int i = 0; i <= N; ++i) {
		sum += f(a + dx * i);
	}

	return sum * dx / 2;
}

/*Использование среднего значения в интегрирование*/
template <typename T>
T integ_sr(T f(T), T dx, T a, T b) {
	T sum = 0;
	T N = (b - a) / dx;
	for (int i = 0; i <= N; ++i) {
		sum += f(a + dx * i);
	}

	return sum * (b - a) / N;
}