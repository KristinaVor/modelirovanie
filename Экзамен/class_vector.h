#pragma once
#include <iostream>
#include <math.h>
#include <cmath>
#include "const.h"


class vector
{
public:

	double x = 0; //проекция вектора на ось x
	double y = 0; //проекция вектора на ось y
	double z = 0; //проекция вектора на ось z

	vector() {	} // конструктор
	
	vector(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	vector operator +(vector v) // оператор для перегрузки функции сумма
	{
		vector tmp = {};
		tmp.x = this->x + v.x;
		tmp.y = this->y + v.y;
		tmp.z = this->z + v.z;
		return tmp;
	}

	vector operator -(vector v)
	{
		vector tmp = {};
		tmp.x = this->x - v.x;
		tmp.y = this->y - v.y;
		tmp.z = this->z - v.z;
		return tmp;
	}

	vector operator /(double c)
	{
		vector tmp = {};
		tmp.x = this->x / c;
		tmp.y = this->y / c;
		tmp.z = this->z / c;
		return tmp;
	}

	vector operator *(double c)
	{
		vector tmp = {};
		tmp.x = this->x * c;
		tmp.y = this->y * c;
		tmp.z = this->z * c;
		return tmp;
	}

	vector operator =(const vector& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	double scal(vector v)
	{
		double c = this->x * v.x + this->y * v.y + this->z * v.z;
		return c;
	}

	double abs()
	{
		double c = sqrt(scal(*this));
		return c;
	}

	vector vec(vector v)
	{
		vector tmp = {};
		tmp.x = this->y * v.z - this->z * v.y;
		tmp.y = this->z * v.x - this->x * v.z;
		tmp.z = this->x * v.y - this->y * v.x;
		return tmp;
	}

	vector norm()
	{
		vector tmp = {};

		tmp.x = this->x / abs();
		tmp.y = this->y / abs();
		tmp.z = this->z / abs();

		return tmp;
	}

	double pscal(vector v)
	{
		double c = scal(v.norm());
		return c;
	}

	vector pvec(vector v)
	{
		vector vec = v.norm() * pscal(v);

		return vec;
	}

	double ugol(vector v)
	{
		double tmp = {};
		tmp = acos(scal(v) / (v.abs() * abs()));
		return tmp;
	}

	vector sum(vector v1, vector v2)
	{
		vector tmp = {};

		tmp.x = v1.x + v2.x;
		tmp.y = v1.y + v2.y;
		tmp.z = v1.z + v2.z;

		return tmp;
	}

	void show()
	{
		std::cout << x << " " << y << " " << z << std::endl;
	}

	~vector(){	} // деструктор
};