#pragma once
#include "class_vector.h"
#include "const.h"

class particle
{
public:
	long double m = {}; // Масса частицы
	long double q = {};
	// long double W = {}; // средняя энергия частиц для задачи с импульсом, когда мы не задаем энергию

	vector r = {}; // радиус-вектор
	vector p = {}; // вектор импульса

	particle(long double m, vector r, vector p, long double q)
	{
		this->m = m;
		this->r = r;
		this->p = p;
		this->q = q;
	}

	long double get_W(void)
	{
		long double gam = sqrt(1 + (p.abs() / (m * c_sgs)) * (p.abs() / (m * c_sgs)));
		return gam * m * c_sgs * c_sgs;
	}

	long double get_W_z(void)
	{
		return sqrt(p.z * p.z + m * m * c_sgs * c_sgs) * c_sgs;
	}

	vector get_E_static( vector r_tn)
	{
		vector R = r_tn - r;
		return R * (q / powl(R.abs(), 3));
	}

	particle() {	}

	long double gamma(void) 
	{
		vector p_mc = p / (m * c_sgs);
		return sqrt(p_mc.scal(p_mc) + 1);
	}

	particle operator =(particle part)
	{
		m = part.m;
		q = part.q;
		r = part.r;
		p = part.p;
		return part;
	}

	~particle() {}
};