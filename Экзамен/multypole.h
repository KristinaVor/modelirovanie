#pragma once
#include "class_vector.h"
#include "const.h"
#include "move.h"

class multy
{
public:

	int n = {}; // ���������� ��� ����������
	int N = {}; // ����� ���������� ������
	int nx = {}; // ���������� ����� �� x
	int ny = {}; // ���������� ����� �� y
	int nz = {}; // ���������� ����� �� z
	long double q = {}; // ����� ������
	

	vector r0 = {}; // ������ � �����
	vector r_ijk = {}; // ������ ��� ������ �����
	vector l = {}; // ������� � ����������� ����� 2 �������
	vector lx = {}; // ������� � ����������� ����� 2 ������� �� x
	vector ly = {}; // ������� � ����������� ����� 2 ������� �� y
	vector lz = {}; // ������� � ����������� ����� 2 ������� �� z
	vector h = {}; // ������ 90 �������� �� l

	particle* mass_toch; // ��������� �� ������

	multy(int n, long double q, vector r0, vector l, vector h)
	{
		this->n = n;
		this->q = q;
		this->r0 = r0;
		this->l = l;
		this->h = h;

		N = 1 << (n - 1);

		mass_toch = new particle[N];

		nx = 1 << ((n + 1) / 3);
		ny = 1 << (n / 3);
		nz = 1 << ((n - 1) / 3);

		lx = l;
		lz = (l.vec(h)).norm() * (l).abs();
		ly = lz.vec(lx) / (l).abs();

		int d = 0; // ������ �������

		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				for (int k = 0; k < nz; k++)
				{
					r_ijk = r0 + lx * (i - ((long double)nx - 1) / 2) + ly * (j - ((long double)ny - 1) / 2) + lz * (k - ((long double)nz - 1) / 2);
					particle p;
					p.q = (((i + j + k) & 1) == 0) ? q : (-q); // ���������� + ��� - ������
					p.r = r_ijk;

					mass_toch[d] = p;

					d += 1;
				}
			}
		}

	}

	vector get_E_static(vector r_tn) // tn - ����� ����������
	{
		vector E(0, 0, 0); // ��������� ����

		for (int i = 0; i < N; i++)
		{
			E = E + mass_toch[i].get_E_static(r_tn);
		}
		return E;
	}

	vector get_F_static(multy* mp) // mp - ����������
	{
		vector F(0, 0, 0); // ��������� ����

		for (int i = 0; i < N; i++)
		{

			F = F + mp->get_E_static(mass_toch[i].r) * mass_toch[i].q;
		}
		return F;
	}

	void set_r0(vector new_r0)
	{
		vector delta_r0 = new_r0 - r0;
		for (int i = 0; i < N; i++)
		{

			 mass_toch[i].r = mass_toch[i].r + delta_r0;
		}

		r0 = new_r0;
	}


	~multy() 
	{
		delete[] mass_toch;
	}
};