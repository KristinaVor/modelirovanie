#pragma once
#include "move.h"

typedef vector(*force_fn_type)(particle *a, particle *b, int index_a, int index_b); // сила котора€ действует на частицу а со стороны б
typedef vector(*static_force_fn_type)(particle *a, int index_a); //функци€ описывающа€ силу со стороны пол€
typedef  bool(*stop_cond_fn_type)(long double t, int N, particle *particles); //ф. услви€ остановки, возвращает тру если нужно остановить


void mega_func(long double dt,
	int N, particle *particles,
	int M, particle *static_particles,
	force_fn_type force_fn,
	static_force_fn_type static_force_fn,
	bool indep,
	long double t_max,
	int *n_iter,
	long double **mass_time,
	vector ***mass_r,
	vector ***mass_p,
	bool modif_euler,
	stop_cond_fn_type stop_cond_fn);