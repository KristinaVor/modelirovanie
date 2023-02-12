#pragma once

template <typename T>
T* createArray(int n_string)
{
	T* array = new T[n_string];

	return array;
}

template <typename T>
void  deleteArray(T* array)
{
	delete[] array;
}

template <typename T>
T** createArray2(int n_string, int n_column)
{
	T** array = new T * [n_string];
	
	for (int i = 0; i < n_string; ++i)
	{
		array[i] = new T[n_column];
	}

	return array;
}

template <typename T>
void deleteArray2(T** array, int n_string)
{
	for (int i = 0; i < n_string; ++i)
	{
		delete[]  array[i];
	}

	delete[] array;
}