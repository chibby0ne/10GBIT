//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Namespace contains mathematic functions in the buffer class
/// \author Timo Lehnigk-Emden
/// \date   2011/05/31
//

#ifndef BUFFER_MATH_H_
#define BUFFER_MATH_H_

#include "buffer.h"

namespace buffer_math {

template<typename T> Buffer<T, 2> Identity_Matrix(unsigned int size)
{
	Buffer<T,2> ident(size, size);

	ident.Clear();

	for(unsigned int i = 0; i < size; i++)
		ident[i][i] = static_cast<T>(1);

	return ident;
}

// Binary, GF(2) operations
namespace bin {

template<typename T> T Mult_Bin(Buffer<T,1>& a, Buffer<T,1>& b)
{
	T res = 0;
	for(unsigned int j = 0; j < a.length(); j++)
	{
		res ^= a[j] & b[j];
	}
	return res;
}

template<typename T> Buffer<T, 2> Mult_Bin(Buffer<T,2>& a, Buffer<T,2>& b)
{
	Buffer<T,2> res(a.dim1(), a.dim2());

	res.Clear();


	for(unsigned int k = 0; k < a.dim1(); k++)
		for(unsigned int i = 0; i < a.dim1(); i++)
			for(unsigned int j = 0; j < a.dim2(); j++)
			{
				res[k][i] ^= a[k][j] & b[j][i];
			}

	return res;
}

template<typename T> Buffer<T, 2> Add_Bin(Buffer<T,2>& a, Buffer<T,2>& b)
{
	Buffer<T,2> res(a.dim1(), a.dim2());
	for(unsigned int i = 0; i < a.dim1(); i++)
		for(unsigned int j = 0; j < a.dim2(); j++)
			res[i][j] = a[i][j] ^ b[i][j];

	return res;
}

} // End namespace bin

/// Find the index of the minimum value in a buffer
template<typename T> unsigned int Min_Index(Buffer<T,1>& buffer)
{
	unsigned int index = 0;

	for(unsigned int i = 1; i < buffer.length(); i++)
	{
		if(buffer[i] < buffer[index])
			index = i;
	}
	return index;
}

} // End namespace buffer_math

#endif

