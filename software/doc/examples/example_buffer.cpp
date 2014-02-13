//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Example how to use the Buffer class
/// \author Timo Lehnigk-Emden
/// \date   2010/12/01
//
#define ITPP_SUPPORT

#include "buffer.h"
#include <itpp/itbase.h>

using namespace itpp;

int main()
{
	// Create new variable and arrays
	int a,b;
	int c[10] = {0};

	// Create itpp Vec<int> object of size 10
	ivec itv(10);
	itv.zeros();

	// Create new buffer of type int buffer_int of size 20
	Buffer<int> buffer_int(20);

	// Create new buffer of type float buffer_float of size 100
	Buffer<float,1> buffer_float(100);

	// Set / Write single element of the buffer
	buffer_int[1] = 11;

	// Print single value
	cout << buffer_int[1] << endl;

	itv[0] = 2;

	// Copy values from itv into the buffer
	buffer_int = itv;

	// Copy values from the buffer to itv
	buffer_int.Read(itv);

	// Returns a new vector which contains the buffer values
	itv = buffer_int.Read();

	// Print the whole buffer
	cout << buffer_int;

	// Write the elements from itv into the buffer
	buffer_int.Write(itv);

	// Create two dimensional buffer
	Buffer<int, 2> dfs(20, 20);

	// Set single value
	dfs[2][2] = 5;

	// Print single value
	cout << dfs[2][2];

	// Create it++ mat object of size 10x10
	Mat<int> matt(10, 10);

	// Copy values of mat object into the buffer
	dfs=matt;

	// Create several it++ data types
	bvec bin;
	ivec itb;
	cvec vtest;

	// Create complex buffer
	Buffer<complex<double> > cbuffer(20);

	// Copy values from it++ objects into the buffer
	cbuffer = vtest;
	buffer_int = bin;
	buffer_int = itb;

	return 0;
}
