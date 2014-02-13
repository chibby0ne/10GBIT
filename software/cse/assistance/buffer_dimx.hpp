//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Implementation of the buffer template class
/// \author Timo Lehnigk-Emden
/// \date   2010/11/02
//



/**
 * \example example_buffer.cpp
 * Example how to use this class
 */

#include "buffer_dimx.h"

template <class T,int DIM> Buffer<T,DIM>::Buffer()
{
	for(int i=0;i < DIM;i++)
		dim_[i] = 0;
	data_ptr_ = 0;
}

// Copy constructor
template <class T,int DIM> Buffer<T,DIM>::Buffer(const Buffer<T,DIM>& buffer)
{
	for(int i=0;i < DIM;i++)
		dim_[i] = 0;
	data_ptr_ = 0;

	*this = buffer;
}

template <class T,int DIM> Buffer<T,DIM>::~Buffer()
{
	Free_Mem();
}

template <class T,int DIM> void Buffer<T,DIM>::Free_Mem()
{
	delete[] data_ptr_;
	data_ptr_ = 0;
}

template <class T,int DIM> void Buffer<T,DIM>::Resize(Buffer<unsigned int,1> size)
{
	if(size.length() != DIM)
	{
		std::cerr << "\nError [Buffer<T,dim>]: Resize function was called with a wrong size "
		             "of parameters!. "
		             "Buffer has the dimension " + hlp_fct::str::Conv_Num_To_Str(DIM); +". "
		             "Resize Buffer has " + hlp_fct::str::Conv_Num_To_Str(size.length()) +
		             " entries.\n";
		throw;
	}
	Resize(size.Data_Ptr());
}

/*template <class T,int dim> void Buffer<T,dim>::Resize(unsigned int size0, ...)
{
	// Extract length information
	unsigned int size[dim];
	va_list vl;
	va_start(vl,size0);
	for (unsigned int i = 1; i < dim; i++)
		size[i] = va_arg(vl, unsigned int);
	va_end(vl);
	Resize(size);
}*/

template <class T,int DIM> void Buffer<T,DIM>::Resize(const unsigned int size[DIM])
{
	// Instantiate new objects with the next smaller dimension
	Buffer<T,DIM-1> *temp_data_ptr=0;
	try
	{
		temp_data_ptr = new Buffer<T,DIM-1>[size[0]];
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "\nError [Buffer<T,dim>]: Could not allocate memory correctly of instance "
		          << instance_name_ << " !\n" << std::endl;
		throw;
	}

	unsigned int temp_dim[DIM-1];

	for(int i=0;i < DIM-1;i++)
		temp_dim[i] = size[i+1];

	// Delete old instances (recursively)
	Free_Mem();

	// Call resize function of all new objects
	for(unsigned int i=0;i < size[0];i++)
		temp_data_ptr[i].Resize(temp_dim);

	data_ptr_ = temp_data_ptr;

	for(int i=0;i < DIM;i++)
		dim_[i]=size[i];

#ifndef DISABLE_EVENT_SUPPORT_FLAG
	Modify_Event();
#endif
}
