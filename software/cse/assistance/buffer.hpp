//
//  Copyright (C) 2010 creonic UG (haftungsbeschränkt)
//
//  This file is part of the creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Implementation of the buffer template class
/// \author Timo Lehnigk-Emden
/// \date   2010/11/02
//


#include <iostream>

/**
 * \example example_buffer.cpp
 * Example how to use this class
 */

template <class T,int dim> Buffer<T,dim>::Buffer()
{
	for(int i=0;i < dim;i++)
		dim_[i] = 0;
	data_ptr_ = 0;
}

// Copy constructor
template <class T,int dim> Buffer<T,dim>::Buffer(const Buffer<T,dim>& buffer)
{
	for(int i=0;i < dim;i++)
		dim_[i] = 0;
	data_ptr_ = 0;

	*this = buffer;
}

template <class T,int dim> Buffer<T,dim>::~Buffer()
{
	Free_Mem();
}

template <class T,int dim> void Buffer<T,dim>::Free_Mem()
{
	delete[] data_ptr_;
	data_ptr_ = 0;
}

template <class T,int dim> void Buffer<T,dim>::Resize(const unsigned int size[dim])
{
	// Instantiate new objects with the next smaller dimension
	Buffer<T,dim-1> *temp_data_ptr=0;
	try
	{
		temp_data_ptr = new Buffer<T,dim-1>[size[0]];
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "\nError [Buffer<T,dim>]: Could not allocate memory correctly of instance "
		          << instance_name_ << " !\n" << std::endl;
		throw;
	}

	unsigned int temp_dim[dim-1];

	for(int i=0;i < dim-1;i++)
		temp_dim[i] = size[i+1];

	// Delete old instances (recursively)
	Free_Mem();

	// Call resize function of all new objects
	for(unsigned int i=0;i < size[0];i++)
		temp_data_ptr[i].Resize(temp_dim);

	data_ptr_ = temp_data_ptr;

	for(int i=0;i < dim;i++)
		dim_[i]=size[i];

#ifndef DISABLE_EVENT_SUPPORT_FLAG
	Modify_Event();
#endif
}
