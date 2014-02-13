//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Specialized implementation of the buffer class with dimension 2
/// \author Timo Lehnigk-Emden
/// \date   2010/11/02
//

#include "buffer_dim2.h"

// Specialized class for dim = 2 (2 dimension array)
template <class T> Buffer<T,2>::Buffer()
{
	dim_[0] = 0;
	dim_[1] = 0;
	data_ptr_ = 0;
	is_view_ = false;
}

template <class T> Buffer<T,2>::Buffer(unsigned int Dim1, unsigned int Dim2)
{
	dim_[0] = Dim1;
	dim_[1] = Dim2;

	data_ptr_ = Alloc_Mem(dim_[0], dim_[1]);
	is_view_ = false;
}

// Copy constructor
template <class T> Buffer<T,2>::Buffer(const Buffer<T,2>& buffer)
{
	dim_[0] = 0;
	dim_[1] = 0;
	data_ptr_ = 0;
	*this = buffer;
	is_view_ = false;
}

template <class T> Buffer<T,1>* Buffer<T,2>::Alloc_Mem(unsigned int Rows, unsigned int Cols)
{
	Buffer<T,1>* data_ptr_temp;

	try
	{
		data_ptr_temp = new Buffer<T,1>[Rows];
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "\nError [Buffer<T,2>]: Could not allocate memory correctly of instance "
		          << instance_name_ << "!\n"
		          "Requested dimension was "<< dim_[0] << "x" << dim_[1] << "!\n\n" << std::endl;
		throw;
	}

	for(unsigned int i = 0; i < Rows; i++)
	{
		data_ptr_temp[i].Resize(Cols);
	}
	return data_ptr_temp;
}

template <class T> void Buffer<T,2>::Free_Mem()
{
	if(is_view_ == false)
	{
		delete[] data_ptr_;
	}
	data_ptr_ = 0;
}


template <class T> Buffer<T,2>::~Buffer()
{
	Free_Mem();
}

template <class T> void Buffer<T,2>::Resize(unsigned int dim1, unsigned int dim2, bool copy)
{
	// Only resizing if any dimension has changed
	if(dim_[0] != dim1 || dim_[1] != dim2)
	{
		Buffer<T,1>* new_mem = Alloc_Mem(dim1, dim2);

		if (copy == true)
			for (unsigned int i = 0; i < dim_[0]; i++)
				for (unsigned int k = 0; k < dim_[1]; k++)
					new_mem[i][k] = data_ptr_[i][k];

		// Delete old memory if exist
		if(data_ptr_ != 0)
			delete[] data_ptr_;

		data_ptr_ = new_mem;
		dim_[0] = dim1;
		dim_[1] = dim2;

#ifndef DISABLE_EVENT_SUPPORT_FLAG
		Modify_Event();
#endif
	}
}


template <class T> Buffer<T,2>& Buffer<T,2>::operator=(const Buffer<T,2> &input)
{
	if(dim_[0] != input.dim1() || dim_[1] != input.dim2())
		Resize(input.dim1(),input.dim2());

	for(unsigned int i = 0; i < dim_[0]; i++)
		for(unsigned int k = 0; k < dim_[1]; k++)
			data_ptr_[i][k] = input[i][k];

	return *this;
}

// Copy the values from an 2 dim array to the buffer, no boundary check of input buffer
template <class T> Buffer<T,2> &Buffer<T,2>::operator=(T **input)
{
	for(unsigned int i = 0; i < dim_[0]; i++)
		for(unsigned int k = 0; k < dim_[1]; k++)
			data_ptr_[i][k] = input[i][k];

	return *this;
}

template <class T> void Buffer<T,2>::Write(T **input)
{
	(*this) = input;
}

template <class T> void Buffer<T,2>::Read(T **output)
{
	for(unsigned int i = 0; i < dim_[0]; i++)
		for(unsigned int k = 0; k < dim_[1]; k++)
			output[i][k] = data_ptr_[i][k];
}

template <class T> void Buffer<T,2>::Dump(std::string filename, bool append)
{
	// If filename is empty (standard value)
	if(filename.size() == 0)
		std::cout << (*this);
	else
	{
		std::ofstream dump_file;
		try
		{
			if(append == false)
				dump_file.open(filename.c_str(), std::ios::trunc);
			else
				dump_file.open(filename.c_str(), std::ios_base::app);
		}
		catch(...)
		{
			std::clog << "\nWarning [buffer_dim2.hpp:" << instance_name_ << "]: Output file: " << filename << "cannot be written !\n" << std::endl;
			throw;
		}

		dump_file << (*this);
		dump_file.close();
	}
}

template <class T> void Buffer<T,2>::Dump_Modulo(int mod_bw, std::string filename, bool append, const char* separator)
{
	// If filename is empty (standard value)
	if(filename.size() == 0)
		std::cout << Stream_Modulo(mod_bw, separator);
	else
	{
		std::ofstream dump_file;
		try
		{
			if(append == false)
				dump_file.open(filename.c_str(),std::ios::trunc);
			else
				dump_file.open(filename.c_str(),std::ios_base::app);
		}
		catch(...)
		{
			std::clog << "\nWarning[buffer_dim2.hpp:" << instance_name_ << "]: Output file: " << filename << "cannot be written !\n" << std::endl;
			throw;
		}

		dump_file << Stream_Modulo(mod_bw, separator);
		dump_file.close();
	}
}

template <class T> std::string Buffer<T,2>::Stream_Modulo(int mod_bw, const char* separator)
{
	std::ostringstream os;
	Buffer<int,2> *b = dynamic_cast<Buffer<int,2>*> (this);
	if (b)
	{
		int mod_min  = -(1<<(mod_bw-1));
		int mod_max  =  (1<<(mod_bw-1))-1;
		int mod_comp =  (1<<(mod_bw-1));
		int val;
		for(unsigned int i = 0; i < dim_[0]; i++)
		{
			for(unsigned int j = 0; j < dim_[1]; j++)
			{
				val = (*b)[i][j];
				if (val & mod_comp) os << (val | mod_min) << separator;
				else os << (val & mod_max) << separator;
			}
			os << "\n";
		}
	}
	else
	{
		os << "ERROR. Modulo streaming only supported for Buffer<int,2>." << separator;
	}
	return os.str();
}

template <class T> void Buffer<T,2>::Check_Internal_Dimensions()
{
	if(dim_[0] > 0)
		dim_[1] = (*this)[0].length();

	for(unsigned int i = 0; i < dim_[0]; i++)
		if((*this)[0].length() != dim_[1])
			std::clog << "\nWarning[buffer_dim2.hpp:" << instance_name_ << "]:"
			          << "One dimensional buffer objects in current buffer instance "
			          << "have not the same length!" << std::endl;
}

// Section for special functions of it++ library types
#ifdef ITPP_SUPPORT
template <class T> Buffer<T,2>& Buffer<T,2>::operator=(itpp::Mat<T> &input)
{
	for(int i = 0; i < input.rows(); i++)
		for(int k = 0; k < input.cols(); k++)
			(*this)[i][k] = input(i,k);

	return *this;
}

template <class T> void Buffer<T,2>::Write(itpp::Mat<T> &input)
{
	*this=input;
}

template <class T> void Buffer<T,2>::Read(itpp::Mat<T> &output)
{
	if(dim_[0] > (unsigned int)output.rows() && dim_[1] > (unsigned int)output.cols())
		output.set_size(dim_[0],dim_[1]);


	for(unsigned int i = 0; i < dim_[0]; i++)
		for(unsigned int k = 0; k < dim_[1]; k++)
			output(i,k) = (*this)[i][k];
}

template <class T> itpp::Mat<T> Buffer<T,2>::Read()
{
	itpp::Mat<T> output;
	Read(output);
	return output;
}
#endif

// Special view functions
template <class T> Buffer<T,2>::Buffer(const Buffer<T,2>& buffer, unsigned int start_pos,
                                       unsigned int length, Buffer_Abstract::BUFFER_TYPE_ENUM buffer_type)
{
	// Check parameter
	if(start_pos + length > buffer.dim_[0])
	{
		std::cerr << "\nError [Buffer<T,2>]: Cannot create view on given buffer "
		             "with the given dimensions. Given Buffer has the dimension "
		             << dim_[0] << " x " << dim_[1] <<". "
		             << "View on that buffer was requested at the starting position " << start_pos
		             << " and a length of " << length << std::endl;
		throw;
	}

	if(length == 0)
		length = buffer.dim1() - start_pos;

	if(buffer_type == Buffer_Abstract::COPY_BUFFER)
	{
		for(unsigned int i = 0; i < length; i++)
			(*this)[i] = buffer[i + start_pos];

		is_view_ = false;
	}
	else if(buffer_type == Buffer_Abstract::VIEW_ON_BUFFER)
	{
		data_ptr_ = buffer.data_ptr_ + start_pos;
		dim_[0] = length;
		dim_[1] = buffer.dim_[1];
		is_view_ = true;
	}
}
