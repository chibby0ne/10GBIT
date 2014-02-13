//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the generic buffer class
/// \author Timo Lehnigk-Emden
/// \date   2011/06/10
//

#ifndef BUFFER_DIMX_H_
#define BUFFER_DIMX_H_

#include "buffer_abstract.h"
#include "buffer_dim1.h"
/// Buffer class with type and optional dimension as template parameter.
/**
 * Buffer class stores data in a multidimensional array. The size of the
 * buffer can be adjusted during runtime. Single values or all data of the
 * buffer can be read or written. You can directly use the buffer at the
 * beginning or/and end of your operations to read/write your values without
 * copying them before or after the operation.
 *  \tparam <class T> Data type of data which can be stored in class.
 *  \tparam <int dim> Number of dimension of the buffer
 *  \ingroup base_class
 *  \ingroup assistance
 */
template <class T, int DIM> class Buffer : public Buffer_Abstract
{
public:

	/// Copy constructor
	explicit Buffer(const Buffer<T,DIM>& buffer);

	/// Create a new buffer with 0 elements
	Buffer();

	// Destructor
	virtual ~Buffer();

	/// Resize the buffer from the current size to a new size, all current data are discarded
	/**
	 * @param size New size of the buffer for each dimension
	 */
	/*void Resize(unsigned int size0, ...);*/

	///\overload
	void Resize(const unsigned int size[DIM]);

	///\overload
	void Resize(Buffer<unsigned int,1> size);

	/// Get or Set a single element of the buffer
	/**
	 * @param i position to set or get, beginning with zero
	 */
	Buffer<T,DIM-1>& operator[](int i) { return data_ptr_[i]; }

	/// Get a single element of the buffer (constant)
	/**
	 * @param i position to set or get, beginning with zero
	 */
	const Buffer<T,DIM-1>& operator[](int i) const { return data_ptr_[i]; }

	/// Return a Buffer with the current dimensions of the buffer
	Buffer<unsigned int, 1> dim()
	{
		Buffer<unsigned int, 1> ret(DIM);
		ret = dim_;
		return ret;
	}

	/// Streaming operator
	inline friend std::ostream & operator<<(std::ostream & os, Buffer<T,DIM> & data)
	{
		for(unsigned int i=0;i<data.dim_[0];i++)
			os << data.data_ptr_[i] << " ";

		return os;
	}

	/// Overload of >> operator
	inline friend std::istream & operator>>(std::istream & os, Buffer<T,DIM> & data)
	{
		for(unsigned int i=0;i<data.dim_[0];i++)
			os >> data.data_ptr_[i];

		return os;
	}


private:

	// Frees the internal memory of the buffer
	void Free_Mem();

	// Internal memory allocation function
	Buffer<T,DIM-1>* Alloc_Mem(int dimension[DIM]);

	// Array to store the dimensions
	unsigned int dim_[DIM];

	Buffer<T,DIM-1>* data_ptr_;
};
#endif
