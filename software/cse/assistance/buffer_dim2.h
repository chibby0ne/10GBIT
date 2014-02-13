//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Specialization of the buffer class for dim = 2
/// \author Timo Lehnigk-Emden
/// \date   2011/06/10
//

#ifndef BUFFER_DIM2_H_
#define BUFFER_DIM2_H_


#include "buffer_abstract.h"

// Enable IT++ Support and library dependency
#ifdef ITPP_SUPPORT
#include <itpp/itbase.h>
#endif

/// Generic data container with a dimensions of two (matrix)
/**
 * Buffer class to store data into a matrix array.
 * Class is a specialization of the buffer<T,dim> class with dim = 2.
 * In addition if the ITPP_SUPPORT
 * definition is set the class member functions support the Mat<T> data type of
 * it++ library.
 *
 * \ingroup base_class
 * \ingroup assistance
 */
template <class T> class Buffer<T,2> : public Buffer_Abstract
{
public:
	/// Create new buffer
	Buffer();

	/// Copy constructor
	Buffer(const Buffer<T,2>& buffer);

	/// Create new 2-dimensional buffer with the given dimensions
	/**
	 * \param dim1 number of entries for the first dimension
	 * \param dim2 number of entries for the second dimension
	 */
	Buffer(unsigned int dim1, unsigned int dim2);

	/// Buffer view constructor (memory is shared with another existing buffer)
	/**
	 * \param buffer Buffer object to create a view on.
	 * \param start_pos starting index of the first dimension for the buffer
	 * \param length length of the view for the first dimension
	 * \param buffer_type determine weather the buffer is a copy or a view of the given buffer.
	 */
	Buffer(const Buffer<T,2>& buffer, unsigned int start_pos, unsigned int length, Buffer_Abstract::BUFFER_TYPE_ENUM buffer_type);

	/// Destructor
	virtual ~Buffer();

	/// Resize the buffer memory
	/**
		 * @param dim1 number of elements of the first dimension (rows)
		 * @param dim2 number of elements of the second dimension (columns)
		 * @param copy if true the buffer keeps its current values
	*/
	void Resize(unsigned int dim1,unsigned int dim2, bool copy=false);

	/// Resize the buffer memory, the two dimension are read from the size array
	void Resize(const unsigned int size[2], bool copy=false)
	{
		Resize(size[0],size[1],copy);
	}

	/// Get or Set single row i, beginning with zero
	Buffer<T,1>& operator[](int i) { return data_ptr_[i]; }

	/// Get single row i, beginning with zero (constant)
	const Buffer<T,1>& operator[](int i) const { return data_ptr_[i]; }

	/// Copy the values from buffer to buffer
	/**
	 * @param input buffer to copy.
	 */
	Buffer<T,2> &operator=(const Buffer<T,2> &input);

	/// Write the date of input into the buffer, input must have equal or more elements than the buffer
	Buffer<T,2> &operator=(T **input);

	/// Same as = operator
	void Write(T **input);

	/// Read the hole buffer into output. Output array must be big enough to store all data.
	void Read(T **output);


	/// Set all elements in the buffer to zero
	void Clear()
	{
		for(unsigned int i = 0; i < dim_[0]; i++)
			(*this)[i].Clear();
	}

	// Functions for it++ Vec<T> data type
#ifdef ITPP_SUPPORT
	/// Copy all elements of the input matrix into the buffer
	Buffer<T,2> &operator=(itpp::Mat<T> &input);

	/// Copy all elements of the input matrix into the buffer
	void Write(itpp::Mat<T> &input);

	/// Copy the data from the buffer to the given output matrix
	void Read(itpp::Mat<T> &output);

	/// Return a new matrix with a copy the data in the buffer
	itpp::Mat<T> Read();
#endif

	/// Overload of << operator
	inline friend std::ostream & operator<<(std::ostream & os,  Buffer<T,2> & data)
	{
		for(unsigned int i = 0; i < data.dim_[0]; i++)
		{
			for(unsigned int j = 0; j < data.dim_[1]; j++)
				os << data.data_ptr_[i][j] << ", ";

			os << "\n";
		}

		return os << std::endl;
	}

	/// Overload of >> operator
	inline friend std::istream & operator>>(std::istream & os,  Buffer<T,2> & data)
	{
		for(unsigned int i = 0; i < data.dim_[0]; i++)
		{
			for(unsigned int j = 0; j < data.dim_[1]; j++)
				os >> data.data_ptr_[i][j];
		}
		return os;
	}


	/// Dump the content of the buffer to stdout or a file, 1dim buffers separated by linebreak, values separated by blank.
	/**
	 * @param filename If a file name is give, the buffer is written to a file
	 * @param append = false, file content is overwritten, = true append the date to the file
	 */
	void Dump(std::string filename="", bool append=false);


	/// Dump the content of the buffer modulo 2^mod_bw to stdout or a file.
	/**
	 * @param filename If a file name is give, the buffer is written to a file
	 * @param append = false, file content is overwritten, = true append the date to the file
	 * @param mod_bw    All values are dumped modulo 2^mod_bw
	 * @param separator Values are separated by this string. Default is a blank. 1dim buffers are separated by linebreak.
	 */
	void Dump_Modulo(int mod_bw, std::string filename="", bool append=false, const char* separator=" ");


	/// Stream the content of the buffer modulo 2^mod_bw
	/**
	 * @param mod_bw    All values are dumped modulo 2^mod_bw
	 * @param separator Values are separated by this string. Default is a blank. 1dim buffers are separated by linebreak.
	 */
	std::string Stream_Modulo(int mod_bw, const char* separator=" ");

	/// Check and update the internal dimensions of the sub buffer objects
	void Check_Internal_Dimensions();

	/// Return the size of the dimension 1 of the buffer
	unsigned int dim1() const {return dim_[0];}

	/// Return the size of the dimension 2 of the buffer
	unsigned int dim2() const {return dim_[1];}

	/// Return a Buffer with the current dimensions of the buffer
	Buffer<unsigned int, 1> dim()
	{
		Buffer<unsigned int, 1> ret(2);
		ret = dim_;
		return ret;
	}

private:
	// Frees the internal memory of the buffer
	void Free_Mem();
	// Name of the Current Buffer instance
	std::string instance_name_;

	// Internal memory allocation function
	Buffer<T,1>* Alloc_Mem(unsigned int Rows, unsigned int Cols);


	// Array to store the dimensions
	unsigned int dim_[2];

	// Pointer to the memory
	Buffer<T,1>* data_ptr_;
};
#endif
