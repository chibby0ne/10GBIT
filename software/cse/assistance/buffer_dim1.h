//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Specialization of the buffer class for dim = 1
/// \author Timo Lehnigk-Emden
/// \date   2011/06/10
//

#ifndef BUFFER_DIM1_H_
#define BUFFER_DIM1_H_

#include "buffer_abstract.h"

// Enable IT++ Support and library dependency
#ifdef ITPP_SUPPORT
#include <itpp/itbase.h>
#endif

/// Buffer class with one dimension (array)
/**
 * Buffer class to store data into a matrix array.
 * Class is a specialization of the buffer<T,dim> class with dim = 1.
 * In addition if the ITPP_SUPPORT definition is set the class member functions support
 * the Vec<T> data type of the it++ library.
 * Class as an explicit it++ support, that mean to operate with the Vec<T>
 * type from it++ instantiate the class with the type T and NOT with the type
 * Vec<T> !
 * In addition it has many functions to read and write data to native C++ arrays.
 *
 * \ingroup base_class
 * \ingroup assistance
 */
template <class T> class Buffer<T,1> : public Buffer_Abstract
{

public:

	/// Buffer View, a new buffer object is created which uses the memory of the given buffer
	Buffer(const Buffer<T,1>& buffer, Buffer_Abstract::BUFFER_TYPE_ENUM buffer_type);

	/// Copy constructor
	Buffer(const Buffer<T,1>& buffer);

	/// Create a new buffer with zero elements
	Buffer();

	/// Create a new buffer and initialize it with the given values
	explicit Buffer(const char* init_value);

	/// Create a new buffer with a given size
	/**
	 * @param size Number of elements which can the buffer store
	 */
	explicit Buffer(unsigned int size);

	/// Creates a partial view or a copy of the given buffer
	Buffer(const Buffer<T,1>& buffer, unsigned int start_pos, unsigned int length, Buffer_Abstract::BUFFER_TYPE_ENUM buffer_type);

	/// Destructor
	virtual ~Buffer();

	/// Append the given value to the buffer and return a buffer reference.
	Buffer<T,1>& operator()(T value)
	{
		this->Append(value);
		return (*this);
	}

	/// Resize the buffer from the current size to a new size
	/**
	 * \param size New size of the buffer
	 * \param copy Holds the data if copy=true, if the new buffer size is smaller than the previous one, the oversized date were discarded.
	 * \param error_out Enable or disable the error messages. An error can occurs if you want to keep the values,
	 *  but the new buffer size is smaller than the old one.
	 */
	void Resize(unsigned int size, bool copy = false, bool error_out = true);

	/// Transform the current buffer instance to a (partial) view of another buffer.
	/**
	 * \param src_buffer Buffer to create the view on.
	 * \param start_pos Starting position of the view of the source buffer
	 * \param length Length of the view (number of elements).
	 *  = 0: Create a view up to the end of the original buffer
	 */
	void Transform_To_View(const Buffer<T,1>& src_buffer, unsigned int start_pos = 0, unsigned int length = 0);


	/// Get a single element of the buffer (constant)
	/**
	 * @param i position to set or get, beginning with zero
	 */
	const T& operator[](int i) const { return mem_ptr_[i]; }
	//T operator[](int i) const { return mem_ptr_[i]; }


	/// Get or Set a single element of the buffer
	/**
	 * @param i position to set or get, beginning with zero
	 */
	T& operator[](int i) { return mem_ptr_[i]; }


	/// Copy the values from buffer to buffer
	/**
	 * @param input buffer to copy.
	 */
	Buffer<T,1> &operator=(const Buffer<T,1>& input);

	/// Set the value in the buffer from the given string
	/**
	 * The format is val1, val2, val3. The operator throws the invalid_argument exception in
	 * case of an conversion error.
	 * The buffer is automatically resized to the number of elements read
	 */
	Buffer<T,1>& operator=(const std::string& input);


	/// Write date into the buffer
	/**
	 * @param input Pointer to data to write into the buffer. Number of
	 * elements in input must be equal or greater than the buffer.
	 */
	Buffer<T,1> &operator=(T *input);

	/// Set all elements to zero (clear the buffer)
	void Clear()
	{
		for(unsigned int i = 0; i < length_; i++)
			(*this)[i] = static_cast<T>(0);
	}

	/// Set all elements to the given value.
	void Set(T input)
	{
		for(unsigned int i = 0; i < length_; i++)
			(*this)[i] = static_cast<T>(input);
	}

	/// Get the pointer of the internal memory of the buffer, only for advanced users
	T* Data_Ptr();

	/// Get the pointer of the internal memory of the buffer, only for advanced users
	const T* Data_Ptr() const;


	/// Write data into the buffer, same as = operator
	/**
	 * @param input Pointer the array to copy the data. Number of elements in input must be equal or greater than the buffer.
	 */
	void Write(T *input);


	/// Write the Sub array/vector into the buffer
	/**
	 * @param[in] input Pointer to the date to write into the buffer
	 * @param start Index of the buffer to start writing
	 * @param end Index of the buffer to stop writing the data, if =-1 copy all data up to the end
	 */
	void Write_Sub_Vector(T* input,int start, int end=-1);


	/// Read the whole buffer
	/**
	 * @param[out] output Pointer to the array to copy the date form the buffer into the array. output array must be great enough to store all data.
	 */
	void Read(T* output);


	/// Read the sub array/vector out of the buffer.
	/**
	 * @param[out] output Pointer to the array to copy the date from the buffer into the array. output array must be great enough to store all data.
	 * @param start Index of the buffer to start reading
	 * @param end Index of the buffer to stop reading, if =-1 copy all data up to the end
	 */
	void Read_Sub_Vector(T* output, int start, int end=-1);

	/// Compare the value of the current buffer with another one and return the number of differences
	int Compare(const Buffer<T,1> &buffer);

	/// Compare the value of the current buffer with another one
	bool operator==(const Buffer<T,1> &input);

	/// Append a single element to the Buffer (Resize included).
	void Append(T input);

	// Functions for it++ Vec<T> data type
#ifdef ITPP_SUPPORT

	/**
	 * For use with the it++ data type Vec<T>, see operator =()
	 * \overload
	 */
	Buffer<T,1> &operator=(itpp::Vec<T> &input);

	/**
	 * For use with the it++ data type bvec, see operator =()
	 * \overload
	 */
	Buffer<T,1> &operator=(itpp::bvec &input);

	/**
	 * For use with the it++ data type itpp::Vec<T>, see Read()
	 * \overload void Read(itpp::Vec<T> &output)
	 */
	void Read(itpp::Vec<T> &output);

	/**
	 * For use with the it++ data type itpp::Vec<T>, see Read_Sub_Vector()
	 * \overload void Read_Sub_Vector(itpp::Vec<T>& output, int start, int end=-1);
	 */
	void Read_Sub_Vector(itpp::Vec<T> &output,int start, int end=-1);


	/// Create a new vector which contains a copy of all data of the buffer
	/**
	* @return New vector of type itpp::Vec<T>
	*/
	itpp::Vec<T> Read();


	/// Return a new vector which contains parts of the buffer values
	/**
	 * @param start First index of data to read from buffer
	 * @param end Last index of data to read
	 * @return New vector of type itpp::Vec<T> which contains the requested date
	 */
	itpp::Vec<T> Read_Sub_Vector(int start, int end=-1);

	/// Read the converted values from buffer and stores them into the given itpp vec container.
	/**
	 * \param destination itpp vec container to write the converted data to type T_IT from the buffer.
	 * The size of the destination container is adapted automatically to the size of the buffer.
	 */
	template<typename T_IT> void Read_Convert(itpp::Vec<T_IT>& destination);

	/**
	 * For use with the it++ data type itpp::Vec<T>, see Write()
	 * \overload void Write(itpp::Vec<T> &input);
	 */
	void Write(itpp::Vec<T> &input);


	/**
	 * For use with the it++ data type itpp::Vec<T>, see Write_Sub_Vector()
	 * \overload void Write_Sub_Vector(itpp::Vec<T> &input, int start, int end=-1);
	 */
	void Write_Sub_Vector(itpp::Vec<T> &input, int start, int end = -1);

	/// Write the converted values from the given itpp vec container into the buffer.
	/**
	 * \param source itpp container to read the data, convert them to the buffer type T and write them into the buffer.
	 * The size of the buffer is adapted automatically to the size of the source container.
	 */
	template<typename T_IT> void Write_Convert(const itpp::Vec<T_IT>& source);

#endif

	/// Overload of << operator
	inline friend std::ostream & operator<<(std::ostream & os, Buffer<T,1> & data)
	{
		if(data.length_ != 0)
		{
			for(unsigned int i=0; i < data.length_ - 1; i++)
				os << data.mem_ptr_[i] << ", ";

			// Print last data without comma and whitespace at the end
			os << data.mem_ptr_[data.length_ - 1];
		}
		return os;
	}

	/// Overload of >> operator
	inline friend std::istream & operator>>(std::istream & os, Buffer<T,1> & data)
	{
		for(unsigned int i = 0; i < data.length_; i++)
			os >> data.mem_ptr_[i];

		return os;
	}

	/// Dump the content of the buffer to stdout or a file.
	/**
	 * @param filename If a file name is give, the buffer is written to a file
	 * @param append = false, file content is overwritten, = true append the date to the file
	 * @param linebreak = false, values separated by blank, = true values separated by linebreak
	 */
	void Dump(std::string filename="", bool append=false, bool linebreak=false);


	/// Dump the content of the buffer modulo 2^mod_bw to stdout or a file.
	/**
	 * @param filename If a file name is give, the buffer is written to a file
	 * @param append = false, file content is overwritten, = true append the date to the file
	 * @param mod_bw    All values are dumped modulo 2^mod_bw
	 * @param separator Values are separated by this string. Default is a blank.
	 */
	void Dump_Modulo(int mod_bw, std::string filename="", bool append=false, const char* separator=" ");


	/// Stream the content of the buffer modulo 2^mod_bw
	/**
	 * @param mod_bw    All values are dumped modulo 2^mod_bw
	 * @param separator Values are separated by this string. Default is a blank.
	 */
	std::string Stream_Modulo(int mod_bw, const char* separator=" ");

	/// Return the number of elements (length) of the buffer
	unsigned int length() const {return length_;}

	/// Return a Buffer with the current dimensions of the buffer
	Buffer<unsigned int, 1> dim()
	{
		Buffer<unsigned int, 1> ret(1);
		ret[0] = length_;
		return ret;
	}

	/// Return the maximum value in the buffer
	T Max();

private:

	/// Set the buffer values from a string
	void Set_From_String(const std::string& input);

	void View_Parent_Changed();

	/// Name of the current buffer instance
	std::string instance_name_;

	/// Internal memory allocation function
	void Alloc_Mem();

	/// Store the length of the current memory
	unsigned int length_;

	/// Pointer to the memory
	T *mem_ptr_;

	/// Buffer View of buffer
	bool view_;

	/// Original view parameter
	class View_Parameter
	{
	public:
		unsigned int length_;   /// length of the view
		unsigned int start_pos_;  /// start position of the view
		Buffer<T, 1>* parent;   /// Pointer to the original buffer object
	} view_parameter_;

};
#endif
