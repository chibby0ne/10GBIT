//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Specialized implementation of the buffer class for dimension = 1
/// \author Timo Lehnigk-Emden
/// \date   2010/11/02
//


#include "buffer_dim1.h"

template <class T>
Buffer<T,1>::Buffer()
{
	length_ = 0;
	mem_ptr_ = 0;
	view_ = false;
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;
}


template <class T>
Buffer<T,1>::Buffer(unsigned int size)
{
	mem_ptr_ = 0;
	length_ = size;
	view_ = false;
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;

	Alloc_Mem();
}


// Copy constructor
template <class T>
Buffer<T,1>::Buffer(const Buffer<T,1>& buffer)
{
	length_ = 0;
	mem_ptr_ = 0;
	view_ = false;
	*this = buffer;
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;
}


template <class T>
Buffer<T,1>::Buffer(const Buffer<T,1>& buffer, Buffer_Abstract::BUFFER_TYPE_ENUM buffer_type)
{
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;

	if(buffer_type == VIEW_ON_BUFFER)
		view_ = true;
	else
		view_ = false;

	if(view_ == true)
	{
		Transform_To_View(buffer, 0, 0);
	}
	else
	{
		length_= 0;
		mem_ptr_ = 0;
		*this = buffer;
	}
}


template <class T>
Buffer<T,1>::Buffer(const char* init_values)
{
	length_ = 0;
	mem_ptr_ = 0;
	view_ = false;
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;

	std::string str_tmp = init_values;
	Set_From_String(str_tmp);
}


template <class T>
void Buffer<T,1>::Alloc_Mem()
{
	try
	{
		mem_ptr_ = new T[length_];
	}
	catch(std::bad_alloc&)
	{
		std::cerr << "\nError[buffer.cpp]:Could not allocate memory correctly of instance" <<
		             instance_name_ << " !\n Requested length was "<< length_ << "!\n\n" << std::endl;
		throw;
	}
	Clear();
}


template <class T>
void Buffer<T,1>::Resize(unsigned int size, bool copy, bool error_out)
{
	if (!view_)
	{
		T* temp_ptr;

		// Make resize only if size has changed
		if (size != length_)
		{
			try
			{
				temp_ptr = new T[size];
			}
			catch (std::bad_alloc&)
			{
				std::cerr << "\nError [Buffer<T,1>]: Could not allocate memory correctly of instance"
				          << instance_name_ << " !\n Requested length was "
				          << length_ << "!\n" << std::endl;
				throw ;
			}

			if(copy == true)
			{
				if(size < length_)
				{
					if(error_out == true)
					{
						std::clog << "\nWarning [Buffer<T,1>]: Resizing of instance " << instance_name_
						          << ", new buffer is smaller than the old one, "
						          << "not all data can be copied!\n" << std::endl;
					}

					length_ = size;
				}

				// Copy the old data into the new memory
				for(unsigned int i = 0; i < length_; i++)
					temp_ptr[i] = mem_ptr_[i];

				// Clear new memory
				for(unsigned int i = length_; i < size; i++)
					temp_ptr[i] = static_cast<T> (0);
			}

			// Delete old memory
			delete [] mem_ptr_;

			// Set pointer to new memory
			mem_ptr_ = temp_ptr;

			// Set correct size
			length_ = size;

			// Clear new memory if values where not copied
			if(!copy)
				Clear();

#ifndef DISABLE_EVENT_SUPPORT_FLAG
			/*
			 * Call event function of connected and registered modules only if
			 * compiled with MODULES_EVENT_SUPPORT only !!!
			 */
			Modify_Event();
#endif
			Trigger_Views();
		}
	}
	else
	{
		std::cerr << "\nError [Buffer<T,1>]: instance "<< instance_name_
		          << "cannot be resized, because it's only a reference object to a buffer!\n"
		          << std::endl;
		throw;
	}
}

template <class T>
Buffer<T,1>::~Buffer()
{
	if(view_ == false)
	{
		if(mem_ptr_ != NULL)
		{
			delete[] mem_ptr_;
			mem_ptr_ = NULL;
		}
	}
	else
	{
		// Unregister this instance if a parent exits
		if(view_parameter_.parent != 0)
			view_parameter_.parent->Unregister_As_View(this);
	}

	// Mark all view children as invalid with a size of zero.
	Buffer_Abstract::View_List_Iterator iter;
	for(iter = view_list_.begin(); iter != view_list_.end(); iter++)
	{
		(static_cast<Buffer<T, 1>* >(*iter))->view_parameter_.parent = 0;
		(static_cast<Buffer<T, 1>* >(*iter))->length_ = 0;
	}
}


template <class T>
Buffer<T,1>& Buffer<T,1>::operator=(const Buffer<T,1>& input)
{
	// Resize the current buffer instance to the source buffer size
	if(length_ != input.length())
		Resize(input.length());

	// Copy the data
	for(unsigned int i = 0; i < length_; i++)
		mem_ptr_[i] = input[i];

	return *this;
}


template <class T>
Buffer<T,1>& Buffer<T,1>::operator=(T* input)
{
	for(unsigned int i = 0; i < length_; i++)
		mem_ptr_[i] = input[i];

	return *this;
}


template <class T>
Buffer<T,1>& Buffer<T,1>::operator=(const std::string& input)
{
	Set_From_String(input);
	return *this;
}


template <class T>
void Buffer<T,1>::Set_From_String(const std::string& input)
{
	std::string result_temp;      // String between two subsequent commas
	size_t pos_first, pos_second; // First and second position of subsequent commas
	unsigned int i = 0;           // Iteration variable for the buffer length

	pos_first = 0;

	do
	{
		// Find following comma.
		pos_second = input.find(',', pos_first);

		// Get String part before the following comma.
		result_temp = input.substr(pos_first, pos_second - pos_first);

		// Update the beginning pointer to the first position after the comma.
		pos_first = pos_second + 1;

		// Resize if needed
		if(i >= length_)
			Resize(i + 1000, true);

		// Read single value if no ":" or "-" is present in the string
		hlp_fct::str::Convert_Strict(result_temp, (*this)[i]);
		i++;
	}
	while (pos_second != std::string::npos);

	Resize(i, true, false);
}


template <class T>
inline void Buffer<T,1>::Write(T *input)
{
	(*this) = input;
}


template <class T>
inline void Buffer<T,1>::Write_Sub_Vector(T* input,int start, int end)
{
	if (end == -1)
		end = length_;

	for(int i=start;i<=end;i++)
		mem_ptr_[i]=input[i-start];
}


template <class T>
inline void Buffer<T,1>::Read_Sub_Vector(T* output,int start, int end)
{
	if (end == -1)
		end = length_ - 1;

	for(int i=start;i<=end;i++)
		output[i-start]=mem_ptr_[i];
}


template <class T>
inline int Buffer<T,1>::Compare(const Buffer<T,1> &comp_buffer)
{
	int difference_counter = 0;

	for(unsigned int i=0; i < length_ ; i++)
		if((*this)[i] != comp_buffer[i])
			difference_counter++;

	return difference_counter;
}


template <class T>
inline bool Buffer<T,1>::operator==(const Buffer<T,1> &input)
{
	if(Compare(input) == 0)
		return true;
	else
		return false;
}


template <class T>
void Buffer<T,1>::Append(T input)
{
	this->Resize(this->length() + 1, true);
	(*this)[this->length() - 1] = input;
}


template <class T>
inline void Buffer<T,1>::Read(T* output)
{
	Read_Sub_Vector(output,0,-1);
}


template <class T>
inline T* Buffer<T,1>::Data_Ptr()
{
	return mem_ptr_;
}


template <class T>
inline const T* Buffer<T,1>::Data_Ptr() const
{
	return mem_ptr_;
}


// Section for special functions of it++ library types
#ifdef ITPP_SUPPORT
template <class T>
inline void Buffer<T,1>::Read(itpp::Vec<T> &output)
{
	if(length_ > (unsigned int)output.size())
		output.set_size(length_);

	for(unsigned int i=0;i<length_;i++)
		output[i]=mem_ptr_[i];
}


template <class T>
inline itpp::Vec<T> Buffer<T,1>::Read()
{
	itpp::Vec<T> output;
	Read(output);
	return output;
}


template <class T>
inline itpp::Vec<T> Buffer<T,1>::Read_Sub_Vector(int start, int end)
{
	itpp::Vec<T> output;
	ReadSubVector(output, start, end);
	return output;
}


template <class T>
inline void Buffer<T,1>::Read_Sub_Vector(itpp::Vec<T> &output,int start, int end)
{
	// end value == -1 needs all values to the end
	if (end == -1)
		end = length_-1;

	if((end-start+1) > output.size())
		output.set_size(end-start+1);

	for(int i=start;i<=end;i++)
		output[i-start]=mem_ptr_[i];
}


template <class T>
inline void Buffer<T,1>::Write(itpp::Vec<T> &input)
{
	*this=input;
}


template <class T>
inline Buffer<T,1>& Buffer<T,1>::operator=(itpp::Vec<T> &input)
{
	if(length_ < (unsigned int) input.size())
		Resize(input.size());
	for(int i=0;i<input.size();i++)
		mem_ptr_[i]=input[i];

	return *this;
}


template <class T>
inline Buffer<T,1>& Buffer<T,1>::operator=(itpp::bvec &input)
{
	if(length_ < (unsigned int)input.size())
		Resize(input.size());
	for(int i=0;i<input.size();i++)
		mem_ptr_[i] = static_cast<bool> (input[i]);

	return *this;
}


template <class T>
inline void Buffer<T,1>::Write_Sub_Vector(itpp::Vec<T> &input, int start, int end)
{
	if (end == -1)
		end = length_-1;

	if((unsigned int) (end-start+1) > length_)
		Resize(end-start+1);

	for(int i=start;i<=end;i++)
		mem_ptr_[i]=input[i-start];
}


template <class T>
template<typename T_IT> void Buffer<T,1>::Read_Convert(itpp::Vec<T_IT>& destination)
{
	// Resize destination
	destination.set_length(length_);

	// Copy and convert the data
	for(unsigned int i = 0; i < length_; i++)
		destination[i] = static_cast<T_IT> (mem_ptr_[i]);
}


template <class T>
template<typename T_IT> void Buffer<T,1>::Write_Convert(const itpp::Vec<T_IT>& source)
{
	// Resize buffer
	Resize(source.length());

	// Copy and convert the data from the IT++ container into the buffer
	for(unsigned int i = 0; i < length_; i++)
		mem_ptr_[i] = static_cast<T_IT> (source[i]);
}

#endif

template <class T>
void Buffer<T,1>::Dump(std::string filename, bool append, bool linebreak)
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
			std::cout << "\nWarning [buffer.cpp:" << instance_name_ << "]: Output file: " << filename << "cannot be written !\n" << std::endl;
			throw;
		}

		if (!linebreak)
			dump_file << (*this);
		else
		{
			for(unsigned int i = 0; i < length_; i++)
				dump_file << mem_ptr_[i] << "\n";
		}
		dump_file.close();
	}
}


template <class T>
void Buffer<T,1>::Dump_Modulo(int mod_bw, std::string filename, bool append, const char* separator)
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
			std::clog << "\nWarning [Buffer<T,1>]: Instance: " << instance_name_ << ". "
			          << "Output file: " << filename << "cannot be written !\n" << std::endl;
			throw;
		}
		dump_file << Stream_Modulo(mod_bw, separator);
		dump_file.close();
	}
}


template <class T>
std::string Buffer<T,1>::Stream_Modulo(int mod_bw, const char* separator)
{
	std::ostringstream os;
	Buffer<int,1> *b = dynamic_cast<Buffer<int,1>*> (this);
	if (b)
	{
		int mod_min  = -(1<<(mod_bw-1));
		int mod_max  =  (1<<(mod_bw-1))-1;
		int mod_comp =  (1<<(mod_bw-1));
		int val;
		for(unsigned int i = 0; i < length_; i++)
		{
			val = (*b)[i];
			if (val & mod_comp)
				os << (val | mod_min) << separator;
			else
				os << (val & mod_max) << separator;
		}
	}
	else
	{
		os << "Error [Buffer<T,1>]: Modulo streaming only supported for Buffer<int>." << separator;
	}
	return os.str();
}


template <class T>
T Buffer<T,1>::Max()
{
	if(length_ > 0)
	{
		T max_val = (*this)[0];
		for(unsigned int i = 1; i < length_; i++)
		{
			if((*this)[i] > max_val)
				max_val = (*this)[i];
		}
		return max_val;
	}
	else
		return 0;
}

// **********
// Section for the view specific functions
// **********

template <class T>
Buffer<T,1>::Buffer(const Buffer<T,1>& src_buffer, unsigned int start_pos,
                    unsigned int length, BUFFER_TYPE_ENUM buffer_type)
{
	view_parameter_.length_ = 0;
	view_parameter_.parent = 0;
	view_parameter_.start_pos_ = 0;

	if (buffer_type == COPY_BUFFER)
	{
		view_ = false;
		this->Resize(length);
		for(unsigned int i = 0; i < length; i++)
			(*this)[i] = src_buffer[start_pos + i];
	}
	else if(buffer_type == VIEW_ON_BUFFER)
	{
		Transform_To_View(src_buffer, start_pos, length);
	}
}


template <class T>
void Buffer<T,1>::Transform_To_View(const Buffer<T,1>& src_buffer, unsigned int start_pos,
                                    unsigned int length)
{
	view_parameter_.length_ = length;
	view_parameter_.start_pos_ = start_pos;

	// Unregister and register current instance if parent has changed
	if(view_)
		if(view_parameter_.parent != &(src_buffer))
		{
			if(view_parameter_.parent != 0)
				view_parameter_.parent->Unregister_As_View(static_cast<Buffer_Abstract*> (this));

			src_buffer.Register_As_View(static_cast<Buffer_Abstract*> (this));
		}

	// Store new parent.
	view_parameter_.parent = &(const_cast<Buffer<T,1>& >(src_buffer));

	// If given length = 0 = default value, view goes up to the end of the source buffer
	if(length == 0)
		length = src_buffer.length() - start_pos;

	// Check Range
	if (length + start_pos > src_buffer.length())
	{
		std::cerr << "\nError [Buffer<T,1>]:  instance "<< instance_name_
		          << ": The size of the buffer view ("<< length << ") "
		          "cannot be larger than the size of the referenced sub buffer ("
		          << src_buffer.length()-start_pos << ")!\n" << std::endl;
		throw;
	}

	// Free memory if needed
	if(view_ == false)
	{
		delete[] mem_ptr_;
		mem_ptr_ = 0;

		// Register this instance as view of the given buffer
		src_buffer.Register_As_View(static_cast<Buffer_Abstract*> (this));
	}

	view_ = true;

	// Call event function of connected and registered modules only if length or memory address have changed!
#ifndef DISABLE_EVENT_SUPPORT_FLAG
	if(mem_ptr_ != const_cast<T*> (src_buffer.Data_Ptr()) + start_pos || length_ != length)
		Modify_Event();
#endif

	if(mem_ptr_ != const_cast<T*> (src_buffer.Data_Ptr()) + start_pos || length_ != length)
		Trigger_Views();

	// Update view
	mem_ptr_ = const_cast<T*> (src_buffer.Data_Ptr()) + start_pos;
	length_ = length;
}

template <class T>
void Buffer<T,1>::View_Parent_Changed()
{
	// Update view
	this->Transform_To_View(*view_parameter_.parent, view_parameter_.start_pos_, view_parameter_.length_);
}
