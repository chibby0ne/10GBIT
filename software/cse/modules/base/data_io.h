//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the data I/O classes of the modules
/// \author Timo Lehnigk-Emden
/// \date   2010/11/05
//
#ifndef DATA_IO_H_
#define DATA_IO_H_

#include "../assistance/buffer.h"
#include "data_io_abstract.h"
#include "data_io_list.h"

namespace cse_lib {
/// Input data port class of the modules
/** 
 * The parameter \tparam T determines the type of the output port.
 * The parameter \tparam DIM the dimension of the output port.
 * The class provides a () operator to get easy access the the data.
 * The data are stored in a Buffer class instance.
 *
 * \ingroup interface
 */
template <class T, unsigned int DIM = 1> class Data_In : public Data_IO_Abstract
{

public:

	Data_In()
	{
		input_port_modified_ = 0;
		input_ = &dummy_buffer;
		input_buffer_valid_ = false;
	}

	virtual ~Data_In()
	{
		// Unregister port if input buffer is still valid (exists).
		if(input_buffer_valid_ == true && input_port_modified_ != 0)
		{
			input_->Unregister_Flag((*input_port_modified_));
		}

	}

	/// Return a reference to the Buffer instance
	Buffer<T, DIM>& operator()() const
	{
		return (*input_);
	}


	/// Connect the object to a Buffer instance
	/**
	 * \param in_data Buffer instance to connect with the input port
	 */
	void operator() (Buffer<T, DIM>& in_data)
	{

		if(input_port_modified_ != 0)
		{
			// Unregister module by the old buffer
			if(input_buffer_valid_ == true)
				input_->Unregister_Flag((*input_port_modified_));

			// Register module at the new buffer
			in_data.Register_Flag((*input_port_modified_), input_buffer_valid_);
			(*input_port_modified_) = true;
		}
		input_ = &in_data;
		input_buffer_valid_ = true;
	}


	/// Connect the object to a buffer instance
	/**
	 * \param in_buffer Pointer to Buffer instance to connect with the input port
	 */
	void Set(Buffer_Abstract* in_buffer)
	{
		input_ = dynamic_cast<Buffer<T, DIM>* > (in_buffer);
		if(input_ == 0)
		{
			if(in_buffer == 0)
				Msg(ERROR, "Data_In::Set", "Given pointer to input buffer pointer has no object instance!");
			else
				Msg(ERROR, "Data_In::Set", "Given buffer (name: " + in_buffer->instance_name() +
				    ") cannot be assigned because it has the wrong type!");

			input_ = &dummy_buffer;
		}
	}

	/// Return the pointer to the connect Buffer instance.
	Buffer_Abstract* Get()
	{
		return dynamic_cast<Buffer_Abstract*> (input_);
	}

	/// Register function
	void Register(string name, Data_IO_List &list)
	{
		this->instance_name(name);
		list.Append(this);
		input_port_modified_= & (list.port_modified());
	}

private:

	/// Pointer to hold the address of the connected buffer.
	Buffer<T, DIM>* input_;

	/// Pointer to a flag to recognize changes of the input buffer port
	bool *input_port_modified_;

	/// Dummy buffer to avoid dereference of pointer to zero
	Buffer<T, DIM> dummy_buffer;

	/// Status flag to indicate if a valid input buffer is present
	bool input_buffer_valid_;
};


/// Output data port class of the modules.
/** 
 * The parameter \tparam T determines the type of the input port.
 * The class provides a () operator to get easy access to the data.
 * The data are stored in a Buffer class instance.
 *
 * \ingroup interface
 */
template <class T, int dim = 1> class Data_Out : public Data_IO_Abstract
{

public:

	Data_Out() {}

	virtual ~Data_Out() {}

	/// Return a reference to the Buffer instance.
	Buffer<T, dim>& operator()()
	{
		return output_;
	}

	/// Return a pointer to the Buffer instance.
	Buffer_Abstract* Get()
	{
		return dynamic_cast<Buffer_Abstract*> (&(output_));
	}

	/// Initialization function to register the port in the list and set the name of it.
	void Register(string name, Data_IO_List &list)
	{
		this->instance_name(name);
		list.Append(this);
	}

private:

	/**
	 * Must be implemented because of virtual abstract base class.  Not
	 * accessible from outside, because setting an output port does not make
	 * sense, because the module manage the instance themselves.
	 */
	void Set(Buffer_Abstract* in){}

	/// Buffer object to store the output data.
	Buffer<T, dim> output_;
};
}
#endif
