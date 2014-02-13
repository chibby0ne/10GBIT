//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Date containers used by the status_out class
/// \author Timo Lehnigk-Emden
/// \date   2011/02/04
//

#ifndef STATUS_OUT_CONTAINER_H_
#define STATUS_OUT_CONTAINER_H_


#include "status_out_abstract.h"

namespace cse_lib {

/// Container to store all plug-in instances of one status_out_port
template<class P, unsigned int DIM> class Status_Out_Plugin_Container
{
public:
	Status_Out_Plugin_Container() : dim_(0) {children_ = 0;};

	virtual ~Status_Out_Plugin_Container() { }

	/// Return the stored plug-in container instance with the index i
	/** If the index exceeds the current dimension, the container is automatically extended
	 * and filled with reseted plug-in instances.
	 */
	Status_Out_Plugin_Container<P, DIM-1>& operator()(unsigned int i)
	{
		if(i >= dim_)
			Extent_Dim(i+1);
		return children_[i];
	}

	/// Return the current dimension of the container
	unsigned int dim()
	{
		return dim_;
	}

private:
	/// Return a pointer to the stored container list
	Status_Out_Plugin_Container<P, DIM-1>* children()
	{
		return children_;
	}

	// Free the memory recursively
	void Delete_Recursive()
	{
		// Recursive call to delete all sub children
		for(unsigned int i = 0; i < dim_; i++)
			children_[i].Delete_Recursive();

		// Delete all own children
		delete[] children_;
		children_ = 0;
	}

	void Extent_Dim(int new_dim)
	{
		// Extent own dimension
		Status_Out_Plugin_Container<P, DIM-1>* tmp;

		try
		{
			tmp = new Status_Out_Plugin_Container<P, DIM-1>[new_dim];
		}
		catch(bad_alloc&)
		{
			Msg(ERROR,"Status_Out_Plugin_Container","Memory allocation failure !"
			    "Requested size was " + Conv_Num_To_Str(new_dim) + ".");
			throw;
		}

		// Create recursively all children and copy the plug-ins instance
		for(unsigned int i = 0; i < dim_; i++)
			Copy_Elements(children_[i],tmp[i]);

		// Delete old instances
		delete[] children_;

		// Update new dimension
		dim_ = new_dim;

		// Update pointer to the children
		children_ = tmp;
	}


	// Recursive copy call
	template<unsigned int DEEPS> void Copy_Elements(Status_Out_Plugin_Container<P, DEEPS>& source, Status_Out_Plugin_Container<P, DEEPS>& dest)
	{
		// New children
		Status_Out_Plugin_Container<P, DEEPS-1>* dest_temp;
		try
		{
			dest_temp = new Status_Out_Plugin_Container<P, DEEPS-1>[source.dim()];
		}
		catch(bad_alloc&)
		{
			Msg(ERROR,"Status_Out_Plugin_Container","Memory allocation failure !"
			    "Requested size was " + Conv_Num_To_Str(source.dim()) + ".");
			throw;
		}

		// Recursive copy and instanziation
		for(unsigned int i = 0; i < source.dim(); i++)
			Copy_Elements(source.children_[i],dest_temp[i]);

		// Delete old children instances
		delete[] source.children_;
		// Copy pointer to children
		dest.children_ = dest_temp;
		// Update dimension
		dest.dim_ = source.dim();
	}

	// Copy the single plug-in
	void Copy_Elements(Status_Out_Plugin_Container<P, 0>& source, Status_Out_Plugin_Container<P, 0>& dest)
	{
		// Copy the pointer to the plugin instance
		dest.plug_inst_ = source.plug_inst_;
	}
private:
	// Pointer for the children
	Status_Out_Plugin_Container<P, DIM-1>* children_;

	// Dimension, number of the child elements
	unsigned int dim_;

	// Instance with higher dimension are friends
	//friend class Status_Out_Plugin_Container<P, DIM+1>;
	//friend class Status_Out_Plugin_Container<P, DIM+2>;
	template <typename PP, unsigned int DIMMM> friend class Status_Out_Plugin_Container;

	// Make all instance of Status_Out friends
	template <class T, unsigned int DIMM, template <typename TT> class F> friend class Status_Out;
};

/// Specialization for DIM = 0
template<class P> class Status_Out_Plugin_Container<P,0>
{
public:
	Status_Out_Plugin_Container()
	{
		plug_inst_ = 0;
		write_into_results_ = true;
	}

	virtual ~Status_Out_Plugin_Container() {}

	/// Return a reference to the plug-in instance
	P& operator()()
	{
		Extend_Dim();
		return *plug_inst_;
	}

	/// Get the write_into_result flag
	bool write_into_results()
	{
		return write_into_results_;
	}

	/// Change the result file behavior of the plugin instance
	void write_into_results(bool value)
	{
		write_into_results_ = value;
	}

private:
	// Create a new plug-in instance if no instance exists
	void Extend_Dim()
	{
		if(plug_inst_ == 0)
			plug_inst_ = new P;
	}

	// Delete the plug-in instance, should be called by frind classes to free the memory
	void Delete_Recursive()
	{
		delete plug_inst_;
		plug_inst_ = 0;
	}

	/// Pointer to the plug-in instance
	P* plug_inst_;

	//friend class Status_Out_Plugin_Container<P, 1>;
	template <typename PP, unsigned int DIMMM> friend class Status_Out_Plugin_Container;
	template <typename T, unsigned int DIMM, template <typename TT> class F> friend class Status_Out;

	// Flag to indicate if the result from the plugin instance should appear in the result file
	bool write_into_results_;
};

}

#endif
