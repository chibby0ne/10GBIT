//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Definition of the status output class for use in the modules
/// \author Timo Lehnigk-Emden
/// \date   2010/12/09
//
#ifndef STATUS_OUT_H_
#define STATUS_OUT_H_

#include "status_out_abstract.h"
#include "status_out_container.h"

// Include all available plug-ins
#include "status_out_plugins.h"
#include "status_out_plugin_mean.h"
#include "status_out_plugin_hist.h"

namespace cse_lib {

/// Output status port class of the modules.
/**
 * \tparam T determines the data type of the status port.
 * \tparam F gives the plug-in type class (e.g, calculate the mean value of all
 * values which will be written. The standard plug-in is the simple data
 * storage.
 * \tparam dim determines the dimension of the port
 *
 * \ingroup interface
 */
template <class T,unsigned int DIM = 0, template <class TT> class F = Status_Out_Plugin_Default>
class Status_Out : public Status_Out_Abstract
{
public:
	Status_Out()
	{
		for(unsigned int i=0; i< DIM; i++)
			dim_name_[i] = string("dim_") + Conv_Num_To_Str(i);
	}

	virtual ~Status_Out()
	{
		// Free Memories of plugin container
		//plugin_inst_.Free_Mem();
		plugin_inst_.Delete_Recursive();
	}

	/// Reset all plug-in values of the port
	void Reset()
	{
		// Reset all instances
		for(unsigned int i=0; i < plugin_inst_.dim(); i++)
			plugin_inst_(i)().Reset();
	}

	/// Return the name of dimension with the index i
	string& dim_name(unsigned int i)
	{
		return dim_name_[i];
	}

	/// Set the names of the port dimensions
	/*
	 * \param i Dimension to set the name
	 * \param str New name of the dimension
	 */
	void dim_name(unsigned int i, string str)
	{
		if(i < DIM)
			dim_name_[i] = str;
	}

	/// Return the plug-in data container, used to address and get access to the plug-in instance
	Status_Out_Plugin_Container<F<T>, DIM - 1>& operator()(int i)
	{
		return plugin_inst_(i);
	}

	/// Read one value into the given data structure, used for automatic processing of the values
	bool Get_Value_Str(Status_Out_Data_Struct& data)
	{
		return Get_Values(plugin_inst_, data);
	}

private:

	// Go thru the dimensions and propagate the values
	template<unsigned int DEEPS> bool Get_Values(Status_Out_Plugin_Container<F<T>,DEEPS > cont, Status_Out_Data_Struct& data)
	{
		bool ret;
		static unsigned int read_cnt = 0;
		do
		{
			data.Clear();

			// Check if children exists (this is the case if the port was written for once.)
			// If no children exist return false to mark that no data are available.
			if(cont.children() != 0)
				ret = Get_Values(cont.children()[read_cnt], data);
			else
				ret = false;

			if(ret == false)
				read_cnt ++;
			else	// Correct address of the data set (Add dimension of the port)
			{
				data.address.Insert(Conv_Num_To_Str(read_cnt),0);
				data.dim_name.Insert(dim_name_[DEEPS-1],0);
			}

		}
		// Stop if current value is valid or all dimension have been read
		while(ret == false && read_cnt < cont.dim());

		// Reset counter if last value from all child instances have been read
		if(ret == false && read_cnt == cont.dim())
			read_cnt = 0;

		return ret;
	}

	bool Get_Values(Status_Out_Plugin_Container<F<T>,0 > cont, Status_Out_Data_Struct& data)
	{
		// Check if plugin instance exist and if values should be written into the result file
		if(cont.plug_inst_ != 0)
		{
			if(cont.write_into_results() == true)
				return cont().Get_Value_Str(data);
		}

		return false;
	}

	// Object instance of the plug-in
	Status_Out_Plugin_Container<F<T>, DIM> plugin_inst_;

	// Name of the dimensions
	string dim_name_[DIM];
};

// Not needed anymore, because of multidimensional template class
// Only for debug reasons !
/// Specialization of Status_out class with template parameter DIM = 1
/**
 * \ingroup interface
 */
//template <class T, template <class T> class F>
//class Status_Out<T,F,1> : public Status_Out_Abstract
//{
//public:
//
//	Status_Out()
//	{
//		dim_name_ = "dim_0";
//		Reset();
//	}
//
//	virtual ~Status_Out()
//	{
//		plugin_inst_.Free_Mem();
//	}
//
//	void Reset()
//	{
//		read_cnt_ = 0;
//		for(unsigned int i=0; i < plugin_inst_.dim(); i++)
//			plugin_inst_(i)().Reset();
//	}
//
//	F<T>& operator()(unsigned int index)
//	{
//		return plugin_inst_(index)();
//	}
//
//	/// Write one value in a data struct including his address
//	/**
//	 * \return true: if value is valid, false: all values have been read
//	 */
//	bool Get_Value_Str(Status_Out_Data_Struct& data)
//	{
//		bool ret = false;
//
//		if(write_into_results_)
//		{
//			do
//			{
//				data.Clear();
//				ret = plugin_inst_(read_cnt_)().Get_Value_Str(data);
//
//				// Go to next port dimension if all value of the current dimension have been read
//				if(ret == false)
//					read_cnt_ ++;
//				else	// Correct address of the data set (Add dimension of the port)
//				{
//					data.address.Insert_Alt(Conv_Num_To_Str(read_cnt_),0);
//					data.dim_name.Insert(dim_name_,0);
//				}
//
//			}
//			// Stop if current value is valid or all dimension have been read
//			while(ret == false && read_cnt_ < plugin_inst_.dim());
//
//			// Find last value read
//			if(ret == false && read_cnt_ == plugin_inst_.dim())
//				read_cnt_ = 0;
//		}
//
//		return ret;
//	}
//
//private:
//	// Object instance of the plug-in
//	Status_Out_Plugin_Container<F<T>, 1> plugin_inst_;
//
//	// Read counter
//	unsigned int read_cnt_;
//
//	// Name of the port dimension
//	string dim_name_;
//};


/// Specialization of Status_out class with template parameter DIM = 0
/**
 * This specialization is needed to adapt the access operator ()
 * \ingroup interface
 */
template <typename T, template <typename TT> class F>
class Status_Out<T, 0, F> : public Status_Out_Abstract
{
public:

	Status_Out() {}

	virtual ~Status_Out()
	{
		//plugin_inst_.Free_Mem();
		plugin_inst_.Delete_Recursive();
	}

	/// Return the stored plug-in instance, used to get access to the plug-in
	F<T>& operator()()
	{
		return plugin_inst_();
	}

	/// Return the value with a given index as string
	bool Get_Value_Str(Status_Out_Data_Struct& data)
	{
		if(write_into_results_)
			return plugin_inst_().Get_Value_Str(data);
		else
		{
			//plugin_inst_().Reset();
			return false;
		}
	}

private:
	// Object instance of the plug-in
	Status_Out_Plugin_Container<F<T>, 0> plugin_inst_;
};
}

#endif


