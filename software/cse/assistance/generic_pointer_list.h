//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Template class definition for a list of pointer in container from the STL library
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef GENERIC_POINTER_LIST_H_
#define GENERIC_POINTER_LIST_H_

#include <iostream>

/// Template class to manage objects (list of object pointers)
/**
 * \tparam Cont can be any class (container) from the STL library, which can store objects pointers
 *  of type T.
 * \tparam T can be any type which can be stored in the STL list types
 *
 *  This class as an internal list to manage the inserted objects.
 *  Only the object pointers are stored. That means no copy will be created,
 *  but the list does not care about stored objects which were deleted!
 *  If you delete an inserted object outside you have to removed it out of the list before.
 *  \ingroup base_class
 *  \ingroup assistance
 */
template <class T, template <class T, typename> class Cont> class Generic_Pointer_List
{

public:

	/// Create a new list
	Generic_Pointer_List()
	{
		obj_list_.clear();
	}

	virtual ~Generic_Pointer_List(){};

	/// Add a new object to the end of the list
	void Append(T* const obj)
	{
		if(obj != 0)
			obj_list_.push_back(obj);
	}

	/// Add a new object to the end of the list
	void Append(T& obj)
	{
		Append(&obj);
	}

	/// Iterator type
	typedef typename Cont< T*, std::allocator < T* > >::iterator iterator;

	/// Iterator type for reverse directions
	typedef typename Cont< T*, std::allocator < T* > >::reverse_iterator reverse_iterator;

	/// Delete all list entries
	void Clear()
	{
		obj_list_.clear();
	}

	/// Insert a new element at the position index
	void Insert(const T &Obj, unsigned int index)
	{
		obj_list_.insert(obj_list_.begin()+index,&Obj);
	}

	/// Return the number of objects pointers in the list
	unsigned int length() const
	{
		return obj_list_.size();
	}

	/// Return the element with the number index
	T* Get(unsigned int index)
	{
		if(index >= length())
			return 0;
		else
			return *(obj_list_.begin()+index);
	}

	/// Return the first iterator object of the list
	iterator begin()
	{
		return obj_list_.begin();
	}

	/// Return the first iterator object of the list for reverse iteration
	reverse_iterator rbegin()
	{
		return obj_list_.rbegin();
	}

	/// Return the last iterator object of the list
	iterator end()
	{
		return obj_list_.end();
	}

	/// Return the last iterator object of the list for reverse iteration
	reverse_iterator rend()
	{
		return obj_list_.rend();
	}

	/// stream all objects in the list
	inline friend std::ostream & operator<<(std::ostream & os, Generic_Pointer_List<T,Cont> & data)
	{
		for(data.iter_ = data.obj_list_.begin();data.iter_ != data.obj_list_.end();data.iter_++)
			os << *(*data.iter_);
		return os;
	}

	/// Public iterator for direct usage
	iterator iter;

	bool operator== (Generic_Pointer_List<T, Cont>& rhs )
	{
		bool ret = true;

		if(rhs.length() != this->length())
			return false;

		rhs.iter = rhs.begin();

		for(iter_ = obj_list_.begin(); iter_ != obj_list_.end(); iter_++)
		{
			if((*iter) != (*rhs.iter))
			{
				ret = false;
				break;
			}
			rhs.iter++;
		}
		return ret;
	}


protected:

	// Internal parameters

	/// List variable
	Cont< T*, std::allocator < T* > > obj_list_;

	/// Iterator over the list
	iterator iter_;

};
#endif
