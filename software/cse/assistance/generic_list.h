//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Template class definition for a list of objects in a container from the STL library
/// \author Timo Lehnigk-Emden
/// \date   2010/11/03
//
#ifndef GENERIC_LIST_H_
#define GENERIC_LIST_H_

/// Template class to manage objects
/**
 * \tparam Cont can be any container class from the STL library, which can store objects of type T
 * (e.g., vector ord dequeue).
 * \tparam T can be any type which can be stored in the STL list types
 *
 *  This class as an internal list to manage the inserted objects.
 *  A copy of the given object is stored.
 *
 *  \ingroup base_class
 *  \ingroup assistance
 */
template <class T, template <class T, typename> class Cont> class Generic_List
{

public:

	/// Create a new list
	Generic_List()
	{
		obj_list_.clear();
	}

	virtual ~Generic_List() {};

	/// Append a new object to the end of the list
	void Append(T* const obj)
	{
		if(obj != 0)
			obj_list_.push_back(*obj);
	}

	/// Append a new object copy to the end of the list
	void Append(const T& obj)
	{
		obj_list_.push_back(obj);
	}

	/// Append a new object to the end of the list
	/*void Add(T obj)
	{
		obj_list_.push_back(obj);
	}*/

	/// Delete all list entries
	void Clear()
	{
		obj_list_.clear();
	}

	/// Insert a new element at the position index
	void Insert(const T &Obj, unsigned int index)
	{
		obj_list_.insert(obj_list_.begin()+index,Obj);
	}

	// Not needed because of const reference in Insert function
	/// Insert a new element at the position index
	/*void Insert_Alt(T Obj, unsigned int index)
	{
		obj_list_.insert(obj_list_.begin()+index,Obj);
	}*/

	/// Return the reference of the object with the index i out of the list
	/// Return a reference to the object with the index i out of the list
	T& Get_Ref(int i)
	{
		return (*(obj_list_.begin() + i));
	}

	/// Return a copy of the object with the index i out of the list
	T operator[](int i) const
	{
		return (*(obj_list_.begin() + i));
	}

	/// Return the index of the object with the given value
	int Search(const T& value)
	{
		// Iterate over the list and search for the value
		/*iter_ = obj_list_.begin();
		int index=0;
		do
		{
			if((*iter_) == value)
				return index;
			iter_++;
			index ++;
		}
		while(iter_ != obj_list_.end());*/

		int index=0;
		for(iter_ = obj_list_.begin(); iter_ != obj_list_.end(); iter_++)
		{
			if((*iter_) == value)
				return index;
			index ++;
		}

		// Return -1 if no object was found
		return -1;
	}

	/// Return the number of objects pointers in the list
	unsigned int length() const
	{
		return obj_list_.size();
	}

	/// stream all objects in the list
	inline friend std::ostream & operator<<(std::ostream& os, Generic_List<T,Cont>& data)
	{
		for(data.iter_ = data.obj_list_.begin(); data.iter_ != data.obj_list_.end() - 1; data.iter_++)
			os << (*data.iter_) << ", ";
		os << (*data.iter_);
		return os;
	}

	/// Iterator type definition
	typedef typename Cont< T, std::allocator < T > >::iterator iterator;

	/// Reverse Iterator type definition
	typedef typename Cont< T, std::allocator < T > >::reverse_iterator reverse_iterator;

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

protected:

	// Internal parameters

	/// List variable
	Cont< T, std::allocator < T > > obj_list_;

	/// Iterator over the list
	iterator iter_;
};
#endif
