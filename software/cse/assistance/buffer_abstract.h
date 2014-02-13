//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief Abstract buffer class interface definition
/// \author Timo Lehnigk-Emden
/// \date 2010/11/03
//
#ifndef BUFFER_ABSTRACT_H_
#define BUFFER_ABSTRACT_H_

#include <string>
#include <deque>

/// Abstract class, Only to provide a common type to store buffer classes in a list
/**
 * \ingroup base_class
 * \ingroup assistance
 */
class Buffer_Abstract
{
public:
	Buffer_Abstract()
	{
		instance_name_.clear();
	}

	virtual ~Buffer_Abstract()
	{
#ifndef DISABLE_EVENT_SUPPORT_FLAG
		Delete_Event();
#endif
	}

	/// Return the name of ther buffer instance
	std::string& instance_name() {return instance_name_;}

	/// Set the name of the buffer instance
	void instance_name(std::string instance_name) { instance_name_ = instance_name;}

#ifndef DISABLE_EVENT_SUPPORT_FLAG
	/// Register a boolean variable (flag) to set if buffer was modified (i.e., resized)
	bool& Register_Flag(bool& modify_flag, bool& valid_flag)
	{
		flag_list_.push_back(&modify_flag);
		valid_list_.push_back(&valid_flag);
		return *valid_list_.back();
	}

	/// Unregister a boolean variable, the variable will not be modified by the buffer object anymore
	void Unregister_Flag(bool& flag)
	{
		for(unsigned int i = 0; i < flag_list_.size(); i++)
		{
			if(flag_list_[i] == &flag)
			{
				flag_list_.erase(flag_list_.begin()+i);
				valid_list_.erase(valid_list_.begin()+i);
			}
		}
	}
#endif

	/// Enum type to determine the internal memory structure in case of creating a new buffer.
	enum BUFFER_TYPE_ENUM {COPY_BUFFER, VIEW_ON_BUFFER};
	/**
	 * \var COPY_BUFFER
	 * The buffer copies the data from the given buffer to its own memory.
	 * \var VIEW_ON_BUFFER
	 * The buffer uses the same memory than the given buffer (view on memory):
	 */

protected:

	/// Type of the list to store the views on the buffer
	typedef std::deque<Buffer_Abstract*> View_List;

	/// Type of the iterator for the View_List
	typedef std::deque<Buffer_Abstract*>::iterator View_List_Iterator;

	/// List to store the views on the instance.
	mutable View_List view_list_;

	/// Register the given buffer object as view of this instance.
	void Register_As_View(Buffer_Abstract* const view) const
	{
		view_list_.push_back(view);
	}

	/// Unregister a view from the current instance.
	void Unregister_As_View(Buffer_Abstract* const view)
	{
		View_List_Iterator iter;
		for(iter = view_list_.begin(); iter != view_list_.end(); iter++)
		{
			if((*iter) == view)
			{
				view_list_.erase(iter);
				break;
			}
		}
	}

	/// Trigger all view, in case of the original buffer has been changed.
	void Trigger_Views()
	{
		Buffer_Abstract::View_List_Iterator iter;
		for(iter = view_list_.begin(); iter != view_list_.end(); iter++)
			(*iter)->View_Parent_Changed();
	}

	/// Virtual Changed function which should be overwritten by the child classe.
	virtual void View_Parent_Changed() {}


#ifndef DISABLE_EVENT_SUPPORT_FLAG
	/// Set all registered objects to true (modified event).
	void Modify_Event()
	{
		for(bool_list_iter_ i = flag_list_.begin(); i != flag_list_.end(); i++)
			(*(*i)) = true;
	}

	/// Notify all connected flags if the buffer object is deleted (delete event).
	void Delete_Event()
	{
		for(bool_list_iter_ i = valid_list_.begin(); i != valid_list_.end(); i++)
			(*(*i)) = false;
	}

	/// Variable to store the flags for the event support
	std::deque<bool*> flag_list_;

	/// Type of the iterator for a deques list of type bool (can be used for the valid_list_ and the flag_list_)
	typedef std::deque<bool*>::iterator bool_list_iter_;

	/// list to store the status of the registered flags (valid or invalid)
	std::deque<bool*> valid_list_;

#endif

	/// Name of the Current Buffer instance
	std::string instance_name_;

	/// Flag to store weather the instance is a view on a Buffer or a Buffer object
	bool is_view_;
};
#endif
