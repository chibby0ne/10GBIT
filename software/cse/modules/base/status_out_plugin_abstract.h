//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Interface definition for the status out plugins
/// \author Timo Lehnigk-Emden
/// \date   2011/02/02
//

#ifndef STATUS_OUT_PLUGIN_ABSTRACT_H_
#define STATUS_OUT_PLUGIN_ABSTRACT_H_

namespace cse_lib
{

/// Abstract class to define the minimal interface of a status out plug-in.
/**
 *
 */
class Status_Out_Plugin_Abstract
{
public:
	virtual ~Status_Out_Plugin_Abstract() {}

	/// Reset the plug-in
	virtual void Reset() = 0;

	/// Function to read all data out of the plug-in. The data are stored in the given data structure
	virtual bool Get_Value_Str(Status_Out_Data_Struct& data) = 0;
};
}
#endif
