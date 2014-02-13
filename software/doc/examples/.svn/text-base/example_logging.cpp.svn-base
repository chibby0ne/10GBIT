//
//  Copyright (C) 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  CSE logging example
/// \author Timo Lehnigk-Emden
/// \date   2011/02/03
//


#include "cse_lib.h"

using namespace cse_lib;

int main()
{
	// Create an instance
	Logging_Management log;

	// Disable the output of the info messages and all data written to cout!
	log.Disable(INFO);

	// Redirect the output of the warning messages and all data written to clog to the file warn.log!
	log.Redirect(WARNING, "warn.log");

	// (Re)Enable the output of the info messages and all data written to cout!
	log.Enable(INFO);

	// Reset the output behavior for all error messages and all data written to cerr!
	log.Reset(ERROR);

	// Normal program
}
