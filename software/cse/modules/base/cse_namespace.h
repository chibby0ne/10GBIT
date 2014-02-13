//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Include file for the cse_lib namespace definition
/// \author Timo Lehnigk-Emden
/// \date   2010/12/13
//
#ifndef CSE_NAMESPACE_H_
#define CSE_NAMESPACE_H_

// Include common header files
#include <complex>
#include <string>
#include <iostream>
#include <vector>
#include "../assistance/buffer.h"
#include "../assistance/helper_functions.h"
#include "../assistance/generic_list.h"

/// Namespace for the Creonic simulation environment (CSE)
/**
 * This namespace contains all modules and their parent classes of the CSE.
 * The standard namespace and the namespace of the helper functions is included in this namespace!
 */
namespace cse_lib
{
	// Make std and helper functions namespaces visible in cse_lib namespace.
	using namespace std;
	using namespace hlp_fct;
	using namespace hlp_fct::logging;
	using namespace hlp_fct::math;
	using namespace hlp_fct::bit_operations;
	using namespace hlp_fct::str;

	// **********
	// Type definitions
	// **********

	/// Abbreviation for the 64-bit unsigned integer type from C++11
	typedef unsigned long long int ull_int;

	/// Abbreviation for the 64-bit signed integer type from C++11
	typedef signed long long int ll_int;

	/// Abbreviation for the 32-bit unsigned integer
	typedef unsigned int u_int;

	/// Enum type for the return values of the Run functions
	/**
	 * This type is used as return type of the Run function of the modules.
	 * The module indicates whether and which stopping event occurs
	 * (e.g., iteration control or all simulation finished).
	 */
	enum RETURN_VALUE_RUN_ENUM {NO_STOPPING_EVENT = 0, STOP_CURRENT_SIMULATION_POINT = -1,
	                            STOP_CURRENT_BLOCK = 1, STOP_COMPLETE_SIMULATION = 2};
	/*! \var NO_STOPPING_EVENT
	 *  No stopping event occurs, normal case
	 *  \var STOP_CURRENT_SIMULATION_POINT
	 *  Stop the current simulation point (e.g., all blocks have been simulated)
	 *  \var STOP_CURRENT_BLOCK
	 *  Stop the current processed block (e.g., iteration control has detected a valid code word)
	 *  \var STOP_COMPLETE_SIMULATION
	 *  Stop the complete simulation (e.g., further simulation point does not make sense)
	 */

	/// Enum type to control the internal state reset behavior of the modules
	enum INTERNAL_STATE_RUN_ENUM {RESET = 0, KEEP = 1};
	/*! \var RESET
	 *  Reset all internal states, default behavior
	 *  \var KEEP
	 *  Keeps the internal states untoched. That measn the module begins with the state wich was
	 *  set from the last call of the Run function
	 */

	// *********
	// Data structs and classed
	// *********

	/// Datatype (class) for unmodified strings (raw) to use as type for the Param class.
	class raw_string : public std::string
	{
	public:
		/// Creates an empty raw_string.
		raw_string() {}

		/// Creates an string with the given characters.
		raw_string(const char* s) : std::string(s) {}

		/// Creates an instance with the given string as content.
		raw_string(const std::string& s) : std::string(s) {}
	};

	/// Template class to hold the pairs of enum values and the corresponding string.
	/**
	 * This data structure can be used to store the associates string to enum types.
	 * \tparam T type of the data to associated with a string.
	 */
	template<typename T> class DS_Enum_String
	{
	public:
		/// Associate the given value to the given string. Both values are added to a list.
		DS_Enum_String& operator()(const T& val1, const std::string& val2)
		{
			val1_.Append(val1);
			val2_.Append(val2);
			return (*this);
		}

		/// List to store the enum value.
		Generic_List<T, std::vector> val1_;

		/// List to store the corresponding string the the enum values.
		Generic_List<std::string, std::vector> val2_;
	};
}
#endif
