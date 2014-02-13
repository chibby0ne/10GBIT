/**
 *
 * \file debug_macros.h
 * \author Timo Lehnigk-Emden (code, documentation), Daniel Schmidt (code, documentation), Christian de Schryver (documentation)
 * \version 0.7
 * \date 2011-10-19
 *
 * \brief A list of useful C++ debug macros as a substitute for cout debugging and logging.
 * It allows control of logging and debug message generation an a per-module level basis.
 *
 * Of nearly all macros here exist two forms, one of which starts with DM_ the other one with LM_.
 * The DM_ macros are for debugging purposes and stream the generated messages to cerr, whereas
 * the LM_ macros are for logging and output to cout by standard. The DM_ macros will also
 * output filename and line-number information for every message.
 *
 * The LM_OUT_* and DM_OUT_*-macros output to their default streams. There exists
 * a second form of all these macros starting with LM_STREAM_* and DM_STREAM_* respectively.
 * These macros take as first argument a stream object to stream the generated messages to.
 * This way you can easily log to a file stream, for example.
 *
 * Logging can be globally turned off by \#define NLOGGING.
 *
 * Debugging can be globally turned of by \#define NDEBUG.
 *
 * Debug and logging messages can be generated on a per-module level basis. Module names can be
 * arbitrarily chosen; all DM_*LEVEL and LM_*LEVEL macros for one module name can be switched on and
 * off by setting the debug level or logging level for this specific module (e.g. in a seperate header 
 * file containing global parameters for debug levels and logging levels):
 *
 * \code
 * \#define DEBUG_LEVEL_<modulename> \<debug_level_set\>
 * \#define LOG_LEVEL_<modulename> \<logging_level_set\>
 * \endcode
 * sets the debug level for the module \<modulename\> to \<debug_level_set\>. Now all messages for \<modulename\> set
 * for debug levels equal or less to \<debug_level_set\> are generated. Same holds for the logging level.
 *
 * An expect macro (DM_EXPECT) is provided as well, which works similar to an
 * assertion without stopping execution of the program, however. 
 */

#ifndef __DEBUG_MACROS_H__
#define __DEBUG_MACROS_H__

/**\brief Standard output stream for DM_* macros */
#define DM_STD_OUT_STREAM std::cerr
/**\brief Standard output stream for LM_* macros */
#define LM_STD_OUT_STREAM std::cout

/// convert x to string for DM_* macros
#define DM_STRINGIFY(x) #x
/// convert x to string for LM_* macros
#define LM_STRINGIFY(x) #x

/// empty message for DM_* macros
#define DM_NULL_MSG ""
/// empty message for LM_* macros
#define LM_NULL_MSG ""

// lots of helper and formatting macros
#define DM_STR_WHERE __FILE__ << ", " << __LINE__

#define DM_STR_LEVEL(module, level) "Module " << #module << "(" << level << "): "
#define DM_FORMAT_LEVEL_START(module, level) DM_STR_WHERE << ": " << DM_STR_LEVEL((module), (level))
#define DM_FORMAT_LEVEL_END(module, level) std::endl

#define DM_STR_COND(cond) #cond << ": " << (cond) << ": "
#define DM_FORMAT_COND_START(cond) DM_STR_WHERE << ": " << DM_STR_COND(cond)
#define DM_FORMAT_COND_END(cond) std::endl

#define DM_STR_VAR(var) #var << " = " << var
#define DM_FORMAT_VAR_START(var) DM_STR_WHERE << ": " << DM_STR_VAR(var) << " : "
#define DM_FORMAT_VAR_END(var) std::endl

#define DM_EXPECTATION_FAILED_MSG(cond) DM_STR_WHERE << " -- Expected: (" << #cond << ") failed!" << std::endl

#define DM_FORMAT_LEVEL(module, level, message) DM_FORMAT_LEVEL_START((module), (level)) << message << DM_FORMAT_LEVEL_END((module), (level))
#define DM_FORMAT_VAR(var, message) DM_FORMAT_VAR_START(var)  << message << DM_FORMAT_VAR_END(var)
#define DM_FORMAT_COND(cond, message) DM_FORMAT_COND_START(cond) << message << DM_FORMAT_COND_END(cond)

#define LM_STR_WHERE __FILE__ << ", " << __LINE__

#define LM_STR_LEVEL(module, level) "Module " << #module << "(" << level << "): "
#define LM_FORMAT_LEVEL_START(module, level) LM_NULL_MSG
#define LM_FORMAT_LEVEL_END(module, level) LM_NULL_MSG

#define LM_STR_COND(cond) #cond << ": " << (cond) << ": "
#define LM_FORMAT_COND_START(cond) LM_NULL_MSG
#define LM_FORMAT_COND_END(cond) LM_NULL_MSG

#define LM_STR_VAR(var) #var << " = " << var
#define LM_FORMAT_VAR_START(var) LM_STR_VAR(var) << " "
#define LM_FORMAT_VAR_END(var) LM_NULL_MSG

#define LM_FORMAT_LEVEL(module, level, message) LM_FORMAT_LEVEL_START((module), (level)) << message << LM_FORMAT_LEVEL_END((module), (level))
#define LM_FORMAT_VAR(var, message) LM_FORMAT_VAR_START(var) << message << LM_FORMAT_VAR_END(var)
#define LM_FORMAT_COND(cond, message) LM_FORMAT_COND_START(cond) << message << LM_FORMAT_COND_END(cond)



// These are the macros intended for use
/// output \<message\> for debugging purposes
#define DM_OUT(message) DM_STREAM(DM_STD_OUT_STREAM, message)
/// output \<message\> for debugging purposes in module \<module\> if debug level for this module is greater or equal to \<level\>
#define DM_OUT_LEVEL(module, level, message) DM_STREAM_LEVEL(DM_STD_OUT_STREAM, module, level, message)
/// output \<var\> = \<value of var\> for debugging purposes
#define DM_OUT_VAR(var) DM_STREAM_VAR(DM_STD_OUT_STREAM, var)
/// output \<var\> = \<value of var\> and \<message\> for debugging purposes
#define DM_OUT_VAR_MSG(var, message) DM_STREAM_VAR_MSG(DM_STD_OUT_STREAM, var, message)
/// output \<message\> for debugging purposes if \<cond\> is true
#define DM_OUT_COND(cond, message) DM_STREAM_COND(DM_STD_OUT_STREAM, cond, message)
/// output \<message\> for debugging purposes, formatted, with filename, position, and line break
#define DM_OUT_POS(message) DM_STREAM(DM_STD_OUT_STREAM, DM_STR_WHERE << ": " << message << std::endl)

/// output \<message\> for debugging purposes to \<stream\> in module \<module\> if debug level for this module is greater or equal to \<level\>
#define DM_STREAM_LEVEL(stream, module, level, message) DM_LEVEL_EXEC(module, level, DM_STREAM(stream, DM_FORMAT_LEVEL(module, level, message)))
/// output \<var\> = \<value of var\> to \<stream\> for debugging purposes
#define DM_STREAM_VAR(stream, var) DM_STREAM_VAR_MSG(stream, var, DM_NULL_MSG)
/// output \<var\> = \<value of var\> and \<message\> to \<stream\> for debugging purposes
#define DM_STREAM_VAR_MSG(stream, var, message) DM_STREAM(stream, DM_FORMAT_VAR(var, message))
/// output \<message\> to \<stream\> for debugging purposes if \<cond\> is true
#define DM_STREAM_COND(stream, cond, message) DM_COND_EXEC(cond, DM_STREAM(stream, DM_FORMAT_COND(cond, message)))

/// like assert (cond), but without terminating the program
#define DM_EXPECT(cond) DM_COND_EXEC(!(cond), DM_OUT(DM_EXPECTATION_FAILED_MSG(cond)))


/// output \<message\> for logging purposes
#define LM_OUT(message) LM_STREAM(LM_STD_OUT_STREAM, message)
/// output \<message\> for logging purposes in module \<module\> if debug level for this module is greater or equal to \<level\>
#define LM_OUT_LEVEL(module, level, message) LM_STREAM_LEVEL(LM_STD_OUT_STREAM, module, level, message)
/// output \<var\> = \<value of var\> for logging purposes
#define LM_OUT_VAR(var) LM_STREAM_VAR(LM_STD_OUT_STREAM, var)
/// output \<var\> = \<value of var\> and \<message\> for logging purposes
#define LM_OUT_VAR_MSG(var, message) LM_STREAM_VAR_MSG(LM_STD_OUT_STREAM, var, message)
/// output \<message\> for logging purposes if \<cond\> is true
#define LM_OUT_COND(cond, message) LM_STREAM_COND(LM_STD_OUT_STREAM, cond, message)

/// output \<message\> for logging purposes to \<stream\> in module \<module\> if debug level for this module is greater or equal to \<level\>
#define LM_STREAM_LEVEL(stream, module, level, message) LM_LEVEL_EXEC(module, level, LM_STREAM(stream, LM_FORMAT_LEVEL(module, level, message)))
/// output \<var\> = \<value of var\> to \<stream\> for logging purposes
#define LM_STREAM_VAR(stream, var) LM_STREAM_VAR_MSG(stream, var, LM_NULL_MSG)
/// output \<var\> = \<value of var\> and \<message\> to \<stream\> for logging purposes
#define LM_STREAM_VAR_MSG(stream, var, message) LM_STREAM(stream, LM_FORMAT_VAR(var, message))
/// output \<message\> to \<stream\> for logging purposes if \<cond\> is true
#define LM_STREAM_COND(stream, cond, message) LM_COND_EXEC(cond, LM_STREAM(stream, LM_FORMAT_COND(cond, message)))



#ifdef NLOGGING
	#define LM_STREAM(stream, message)
	
	#define LM_EXEC(statement)
	#define LM_COND_EXEC(cond, statement)
	#define LM_LEVEL_EXEC(module, level, statement)
	#define LM_EXEC_NOBLOCK(statement)
#else
	#include <iostream>
	#define LM_STREAM(stream, message) do { stream << message; } while (0)
	
	#define LM_EXEC(statement) do { statement; } while (0)
	#define LM_COND_EXEC(cond, statement) do { if ((cond)) do { statement; } while (0); } while (0)
	#define LM_LEVEL_EXEC(module, level, statement) do { if ((LOG_LEVEL_ ## module) >= level) do { statement; } while (0); } while (0)
	/// unsafe version of LM_EXEC, statement is not in a block
	#define LM_EXEC_NOBLOCK(statement) statement
#endif //def NLOGGING


#ifdef NDEBUG
	#define DM_STREAM(stream, message)
	
	#define DM_COND_EXEC(cond, statement)
	#define DM_LEVEL_EXEC(module, level, statement)
#else
	#include <iostream>
	#define DM_STREAM(stream, message) do { stream << message; } while(0)
	
	#define DM_COND_EXEC(cond, statement) do { if ((cond)) do { statement; } while (0); } while (0)
	#define DM_LEVEL_EXEC(module, level, statement) do { if ((DEBUG_LEVEL_ ## module) >= level) do { statement; } while (0); } while (0)
#endif //def NDEBUG

#endif //def __DEBUG_MACROS_H__

