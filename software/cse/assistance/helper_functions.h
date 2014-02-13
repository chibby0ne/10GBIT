//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Helper functions declarations / implementation
/// \author Timo Lehnigk-Emden
/// \date   2010/11/09
//

#ifndef HELPER_FUNCTIONS_H_
#define HELPER_FUNCTIONS_H_

#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <fstream>
#include <stdexcept>
#include <typeinfo>

/// Name space of the helper functions
/**
 * \ingroup assistance
 */
namespace hlp_fct
{

/// Name space for C++ types helper functions
namespace types
{
	/// Generates the type name of a C++ type as string
	/**
	 * \tparam T C++ type to generate the type name.
	 */
	template <typename T>
	class TypeName
	{
	public:
		/// Return the type name as string
		static const char* Get() {return typeid(T).name();}
	};

	// \def ENABLE_TYPENAME(A) \internal Specialized classes for a common string of the used types
	#define ENABLE_TYPENAME(A) template<> struct TypeName<A> { static const char* Get() { return #A; }};
	ENABLE_TYPENAME(int)
	ENABLE_TYPENAME(unsigned int)
	ENABLE_TYPENAME(long long int)
	ENABLE_TYPENAME(unsigned long long int)
	ENABLE_TYPENAME(bool)
	ENABLE_TYPENAME(std::string)
	ENABLE_TYPENAME(float)
	ENABLE_TYPENAME(double)
	#undef ENABLE_TYPENAME
}

/// Name space for the string helper function
namespace str
{
	/// Remove all white characters at the beginning and at the end of the string and all newline characters
	/**
	 * \param[in] input string to process
	 * \return return a new string without any newline character in it
	 */
	::std::string Remove_White_Newline(const ::std::string& input);

	/// Convert a string into another type with the stringstream object
	/**
	 * \param[in] input string to convert
	 * \param[out] output object to store the converted value
	 * \return return the bitmask variable which contains stream error state flags.
	 */
	template <class T>
	inline ::std::ios_base::iostate Convert(const ::std::string& input, T& output)
	{
		::std::stringstream Str;
		Str.clear();
		Str << (input);
		Str >> output;

		return Str.rdstate();
	}

	/// Check the IO state from a stream operation an throw an exceptions in case of an error.
	void Check_Io_State(const ::std::ios_base::iostate& conv_state,
	                    const ::std::string& input,
	                    const ::std::string& type_name);

	/// Convert a string into another type with a very strict error handling
	/**
	 * In case of any error during the conversion or in case not all characters can be converted
	 * the function throws the invalid_argument exception from the C++ standard library.
	 */
	template <class T>
	inline void Convert_Strict(const ::std::string& input, T& output)
	{
		::std::ios_base::iostate conv_state;

		// Remove whitespace or newlines and make conversion
		conv_state = Convert(Remove_White_Newline(input), output);

		Check_Io_State(conv_state, input, types::TypeName<T>().Get());
	}

		/// Convert each character in a string to upper case
	inline ::std::string Convert_To_Uppper(const ::std::string& in)
	{
		::std::string out;
		for( unsigned int i=0;i< in.size();i++)
			out += toupper(in[i]);
		return out;
	}

	/// Convert a string contains hexadecimal or octal or decimal number as as string into another type
	/**
	 * A hexadecimal number begins with the "#" sign or the "0x" signs.
	 * An octal numbers begins with the "0o".
	 * \param[in] raw_string string to convert into a number
	 * \param[out] value variable to write the value
	 * \return return the bitmask variable which contains stream error state flags.
	 */
	template <class T>
	::std::ios_base::iostate Convert_Number_String(::std::string raw_string, T& value)
	{
		enum number_format_type {deci,hexa,octa};
		number_format_type number_format = deci;

		// Detect hex values and remove identification
		if(raw_string.length() >= 2)
		{
			if(raw_string.compare(0,2,"0x") == 0)
			{
				raw_string.erase(0,2);
				number_format = hexa;
			}
			else if(raw_string.compare(0,1,"#") == 0)
			{
				raw_string.erase(0,1);
				number_format = hexa;
			}
			else if(raw_string.compare(0,2,"0o") == 0)
			{
				raw_string.erase(0,2);
				number_format = octa;
			}
		}

		::std::stringstream Str;
		Str.clear();
		switch(number_format)
		{
		case hexa:Str << ::std::hex << (raw_string); break;
		case octa:Str << ::std::oct << (raw_string); break;
		default: Str << (raw_string); break;
		}
		Str >> value;

		return Str.rdstate();
	}

	/// Strict version of hlp_fct::str::Convert_Number_String
	/**
	 * The function throws an exception in case of any conversion error.
	 */
	template <class T>
	inline void Convert_Number_String_Strict(const ::std::string& raw_string, T& value)
	{
		::std::ios_base::iostate conv_state;

		// Conversion with out checking the conversion process
		conv_state = Convert_Number_String(Remove_White_Newline(raw_string), value);

		// Check the conversion process. In case of any error the function throws an exception.
		Check_Io_State(conv_state, raw_string, types::TypeName<T>().Get());
	}

	/// Convert the given string to a number
	template<class T>
	inline void Conv_Num_To_Str(T number, std::string& str)
	{
		std::stringstream temp;
		temp.clear();
		str.clear();
		temp << number;
		temp >> str;
	}

	/// Convert the given string to a number and return the string object
	template<class T>
	inline std::string Conv_Num_To_Str(T number)
	{
		std::string str_tmp;
		Conv_Num_To_Str(number, str_tmp);
		return str_tmp;
	}

	/// Convert Number to a string according to the given base.
	/**
	 * \param value Input value.
	 * \param base Base for conversion of number (has to be between 2 and 16).
	 * \param num_digits number of output digits (fill with "0")
	 */
	inline std::string Conv_Num_To_Str_Base(int value, unsigned int base, unsigned int num_digits)
	{
		std::string digits = "0123456789abcdef";
		std::string result;
		if((base < 2) || (base > 16))
		{
			result = "Error [Conv_Num_To_Str_Base]: Base out of range.";
			throw;
		}
		else
		{
			do
			{
				result = digits[value % base] + result;
				value /= base;
				num_digits--;
			}
			while(num_digits);
		}
		return result;
	}


	/// Return the given string with embedded in quotation markers
	inline std::string Quote_Str(const std::string& input)
	{
		return (std::string("\"") + input + std::string("\""));
	}

	/// Give the position of the n-th given sequence in the string
	/**
	 * \param str string to search in
	 * \param search string to search for in the string str
	 * \param num occurrence from 1 to N
	 * Example: str = "abcdabcdabcd"; search = "a"; num = 2; return = 5
	 */
	inline size_t Find_Pos_Of_Nth(const std::string &str, const std::string& search, unsigned int num)
	{
		size_t pos = 0;
		unsigned int num_cnt = 0;
		do
		{
			pos = str.find(search, pos + 1);
			// Check if found
			if(pos != std::string::npos)
				num_cnt++;
		}
		while(pos != std::string::npos && num_cnt < num);

		return pos;
	}

	/// Check if the the given string contains the search string
	inline bool Check_Str_In_Str(const std::string &str, const std::string& search)
	{
		size_t pos = str.find(search);
		if(pos != std::string::npos)
			return true;
		else
			return false;
	}
}
	/// Namespace contains comfort function to print out information, warning and error messages
	namespace logging
	{
		/// Enum for the type of the message to give out
		enum MSG_TYPE_ENUM {INFO, WARNING, ERROR};
		/*!\var INFO
		 * Indicates information messages (print the message to cout)
		 * \var WARNING
		 * Indicates warning messages (print the message to clog)
		 * \var ERROR
		 * Indicates error messages (print the message to cerr)
		 */

		/// Print a formated message with a given type
		/**
		 * \param msg_type Chose the output stream and the message prefix
		 * \param class_name Give the location of the message
		 * \param message Message to print out
		 * \param print_meta if true (=default) the class_name string is added to the output message
		 */
		inline void Msg(MSG_TYPE_ENUM msg_type, std::string class_name, std::string message, bool print_meta=true)
		{
			if(print_meta)
			{
				switch(msg_type)
				{
				case INFO: std::cout << "Info [" << class_name <<"]: " << message << std::endl; break;
				case WARNING: std::clog << "Warning [" << class_name <<"]: " << message << std::endl; break;
				case ERROR: std::cerr << "Error [" << class_name <<"]: " << message << std::endl; break;
				}
			}
			else
			{
				switch(msg_type)
				{
				case INFO: std::cout << message << std::endl; break;
				case WARNING: std::clog << message << std::endl; break;
				case ERROR: std::cerr << message << std::endl; break;
				}
			}
		}

		///\overload
		inline void Msg(MSG_TYPE_ENUM msg_type, std::string class_name, std::stringstream& message, bool print_meta=true)
		{
			Msg(msg_type, class_name, message.str(), print_meta);
		}

		/// Class to manage the output of the standard streams cout, clog and cerr
		/**
		 * This class manage the output behavior for the tree standard output streams
		 * (cout, clog and cerr). Each stream can be disabled, enabled or redirected to a file.
		 *
		 * \ingroup assistance
		 */
		class Logging_Management
		{
		public:
			/// Create a new instance to control the logging levels
			Logging_Management();

			virtual ~Logging_Management();

			/// Enable the output of the message with the given type
			void Enable(MSG_TYPE_ENUM msg_type);

			/// Disable the output of the message with the given type
			void Disable(MSG_TYPE_ENUM msg_type);

			/// Redirect the output of the message with the given type to a file
			void Redirect(MSG_TYPE_ENUM msg_type, std::string file_name, bool append = false);

			/// Reset the output behavior for the given type
			void Reset(MSG_TYPE_ENUM msg_type);

		private:
			/// Defines a "null" ostream, a stream with no output
			class nullstream: public std::ostream
			{
			public:
				nullstream(): std::ios(&null_sbuf_), std::ostream(&null_sbuf_) {}

				virtual ~nullstream() {}
			private:
				// Defines a "null" streambuf used for the ostream
				class nullbuf: public std::streambuf
				{
				protected:
					int overflow(int c) { return traits_type::not_eof(c); }
				}null_sbuf_;
			};

			// Returns a reference to the associated stream object
			std::ostream& Chose_Stream(MSG_TYPE_ENUM msg_type);

			// Returns a number which is associated with the type
			int Get_Type_Number(MSG_TYPE_ENUM msg_type);

			// Null ostream instance
			nullstream null_ostream_;

			// Pointer to store the original streambuf addresses
			std::streambuf* stream_orig_[3];
			// Pointer to store the last streambuf addresses
			std::streambuf* stream_last_[3];
			// Filestream objects
			std::ofstream file_stream_[3];
		};

		/// Print out the plain string onto the given output stream
		/*inline void Msg_Plain(msg_type_enum msg_type, std::string message)
		{
			switch(msg_type)
			{
			case INFO: std::cout << message; break;
			case WARNING: std::clog << message; break;
			case ERROR: std::cerr << message; break;
			}
		}

		inline void Msg_Plain(msg_type_enum msg_type, std::iostream& stream)
		{
			switch(msg_type)
			{
			case INFO: std::cout << stream; break;
			case WARNING: std::clog << stream; break;
			case ERROR: std::cerr << stream; break;
			}
		}*/
	}

	/// Namespace contains mathematical functions
	namespace math
	{
		/// Log function without the singularities, saturates the vales with the min and max values
		/**
		 *
		 * \param x argument of the log function
		 * \return the logarithm of x
		 */
		inline double Trunc_Log(double x)
		{
			const double log_double_max = std::log(std::numeric_limits<double>::max());
			const double log_double_min = std::log(std::numeric_limits<double>::min());

			if (std::numeric_limits<double>::is_iec559)
			{
				if (x == std::numeric_limits<double>::infinity())
					return log_double_max;
				if (x <= 0)
					return log_double_min;
			}
			return std::log(x);
		}


		/// Saturate the given unsigned value in place to the given magnitude.
		inline void Saturate_Value(unsigned int &input, unsigned int max_magnitude)
		{
			if (input > max_magnitude)
				input = max_magnitude;
		}


		/// Saturate the given signed value in place to the given magnitude.
		inline void Saturate_Value(int &input, unsigned int max_magnitude)
		{
			if (input > static_cast<int> (max_magnitude))
				input = static_cast<int> (max_magnitude);
			else if (input < static_cast<int> (- max_magnitude))
				input = static_cast<int> (- max_magnitude);
		}


		/// Saturate the given signed value in place to the given maximum and minimum values.
		inline void Saturate_Value(int &input, int max_value, int min_value)
		{
			if (input > max_value)
				input = max_value;
			else if (input < min_value)
				input = min_value;
		}

		/// Calculate 2^exponent.
		inline unsigned int Pow_2(unsigned int exponent)
		{
			return 1 << exponent;
		}

	}


	namespace bit_operations
	{
		/// Convert an unsigned decimal number to the corresponding value in Gray representation
		inline unsigned int Unsigned_2_Gray(unsigned int unsigned_value)
		{
			unsigned int unsigned_value_shifted = unsigned_value << 1;

			return (unsigned_value ^ unsigned_value_shifted) >> 1;
		}


		/// Convert a value in Gray representation to an unsigned decimal number.
		inline unsigned int Gray_2_Unsigned(unsigned int gray_value)
		{
			unsigned int result = 0;
			unsigned int last_bit = 0;

			for (int i = 31; i >= 0; i--)
			{
				unsigned int gray_msb = (gray_value >> i) & 1;
				unsigned int new_bit = gray_msb ^ last_bit;

				result = (result << 1) | new_bit;
				last_bit = new_bit;
			}

			return result;
		}
	}


	/// Helper function and type for template based meta programming.
	/**
	 * The namespace hold functions and data types for C++ meta programming.
	 * Meta programming allows for a compile time resolution of expressions.
	 * This can be used to generate code or to speed up the run-time.
	 */
	namespace meta_prog
	{
		// usage if(IsSame<T1,T>::Result::value)

		/// Type "False" for meta programming
		struct FalseType { enum { value = false }; };

		/// Type "True" for meta programming
		struct TrueType { enum { value = true }; };

		/// Template for a equal comparison of types in meta programming (compile time resolution), false result.
		template <typename T1, typename T2>
		struct IsSame
		{
			/// Type for the result false, if the two types are different.
			typedef FalseType Result;
		};

		/// Template for a equal comparison of types in meta programming (compile time resolution), true result.
		template <typename T>
		struct IsSame<T,T>
		{
			/// Type for the result true, if the two types are the same.
			typedef TrueType Result;
		};
	}
}
#endif /* HELPER_FUNCTIONS_H_ */
