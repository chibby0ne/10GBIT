//
//  Copyright (C) 2010 - 2011 Creonic GmbH
//
//  This file is part of the Creonic simulation environment (CSE)
//  for communication systems.
//
/// \file
/// \brief  Wrapper header file for buffer class class and its specializations.
/// \author Timo Lehnigk-Emden
/// \date   2010/11/02
//
#ifndef BUFFER_H_
#define BUFFER_H_

//for template safety check if template parameter are already defined!
#if (defined(T) || defined(DIM))
#error One or more of the following is defined: T, DIM. Definition conflicts with their usage
       as template parameters.
#endif

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include "helper_functions.h"
// Include the implementation of the functions of the template classes to prevent linker problems

//
// Forward declaration of generic Buffer class. The generic class uses the specialized classes,
// therefore the definition of the specialized classes
// are done before the definition of the generic class.
template <class T, int DIM = 1> class Buffer;
template <class T> class Buffer<T,1>;  // Forward declaration of the specialized versions.
template <class T> class Buffer<T,2>;
#include "buffer_dim1.hpp"
#include "buffer_dim2.hpp"
#include "buffer_dimx.hpp"

#endif

