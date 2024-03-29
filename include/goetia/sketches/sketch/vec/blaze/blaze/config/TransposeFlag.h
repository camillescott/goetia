//=================================================================================================
/*!
//  \file blaze/config/TransposeFlag.h
//  \brief Configuration of the default transpose flag for all vectors of the Blaze library
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================


//*************************************************************************************************
/*!\brief The default transpose flag for all vectors of the Blaze library.
// \ingroup config
//
// This value specifies the default transpose flag for all vectors of the Blaze library.
// In case no explicit transpose flag is specified with the according vector type, this
// setting is used.

   \code
   // Explicit specification of the transpose flag => column vector
   StaticVector<double,3UL,columnVector> a;

   // No explicit specification of the transpose flag => use of the default transpose flag
   StaticVector<double,3UL> b;
   \endcode

// Valid settings for the BLAZE_DEFAULT_TRANSPOSE_FLAG are blaze::rowVector and blaze::columnVector.
//
// \note It is possible to specify the default transpose flag via command line by defining this
// symbol manually before including any Blaze header file:

   \code
   #define BLAZE_DEFAULT_TRANSPOSE_FLAG blaze::columnVector
   #include <blaze/Blaze.h>
   \endcode
*/
#ifndef BLAZE_DEFAULT_TRANSPOSE_FLAG
#define BLAZE_DEFAULT_TRANSPOSE_FLAG blaze::columnVector
#endif
//*************************************************************************************************
