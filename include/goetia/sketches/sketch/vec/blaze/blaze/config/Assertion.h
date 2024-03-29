//=================================================================================================
/*!
//  \file blaze/config/Assertion.h
//  \brief Configuration of the run time assertion macros
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
/*!\brief Compilation switch for internal assertions.
// \ingroup config
//
// This compilation switch triggers internal assertions, which are used to verify the program
// itself. The internal assertions can also be deactivated by defining \a NDEBUG during the
// compilation.
//
// Possible settings for the internal assertion switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1
//
// \note It is possible to (de-)activate internal assertions via command line by defining this
// symbol manually before including any Blaze header file:

   \code
   #define BLAZE_INTERNAL_ASSERTION 1
   #include <blaze/Blaze.h>
   \endcode
*/
#ifndef BLAZE_INTERNAL_ASSERTION
#define BLAZE_INTERNAL_ASSERTION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for user assertions.
// \ingroup config
//
// This compilation switch triggers user assertions, which are used to check user specified
// function parameters and values. The user assertions can also be deactivated by defining
// \a NDEBUG during the compilation.
//
// Possible settings for the user assertion switch:
//  - Deactivated: \b 0
//  - Activated  : \b 1
//
// \note It is possible to (de-)activate user assertions via command line by defining this
// symbol manually before including any Blaze header file:

   \code
   #define BLAZE_USER_ASSERTION 1
   #include <blaze/Blaze.h>
   \endcode
*/
#ifndef BLAZE_USER_ASSERTION
#define BLAZE_USER_ASSERTION 0
#endif
//*************************************************************************************************
