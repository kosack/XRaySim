/*
	File:			VLConfig.h

	Function:		Contains configuration options for compiling the VL 
					library.
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
*/

//
//	Current options are as follows:
//
//	VL_NO_BOOL		- Compiler has no bool type.
//	VL_NO_TF		- true and false are not predefined.
//	VL_HAS_ABSF		- has the absf() call. 
//	VL_NO_INST_TMPL - compiler doesn't accept ansi syntax for template 
//                    instantiation.
//	VL_SGI_INST		- use sgi's weirdo instantiation syntax.
//	VL_ROW_ORIENT	- Use row-oriented transforms, so you can swap 'em with 
//					  OpenGL. If this is defined, transformations are 
//                    v = v * Rot3x(u, 0.5), rather than v = Rot3x(u, 0.5) * v.
//                    This is off by default.
//

// --- SGI -------------------------------

#ifdef __SGI__
#define VL_NO_BOOL
#define VL_NO_TF
#define VL_HAS_ABSF
#define VL_SGI_INST
#endif

// --- GCC -------------------------------

#ifdef __GCC__
#endif

// --- Metrowerks, MacOS -----------------

#ifdef __MACOS__
#endif

