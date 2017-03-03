/*
	File:			VLMath.h

	Function:		Various math definitions for VL
					
	Author(s):		Andrew Willmott

	Copyright:		Copyright (c) 1995-1996, Andrew Willmott
 */

#ifndef __VL_MATH__
#define __VL_MATH__

// --- Inlines -------------------------------------------------------

#ifdef VL_HAS_ABSF
inline Float abs(Float x)
{ return (fabsf(x)); }
inline Float len(Float x)
{ return (fabsf(x)); }
#endif

inline Double abs(Double x)
{ return (fabs(x)); }
inline Double len(Double x)
{ return (fabs(x)); }

inline Float sqrlen(Float r)
{ return(sqr(r)); }
inline Double sqrlen(Double r)
{ return(sqr(r)); }

inline Float mix(Float a, Float b, Float s)
{ return((1.0 - s) * a + s * b); }
inline Float mix(Double a, Double b, Double s)
{ return((1.0 - s) * a + s * b); }

inline Void SetReal(Float &a, Double b)
{ a = b; }
inline Void SetReal(Double &a, Double b)
{ a = b; }

#endif
