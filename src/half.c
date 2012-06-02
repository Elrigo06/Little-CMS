/******************************************************************************
 *
 * Filename:    ieeehalfprecision.c
 * Programmer:  James Tursa
 * Version:     1.0
 * Date:        March 3, 2009
 * Copyright:   (c) 2009 by James Tursa, All Rights Reserved
 *
 * Edited by:   James Cloos
 * Edited on:   June 2, 2012
 * Copyright:   (c) 2012 by James Cloos, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * This file contains C code to convert between IEEE double, single, and half
 * precision floating point formats. The intended use is for standalone C code
 * that does not rely on MATLAB mex.h. The bit pattern for the half precision
 * floating point format is stored in a 16-bit unsigned int variable. The half
 * precision bit pattern definition is:
 *
 * 1 bit sign bit
 * 5 bits exponent, biased by 15
 * 10 bits mantissa, hidden leading bit, normalized to 1.0
 *
 * Special floating point bit patterns recognized and supported:
 *
 * All exponent bits zero:
 * - If all mantissa bits are zero, then number is zero (possibly signed)
 * - Otherwise, number is a denormalized bit pattern
 *
 * All exponent bits set to 1:
 * - If all mantissa bits are zero, then number is +Infinity or -Infinity
 * - Otherwise, number is NaN (Not a Number)
 *
 * For the denormalized cases, note that 2^(-24) is the smallest number that can
 * be represented in half precision exactly. 2^(-25) will convert to 2^(-24)
 * because of the rounding algorithm used, and 2^(-26) is too small and underflows
 * to zero.
 *
 ********************************************************************************/

// Includes -------------------------------------------------------------------

#include <string.h>
#include "lcms2_internal.h"

//-----------------------------------------------------------------------------
//
// Routine:  singles2halfp
//
// Input:  source = Address of 32-bit floating point data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 16-bit data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

cmsBool singles2halfp(void *target, void *source, cmsUInt32Number numel)
{
    cmsUInt16Number *hp = (cmsUInt16Number *) target; // Type pun output as an unsigned 16-bit int
    cmsUInt32Number *xp = (cmsUInt32Number *) source; // Type pun input as an unsigned 32-bit int
    cmsUInt16Number    hs, he, hm;
    cmsUInt32Number x, xs, xe, xm;
    cmsUInt32Number hes;
    static cmsUInt32Number next;  // Little Endian adjustment
    static cmsBool checkieee = TRUE;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    cmsUInt32Number *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (cmsUInt32Number *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return FALSE;  // Floating point bit pattern is not IEEE 754
        }
        checkieee = FLASE; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return TRUE;
    }
    
    while( numel-- ) {
        x = *xp++;
        if( (x & 0x7FFFFFFFu) == 0 ) {  // Signed zero
            *hp++ = (cmsUInt16Number) (x >> 16);  // Return the signed zero
        } else { // Not zero
            xs = x & 0x80000000u;  // Pick off sign bit
            xe = x & 0x7F800000u;  // Pick off exponent bits
            xm = x & 0x007FFFFFu;  // Pick off mantissa bits
            if( xe == 0 ) {  // Denormal will underflow, return a signed zero
                *hp++ = (cmsUInt16Number) (xs >> 16);
            } else if( xe == 0x7F800000u ) {  // Inf or NaN (all the exponent bits are set)
                if( xm == 0 ) { // If mantissa is zero ...
                    *hp++ = (cmsUInt16Number) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else {
                    *hp++ = (cmsUInt16Number) 0xFE00u; // NaN, only 1st mantissa bit set
                }
            } else { // Normalized number
                hs = (cmsUInt16Number) (xs >> 16); // Sign bit
                hes = ((int)(xe >> 23)) - 127 + 15; // Exponent unbias the single, then bias the halfp
                if( hes >= 0x1F ) {  // Overflow
                    *hp++ = (cmsUInt16Number) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else if( hes <= 0 ) {  // Underflow
                    if( (14 - hes) > 24 ) {  // Mantissa shifted all the way off & no rounding possibility
                        hm = (cmsUInt16Number) 0u;  // Set mantissa to zero
                    } else {
                        xm |= 0x00800000u;  // Add the hidden leading bit
                        hm = (cmsUInt16Number) (xm >> (14 - hes)); // Mantissa
                        if( (xm >> (13 - hes)) & 0x00000001u ) // Check for rounding
                            hm += (cmsUInt16Number) 1u; // Round, might overflow into exp bit, but this is OK
                    }
                    *hp++ = (hs | hm); // Combine sign bit and mantissa bits, biased exponent is zero
                } else {
                    he = (cmsUInt16Number) (hes << 10); // Exponent
                    hm = (cmsUInt16Number) (xm >> 13); // Mantissa
                    if( xm & 0x00001000u ) // Check for rounding
                        *hp++ = (hs | he | hm) + (cmsUInt16Number) 1u; // Round, might overflow to inf, this is OK
                    else
                        *hp++ = (hs | he | hm);  // No rounding
                }
            }
        }
    }
    return TRUE;
}

//-----------------------------------------------------------------------------
//
// Routine:  doubles2halfp
//
// Input:  source = Address of 64-bit floating point data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 16-bit data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

cmsBool doubles2halfp(void *target, void *source, cmsUInt32Number numel)
{
    cmsUInt16Number *hp = (cmsUInt16Number *) target; // Type pun output as an unsigned 16-bit int
    cmsUInt32Number *xp = (cmsUInt32Number *) source; // Type pun input as an unsigned 32-bit int
    cmsUInt16Number    hs, he, hm;
    cmsUInt32Number x, xs, xe, xm;
    cmsUInt32Number hes;
    static cmsUInt32Number next;  // Little Endian adjustment
    static cmsBool checkieee = TRUE;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    cmsUInt32Number *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (cmsUInt32Number *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return FALSE;  // Floating point bit pattern is not IEEE 754
        }
        checkieee = FALSE; // Everything checks out OK
    }

    xp += next;  // Little Endian adjustment if necessary
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return TRUE;
    }
    
    while( numel-- ) {
        x = *xp++; xp++; // The extra xp++ is to skip over the remaining 32 bits of the mantissa
        if( (x & 0x7FFFFFFFu) == 0 ) {  // Signed zero
            *hp++ = (cmsUInt16Number) (x >> 16);  // Return the signed zero
        } else { // Not zero
            xs = x & 0x80000000u;  // Pick off sign bit
            xe = x & 0x7FF00000u;  // Pick off exponent bits
            xm = x & 0x000FFFFFu;  // Pick off mantissa bits
            if( xe == 0 ) {  // Denormal will underflow, return a signed zero
                *hp++ = (cmsUInt16Number) (xs >> 16);
            } else if( xe == 0x7FF00000u ) {  // Inf or NaN (all the exponent bits are set)
                if( xm == 0 ) { // If mantissa is zero ...
                    *hp++ = (cmsUInt16Number) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else {
                    *hp++ = (cmsUInt16Number) 0xFE00u; // NaN, only 1st mantissa bit set
                }
            } else { // Normalized number
                hs = (cmsUInt16Number) (xs >> 16); // Sign bit
                hes = ((int)(xe >> 20)) - 1023 + 15; // Exponent unbias the double, then bias the halfp
                if( hes >= 0x1F ) {  // Overflow
                    *hp++ = (cmsUInt16Number) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else if( hes <= 0 ) {  // Underflow
                    if( (10 - hes) > 21 ) {  // Mantissa shifted all the way off & no rounding possibility
                        hm = (cmsUInt16Number) 0u;  // Set mantissa to zero
                    } else {
                        xm |= 0x00100000u;  // Add the hidden leading bit
                        hm = (cmsUInt16Number) (xm >> (11 - hes)); // Mantissa
                        if( (xm >> (10 - hes)) & 0x00000001u ) // Check for rounding
                            hm += (cmsUInt16Number) 1u; // Round, might overflow into exp bit, but this is OK
                    }
                    *hp++ = (hs | hm); // Combine sign bit and mantissa bits, biased exponent is zero
                } else {
                    he = (cmsUInt16Number) (hes << 10); // Exponent
                    hm = (cmsUInt16Number) (xm >> 10); // Mantissa
                    if( xm & 0x00000200u ) // Check for rounding
                        *hp++ = (hs | he | hm) + (cmsUInt16Number) 1u; // Round, might overflow to inf, this is OK
                    else
                        *hp++ = (hs | he | hm);  // No rounding
                }
            }
        }
    }
    return TRUE;
}

//-----------------------------------------------------------------------------
//
// Routine:  halfp2singles
//
// Input:  source = address of 16-bit data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 32-bit floating point data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

cmsBool halfp2singles(void *target, void *source, cmsUInt32Number numel)
{
    cmsUInt16Number *hp = (cmsUInt16Number *) source; // Type pun input as an unsigned 16-bit int
    cmsUInt32Number *xp = (cmsUInt32Number *) target; // Type pun output as an unsigned 32-bit int
    cmsUInt16Number h, hs, he, hm;
    cmsUInt32Number xs, xe, xm;
    cmsInt32Number xes;
    cmsUInt32Number e;
    static cmsUInt32Number next;  // Little Endian adjustment
    static cmsBool checkieee = TRUE;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    cmsUInt32Number *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (cmsUInt32Number *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return FALSE;  // Floating point bit pattern is not IEEE 754
        }
        checkieee = FALSE; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return TRUE;
    
    while( numel-- ) {
        h = *hp++;
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            *xp++ = ((cmsUInt32Number) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((cmsUInt32Number) hs) << 16; // Sign bit
                xes = ((cmsInt32Number) (he >> 10)) - 15 + 127 - e; // Exponent unbias the halfp, then bias the single
                xe = (cmsUInt32Number) (xes << 23); // Exponent
                xm = ((cmsUInt32Number) (hm & 0x03FFu)) << 13; // Mantissa
                *xp++ = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                if( hm == 0 ) { // If mantissa is zero ...
                    *xp++ = (((cmsUInt32Number) hs) << 16) | ((cmsUInt32Number) 0x7F800000u); // Signed Inf
                } else {
                    *xp++ = (cmsUInt32Number) 0xFFC00000u; // NaN, only 1st mantissa bit set
                }
            } else { // Normalized number
                xs = ((cmsUInt32Number) hs) << 16; // Sign bit
                xes = ((cmsInt32Number) (he >> 10)) - 15 + 127; // Exponent unbias the halfp, then bias the single
                xe = (cmsUInt32Number) (xes << 23); // Exponent
                xm = ((cmsUInt32Number) hm) << 13; // Mantissa
                *xp++ = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
    }
    return TRUE;
}

//-----------------------------------------------------------------------------
//
// Routine:  halfp2singles
//
// Input:  source = address of 16-bit data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 32-bit floating point data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

cmsBool halfp2doubles(void *target, void *source, cmsUInt32Number numel)
{
    cmsUInt16Number *hp = (cmsUInt16Number *) source; // Type pun input as an unsigned 16-bit int
    cmsUInt32Number *xp = (cmsUInt32Number *) target; // Type pun output as an unsigned 32-bit int
    cmsUInt16Number h, hs, he, hm;
    cmsUInt32Number xs, xe, xm;
    cmsInt32Number xes;
    cmsUInt32Number e;
    static cmsUInt32Number next;  // Little Endian adjustment
    static cmsBool checkieee = TRUE;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    cmsUInt32Number *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (cmsUInt32Number *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return FALSE;  // Floating point bit pattern is not IEEE 754
        }
        checkieee = FALSE; // Everything checks out OK
    }

    xp += next;  // Little Endian adjustment if necessary
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return TRUE;
    
    while( numel-- ) {
        h = *hp++;
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            *xp++ = ((cmsUInt32Number) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((cmsUInt32Number) hs) << 16; // Sign bit
                xes = ((cmsInt32Number) (he >> 10)) - 15 + 1023 - e; // Exponent unbias the halfp, then bias the double
                xe = (cmsUInt32Number) (xes << 20); // Exponent
                xm = ((cmsUInt32Number) (hm & 0x03FFu)) << 10; // Mantissa
                *xp++ = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                if( hm == 0 ) { // If mantissa is zero ...
                    *xp++ = (((cmsUInt32Number) hs) << 16) | ((cmsUInt32Number) 0x7FF00000u); // Signed Inf
                } else {
                    *xp++ = (cmsUInt32Number) 0xFFF80000u; // NaN, only the 1st mantissa bit set
                }
            } else {
                xs = ((cmsUInt32Number) hs) << 16; // Sign bit
                xes = ((cmsInt32Number) (he >> 10)) - 15 + 1023; // Exponent unbias the halfp, then bias the double
                xe = (cmsUInt32Number) (xes << 20); // Exponent
                xm = ((cmsUInt32Number) hm) << 10; // Mantissa
                *xp++ = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
        xp++; // Skip over the remaining 32 bits of the mantissa
    }
    return TRUE;
}
