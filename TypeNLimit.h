/*

   TypeNLimit.h		Miscellaneous Constants

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef __TYPENLIMIT_H__
#define __TYPENLIMIT_H__

#include <limits.h>

#define BITS_IN_WORD 32
#define BITS_IN_WORD_MINUS_1 31
#define BITS_IN_WORD_MASK 0x0000001F
#define BITS_IN_WORD_SHIFT 5
#define BITS_IN_HALF_WORD 16
#define BITS_IN_4_WORD 128
#define BITS_IN_4_WORD_MINUS_1 127
#define BITS_IN_4_WORD_SHIFT 7
#define FIRST_BIT_MASK 0x80000000
#define ALL_BUT_FIRST_BIT_MASK 0x7FFFFFFF
#define ALL_ONE_MASK 0xFFFFFFFF
#define FOUR_MULTIPLE_MASK 0xFFFFFFFC
#define BITS_IN_BYTE 8
#define BITS_IN_BYTE_SHIFT 3
#define BYTES_IN_WORD 4

#define TRUE    1
#define FALSE   0

// Compatibilities
//#define fopen64		fopen			// For compatibility of environment without large file support
//#define ftello64	ftell			// For compatibility of environment without large file support
#define INLINE		inline		// For compatibility

#define MAX_FILENAME_LEN 256

#endif
