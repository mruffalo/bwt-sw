/*

   BWTFormatdb.h		Build index for FASTA database

   This program builds index for FASTA database for use of BWTBlastn.

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

#ifndef __BWT_FORMAT_DB_H__
#define __BWT_FORMAT_DB_H__

#include "TypeNLimit.h"

// Database and ini
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void ProcessIni();
void ValidateIni();
void PrintIni();

void ProcessFileName(char *outputFileName, const char *inputFileName, const char *databaseName);


#endif
