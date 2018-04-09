//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup PH1_group
 * @brief   <b>LOCAL</b> I/O for four_point<PH1>
 * @author  N.ishii
 * @since   Mon May 24 16:29:27 JST 2010
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_PH1_4PT_H
#define IS_INCLUdED_PH1_4PT_H

#include <yukawa/yukawa.h>
#include <yukawa/PH1.h>
#include <yukawa/4pt.h>

Start_Name_Space_Yukawa

extern void Read(four_point<PH1>& four,
		 FILE* fp,
		 bool is_change_endian=false);

extern void Read(four_point<PH1>& four,
		 const char* fname,
		 bool is_change_endian=false);

extern void Write(const four_point<PH1>& four,
		  FILE* fp,
		  bool is_change_endian=false);

extern void Write(const four_point<PH1>& four,
		  const char* fname,
		  bool is_change_endian=false);

End_Name_Space_Yukawa

#endif
