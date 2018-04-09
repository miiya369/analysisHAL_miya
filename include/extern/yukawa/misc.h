//--------------------------------------------------------------------------
/**
 * @file
 * @brief   miscelleanous things
 * @ingroup General_Util
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Sat May 22 20:25:57 JST 2010
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_MISC_H_YUKAWA
#define IS_INCLUDED_MISC_H_YUKAWA

#include <stdint.h>

#include <yukawa/yukawa.h>

Start_Name_Space_Yukawa

extern void change_endian(uint32_t buf[], int len);
extern void change_endian( int32_t buf[], int len);
extern void change_endian(double   buf[], int len);
extern bool is_bigendian();

extern const char* LocalTime();

End_Name_Space_Yukawa

#endif
