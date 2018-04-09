//--------------------------------------------------------------------------
/**
 * @file
 * @brief   miscelleneous things
 * @ingroup General_Util
 * @author  N.Ishii
 * @since   Sat May 22 20:31:03 JST 2010
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>


#include <yukawa/yukawa.h>
#include <yukawa/misc.h>

Start_Name_Space_Yukawa

//---------------------------------------------------------------------
/**
 * @brief   changes endian for an array of @a uint32_t of length @a len
 */
//---------------------------------------------------------------------

void change_endian(uint32_t buf[], int len)
{
  for(int i=0; i<len; i++){
    char c[sizeof(uint32_t)];
    for(int j = 0; j < (int)sizeof(uint32_t); j++)
      c[j] = ((char*)buf)[j];

    for(int j = 0; j < (int)sizeof(uint32_t); j++)
      ((char*)buf)[j] = c[sizeof(uint32_t) - j - 1];

    buf += 1;
  }
}

//---------------------------------------------------------------------
/**
 * @brief   changes endian for an array of @a int32_t of length @a len
 */
//---------------------------------------------------------------------

void change_endian(int32_t buf[], int len)
{
  for(int i=0; i<len; i++){
    char c[sizeof(int32_t)];
    for(int j = 0; j < (int)sizeof(int32_t); j++)
      c[j] = ((char*)buf)[j];

    for(int j = 0; j < (int)sizeof(int32_t); j++)
      ((char*)buf)[j] = c[sizeof(int32_t) - j - 1];

    buf += 1;
  }
}

//---------------------------------------------------------------------
/**
 * @brief   changes endian for an array of @a double of length @a len
 */
//---------------------------------------------------------------------

void change_endian(double buf[], int len)
{
  for(int i=0; i<len; i++){
    char c[8];
    for(int j=0; j<8; j++) c[j] = ((char*)buf)[j];
    for(int j=0; j<8; j++) ((char*)buf)[j] = c[7-j];
    buf += 1;
  }
}

//---------------------------------------------------------------------
/**
 * @brief   returns true if the system adopt bigendian
 */
//---------------------------------------------------------------------

bool is_bigendian()
{
  uint32_t num;
  num = 1 + 2*256 + 3*256*256 + 4*256*256*256;

  // big endian
  if((int)((char*)&num)[0] == 4 &&
     (int)((char*)&num)[1] == 3 &&
     (int)((char*)&num)[2] == 2 &&
     (int)((char*)&num)[3] == 1) return true;

  // little endian
  if((int)((char*)&num)[0] == 1 &&
     (int)((char*)&num)[1] == 2 &&
     (int)((char*)&num)[2] == 3 &&
     (int)((char*)&num)[3] == 4) return false;

  fprintf(stderr,
          "ERROR(is_bigendian()): unknown endian: %d,%d,%d,%d\n",
          ((int)((char*)&num)[0]),
          ((int)((char*)&num)[1]),
          ((int)((char*)&num)[2]),
          ((int)((char*)&num)[3]));
  exit(1);
  return false;
}

//---------------------------------------------------------------------
/**
 * @brief   returns the pointer to the character string
 *          containing the information of the local time
 *
 */
//---------------------------------------------------------------------

const char* LocalTime()
{
  struct timeval tp;
  time_t         ptm;
  static char    str[128];

  gettimeofday(&tp,NULL);
  ptm = tp.tv_sec;
  strcpy(str,asctime(localtime(&ptm)));
  str[strlen(str)-1]='\0';

  return str;
}

End_Name_Space_Yukawa
