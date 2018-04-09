//--------------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup General_Util
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Sat May 22 19:28:13 JST 2010
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_YUKAWA_H
#define IS_INCLUDED_YUKAWA_H


#define Start_Name_Space_Yukawa namespace yukawa {

#define End_Name_Space_Yukawa }


//---------------------------------------
// added by T. Doi

#define YUKAWA_NOARCH 0
#define YUKAWA_QCDOC  1
#define YUKAWA_QCDSP  2
#define YUKAWA_x86    3
#define YUKAWA_BGL    4
#define YUKAWA_BGQ    5
#define YUKAWA_SR16K  6
#define YUKAWA_KEI    7

#define YUKAWA_TARGET @ARCH@

// The configure procedure should make this unnecessary, but just in case...
#ifndef YUKAWA_TARGET
#define YUKAWA_TARGET YUKAWA_NOARCH
#endif
//---------------------------------------


#endif
