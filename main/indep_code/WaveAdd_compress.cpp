//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup none
 * @brief   The code for averaging of 1/48 compressed NBSwave data
 * @author  Takaya Miyamoto
 * @since   Tue Jun 21 21:39:44 JST 2016
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

static string *ipath = NULL;
static double *fac   = NULL;
static string  opath;
static int     Nadd;

static int32_t tmp0   = 0;
static int32_t tmp16  = 16;
static int32_t tmp20  = 20;
static int32_t tmp256 = 256;

static int32_t *RepLst = NULL;
static int32_t  Lsize[3];
static int32_t  length;
static int32_t  WaveSize;
static int32_t  dummy;

static long skip_len;

static bool endian_flip;

inline int set_args(int, char**);

inline int read_header();

inline int input (double*, int);
inline int output(double*);

inline bool machine_is_little_endian();
inline void endian_convert(double*, size_t);
inline void endian_convert(int32_t*, size_t);

//========================================================================//
int main(int argc, char **argv)
{
   if (set_args(argc, argv) == 1) return 1;
   
   if (read_header() == 1) return 1;
   
   skip_len = (16 + length) * sizeof(int32_t);
   WaveSize = 32 * length;
   
   double *Wave = new double[WaveSize];
   for (int i=0; i<WaveSize; i++) Wave[i] = 0.0;
   
   for (int i=0; i<Nadd; i++) if (input(Wave, i) == 1) return 1;
   
   if (output(Wave) == 1) return 1;
   
   delete [] Wave;
   if (ipath  != NULL) delete [] ipath;
   if (RepLst != NULL) delete [] RepLst;
   return 0;
}
//========================================================================//

inline int set_args(int argc, char** argv)
{
   if (argc < 3) {
      printf("usage: %s [ofile] [factor1] [ifile1] [factor2] [ifile2] ...\n", argv[0]);
      return 1;
   }
   
   if ((argc - 2) % 2 != 0) {
      printf("Arguments error: Should be put each '[factor] [ifile]'\n");
      return 1;
   }
   
   Nadd = (argc - 2) / 2;
   ipath = new string [Nadd];
   fac   = new double [Nadd];
   
   for (int i=0; i<Nadd; i++) {
      fac  [i] = atof(argv[2 + 0 + 2*i]);
      ipath[i] =      argv[2 + 1 + 2*i];
      
//      printf("%lf\n", fac[i]); // For Debug
   }
   opath = argv[1];
   
   endian_flip = !machine_is_little_endian();
   
//   printf("Nave = %d\n", Naverage); // For Debug
//   return 1;
   
   return 0;
}

inline int read_header()
{
   FILE* fp;
   if ((fp = fopen(ipath[0].c_str(), "rb")) == NULL) {
      printf("The file '%s' is not found.\n", ipath[0].c_str());
      return 1;
   }
   fread((void*)&dummy, sizeof(int32_t), 1, fp);
   if (endian_flip) endian_convert(&dummy, 1);
   if (dummy != 16) {
      printf("The file '%s' may be not Ishii-san's compressed NBS data.", ipath[0].c_str());
      return 1;
   }
   fread((void*) Lsize,  sizeof(int32_t), 3, fp);
   if (endian_flip) endian_convert(Lsize, 3);
   fread((void*)&length, sizeof(int32_t), 1, fp);
   if (endian_flip) endian_convert(&length, 1);
   
   RepLst = new int32_t[length];
   
   fread((void*)&dummy, sizeof(int32_t), 1, fp);
   fread((void*)&dummy, sizeof(int32_t), 1, fp);
   if (endian_flip) endian_convert(&dummy, 1);
   if (dummy != length * sizeof(int32_t)) {
     printf("The file '%s' may be not Ishii-san's compressed NBS data.", ipath[0].c_str());
     return 1;
   }
   fread((void*)RepLst, sizeof(int32_t), length, fp);
   if (endian_flip) endian_convert(RepLst, length);
   
   fread((void*)&dummy, sizeof(int32_t), 1, fp);
   if (endian_flip) endian_convert(&dummy, 1);
   if (dummy != length * sizeof(int32_t)) {
     printf("The file '%s' may be not Ishii-san's compressed NBS data.", ipath[0].c_str());
     return 1;
   }
   fclose(fp);

   return 0;
}

inline int input(double* wave, int data_num)
{
   double *tmp_wave = new double[WaveSize];
   
   FILE* fp;
   if ((fp = fopen(ipath[data_num].c_str(), "rb")) == NULL) {
      printf("The file '%s' is not found.\n", ipath[data_num].c_str());
      return 1;
   }
   fseek(fp, skip_len, SEEK_SET);
   
   fread((void*)tmp_wave, sizeof(double), WaveSize, fp);
   if (endian_flip) endian_convert(tmp_wave, WaveSize);

   fread((void*)&dummy, sizeof(int32_t), 1, fp);
   if (endian_flip) endian_convert(&dummy, 1);
   if (dummy != length * 256) {
     printf("The file '%s' may be not Ishii-san's compressed NBS data.", ipath[0].c_str());
     return 1;
   }
   fclose(fp);
   
   for (int i=0; i<WaveSize; i++) wave[i] += (fac[data_num] * tmp_wave[i]);
   
   delete [] tmp_wave;
   
   return 0;
}

inline int output(double* wave)
{
   FILE* fp;
   if ((fp = fopen(opath.c_str(), "wb")) == NULL) {
      printf("The file '%s' can not be open.\n", opath.c_str());
      return 1;
   }
   
   int32_t tmp_length0 = length;
   int32_t tmp_length1 = length * sizeof(int32_t);
   int32_t tmp_length2 = length * 256;
   
   if (endian_flip) {
     endian_convert(&tmp16,  1);
     endian_convert( Lsize,  3);
     endian_convert(&length, 1);
     
     endian_convert(&tmp_length1, 1);
     endian_convert( RepLst,      tmp_length0);
     
     endian_convert(&tmp20,  1);
     endian_convert(&tmp0,   1);
     endian_convert(&tmp256, 1);
     
     endian_convert(&tmp_length2, 1);
     endian_convert( wave,        WaveSize);
   }
   
   fwrite((void*)&tmp16,  sizeof(int32_t), 1, fp);
   fwrite((void*) Lsize,  sizeof(int32_t), 3, fp);
   fwrite((void*)&length, sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp16,  sizeof(int32_t), 1, fp);
   
   fwrite((void*)&tmp_length1, sizeof(int32_t), 1,           fp);
   fwrite((void*) RepLst,      sizeof(int32_t), tmp_length0, fp);
   fwrite((void*)&tmp_length1, sizeof(int32_t), 1,           fp);

   fwrite((void*)&tmp20,  sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp0,   sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp0,   sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp0,   sizeof(int32_t), 1, fp);
   fwrite((void*)&length, sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp256, sizeof(int32_t), 1, fp);
   fwrite((void*)&tmp20,  sizeof(int32_t), 1, fp);

   fwrite((void*)&tmp_length2, sizeof(int32_t), 1,        fp);
   fwrite((void*) wave,        sizeof(double),  WaveSize, fp);
   fwrite((void*)&tmp_length2, sizeof(int32_t), 1,        fp);
   
   fclose(fp);
   
   return 0;
}

bool machine_is_little_endian() 
{
  int endianTEST = 1;
  if (*(char*)&endianTEST) return true;
  else                     return false;
}

void endian_convert(double *DATA, size_t DATA_size)
{
  char dummy[8];
   
  for (size_t k=0; k<DATA_size; k++) {
    for (int j=0; j<8; j++) dummy[j] = ((char*)&DATA[k])[j];
    for (int j=0; j<8; j++) ((char*)&DATA[k])[j] = dummy[7-j];
  }
}

void endian_convert(int* DATA, size_t DATA_size)
{
  char dummy[4];
   
  for (size_t k=0; k<DATA_size; k++ ) {
    for(int j=0; j<4; j++) dummy[j] = ((char*)&DATA[k])[j];
    for(int j=0; j<4; j++) ((char*)&DATA[k])[j] = dummy[3-j];
  }
}
