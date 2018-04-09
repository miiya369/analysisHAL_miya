//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup none
 * @brief   Make Jack-Knife bin samples for 1/48 Not-compressed NBSwave data
 * @author  Takaya Miyamoto
 * @since   Wed Sep 28 03:45:59 JST 2016
 */
//--------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#define MAX_N_DATA   1024
#define MAX_LEN_PATH 1024

#define idx(n,i) ((n) + DataSize * (i))

static string ilist;
static string olist;
static int    Nbin, binSize;
static string ifpath;
static string ofpath;
static string repstr;

static char iconf_list[MAX_N_DATA][MAX_LEN_PATH];
static char oconf_list[MAX_N_DATA][MAX_LEN_PATH];

static int Niconf, Noconf;

static size_t DataSize;

static bool endian_flip;

inline int set_args(const int, const char**);

inline int input (double*, const int);
inline int output(const double*, const int);

inline bool machine_is_little_endian();
inline void endian_convert(double*, const size_t);

inline int    set_conf_list(const char*, char conf_list[MAX_N_DATA][MAX_LEN_PATH]);
inline void   MakeJKsample(double*, double*);
inline string StrReplace(string, const string, const string);

//========================================================================//
int main(const int argc, const char **argv)
{
   if (set_args(argc, argv) == 1) return 1;
   
   double *iData = new double[DataSize * Niconf];
   double *oData = new double[DataSize * Noconf];
   
   for (int i=0; i<Niconf; i++) if (input(iData, i) == 1) return 1;

   MakeJKsample(iData, oData);
   
   for (int i=0; i<Noconf; i++) if (output(oData, i) == 1) return 1;
   
   delete [] iData;
   delete [] oData;
   
   return 0;
}
//========================================================================//

inline int set_args(const int argc, const char** argv)
{
   if (argc < 6) {
      printf("usage: %s [Org conf list] [JK.bin conf list] [Bin size] [input path] [output path]\n"
	     "     : \"REPCONF\" at the [File path] to be replaced to the each line of conf list.\n"
	     , argv[0]);
      return 1;
   }
   
   ilist   = argv[1];
   olist   = argv[2];
   binSize = atoi(argv[3]);
   ifpath   = argv[4];
   ofpath   = argv[5];
   
   repstr = "REPCONF";
   
   Niconf = set_conf_list(ilist.c_str(), iconf_list);
   Noconf = set_conf_list(olist.c_str(), oconf_list);
   if (Niconf == 0 || Noconf == 0) return 1;
   
   Nbin = Niconf / binSize;
   
   if (Noconf != Nbin) {
     printf("ERROR: #.JKbin conf != #.bin\n");
     return 1;
   }

   string ifname = StrReplace(ifpath, repstr, iconf_list[0]);
   FILE* fp;   
   if ((fp = fopen(ifname.c_str(), "rb")) == NULL) {
      printf("The file '%s' is not found.\n", ifname.c_str());
      return 1;
   }
   
   fseek (fp, 0, SEEK_END);
   DataSize = ftell(fp);
   fseek (fp, 0, SEEK_SET);
   fclose(fp);
   
   if (DataSize % 8 != 0) {
     printf("Invalid file size: fsize %% sizeof(double) != 0\n");
   }
   DataSize /= 8;
   
   endian_flip = machine_is_little_endian();

//   printf("%d %d\n", Niconf, Noconf); // Debug
   
   return 0;
}

inline int input(double* data, const int iconf)
{
   double *tmp_data = new double[DataSize];

   string ifname = StrReplace(ifpath, repstr, iconf_list[iconf]);
   FILE* fp;
   if ((fp = fopen(ifname.c_str(), "rb")) == NULL) {
      printf("The file '%s' is not found.\n", ifname.c_str());
      return 1;
   }

   fseek  (fp, 0, SEEK_END);
   size_t tmp_DataSize = ftell(fp);
   fseek  (fp, 0, SEEK_SET);
   if (tmp_DataSize != DataSize * 8) {
      printf("ERROR: Defferent file size, %s", ifname.c_str());
      return 1;
   }
   
   fread((void*)tmp_data, sizeof(double), DataSize, fp);
   if (endian_flip) endian_convert(tmp_data, DataSize);
   
   fclose(fp);
   
   for (int n=0; n<DataSize; n++) data[idx(n,iconf)] = tmp_data[n];
   
   delete [] tmp_data;
   
   return 0;
}

inline int output(const double* data, const int iconf)
{
   double *tmp_data = new double[DataSize];
   for (int n=0; n<DataSize; n++) tmp_data[n] = data[idx(n,iconf)];
   
   string ofname = StrReplace(ofpath, repstr, oconf_list[iconf]);
   FILE* fp;
   if ((fp = fopen(ofname.c_str(), "wb")) == NULL) {
      printf("The file '%s' can not be open.\n", ofname.c_str());
      return 1;
   }
   
   if (endian_flip) endian_convert(tmp_data, DataSize);
   
   fwrite((void*)tmp_data, sizeof(double), DataSize, fp);
   
   fclose(fp);
   
   delete [] tmp_data;

   return 0;
}

inline bool machine_is_little_endian()
{
  int endianTEST = 1;
  if (*(char*)&endianTEST) return true;
  else                     return false;
}

inline void endian_convert(double *DATA, const size_t DATA_size)
{
  char dummy[8];

  for (size_t k=0; k<DATA_size; k++) {
    for (int j=0; j<8; j++) dummy[j] = ((char*)&DATA[k])[j];
    for (int j=0; j<8; j++) ((char*)&DATA[k])[j] = dummy[7-j];
  }
}

inline void endian_convert(int32_t* DATA, const size_t DATA_size)
{
  char dummy[4];

  for (size_t k=0; k<DATA_size; k++ ) {
    for(int j=0; j<4; j++) dummy[j] = ((char*)&DATA[k])[j];
    for(int j=0; j<4; j++) ((char*)&DATA[k])[j] = dummy[3-j];
  }
}

inline int set_conf_list(const char* path_conf_file_list, char conf_list[MAX_N_DATA][MAX_LEN_PATH])
{
  int i = 0;
  
  ifstream ifs(path_conf_file_list, ios::in);
  if (!ifs) {
    printf("The file '%s' is not found.\n", path_conf_file_list);
    return 0;
  }

  while (ifs.getline(conf_list[i], sizeof(conf_list[i]))) i++;
  ifs.close();

  return i;
}

inline void MakeJKsample(double *iDATA, double *oDATA)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (int n=0; n<DataSize; n++) {
      double sum = 0.0, tmp_sum;
      for (int i=0; i<Niconf; i++) sum += iDATA[idx(n,i)];
      
      for (int ibin=0; ibin<Nbin; ibin++) {
         tmp_sum = sum;
	 for (int ioconf=0; ioconf<binSize; ioconf++) 
	   tmp_sum -= iDATA[idx(n, ioconf + binSize*ibin)];
	 
	 oDATA[idx(n, ibin)] = tmp_sum / double(Niconf-binSize);
      }
   }
}

inline string StrReplace(string Str, const string Str_from, const string Str_to)
{
  string::size_type  Pos(Str.find(Str_from));

  while(Pos != string::npos) {
    Str.replace(Pos, Str_from.length(), Str_to);
    Pos = Str.find(Str_from, Pos + Str_to.length());
  }
  return Str;
}
