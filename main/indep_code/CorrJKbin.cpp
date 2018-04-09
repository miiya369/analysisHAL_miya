//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup none
 * @brief   Make Jack-Knife bin samples for Text 2pt-correlator data
 * @author  Takaya Miyamoto
 * @since   Mon Sep 19 04:17:44 JST 2016
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

#define NCOL 10 // #.column of correlator files

#define MAX_N_DATA   1024
#define MAX_LEN_PATH 1024

#define idx(n,i) ((n) + DataSize * (i))

static string ilist;
static string olist;
static int    Tsize, Nbin, binSize;
static string ifpath;
static string ofpath;
static string repstr;

static char iconf_list[MAX_N_DATA][MAX_LEN_PATH];
static char oconf_list[MAX_N_DATA][MAX_LEN_PATH];

static int Niconf, Noconf;

static size_t DataSize;

inline int set_args(const int, const char**);

inline int input (double*, const int);
inline int output(const double*, const int);

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
   
   FILE* fp;
   if ((fp = fopen(StrReplace(ifpath, repstr, iconf_list[0]).c_str(), "r")) == NULL) {
     printf("The file '%s' is not found.\n", StrReplace(ifpath, repstr, iconf_list[0]).c_str());
     return 1;
   }
   int count = 0;
   char tmp[1000];
   while(fgets(tmp, 1000, fp) != NULL) count++;
   fclose(fp);

   Tsize    = count;
   DataSize = Tsize * NCOL;

//   printf("%d %d %d\n", DataSize, Niconf, Noconf); // Debug
   
   return 0;
}

inline int input(double* data, const int iconf)
{
   double *tmp_data = new double[DataSize];

   string ifname = StrReplace(ifpath, repstr, iconf_list[iconf]);
   FILE* fp;
   if ((fp = fopen(ifname.c_str(), "r")) == NULL) {
     printf("The file '%s' is not found.\n", ifname.c_str());
     return 1;
   }

   int    t_tmp;
   double tmp[10];
   for (int it=0; it<Tsize; it++) {
     if (fscanf(fp, " %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &t_tmp
		, &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4]
		, &tmp[5], &tmp[6], &tmp[7], &tmp[8], &tmp[9] ) != 11) {
       printf("ERROR: Invalid #.column '%s'\n", ifname.c_str());
       return 1;
     }
     if (it != t_tmp) {
       printf("ERROR: Invalid time size '%s'\n", ifname.c_str());
       return 1;
     }
     for (int loop=0; loop<10; loop++) data[idx(loop+10*it, iconf)] = tmp[loop];
   }
   fclose(fp);

   delete [] tmp_data;
   
   return 0;
}

inline int output(const double* data, const int iconf)
{
   double *tmp_data = new double[DataSize];
   for (int n=0; n<DataSize; n++) tmp_data[n] = data[idx(n,iconf)];
   
   string ofname = StrReplace(ofpath, repstr, oconf_list[iconf]);
   FILE* fp;
   if ((fp = fopen(ofname.c_str(), "w")) == NULL) {
     printf("The file '%s' can not be open.\n", ofname.c_str());
     return 1;
   }

   for (int it=0; it<Tsize; it++)
     if ( fprintf(  fp, "%d %1.16e %1.16e %1.16e %1.16e "
                    "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", it
                    , tmp_data[0+10*it], tmp_data[1+10*it], tmp_data[2+10*it]
                    , tmp_data[3+10*it], tmp_data[4+10*it], tmp_data[5+10*it]
                    , tmp_data[6+10*it], tmp_data[7+10*it]
                    , tmp_data[8+10*it], tmp_data[9+10*it] ) < 0 ) {
       printf("invalid writing for file '%s'\n", ofname.c_str());
       return 1;
     }
   fclose(fp);

   delete [] tmp_data;

   return 0;
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
