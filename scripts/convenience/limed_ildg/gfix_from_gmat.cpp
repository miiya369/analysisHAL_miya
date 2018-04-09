extern "C"
{
#include <lime_config.h>
#include <lime.h>
#include <lime_fixed_types.h>
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gauge_link.h>

#include <iostream>
#include <fstream>

using namespace std;

int ssize;
int tsize;

static double plaquette_criterion = 10e-12;

bool machine_is_little_endian();
void change_endian(double*, size_t);
void read_conf(const char*,double*);

int main(int argc, char **argv)
{
   if (argc != 6) {
      printf(" usage : %s [space size] [time size] [iconf] [Gmat] [oconf]\n", argv[0]);
      printf(" usage : The gauge configuration is must be ildg-lime format.\n");
      return -1;
   }
   ssize = atoi(argv[1]);
   tsize = atoi(argv[2]);
   size_t conf_size = ssize * ssize * ssize * tsize * 4 * 3 * 3 * 2;
   size_t gmat_size = ssize * ssize * ssize * tsize * 3 * 3 * 2;
   
   double *conf = new double[conf_size];
   read_conf(argv[3], conf);
   
   ifstream ifs(argv[4], ios::in | ios::binary);
   if (!ifs) {
      printf(" ERROR : The file \"%s\" is  not found.\n", argv[4]);
      return -1;
   }
   double *gmat = new double[gmat_size];
   ifs.read((char*)&gmat[0], sizeof(double)*gmat_size);   ifs.close();
   
   double before_plaq = plaquette(conf, ssize, tsize);
   
   gfix_by_gmat(conf, gmat, ssize, tsize);
   
   double after_plaq = plaquette(conf, ssize, tsize);
   
   if (abs(before_plaq - after_plaq) > plaquette_criterion) {
     printf("@@@@@@ WARNING @@@@@@ Differ value of plaquette, %1.10lf != %1.10lf, Doesn't Output.\n",
            before_plaq, after_plaq);
     
     delete [] conf;
     delete [] gmat;
     
     return -1;
   }
   
   if (machine_is_little_endian()) change_endian(conf, conf_size);
   
   ofstream ofs(argv[5], ios::out | ios::binary);
   ofs.write((char*)&conf[0], sizeof(double)*conf_size);   ofs.close();
   
   delete [] conf;
   delete [] gmat;
   
   return 0;
}

bool machine_is_little_endian() {

  int endianTEST = 1;
  if (*(char*)&endianTEST) return true;
  else                     return false;
}

void change_endian(double* DATA, size_t DATA_size)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t k=0; k<DATA_size; k++)
    {
      char dummy[8];
      for (int j=0; j<8; j++) dummy[j] = ((char*)&DATA[k])[j];
      for (int j=0; j<8; j++) ((char*)&DATA[k])[j] = dummy[7-j];
    }
}

void read_conf(const char* fname, double* conf)
{
  LimeReader *reader;
  int        status;
  n_uint64_t nbytes, read_bytes;
  char*      lime_type;

  FILE* fp;
  if((fp = fopen(fname, "r"))==NULL){
    fprintf(stderr,
            "cannot fopen '%s' for reading\n", fname);
    exit(1);
  }

  if((reader = limeCreateReader(fp)) == NULL){
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(1);
  }
  while((status = limeReaderNextRecord(reader)) != LIME_EOF){
    if(LIME_SUCCESS != status){
      fprintf(stderr, "limeReaderNextRecord returns status = %d\n", status);
      exit(1);
    }

    nbytes    = limeReaderBytes   (reader);
    lime_type = limeReaderType    (reader);

    if(strcmp(lime_type,"ildg-binary-data")==0)
      {
        if(2 *3*3*4 *ssize *ssize *ssize *tsize *sizeof(double) != nbytes){
          fprintf(stderr,
                  "ERROR(misc::read_ildg_config): inconsistent data size ! %ld != %ld\n",
                  nbytes,
                  2*3*3*4 *ssize *ssize *ssize *tsize *sizeof(double));
          exit(1);
        }

        read_bytes = nbytes;
        status = limeReaderReadData((void*)conf, &read_bytes, reader);

        if(status < 0){
          if(status != LIME_EOR){
            fprintf(stderr,
                    "LIME read error ccurred in reading ildg-binary-data: "
                    "status = %d %llu bytes wanted, %llu read\n",
                    status,
                    (unsigned long long)nbytes,
                    (unsigned long long)read_bytes);
            exit(1);
          }
        }

        if(machine_is_little_endian()){
          change_endian((double*)conf, read_bytes/sizeof(double));
        }
      }
  }
  limeDestroyReader(reader);
  fclose(fp);
}
