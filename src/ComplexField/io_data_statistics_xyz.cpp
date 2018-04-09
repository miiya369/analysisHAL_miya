//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Function for output the 3-demensional complex field data with error
 * @author  Takaya Miyamoto
 * @since   Thu Sep 15 00:10:27 JST 2016
 */
//--------------------------------------------------------------------------

#include <StatisticsTemplate.h>
#include <ComplexField_Sub.h>

//--------------------------------------------------------------------------
/**
 * @brief Function to output data with error
 * @brief with 1/48 data size reduction for 3D-Field type data
 */
//--------------------------------------------------------------------------
template <> void STATISTICS<ComplexField_XYZ> ::
output_data_err(const string outfile_name,
                const double lattice_spacing,
                const bool   is_jack_knife_data) {
   DEBUG_LOG
   FILE *fp;
   if ((fp = fopen(outfile_name.c_str(), "w")) == NULL)
      ERROR_FOPEN (outfile_name.c_str());
   
   double MeVfactor = 1.0, lat_sp = 1.0;
   char   Unit[16]  = "[Lattice Unit]";
   if (lattice_spacing != 0.0) {
      MeVfactor = hbar_c / lattice_spacing;
      lat_sp    = lattice_spacing;
      snprintf(Unit, sizeof(Unit), "[MeV]");
   }
   
   cdouble *tmp_mean = new cdouble[(*this)(0).data_size()];
   cdouble *tmp_err  = new cdouble[(*this)(0).data_size()];
   
   fprintf(fp, "# r, mean.real, err.real, mean.imag, err.imag %s\n", Unit);
   (*this).make_mean_err(tmp_mean, tmp_err, is_jack_knife_data);
   
   int xyz;
   int l_xSIZE = (*this)(0).get_xSIZE();
   int l_ySIZE = (*this)(0).get_ySIZE();
   int l_zSIZE = (*this)(0).get_zSIZE();
   double R;
   for (      int z=0; z<l_zSIZE; z++)
      for (   int y=0; y<l_ySIZE; y++)
         for (int x=0; x<l_xSIZE; x++) {
            if (x>l_xSIZE/2.0 || y>x || z>y) continue;
            R   = sqrt(x*x + y*y + z*z) * lat_sp;
            xyz = x + l_xSIZE*(y + l_ySIZE*z);
            
            fprintf(fp, "%1.16e %1.16e %1.16e %1.16e %1.16e\n", R,
                    tmp_mean[xyz].real() * MeVfactor,
                    tmp_err [xyz].real() * MeVfactor,
                    tmp_mean[xyz].imag() * MeVfactor,
                    tmp_err [xyz].imag() * MeVfactor);
         }
   
   delete [] tmp_mean;
   delete [] tmp_err;
   fclose(fp);
}

//--------------------------------------------------------------------------
/**
 * @brief Function for binary data output
 * @brief with 1/48 data size reduction for 3D-field type data
 */
//--------------------------------------------------------------------------

namespace { int32_t MgNum_multi_reduced_complex_field_3D = 19900518; }

//======================== miyamoto-format notation =======================//
//
//                        !! ALWAYS LITTLE ENDIAN !!
//
//      1) Magic Number (19900518)  (int)
//      2) #.conf                   (int)
//      3) #.data/conf              (int)
//      4) bytes of data            (int)
//
//      5) data coordinates         ((3) * 8 bytes)
//         -> for n = 0 to #.data
//               crd[n]
//
//      6) data                     ((2)*(3)*(4) bytes)
//         -> for    i = 0 to #.conf
//               for n = 0 to #.data
//                  data[n+#.data*i]
//
//=========================================================================//

template <> void STATISTICS<ComplexField_XYZ> ::
output_data_bin_reduce(const string outfile_name,
                       const double lattice_spacing,
                       const bool   is_complex_data) {
   DEBUG_LOG
   ofstream ofs(outfile_name.c_str(), ios::out | ios::binary);
   if (!ofs) ERROR_FOPEN(outfile_name.c_str());
   
   int Magic_Num = MgNum_multi_reduced_complex_field_3D;
   int Conf__Num = (*this).Ndata();
   int Data__Num = (*this)(0).data_size();
   int Data_Byte;
   if (is_complex_data) Data_Byte = sizeof(cdouble);
   else                 Data_Byte = sizeof( double);
   
   Data__Num = anaHAL::reduced_Ndata((*this)(0).get_xSIZE(),
                                     (*this)(0).get_ySIZE(),
                                     (*this)(0).get_zSIZE());
   
   if (!anaHAL::machine_is_little_endian()) {
      anaHAL::endian_convert(&Magic_Num, 1);
      anaHAL::endian_convert(&Conf__Num, 1);
      anaHAL::endian_convert(&Data__Num, 1);
      anaHAL::endian_convert(&Data_Byte, 1);
   }
   ofs.write((char*)&Magic_Num, sizeof(int));
   ofs.write((char*)&Conf__Num, sizeof(int));
   ofs.write((char*)&Data__Num, sizeof(int));
   ofs.write((char*)&Data_Byte, sizeof(int));
   
   double MeVfactor = 1.0, lat_sp = 1.0;
   if (lattice_spacing != 0.0) {
      MeVfactor = hbar_c / lattice_spacing;
      lat_sp    = lattice_spacing;
   }
   double  dtmp;
   cdouble ctmp;
   int l_xSIZE = (*this)(0).get_xSIZE();
   int l_ySIZE = (*this)(0).get_ySIZE();
   int l_zSIZE = (*this)(0).get_zSIZE();
   for (      int z=0; z<l_zSIZE; z++)
      for (   int y=0; y<l_ySIZE; y++)
         for (int x=0; x<l_xSIZE; x++) {
            if (x>l_xSIZE/2.0 || y>x || z>y) continue;
            dtmp = sqrt(x*x + y*y + z*z) * lat_sp;
            if (!anaHAL::machine_is_little_endian())
               anaHAL::endian_convert(&dtmp, 1);
            ofs.write((char*)&dtmp, sizeof(double));
         }
   for(          int i=0; i<(*this).Ndata(); i++)
      for (      int z=0; z<l_zSIZE; z++)
         for (   int y=0; y<l_ySIZE; y++)
            for (int x=0; x<l_xSIZE; x++) {
               if (x>l_xSIZE/2.0 || y>x || z>y) continue;
               if (is_complex_data) {
                  ctmp = (*this)(i)(x,y,z) * MeVfactor;
                  if (!anaHAL::machine_is_little_endian())
                     anaHAL::endian_convert(&ctmp, 1);
                  ofs.write((char*)&ctmp, sizeof(cdouble));
               }
               else {
                  dtmp = (*this)(i)(x,y,z).real() * MeVfactor;
                  if (!anaHAL::machine_is_little_endian())
                     anaHAL::endian_convert(&dtmp, 1);
                  ofs.write((char*)&dtmp, sizeof(double));
               }
            }
   ofs.close();
}
