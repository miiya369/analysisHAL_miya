//--------------------------------------------------------------------------
/**
 * @file
 * @ingroup ComplexField
 * @brief   Implementations of input/output complex field for binary data
 * @brief   The binary data will be stored as BIG endian
 * @author  Takaya Miyamoto
 * @since   Thu Sep 15 00:10:27 JST 2016
 */
//--------------------------------------------------------------------------

#include <ComplexField_Base.h>

#include <yukawa/compress48.h>
#include <yukawa/PH1.h>

using namespace yukawa;
typedef four_point_base::BC BC;
typedef PH1                 Core;
typedef four_point<Core>    Four_Point;
typedef compress48<Core>    Compress48;

namespace {
   int32_t MgNum_single_complex_field = 9190517;
}
//======================== miyamoto-format notation =======================//
//
//                        !! ALWAYS LITTLE ENDIAN !!
//
//      (1) Magic Number (9190517)  (32-bit integer)
//      (2) aSIZE (inner  DoF)      (32-bit integer)
//      (3) xSIZE (x-coordination)  (32-bit integer)
//      (4) ySIZE (y-coordination)  (32-bit integer)
//      (5) zSIZE (z-coordination)  (32-bit integer)
//      (6) tSIZE (t-coordination)  (32-bit integer)
//      (7) bSIZE (outner DoF)      (32-bit integer)
//
//      (8) data ((2)*(3)*(4)*(5)*(6)*(7)* 2 * 64-bit double float)
//
//         -> data[re/im+2*(a+aS*(x+xS*(y+yS*(z+zS*(t+tS*b)))))]
//
//=========================================================================//

void ComplexField_BASE::input_data_miya(const string ifile_name) {
   DEBUG_LOG
   ifstream ifs(ifile_name.c_str(), ios::in | ios::binary);
   if (!ifs) ERROR_FOPEN(ifile_name.c_str());
   
   int32_t MgNum;
   ifs.read((char*)&MgNum, sizeof(int32_t));
   if (!anaHAL::machine_is_little_endian()) anaHAL::endian_convert(&MgNum, 1);
   if (MgNum != MgNum_single_complex_field)
      ERROR_COMMENTS("This file may not be Miyamoto-format binary file.");
   
   int32_t tmp_a, tmp_x, tmp_y, tmp_z, tmp_t, tmp_b;
   ifs.read((char*)&tmp_a, sizeof(int32_t));
   ifs.read((char*)&tmp_x, sizeof(int32_t));
   ifs.read((char*)&tmp_y, sizeof(int32_t));
   ifs.read((char*)&tmp_z, sizeof(int32_t));
   ifs.read((char*)&tmp_t, sizeof(int32_t));
   ifs.read((char*)&tmp_b, sizeof(int32_t));
   
   if (!anaHAL::machine_is_little_endian()) {
      anaHAL::endian_convert(&tmp_a, 1);
      anaHAL::endian_convert(&tmp_x, 1);
      anaHAL::endian_convert(&tmp_y, 1);
      anaHAL::endian_convert(&tmp_z, 1);
      anaHAL::endian_convert(&tmp_t, 1);
      anaHAL::endian_convert(&tmp_b, 1);
   }
   mem_alloc(tmp_a, tmp_x, tmp_y, tmp_z, tmp_t, tmp_b);
   
   ifs.read ((char*)&m_field[0], sizeof(cdouble) * (*this).data_size());
   ifs.close();
   
   if (!anaHAL::machine_is_little_endian())
      anaHAL::endian_convert(m_field, (*this).data_size());
}
void ComplexField_BASE::output_data_miya(const string ofile_name) {
   DEBUG_LOG
   if ((*this).data_size() == 0)
      ERROR_COMMENTS("The data size has not been initialized.");
   
   ofstream ofs(ofile_name.c_str(), ios::out | ios::binary);
   if (!ofs) ERROR_FOPEN(ofile_name.c_str());
   
   int32_t MgNum = MgNum_single_complex_field;
   int32_t tmp_a = (*this).get_aSIZE();
   int32_t tmp_x = (*this).get_xSIZE();
   int32_t tmp_y = (*this).get_ySIZE();
   int32_t tmp_z = (*this).get_zSIZE();
   int32_t tmp_t = (*this).get_tSIZE();
   int32_t tmp_b = (*this).get_bSIZE();
   
   if (!anaHAL::machine_is_little_endian()) {
      anaHAL::endian_convert( &MgNum, 1);
      anaHAL::endian_convert( &tmp_a, 1);
      anaHAL::endian_convert( &tmp_x, 1);
      anaHAL::endian_convert( &tmp_y, 1);
      anaHAL::endian_convert( &tmp_z, 1);
      anaHAL::endian_convert( &tmp_t, 1);
      anaHAL::endian_convert( &tmp_b, 1);
      anaHAL::endian_convert(m_field, (*this).data_size());
   }
   ofs.write((char*)&MgNum, sizeof(int32_t));
   ofs.write((char*)&tmp_a, sizeof(int32_t));
   ofs.write((char*)&tmp_x, sizeof(int32_t));
   ofs.write((char*)&tmp_y, sizeof(int32_t));
   ofs.write((char*)&tmp_z, sizeof(int32_t));
   ofs.write((char*)&tmp_t, sizeof(int32_t));
   ofs.write((char*)&tmp_b, sizeof(int32_t));
   ofs.write((char*)&m_field[0], sizeof(cdouble) * (*this).data_size());
   ofs.close();
   
   if (!anaHAL::machine_is_little_endian())
      anaHAL::endian_convert(m_field, (*this).data_size());
}

//--------------------------------------------------------------------------
/**
 * @brief Input/Output bare binary data
 * @brief This function assumes that the file size equal to data size
 */
//--------------------------------------------------------------------------
void ComplexField_BASE::input_data_bin(const string ifile_name) {
   DEBUG_LOG
   ifstream ifs(ifile_name.c_str(), ios::in | ios::binary);
   if (!ifs) ERROR_FOPEN(ifile_name.c_str());
   
   size_t fSIZE = ifs.seekg(0, ifs.end).tellg(); ifs.seekg(0, ifs.beg);
   if (fSIZE != sizeof(cdouble) * (*this).data_size())
      ERROR_COMMENTS("Unexpected file size.");
   
   ifs.read ((char*)&m_field[0], sizeof(cdouble) * (*this).data_size());
   ifs.close();
   
   if (anaHAL::machine_is_little_endian())
      anaHAL::endian_convert(m_field, (*this).data_size());
}
void ComplexField_BASE::output_data_bin(const string ofile_name) {
   DEBUG_LOG
   if ((*this).data_size() == 0)
      ERROR_COMMENTS("The data size has not been initialized.");
   
   ofstream ofs(ofile_name.c_str(), ios::out | ios::binary);
   if (!ofs) ERROR_FOPEN(ofile_name.c_str());
   
   if (anaHAL::machine_is_little_endian())
      anaHAL::endian_convert(m_field, (*this).data_size());
   
   ofs.write((char*)&m_field[0], sizeof(cdouble) * (*this).data_size());
   ofs.close();
   
   if (anaHAL::machine_is_little_endian())
      anaHAL::endian_convert(m_field, (*this).data_size());
}

//--------------------------------------------------------------------------
/**
 * @brief Input/Output 1/48-compressed NBS wave function (ishii-san's format)
 * @brief This function is using yuakwa library
 */
//--------------------------------------------------------------------------
void ComplexField_BASE::input_data_comp(const string ifile_name, const int it) {
   DEBUG_LOG
   if ((*this).get_aSIZE() != 4 || (*this).get_bSIZE() != 4)
      ERROR_COMMENTS("Only #.external/internal index=4 is allowed "
                     "for reading 1/48-complex NBS wave function.");
   if ((*this).get_tSIZE() < it+1)
      ERROR_COMMENTS("Index of time Overflow");
   
   mapping48 map;
   map.read (ifile_name.c_str());
   
   Compress48 comp;
   comp.read(ifile_name.c_str(), map);
   
   BC Xbc = comp.Xbc;
   BC Ybc = comp.Ybc;
   BC Zbc = comp.Zbc;
   
   map.construct_map_and_phase(Xbc, Ybc, Zbc);
   
   Four_Point four = comp.decompress48();
   
   int Xsites = four.Xsites;
   int Ysites = four.Ysites;
   int Zsites = four.Zsites;
   if ((*this).get_xSIZE() != Xsites || (*this).get_ySIZE() != Ysites ||
       (*this).get_zSIZE() != Zsites) ERROR_COMMENTS("Unexpected x,y,z size.");
   
   int sign = 1;
   for(                  int betaP =0; betaP <2;       betaP++)
      for(               int alphaP=0; alphaP<2;       alphaP++)
         for(            int iz    =0; iz    <Zsites;  iz++)
            for(         int iy    =0; iy    <Ysites;  iy++)
               for(      int ix    =0; ix    <Xsites;  ix++)
                  for(   int beta  =0; beta  <2;       beta++)
                     for(int alpha =0; alpha <2;       alpha++)
                        (*this)(alpha+2*beta,ix,iy,iz,it,alphaP+2*betaP)
                        = four(ix,iy,iz,sign)(alpha,beta,alphaP,betaP);
}
void ComplexField_BASE::output_data_comp(const string   ofile_name,
                                         const string mapfile_name, const int it) {
   DEBUG_LOG
   if ((*this).get_aSIZE() != 4 || (*this).get_bSIZE() != 4)
      ERROR_COMMENTS("Only #.external/internal index=4 is allowed "
                     "for writing 1/48-complex NBS wave function.");
   if ((*this).get_tSIZE() < it+1)
      ERROR_COMMENTS("Index of time Overflow");
   
   mapping48 map;
   map.read(mapfile_name.c_str());
   
   int Xsites = map.Xsites;
   int Ysites = map.Ysites;
   int Zsites = map.Zsites;
   
   if ((*this).get_xSIZE() != Xsites || (*this).get_ySIZE() != Ysites ||
       (*this).get_zSIZE() != Zsites) ERROR_COMMENTS("Unexpected x,y,z size.");
   
   BC Xbc = four_point_base::PBC;
   BC Ybc = four_point_base::PBC;
   BC Zbc = four_point_base::PBC;
   
   Four_Point four(Xsites, Ysites, Zsites, Xbc, Ybc, Zbc);
   int sign = 1;
   for(                  int betaP =0; betaP <2;       betaP++)
      for(               int alphaP=0; alphaP<2;       alphaP++)
         for(            int iz    =0; iz    <Zsites;  iz++)
            for(         int iy    =0; iy    <Ysites;  iy++)
               for(      int ix    =0; ix    <Xsites;  ix++)
                  for(   int beta  =0; beta  <2;       beta++)
                     for(int alpha =0; alpha <2;       alpha++)
                        four(ix,iy,iz,sign)(alpha,beta,alphaP,betaP)
                        = (*this)(alpha+2*beta,ix,iy,iz,it,alphaP+2*betaP);
   
   four = noise_reduction_reflection(four);
   four = noise_reduction_rotation  (four);
   
   Compress48 comp(map, four);
   
   comp.write(ofile_name.c_str());
}
