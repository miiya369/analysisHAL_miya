//--------------------------------------------------------------------------
/**
 * @file
 * @brief
 * @ingroup Compression
 * @ingroup Header_File
 * @author  N.Ishii
 * @since   Sat May 22 01:52:19 JST 2010
 *
 * keep only those points (ix,iy,iz), which satisfy iz >= ix >= iy.
 */
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
/**
 * @defgroup Compression 1/48 compression related topics
 * @ingroup  General_Util
 */
//--------------------------------------------------------------------------

#ifndef IS_INCLUDED_COMPRESS48_H
#define IS_INCLUDED_COMPRESS48_H

#include <stdio.h>
#include <complex>

#include <yukawa/yukawa.h>
#include <yukawa/4pt.h>
#include <yukawa/mapping48.h>
#include <yukawa/misc.h>

#ifndef NDEBUG
#define NDEBUG 0
#endif

Start_Name_Space_Yukawa

template <class T>
struct compress48 {
  typedef std::complex<double> Complex;
  typedef four_point_base::BC  BC;

  const
  mapping48* map_ptr;

  int   length;
  T*    a;
  BC    Xbc, Ybc, Zbc;

  ~compress48();
  compress48();
  compress48(const mapping48& mapping,
	     const four_point<T>& four);

  compress48(const mapping48& mapping,
	     BC xbc = four_point_base::PBC,
	     BC ybc = four_point_base::PBC,
	     BC zbc = four_point_base::PBC);

  compress48(const compress48<T>& other);

  compress48<T>& operator=(const compress48<T>& other);

  const
  T& operator[](int i) const { return a[range(i)]; }
  T& operator[](int i)       { return a[range(i)]; }

  compress48<T>& operator+=(const compress48<T>& other);
  compress48<T>& operator-=(const compress48<T>& other);
  
  compress48<T>& unary_minus();

  compress48<T>& operator*=(double x);
  compress48<T>& operator*=(Complex x);
  compress48<T>& operator*=(const compress48<T>& other);

  compress48<T>& operator/=(double x);
  compress48<T>& operator/=(Complex x);

  four_point<T> decompress48() const;

  const
  compress48<T>& write(FILE* fp) const;
  compress48<T>& read( FILE* fp,const mapping48& map2);

  const
  compress48<T>& write(const char* fname) const;
  compress48<T>& read( const char* fname,const mapping48& map2);

  void is_same_map(const compress48<T>& other) const;
  int  range(int i) const;
};

///
template <class T> inline
void compress48<T>::is_same_map(const compress48<T>& other) const
{
  if(!NDEBUG){
    if (map_ptr != other.map_ptr) {
      fprintf(stderr,
	      "ERROR(compress48::is_same_map): inconsistent maps\n");
      exit(1);
    }
  }
}

/// range checker
template <class T> inline
int compress48<T>::range(int i) const
{
  if(!NDEBUG){
    if(i < 0 || length <= i){
      fprintf(stderr,
	      "ERROR(compress48::range(int)): "
	      "index out of range: i=%d, length = %d\n",
	      i, length);
      exit(1);
    }
  }
  return i;
}

/// destructor
template<class T> inline
compress48<T>::~compress48()
{
  if (a) delete[] a; a = 0;
}

/// the default constructor
template <class T> inline
compress48<T>::compress48()
  : map_ptr(0),
    length(0),
    a(0),
    Xbc(four_point_base::PBC),
    Ybc(four_point_base::PBC),
    Zbc(four_point_base::PBC)
{}

//--------------------------------------------------------------------------
/**
 * @brief a constructor for compression
 */
//--------------------------------------------------------------------------

template <class T> inline
compress48<T>::compress48(const mapping48& mapping, const four_point<T>& four)
  : map_ptr(&mapping),
    length(0),
    a(0),
    Xbc(four.Xbc),
    Ybc(four.Ybc),
    Zbc(four.Zbc)
{
  length = map_ptr->length;
  a      = new T[length];

  int sign = 1;
  int Xsites = map_ptr->Xsites;
  int Ysites = map_ptr->Ysites;
  int Zsites = map_ptr->Zsites;

  for(int index = 0; index < length; index++){
    int ixyz = map_ptr->rep[index];
    int ix = ixyz % Xsites; ixyz /= Xsites;
    int iy = ixyz % Ysites; ixyz /= Ysites;
    int iz = ixyz % Zsites;

    a[index] = four(ix,iy,iz,sign);
  }
}

//--------------------------------------------------------------------------
/**
 * @brief a constructor
 */
//--------------------------------------------------------------------------

template <class T> inline
compress48<T>::compress48(const mapping48& mapping, BC xbc, BC ybc, BC zbc)
  : map_ptr(&mapping),
    length(0),
    a(0),
    Xbc(xbc),
    Ybc(ybc),
    Zbc(zbc)
{
  length = map_ptr->length;
  a      = new T[length];
}


/// the copy constructor
template <class T> inline
compress48<T>::compress48(const compress48<T>& other)
  : map_ptr(other.map_ptr),
    length(other.length),
    a(0),
    Xbc(other.Xbc),
    Ybc(other.Ybc),
    Zbc(other.Zbc)
{
  a   = new T[length];
  for(int i = 0; i < length; i++){
    a[i] = other.a[i];
  }
}

/// the assignment operator
template <class T> inline
compress48<T>& compress48<T>::operator=(const compress48<T>& other)
{
  map_ptr = other.map_ptr;
  length  = other.length;

  if(a) delete[] a;
  a = new T[length];

  for(int i = 0; i < length; i++){
    a[i] = other.a[i];
  }

  Xbc = other.Xbc;
  Ybc = other.Ybc;
  Zbc = other.Zbc;

  return *this;
}

/// += operation
template <class T> inline
compress48<T>& compress48<T>::operator+=(const compress48<T>& other)
{
  is_same_map(other);
  for(int i = 0; i < length; i++)
    a[i] += other.a[i];
  return *this;
}

/// -= operation
template <class T> inline
compress48<T>& compress48<T>::operator-=(const compress48<T>& other)
{
  is_same_map(other);
  for(int i = 0; i < length; i++)
    a[i] -= other.a[i];
  return *this;
}

/// unary minus
template <class T> inline
compress48<T>& compress48<T>::unary_minus()
{
  for(int i = 0; i < length; i++)
    a[i].unary_minus();
  return *this;
}

/// scalar multiplication (double)
template <class T> inline
compress48<T>& compress48<T>::operator*=(double x)
{
  for(int i = 0; i < length; i++)
    a[i] *= x;
  return *this;
}

/// scalar multiplication (Complex)
template <class T> inline
compress48<T>& compress48<T>::operator*=(Complex x)
{
  for(int i = 0; i < length; i++)
    a[i] *= x;
  return *this;
}

/// component-wise multiplication (double)
template <class T> inline
compress48<T>& compress48<T>::operator*=(const compress48<T>& other)
{
  is_same_map(other);
  for(int i = 0; i < length; i++)
    a[i] *= other.a[i];
  return *this;
}

/// scalar division (double)
template <class T> inline
compress48<T>& compress48<T>::operator/=(double x)
{
  for(int i = 0; i < length; i++)
    a[i] /= x;
  return *this;
}

/// scalar division (Complex)
template <class T> inline
compress48<T>& compress48<T>::operator/=(Complex x)
{
  for(int i = 0; i < length; i++)
    a[i] /= x;
  return *this;
}

//--------------------------------------------------------------------------
/**
 * @brief decompression
 */
//--------------------------------------------------------------------------

template <class T> inline
four_point<T> compress48<T>::decompress48() const
{
  typedef mapping48::Intger Intger;

  int Xsites   = map_ptr->Xsites;
  int Ysites   = map_ptr->Ysites;
  int Zsites   = map_ptr->Zsites;
  int XYZsites = map_ptr->XYZsites;

  four_point<T> four(Xsites,Ysites,Zsites, Xbc,Ybc,Zbc);

  Complex* phase = map_ptr->phase;
  Intger*  map   = map_ptr->map;

  for(int ixyz = 0; ixyz < XYZsites; ixyz++){
    int index  = map[ixyz].index;
    int ith    = map[ixyz].ith;
    int parity = map[ixyz].parity;
    
    if(parity == 1){
      four[ixyz] = phase[ixyz] * rotate(a[index], ith);
    }
    else if(parity == -1){
      four[ixyz] = phase[ixyz] * reflection(rotate(a[index], ith));
    }
  }

  return four;
}

//--------------------------------------------------------------------------
/**
 * @brief I/O (write)
 *
 * The format is arranged as<br>
 * mapping48 map;
 * int32_t   size;                  // size = (2 + 3) * sizeof(int32_t)
 * int32_t   Xbc, Ybc, Zbc;         // PBC = 0, APBC = 1
 * int32_t   length, size_of_block; // size_of_block = sizeof(T)
 * int32_t   size;                  // size = (2 + 3) * sizeof(int32_t)
 *
 * int32_t   size;                  // size = length * size_of_block * sizeof(int32_t)
 * T         a[length];
 * int32_t   size;                  // size = length * size_of_block * sizeof(int32_t)
 *
 * Data is written adopting little endian.
 */
//--------------------------------------------------------------------------

template <class T> inline
const
compress48<T>& compress48<T>::write(FILE *fp) const
{
  bool is_change_endian = is_bigendian();

  map_ptr->write(fp);

  int32_t tmp[4];

  //
  // size = (2 + 3) * sizeof(int32_t)
  //
  tmp[0] = (2+3) * sizeof(int32_t);
  if(is_change_endian == true) change_endian(tmp, 1);
  if(fwrite(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr,
	    "ERROR(compress48::write): 0\n");
    exit(1);
  }
  
  //
  // int32_t Xbc, Ybc, Zbc;  (PBC = 0, APBC = 1)
  //
  tmp[0] = Xbc;
  tmp[1] = Ybc;
  tmp[2] = Zbc;
  if(is_change_endian == true) change_endian(tmp, 3);
  if(fwrite(tmp, sizeof(int32_t), 3, fp) != 3){
    fprintf(stderr,
	    "ERROR(compress48::write): 1\n");
    exit(1);
  }

  //
  // int32_t length, size_of_block;
  //
  tmp[0] = length;
  tmp[1] = sizeof(T);
  if(is_change_endian == true) change_endian(tmp, 2);
  if(fwrite(tmp, sizeof(int32_t), 2, fp) != 2){
    fprintf(stderr,
	    "ERROR(compress48::write): 2\n");
    exit(1);
  }

  //
  // size = (2 + 3) * sizeof(int32_t)
  //
  tmp[0] = (2+3) * sizeof(int32_t);
  if(is_change_endian == true) change_endian(tmp, 1);
  if(fwrite(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr,
	    "ERROR(compress48::write): 3\n");
    exit(1);
  }
  


  //
  // size = length * size_of_block
  //
  tmp[0] = length * sizeof(T);
  if(is_change_endian == true) change_endian(tmp, 1);
  if(fwrite(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr,
	    "ERROR(compress48::write): 4\n");
    exit(1);
  }
  
  //
  // T a[length];
  //
  for(int i = 0; i < length; i++){
    a[i].write(fp, is_change_endian);
  }

  //
  // size = length * size_of_block
  //
  tmp[0] = length * sizeof(T);
  if(is_change_endian == true) change_endian(tmp, 1);
  if(fwrite(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr,
	    "ERROR(compress48::write): 4\n");
    exit(1);
  }

  return *this;
}

//--------------------------------------------------------------------------
/**
 * @brief I/O (read)
 *
 * Before reading a contents, read a common map for compression.
 *
 * Data is assumed to be stored adopting little endian.
 */
//--------------------------------------------------------------------------

template <class T> inline
compress48<T>& compress48<T>::read(FILE *fp,const mapping48& map2)
{
  bool is_change_endian = is_bigendian();

  mapping48 map3;
  map3.read(fp);

  if(map2.is_same_mapping(map3) == false){
    fprintf(stderr,
	    "ERROR(compress48::read): inconsistent mapping\n");
    exit(1);
  }
  map_ptr = &map2;

  int32_t tmp[4];


  //
  // size = (2 + 3) * sizeof(int32_t);
  //
  if(fread(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr, "ERROR(compress48::read): 0\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 1);
  if(tmp[0] != 5 * sizeof(int32_t)){
    fprintf(stderr, "ERROR(compress48::read): 1\n");
    exit(1);
  }

  //
  // int32_t Xbc, Ybc, Zbc;  (PBC = 0, APBC = 1)
  //
  if(fread(tmp, sizeof(int32_t), 3, fp) != 3){
    fprintf(stderr, "ERROR(compress48::read): 2\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 3);
  switch(tmp[0]){
  case four_point_base::PBC:  Xbc = four_point_base::PBC;  break;
  case four_point_base::APBC: Xbc = four_point_base::APBC; break;
  default:
    break;
  }
  switch(tmp[1]){
  case four_point_base::PBC:  Ybc = four_point_base::PBC;  break;
  case four_point_base::APBC: Ybc = four_point_base::APBC; break;
  default:
    break;
  }
  switch(tmp[2]){
  case four_point_base::PBC:  Zbc = four_point_base::PBC;  break;
  case four_point_base::APBC: Zbc = four_point_base::APBC; break;
  default:
    fprintf(stderr, "ERROR(compress48::read): 3\n");
    exit(1);
    break;
  }

  //
  // int32_t length, size_of_block;
  //
  if(fread(tmp, sizeof(int32_t), 2, fp) != 2){
    fprintf(stderr, "ERROR(compress48::read): 4\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 2);
  length = (int)tmp[0];
  if(tmp[1] != sizeof(T)){
    fprintf(stderr, "ERROR(compress48::read): 4\n");
    exit(1);
  }
  if(a) delete[] a;
  a = new T[length];

  //
  // size = (2 + 3) * sizeof(int32_t);
  //
  if(fread(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr, "ERROR(compress48::read): 0\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 1);
  if(tmp[0] != 5 * sizeof(int32_t)){
    fprintf(stderr, "ERROR(compress48::read): 5\n");
    exit(1);
  }



  //
  // size = length * sizeof(T)
  //
  if(fread(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr, "ERROR(compress48::read): 6\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 1);
  if((int)tmp[0] != length * (int)sizeof(T)){
    fprintf(stderr, "ERROR(compress48::read): 7\n");
    exit(1);
  }

  //
  // T a[length];
  //
  for(int i = 0; i < length; i++){
    a[i].read(fp,is_change_endian);
  }

  //
  // size = length * sizeof(T)
  //
  if(fread(tmp, sizeof(int32_t), 1, fp) != 1){
    fprintf(stderr, "ERROR(compress48::read): 8\n");
    exit(1);
  }
  if(is_change_endian == true) change_endian((int32_t*)tmp, 1);
  if((int)tmp[0] != length * (int)sizeof(T)){
    fprintf(stderr, "ERROR(compress48::read): 9\n");
    exit(1);
  }

  return *this;
}

/// @overload
template <class T> inline
const compress48<T>&
compress48<T>::write(const char* fname) const
{
  bool is_popen = false;

  if(strcmp(&fname[strlen(fname)-3],".gz")==0){
    is_popen = true;
  }

  FILE* fp;
  if(is_popen == true){
    char str[2048];
    sprintf(str, "gzip -c > %s", fname);
    if((fp = popen(str, "w"))==NULL){
      fprintf(stderr, "cannot popen '%s'\n", str);
      exit(1);
    }
  }
  else{
    if((fp = fopen(fname, "w"))==NULL){
      fprintf(stderr, "cannot fopen '%s'\n", fname);
      exit(1);
    }
  }
  write(fp);
  if(is_popen == true){
    pclose(fp);
  }
  else{
    fclose(fp);
  }
  return *this;
}

/// @overload
template <class T> inline
compress48<T>& compress48<T>::read(const char* fname,
				   const mapping48& map2)
{
  bool is_popen = false;
  
  if(strcmp(&fname[strlen(fname)-3],".gz")==0){
    is_popen = true;
  }

  FILE* fp;
  if(is_popen == true){
    char str[2048];
    sprintf(str, "zcat %s", fname);
    if((fp = popen(str, "r"))==NULL){
      fprintf(stderr, "cannot popen '%s'\n", str);
      exit(1);
    }
  }
  else{
    if((fp = fopen(fname, "r"))==NULL){
      fprintf(stderr, "cannot fopen '%s'\n", fname);
      exit(1);
    }
  }
  read(fp,map2);
  if(is_popen == true){
    pclose(fp);
  }
  else{
    fclose(fp);
  }
  return *this;
}

//--------------------------------------------------------------------------
// helper functions
//--------------------------------------------------------------------------

/// addition
template <class T> inline
compress48<T> operator+(compress48<T> lhs,const compress48<T>& rhs)
{ return lhs += rhs; }

/// subtraction
template <class T> inline
compress48<T> operator-(compress48<T> lhs,const compress48<T>& rhs)
{ return lhs -= rhs; }

/// unary minus
template <class T> inline
compress48<T> operator-(compress48<T> comp)
{ return comp.unary_minus(); }


/// scalar multiplication (double)
template <class T> inline
compress48<T> operator*(compress48<T> comp, double x)
{ return comp *= x; }

/// scalar multiplication (double) 2
template <class T> inline
compress48<T> operator*(double x, compress48<T> comp)
{ return comp *= x; }

/// scalar multiplication (Complex)
template <class T> inline
compress48<T> operator*(compress48<T> comp, std::complex<double> x)
{ return comp *= x; }

/// scalar multiplication (Complex) 2
template <class T> inline
compress48<T> operator*(std::complex<double> x, compress48<T> comp)
{ return comp *= x; }

/// component-wise multiplication
template <class T> inline
compress48<T> operator*(compress48<T> lhs, const compress48<T>& rhs)
{ return lhs *= rhs; }


/// scalar division (double)
template <class T> inline
compress48<T> operator/(compress48<T> comp, double x)
{ return comp /= x; }

/// scalar multiplication (Complex)
template <class T> inline
compress48<T> operator/(compress48<T> comp, std::complex<double> x)
{ return comp /= x; }



End_Name_Space_Yukawa

#endif
