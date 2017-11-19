
#ifndef CVecH
#define CVecH


#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string>


#include <opencv2/opencv.hpp>

using namespace std;

typedef long HRESULT;

template <class T> inline T armsMagnitudeSq(const T &v) { return v*v; }
template <class T> inline T armsReScalar(const T &v) { return v; }
template <class T> inline T armsMagnitude(const T &v) { return fabs(v); }


//template <class ET> class CVec;
template <class ET> class CVec2;
template <class ET> class CVec3;
template <class ET> class CVec4;
//
template <class ET> class CMtx;
template <class ET> class CMtx2x2;
template <class ET> class CMtx3x3;
template <class ET> class CMtx4x4;

////////////////////////////////////CVec/////////////////////////////////
template <class ET>
class CVec 
{
public:
    CVec() : m_iSize(0), m_p(NULL), m_bWrap(false) {}
    CVec(const CVec<ET> &v) : m_iSize(0), m_p(NULL) { operator=(v); }
    CVec(int iSize) : m_iSize(0), m_p(NULL) { Create(iSize); }
    CVec(const ET *p, int iSize) : m_iSize(0), m_p(NULL) { Create(iSize); Init(p, iSize); }
    CVec(const CVec2<ET> &v) : m_iSize(0), m_p(NULL) { operator=(v); }
    CVec(const CVec3<ET> &v) : m_iSize(0), m_p(NULL) { operator=(v); }
    CVec(const CVec4<ET> &v) : m_iSize(0), m_p(NULL) { operator=(v); }

    virtual ~CVec();

    HRESULT Create(int iSize);
    void Free(); // deallocates memory and sets size to 0x0.

    void Wrap(ET *p, int iSize); // wrap around existing memory

    CVec<ET> & operator= (const CVec<ET> &v);
    CVec<ET> & operator= (const CVec2<ET> &v);
    CVec<ET> & operator= (const CVec3<ET> &v);
    CVec<ET> & operator= (const CVec4<ET> &v);

    int Size() const { return m_iSize; }

    const ET *Ptr() const { return m_p; }
    ET *Ptr() { return m_p; }

    const ET &El(int i) const { return m_p[i]; }
    ET &El(int i) { return m_p[i]; }
    const ET &operator[] (int i) const { return m_p[i]; }
    ET &operator[] (int i) { return m_p[i]; }
    const ET &operator() (int i) const { return m_p[i]; }
    ET &operator() (int i) { return m_p[i]; }

    CVec<ET> operator+ () const { return *this; }
    CVec<ET> operator- () const;
    CVec<ET> operator+ (const CVec<ET> &v) const;
    CVec<ET> operator- (const CVec<ET> &v) const;
    CMtx<ET> operator^ (const CVec<ET> &v) const;
    ET operator* (const CVec<ET> &v) const;
    CVec<ET> operator* (const ET &f) const;
    CVec<ET> operator* (const CMtx<ET> &m) const;
    CVec<ET> operator/ (const ET &f) const;
    bool operator== (const CVec<ET> &v) const;
    bool operator!= (const CVec<ET> &v) const { return !(operator==(v)); }

    ET MagnitudeSq() const;
    ET Magnitude() const { return (ET)sqrt(armsReScalar(MagnitudeSq())); }
    CVec<ET> Unit() const { return *this/Magnitude(); }
    //CVec<ET> Dehom() const; // convert to non-homogenous system
    //CVec<ET> Hom() const; // convert to homogenous system
    //void SortAbs(int *piMapTable) const;
    //int MaxAbs() const;
    //CVec<ET> Conj() const;
    //CVec<ET> Extract(int iSize, int iStart) const;
    //CVec<ET> & Update(const CVec<ET> &vec, int iStart);

    CVec<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CVec<ET> & operator= (const ET &f);
    CVec<ET> & operator+= (const CVec<ET> &v);
    CVec<ET> & operator-= (const CVec<ET> &v);
    CVec<ET> & operator*= (const CMtx<ET> &m) { *this = *this * m; return *this; }
    CVec<ET> & operator*= (const ET &f);
    CVec<ET> & operator/= (const ET &f);



private:
    int m_iSize;
    ET *m_p;
    bool m_bWrap;
};

template <class ET>
class CMtx
{
public:
    CMtx() : m_iRows(0), m_iCols(0), m_p(NULL), m_bWrap(false) {}

    CMtx(const CMtx<ET> &m) : m_iRows(0), m_iCols(0), m_p(NULL) { operator=(m); }
    CMtx(int iRows, int iCols) : m_iRows(0), m_iCols(0), m_p(NULL) { Create(iRows, iCols); }
    CMtx(const ET *p, int iRows, int iCols) : m_iRows(0), m_iCols(0), m_p(NULL) 
        { Create(iRows, iCols); Init(p, iRows*iCols); }
    CMtx(const CMtx2x2<ET> &m) : m_iRows(0), m_iCols(0), m_p(NULL) { operator=(m); }
    CMtx(const CMtx3x3<ET> &m) : m_iRows(0), m_iCols(0), m_p(NULL) { operator=(m); }
    CMtx(const CMtx4x4<ET> &m) : m_iRows(0), m_iCols(0), m_p(NULL) { operator=(m); }

    virtual ~CMtx();

    HRESULT Create(int iRows, int iCols);
    void Free(); // deallocates memory and sets size to 0x0.

    void Wrap(ET *p, int iRows, int iCols); // wrap around existing memory

    CMtx<ET> & operator= (const CMtx<ET> &m);
    CMtx<ET> & operator= (const CMtx2x2<ET> &m);
    CMtx<ET> & operator= (const CMtx3x3<ET> &m);
    CMtx<ET> & operator= (const CMtx4x4<ET> &m);

    int Rows() const { return m_iRows; }
    int Cols() const { return m_iCols; }
    int Size() const { return Rows() * Cols(); }

    const ET *Ptr() const { return m_p; }
    ET *Ptr() { return m_p; }

    CVec<ET> GetRow(int iRow) const;
    CVec<ET> GetCol(int iCol) const;
    const ET &El(int iRow, int iCol) const { return m_p[iRow * Cols() + iCol]; }
    ET &El(int iRow, int iCol) { return m_p[iRow * Cols() + iCol]; }
    const ET *operator[] (int i) const { return m_p + i * Cols(); }
    const ET &operator() (int iRow, int iCol) const { return m_p[iRow * Cols() + iCol]; }
    ET *operator[] (int i) { return m_p + i * Cols(); }
    ET &operator() (int iRow, int iCol) { return El(iRow, iCol); }

    //CMtx<ET> & MakeI();
    //CMtx<ET> & MakeDiag(const CVec<ET> &v);
    CMtx<ET> T() const;
    //CMtx<ET> H() const; // conjugate transpose
    //CMtx<ET> Conj() const; // conjugate
    //ET Tr() const;
    //ET Det() const;
    //CMtx<ET> Inv() const; // gaussian elimination with partial pivot
    //void MaxAbs(int &row, int &col) const; // get maximum element
    //CMtx<ET> & MakeSymmetric(); // copies below diagonal entries from above diagonal

    //CMtx<ET> & SetRow(int iRow, const CVec<ET> &v);
    //CMtx<ET> & SetCol(int iCol, const CVec<ET> &v);
    //CVec<ET> GetDiag() const;
    //CMtx<ET> DeleteRowCol(int iRow, int iCol) const;

    CMtx<ET> GetCols(const vector<int> &vCols) const;
    CMtx<ET> GetRows(const vector<int> &vRows) const;
    CMtx<ET> & SetRows(const vector<int> &vRows, const CMtx<ET> &m);
    CMtx<ET> & SetCols(const vector<int> &vCols, const CMtx<ET> &m);

    //CMtx<ET> & Update(int iRow, int iCol, const CMtx<ET> &m); // update sub-matrix
    //CMtx<ET> Extract(int iRow, int iCol, int iRowCount = -1, int iColCount = -1) const; // extract matrix

    CMtx<ET> operator+ () const { return *this; }
    CMtx<ET> operator- () const;
    CMtx<ET> operator+ (const CMtx<ET> &m) const;
    CMtx<ET> operator- (const CMtx<ET> &m) const;
    CMtx<ET> operator* (const CMtx<ET> &m) const;
    CMtx<ET> operator* (const ET &f) const;
    CVec<ET> operator* (const CVec<ET> &v) const;
    CMtx<ET> operator/ (const ET &f) const;
    bool operator== (const CMtx<ET> &m) const;
    bool operator!= (const CMtx<ET> &m) const { return !(operator==(m)); }

    CMtx<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CMtx<ET> & operator= (const ET &f);
    CMtx<ET> & operator+= (const CMtx<ET> &m);
    CMtx<ET> & operator-= (const CMtx<ET> &m);
    CMtx<ET> & operator*= (const CMtx<ET> &m) { *this = *this * m; return *this; }
    CMtx<ET> & operator*= (const ET &f);
    CMtx<ET> & operator/= (const ET &f);



private:
    int m_iRows, m_iCols;
    ET *m_p;
    bool m_bWrap;
};

typedef CVec<float> CVecf;
typedef CVec<double> CVecd;
typedef CMtx<float> CMtxf;
typedef CMtx<double> CMtxd;


template <class ET>
class CVec2
{
public:
    CVec2() {}
    CVec2(const CVec<ET> &v) { operator=(v); }
    CVec2(const CVec2<float> &v) { x=(ET)v.x; y=(ET)v.y; }
    CVec2(const CVec2<double> &v) { x=(ET)v.x; y=(ET)v.y; }
    CVec2(const ET &f) { operator=(f); }
    CVec2(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CVec2(const ET &m0, const ET &m1) { x = m0; y = m1; }

    int Size() const { return 2; }

    const ET *Ptr() const { return &x; }
    ET *Ptr() { return &x; }

    const ET &El(int i) const { return (&x)[i]; }
    ET &El(int i) { return (&x)[i]; }
    const ET &operator[] (int i) const { return (&x)[i]; }
    ET &operator[] (int i) { return (&x)[i]; }
    const ET &operator() (int i) const { return (&x)[i]; }
    ET &operator() (int i) { return (&x)[i]; }

    CVec2<ET> operator+ () const { return *this; }
    CVec2<ET> operator- () const { return CVec2<ET>(-El(0), -El(1)); }
    CVec2<ET> operator+ (const CVec2<ET> &v) const { return CVec2<ET>(El(0) + v[0], El(1) + v[1]); }
    CVec2<ET> operator- (const CVec2<ET> &v) const { return CVec2<ET>(El(0) - v[0], El(1) - v[1]); }
    CMtx2x2<ET> operator^ (const CVec2<ET> &v) const;
    ET operator* (const CVec2<ET> &v) const { return v[0] * El(0) + v[1] * El(1); }
    CVec2<ET> operator* (const ET &f) const { return CVec2<ET>(El(0) * f, El(1) * f); }
    CVec2<ET> operator* (const CMtx2x2<ET> &m) const;
    CVec2<ET> operator/ (const ET &f) const { return CVec2<ET>(El(0) / f, El(1) / f); }
    bool operator== (const CVec2<ET> &v) const { return El(0)==v[0] && El(1)==v[1]; }
    bool operator!= (const CVec2<ET> &v) const { return !(operator==(v)); }

    ET MagnitudeSq() const { return armsMagnitudeSq(El(0)) + armsMagnitudeSq(El(1)); }
    ET Magnitude() const { return sqrt(armsMagnitudeSq(El(0)) + armsMagnitudeSq(El(1))); }
    CVec2<ET> Unit() const { return *this/Magnitude(); }
    CVec3<ET> Hom() const { return CVec3<ET>(*this); }
  //  void SortAbs(int *piMapTable) const;
   // int MaxAbs() const;
    ET Arg() const { return ET(atan2(armsReScalar(El(1)), armsReScalar(El(0)))); }
    //CVec2<ET> Conj() const { return CVec2<ET>(armsConj(El(0)), armsConj(El(1))); }

    CVec2<ET> & operator= (const CVec<ET> &v);
    CVec2<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CVec2<ET> & operator= (const ET &f) { El(0) = f; El(1) = f; return *this; }
    CVec2<ET> & operator+= (const CVec2<ET> &v) { El(0) += v[0]; El(1) += v[1]; return *this; }
    CVec2<ET> & operator-= (const CVec2<ET> &v) { El(0) -= v[0]; El(1) -= v[1]; return *this; }
    CVec2<ET> & operator*= (const CMtx2x2<ET> &m);
    CVec2<ET> & operator*= (const ET &f) { El(0) *= f; El(1) *= f; return *this; }
    CVec2<ET> & operator/= (const ET &f) { El(0) /= f; El(1) /= f; return *this; }


public:
    ET x, y;
};


template <class ET>
class CMtx2x2
{
public:
    CMtx2x2() {}
    CMtx2x2(const CMtx2x2<float> &m) 
        { m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(1,0); m_m[3] = (ET)m(1,1); }
    CMtx2x2(const CMtx2x2<double> &m) 
        { m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(1,0); m_m[3] = (ET)m(1,1); }
  //  CMtx2x2(const CMtx<ET> &m) { operator=(m); }
    CMtx2x2(const ET &f) { operator=(f); }
    CMtx2x2(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CMtx2x2(ET _00, ET _01,
            ET _10, ET _11 ) 
    { m_m[0] = _00; m_m[1] = _01;
      m_m[2] = _10; m_m[3] = _11; }

    int Rows() const { return 2; }
    int Cols() const { return 2; }
    int Size() const { return Rows() * Cols(); }

    const ET *Ptr() const { return m_m; }
    ET *Ptr() { return m_m; }

    CVec2<ET> GetRow(int i) const { return CVec2<ET>((*this)[i], 2); }
    CVec2<ET> GetCol(int i) const { return CVec2<ET>(m_m[i], m_m[i+Cols()]); }
    const ET &El(int iRow, int iCol) const { return m_m[iRow * Cols() + iCol]; }
    ET &El(int iRow, int iCol) { return m_m[iRow * Cols() + iCol]; }
    const ET *operator[] (int i) const { return m_m + i * Cols(); }
    const ET &operator() (int iRow, int iCol) const { return m_m[iRow * Cols() + iCol]; }
    ET *operator[] (int i) { return m_m + i * Cols(); }
    ET &operator() (int iRow, int iCol) { return El(iRow, iCol); }

    //CMtx2x2<ET> & MakeI();
    //CMtx2x2<ET> & MakeDiag(const CVec2<ET> &v) { El(0,0) = v[0]; El(0,1) = 0; El(1,0) = 0; El(1,1) = v[1]; return *this; }
    //CMtx2x2<ET> & MakeScale(ET sx, ET sy) {return MakeDiag(CVec2<ET>(sx, sy)); }
    //CMtx2x2<ET> & MakeRotation(ET fTheta);
    //CMtx2x2<ET> & MakeSymmetric() { m_m[2] = m_m[1]; return *this; } // copies below diagonal entries from above diagonal
    //CMtx2x2<ET> Inv() const;
    //CMtx2x2<ET> T() const;
    //CMtx2x2<ET> H() const; // conjugate transpose
    //CMtx2x2<ET> Conj() const; // conjugate
    ET Tr() const { return m_m[0] + m_m[3]; }
    ET Det() const { return m_m[0] * m_m[3] - m_m[1] * m_m[2]; }
   // void MaxAbs(int &row, int &col) const; // get maximum element

    CMtx2x2<ET> & SetRow(int i, const CVec2<ET> &v) { El(i,0) = v[0]; El(i,1) = v[1]; return *this; }
    CMtx2x2<ET> & SetCol(int i, const CVec2<ET> &v) { El(0,i) = v[0]; El(1,i) = v[1]; return *this; }
    CVec2<ET> GetDiag() const { return CVec2<ET>(El(0,0), El(1,1)); }

    CMtx2x2<ET> operator+ () const { return *this; }
    CMtx2x2<ET> operator- () const;
    CMtx2x2<ET> operator+ (const CMtx2x2<ET> &m) const;
    CMtx2x2<ET> operator- (const CMtx2x2<ET> &m) const;
    CMtx2x2<ET> operator* (const CMtx2x2<ET> &m) const;
    CMtx2x2<ET> operator* (const ET &f) const;
    CVec2<ET> operator* (const CVec2<ET> &v) const;
    CMtx2x2<ET> operator/ (const ET &f) const;
    bool operator== (const CMtx2x2<ET> &m) const;
    bool operator!= (const CMtx2x2<ET> &m) const { return !(operator==(m)); }

    //CMtx2x2<ET> & operator= (const CMtx<ET> &m);
    CMtx2x2<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CMtx2x2<ET> & operator= (const ET &f);
    CMtx2x2<ET> & operator+= (const CMtx2x2<ET> &m);
    CMtx2x2<ET> & operator-= (const CMtx2x2<ET> &m);
    CMtx2x2<ET> & operator*= (const CMtx2x2<ET> &m);
    CMtx2x2<ET> & operator*= (const ET &f);
    CMtx2x2<ET> & operator/= (const ET &f);



private:
    ET m_m[4];
};

typedef CVec2<int> CVec2i;
typedef CVec2<float> CVec2f;
typedef CVec2<double> CVec2d;
typedef CMtx2x2<float> CMtx2x2f;
typedef CMtx2x2<double> CMtx2x2d;

template <class ET>
class CVec3
{
public:
    CVec3() {}
    CVec3(const CVec3<float> &v) { x=(ET)v.x; y=(ET)v.y; z=(ET)v.z; }
    CVec3(const CVec3<double> &v) { x=(ET)v.x; y=(ET)v.y; z=(ET)v.z; }
    CVec3(const CVec2<ET> &v, ET zz=(ET)1) { x=v.x; y=v.y; z=zz; }
    CVec3(const CVec<ET> &v) { operator=(v); }
    CVec3(const ET &f) { operator=(f); }
    CVec3(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CVec3(const ET &m0, const ET &m1, const ET &m2) { El(0) = m0; El(1) = m1; El(2) = m2; }

    int Size() const { return 3; }

    const ET *Ptr() const { return &x; }
    ET *Ptr() { return &x; }

    const ET &El(int i) const { return (&x)[i]; }
    ET &El(int i) { return (&x)[i]; }
    const ET &operator[] (int i) const { return (&x)[i]; }
    ET &operator[] (int i) { return (&x)[i]; }
    const ET &operator() (int i) const { return (&x)[i]; }
    ET &operator() (int i) { return (&x)[i]; }

    CVec3<ET> operator+ () const { return *this; }
    CVec3<ET> operator- () const { return CVec3<ET>(-El(0), -El(1), -El(2)); }
    CVec3<ET> operator+ (const CVec3<ET> &v) const { return CVec3<ET>(El(0) + v[0], El(1) + v[1], El(2) + v[2]); }
    CVec3<ET> operator- (const CVec3<ET> &v) const { return CVec3<ET>(El(0) - v[0], El(1) - v[1], El(2) - v[2]); }
    CMtx3x3<ET> operator^ (const CVec3<ET> &v) const;
    ET operator* (const CVec3<ET> &v) const { return v[0] * El(0) + v[1] * El(1)+ v[2] * El(2); }
    CVec3<ET> operator* (const ET &f) const { return CVec3<ET>(El(0) * f, El(1) * f, El(2) * f); }
    CVec3<ET> operator* (const CMtx3x3<ET> &m) const;
    CVec3<ET> operator/ (const ET &f) const { return CVec3<ET>(El(0) / f, El(1) / f, El(2) / f); }
    bool operator== (const CVec3<ET> &v) const { return El(0)==v[0] && El(1)==v[1] && El(2)==v[2]; }
    bool operator!= (const CVec3<ET> &v) const { return !(operator==(v)); }

    ET MagnitudeSq() const { return armsMagnitudeSq(El(0)) + armsMagnitudeSq(El(1)) + armsMagnitudeSq(El(2)); }
    ET Magnitude() const { return (ET)sqrt(armsReScalar(MagnitudeSq())); }
    CVec3<ET> Unit() const { return *this/Magnitude(); }
    CVec2<ET> Dehom() const { ET sc = (ET)1/z; return CVec2<ET>(sc*x, sc*y); }
    CVec4<ET> Hom() const { return CVec4<ET>(*this); }
    //void SortAbs(int *piMapTable) const;
    //int MaxAbs() const;
    CVec3<ET> Cross(const CVec3<ET>& v) const
    { return CVec3<ET>(y * v.z - v.y * z,
                       z * v.x - v.z * x,
                       x * v.y - v.x * y); }
    //CVec3<ET> Conj() const { return CVec3<ET>(armsConj(El(0)), armsConj(El(1)), armsConj(El(2))); }

    CVec3<ET> & operator= (const CVec<ET> &v);
    CVec3<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CVec3<ET> & operator= (const ET &f) { El(0) = f; El(1) = f; El(2) = f; return *this; }
    CVec3<ET> & operator+= (const CVec3<ET> &v) { El(0) += v[0]; El(1) += v[1]; El(2) += v[2]; return *this; }
    CVec3<ET> & operator-= (const CVec3<ET> &v) { El(0) -= v[0]; El(1) -= v[1]; El(2) -= v[2]; return *this; }
    CVec3<ET> & operator*= (const CMtx3x3<ET> &m);
    CVec3<ET> & operator*= (const ET &f) { El(0) *= f; El(1) *= f; El(2) *= f; return *this; }
    CVec3<ET> & operator/= (const ET &f) { El(0) /= f; El(1) /= f; El(2) /= f; return *this; }

	CVec3<ET> hsv();
	CVec3<ET> srgb2hsv();

	CMtx3x3<ET> CrossMatrix () ;
	CMtx4x4<ET> CrossMatrix (const CVec3<ET> v) ;





public:
    ET x, y, z;
};



template <class ET>
class CMtx3x3
{
public:
    CMtx3x3() {}
    CMtx3x3(const CMtx3x3<float> &m) 
        { 
          m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(0,2);
          m_m[3] = (ET)m(1,0); m_m[4] = (ET)m(1,1); m_m[5] = (ET)m(1,2);
          m_m[6] = (ET)m(2,0); m_m[7] = (ET)m(2,1); m_m[8] = (ET)m(2,2);
        }
    CMtx3x3(const CMtx3x3<double> &m) 
        { 
          m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(0,2);
          m_m[3] = (ET)m(1,0); m_m[4] = (ET)m(1,1); m_m[5] = (ET)m(1,2);
          m_m[6] = (ET)m(2,0); m_m[7] = (ET)m(2,1); m_m[8] = (ET)m(2,2);
        }
	CMtx3x3(const CMtx<ET> &m) { operator=(m); }
    CMtx3x3(const ET &f) { operator=(f); }
    CMtx3x3(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CMtx3x3(ET _00, ET _01, ET _02,
            ET _10, ET _11, ET _12,
            ET _20, ET _21, ET _22 ) 
    { m_m[0] = _00; m_m[1] = _01; m_m[2] = _02;
      m_m[3] = _10; m_m[4] = _11; m_m[5] = _12;
      m_m[6] = _20; m_m[7] = _21; m_m[8] = _22; }

    int Rows() const { return 3; }
    int Cols() const { return 3; }
    int Size() const { return Rows() * Cols(); }

    const ET *Ptr() const { return m_m; }
    ET *Ptr() { return m_m; }
	

    CVec3<ET> GetRow(int i) const { return CVec3<ET>((*this)[i], 3); }
    CVec3<ET> GetCol(int i) const { return CVec3<ET>(m_m[i], m_m[i+Cols()], m_m[i+2*Cols()]); }

    const ET &El(int iRow, int iCol) const { return m_m[iRow * Cols() + iCol]; }
    ET &El(int iRow, int iCol) { return m_m[iRow * Cols() + iCol]; }
    const ET *operator[] (int i) const { return m_m + i * Cols(); }
    ET *operator[] (int i) { return m_m + i * Cols(); }
    const ET &operator() (int iRow, int iCol) const { return El(iRow, iCol); }
    ET &operator() (int iRow, int iCol) { return El(iRow, iCol); }

    CMtx3x3<ET> & MakeI();
    CMtx3x3<ET> & MakeDiag(const CVec3<ET> &v) { memset(m_m, 0, sizeof(ET) * Size());
        El(0,0) = v[0]; El(1,1) = v[1]; El(2,2) = v[2]; return *this; }
    CMtx3x3<ET> & MakeScale(ET sx, ET sy, ET sz=ET(1)) 
        { return MakeDiag(CVec3<ET>(sx, sy, sz)); }
    CMtx3x3<ET> & MakeTranslate(ET tx, ET ty) 
        { *this = CMtx3x3(ET(1), 0, tx, 0, ET(1), ty, 0, 0, ET(1)); return *this; }
    CMtx3x3<ET> & MakeRotation(ET fTheta);
    CMtx3x3<ET> & MakeSimilarity(ET fScale, ET fRotation, ET fDX, ET fDY);
    CMtx3x3<ET> & MakeSymmetric()
        { m_m[3] = m_m[1]; m_m[6] = m_m[2]; m_m[7] = m_m[5]; return *this; } // copies below diagonal entries from above diagonal
    CMtx3x3<ET> Inv() const;
    CMtx3x3<ET> T() const;
    //CMtx3x3<ET> H() const; // conjugate transpose
    //CMtx3x3<ET> Conj() const; // conjugate
    ET Tr() const { return m_m[0] + m_m[4] + m_m[8]; }
    ET Det() const { return m_m[0] * (m_m[4]*m_m[8]-m_m[5]*m_m[7]) 
        - m_m[1] * (m_m[3]*m_m[8]-m_m[6]*m_m[5]) 
        + m_m[2] * (m_m[3]*m_m[7]-m_m[4]*m_m[6]); }

    CMtx2x2<ET> DeleteRowCol(int iRow, int iCol) const;
    CMtx3x3<ET> & SetRow(int i, const CVec3<ET> &v) { El(i,0) = v[0]; El(i,1) = v[1]; El(i,2) = v[2]; return *this; }
    CMtx3x3<ET> & SetCol(int i, const CVec3<ET> &v) { El(0,i) = v[0]; El(1,i) = v[1]; El(2,i) = v[2]; return *this; }
    CVec3<ET> GetDiag() const { return CVec3<ET>(El(0,0), El(1,1), El(2,2)); }
    //void MaxAbs(int &row, int &col) const; // get maximum element

    CMtx3x3<ET> operator+ () const { return *this; }
    CMtx3x3<ET> operator- () const;
    CMtx3x3<ET> operator+ (const CMtx3x3<ET> &m) const;
    CMtx3x3<ET> operator- (const CMtx3x3<ET> &m) const;
    CMtx3x3<ET> operator* (const CMtx3x3<ET> &m) const;
    CMtx3x3<ET> operator* (const ET &f) const;
    CVec3<ET> operator* (const CVec3<ET> &v) const;
    CMtx3x3<ET> operator/ (const ET &f) const;
    bool operator== (const CMtx3x3<ET> &m) const;
    bool operator!= (const CMtx3x3<ET> &m) const { return !(operator==(m)); }

    CMtx3x3<ET> & operator= (const CMtx<ET> &m);
    CMtx3x3<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CMtx3x3<ET> & operator= (const ET &f);
    CMtx3x3<ET> & operator+= (const CMtx3x3<ET> &m);
    CMtx3x3<ET> & operator-= (const CMtx3x3<ET> &m);
    CMtx3x3<ET> & operator*= (const CMtx3x3<ET> &m);
    CMtx3x3<ET> & operator*= (const ET &f);
    CMtx3x3<ET> & operator/= (const ET &f);







private:
    ET m_m[9];
};

typedef CVec3<int> CVec3i;
typedef CVec3<float> CVec3f;
typedef CVec3<double> CVec3d;
typedef CMtx3x3<float> CMtx3x3f;
typedef CMtx3x3<double> CMtx3x3d;

template <class ET>
class CVec4
{
public:
    CVec4() {}
    CVec4(const CVec4<float> &v) { x=(ET)v.x; y=(ET)v.y; z=(ET)v.z; w=(ET)v.w; }
    CVec4(const CVec4<double> &v) { x=(ET)v.x; y=(ET)v.y; z=(ET)v.z; w=(ET)v.w; }
    CVec4(const CVec3<ET> &v, ET ww=(ET)1) { x=v.x; y=v.y; z=v.z; w=ww; }
    CVec4(const CVec<ET> &v) { operator=(v); }
    CVec4(const ET &f) { operator=(f); }
    CVec4(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CVec4(const ET &m0, const ET &m1, const ET &m2, const ET &m3)
        { El(0) = m0; El(1) = m1; El(2) = m2; El(3) = m3; }

    int Size() const { return 4; }

    const ET *Ptr() const { return &x; }
    ET *Ptr() { return &x; }

    const ET &El(int i) const { return (&x)[i]; }
    ET &El(int i) { return (&x)[i]; }
    const ET &operator[] (int i) const { return (&x)[i]; }
    ET &operator[] (int i) { return (&x)[i]; }
    const ET &operator() (int i) const { return (&x)[i]; }
    ET &operator() (int i) { return (&x)[i]; }

    CVec4<ET> operator+ () const { return *this; }
    CVec4<ET> operator- () const { return CVec4<ET>(-El(0), -El(1), -El(2), -El(3)); }
    CVec4<ET> operator+ (const CVec4<ET> &v) const 
        { return CVec4<ET>(El(0) + v[0], El(1) + v[1], El(2) + v[2], El(3) + v[3]); }
    CVec4<ET> operator- (const CVec4<ET> &v) const 
        { return CVec4<ET>(El(0) - v[0], El(1) - v[1], El(2) - v[2], El(3) - v[3]); }
    CMtx4x4<ET> operator^ (const CVec4<ET> &v) const;
    ET operator* (const CVec4<ET> &v) const { return v[0] * El(0) + v[1] * El(1) + v[2] * El(2) + v[3] * El(3); }
    CVec4<ET> operator* (const ET &f) const { return CVec4<ET>(El(0) * f, El(1) * f, El(2) * f, El(3) * f); }
    CVec4<ET> operator* (const CMtx4x4<ET> &m) const;
    CVec4<ET> operator/ (const ET &f) const { return CVec4<ET>(El(0) / f, El(1) / f, El(2) / f, El(3) / f); }
    bool operator== (const CVec4<ET> &v) const { return El(0)==v[0] && El(1)==v[1] && El(2)==v[2] && El(3)==v[3]; }
    bool operator!= (const CVec4<ET> &v) const { return !(operator==(v)); }

    ET MagnitudeSq() const { return armsMagnitudeSq(El(0)) + armsMagnitudeSq(El(1)) 
        + armsMagnitudeSq(El(2)) + armsMagnitudeSq(El(3)); }
    ET Magnitude() const { return (ET)sqrt(armsReScalar(MagnitudeSq())); }
    CVec4<ET> Unit() const { return *this/Magnitude(); }
    CVec3<ET> Dehom() const { ET sc = (ET)1/w; return CVec3<ET>(sc*x, sc*y, sc*z); }
    //void SortAbs(int *piMapTable) const;
    //int MaxAbs() const;
    //CVec4<ET> Conj() const { return CVec4<ET>(armsConj(El(0)), armsConj(El(1)), armsConj(El(2)), armsConj(El(3))); }

    CVec4<ET> & operator= (const CVec<ET> &v);
    CVec4<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CVec4<ET> & operator= (const ET &f) { El(0) = f; El(1) = f; El(2) = f; El(3) = f; return *this; }
    CVec4<ET> & operator+= (const CVec4<ET> &v) { El(0) += v[0]; El(1) += v[1]; El(2) += v[2]; El(3) += v[3]; return *this; }
    CVec4<ET> & operator-= (const CVec4<ET> &v) { El(0) -= v[0]; El(1) -= v[1]; El(2) -= v[2]; El(3) -= v[3]; return *this; }
    CVec4<ET> & operator*= (const CMtx4x4<ET> &m);
    CVec4<ET> & operator*= (const ET &f) { El(0) *= f; El(1) *= f; El(2) *= f; El(3) *= f; return *this; }
    CVec4<ET> & operator/= (const ET &f) { El(0) /= f; El(1) /= f; El(2) /= f; El(3) /= f; return *this; }



public:
    ET x, y, z, w;
};

template <class ET>
class CMtx4x4
{
public:
    CMtx4x4() {}
    CMtx4x4(const CMtx4x4<float> &m) 
        { 
          m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(0,2); m_m[3] = (ET)m(0,3);
          m_m[4] = (ET)m(1,0); m_m[5] = (ET)m(1,1); m_m[6] = (ET)m(1,2); m_m[7] = (ET)m(1,3);
          m_m[8] = (ET)m(2,0); m_m[9] = (ET)m(2,1); m_m[10] = (ET)m(2,2); m_m[11] = (ET)m(2,3);
          m_m[12] = (ET)m(3,0); m_m[13] = (ET)m(3,1); m_m[14] = (ET)m(3,2); m_m[15] = (ET)m(3,3);
        }
		
	CMtx4x4(const CMtx3x3<float> &m) 
        { 
          m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(0,2); m_m[3] = (ET)0;
          m_m[4] = (ET)m(1,0); m_m[5] = (ET)m(1,1); m_m[6] = (ET)m(1,2); m_m[7] = (ET)0;
          m_m[8] = (ET)m(2,0); m_m[9] = (ET)m(2,1); m_m[10] = (ET)m(2,2); m_m[11] = (ET)0;
          m_m[12] = (ET)0; m_m[13] = (ET)0; m_m[14] = (ET)0; m_m[15] = (ET)0;
        }
    CMtx4x4(const CMtx4x4<double> &m) 
        { 
          m_m[0] = (ET)m(0,0); m_m[1] = (ET)m(0,1); m_m[2] = (ET)m(0,2); m_m[3] = (ET)m(0,3);
          m_m[4] = (ET)m(1,0); m_m[5] = (ET)m(1,1); m_m[6] = (ET)m(1,2); m_m[7] = (ET)m(1,3);
          m_m[8] = (ET)m(2,0); m_m[9] = (ET)m(2,1); m_m[10] = (ET)m(2,2); m_m[11] = (ET)m(2,3);
          m_m[12] = (ET)m(3,0); m_m[13] = (ET)m(3,1); m_m[14] = (ET)m(3,2); m_m[15] = (ET)m(3,3);
        }
	CMtx4x4(const CMtx<ET> &m) { operator=(m); }
    CMtx4x4(const ET &f) { operator=(f); }
    CMtx4x4(const ET *p, int iSize) { memcpy(Ptr(), p, sizeof(ET) * min(Size(), iSize)); }
    CMtx4x4(ET _00, ET _01, ET _02, ET _03,
            ET _10, ET _11, ET _12, ET _13,
            ET _20, ET _21, ET _22, ET _23,
            ET _30, ET _31, ET _32, ET _33 ) 
    { m_m[0]  = _00; m_m[1]  = _01; m_m[2]  = _02; m_m[3]  = _03;
      m_m[4]  = _10; m_m[5]  = _11; m_m[6]  = _12; m_m[7]  = _13;
      m_m[8]  = _20; m_m[9]  = _21; m_m[10] = _22; m_m[11] = _23;
      m_m[12] = _30; m_m[13] = _31; m_m[14] = _32; m_m[15] = _33; }

    int Rows() const { return 4; }
    int Cols() const { return 4; }
    int Size() const { return Rows() * Cols(); }

    const ET *Ptr() const { return m_m; }
    ET *Ptr() { return m_m; }

    CVec4<ET> GetRow(int i) const { return CVec4<ET>((*this)[i], 4); }
    CVec4<ET> GetCol(int i) const { return CVec4<ET>(m_m[i], m_m[i+Cols()], m_m[i+2*Cols()], m_m[i+3*Cols()]); }

    const ET &El(int iRow, int iCol) const { return m_m[iRow * Cols() + iCol]; }
    ET &El(int iRow, int iCol) { return m_m[iRow * Cols() + iCol]; }
    const ET *operator[] (int i) const { return m_m + i * Cols(); }
    ET *operator[] (int i) { return m_m + i * Cols(); }
    const ET &operator() (int iRow, int iCol) const { return El(iRow, iCol); }
    ET &operator() (int iRow, int iCol) { return El(iRow, iCol); }


	CVec3<ET> transform_point(const CVec3<ET> &v) const
	{
		CVec3<ET> vr;

		vr(0) = El(0, 0) * v(0) + El(0, 1) * v(1) + El(0, 2) * v(2) + El(0, 3);
		vr(1) = El(1, 0) * v(0) + El(1, 1) * v(1) + El(1, 2) * v(2) + El(1, 3);
		vr(2) = El(2, 0) * v(0) + El(2, 1) * v(1) + El(2, 2) * v(2) + El(2, 3);

		return vr;

	}

	CVec3<ET> transform_vector(const CVec3<ET> &v) const
	{
		CVec3<ET> vr;

		vr(0) = El(0, 0) * v(0) + El(0, 1) * v(1) + El(0, 2) * v(2);
		vr(1) = El(1, 0) * v(0) + El(1, 1) * v(1) + El(1, 2) * v(2);
		vr(2) = El(2, 0) * v(0) + El(2, 1) * v(1) + El(2, 2) * v(2);

		return vr;

	}

	CMtx3x3<ET> Submat3x3() const{
		CMtx3x3<ET> retmat;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				retmat(i, j) = (*this)(i, j);
		}
		return retmat;
	
	}

    CMtx4x4<ET> & MakeI();
    CMtx4x4<ET> & MakeDiag(const CVec4<ET> &v) { memset(m_m, 0, Size() * sizeof(ET));
        El(0,0) = v[0]; El(1,1) = v[1]; El(2,2) = v[2]; El(3,3) = v[3]; return *this; }
    CMtx4x4<ET> & MakeScale(ET sx, ET sy, ET sz, ET sw = ET(1)) 
        { return MakeDiag(CVec4<ET>(sx, sy, sz, sw)); }
    CMtx4x4<ET> & MakeTranslate(ET tx, ET ty, ET tz)
        { MakeI(); El(0, 3) = tx; El(1, 3) = ty; El(2, 3) = tz; return *this; }
    CMtx4x4<ET> & MakeRotationX(ET fTheta);
    CMtx4x4<ET> & MakeRotationY(ET fTheta);
    CMtx4x4<ET> & MakeRotationZ(ET fTheta);
    CMtx4x4<ET> & MakeSymmetric() // copies below diagonal entries from above diagonal
        { m_m[4] = m_m[1]; m_m[8] = m_m[2]; m_m[9] = m_m[6];
          m_m[12] = m_m[3]; m_m[13] = m_m[7]; m_m[14] = m_m[11]; return *this; } 
    CMtx4x4<ET> Inv() const;
    CMtx4x4<ET> T() const;
    //CMtx4x4<ET> H() const; // conjugate transpose
    //CMtx4x4<ET> Conj() const; // conjugate
    ET Tr() const { return m_m[0] + m_m[5] + m_m[10] + m_m[15]; }
    //ET Det() const;

   // CMtx3x3<ET> DeleteRowCol(int iRow, int iCol) const;
    CMtx4x4<ET> & SetRow(int i, const CVec4<ET> &v) { El(i,0) = v[0]; El(i,1) = v[1]; El(i,2) = v[2]; El(i,3) = v[3]; return *this; }
    CMtx4x4<ET> & SetCol(int i, const CVec4<ET> &v) { El(0,i) = v[0]; El(1,i) = v[1]; El(2,i) = v[2]; El(3,i) = v[3]; return *this; }
    CVec4<ET> GetDiag() const { return CVec4<ET>(El(0,0), El(1,1), El(2,2), El(3,3)); }
   // void MaxAbs(int &row, int &col) const; // get maximum element

    CMtx4x4<ET> operator+ () const { return *this; }
    CMtx4x4<ET> operator- () const;
    CMtx4x4<ET> operator+ (const CMtx4x4<ET> &m) const;
    CMtx4x4<ET> operator- (const CMtx4x4<ET> &m) const;
    CMtx4x4<ET> operator* (const CMtx4x4<ET> &m) const;
    CMtx4x4<ET> operator* (const ET &f) const;
    CVec4<ET> operator* (const CVec4<ET> &v) const;
    CMtx4x4<ET> operator/ (const ET &f) const;
    bool operator== (const CMtx4x4<ET> &m) const;
    bool operator!= (const CMtx4x4<ET> &m) const { return !(operator==(m)); }

    CMtx4x4<ET> & operator= (const CMtx<ET> &m);
    CMtx4x4<ET> & Init(const ET *p, int iSize); // iSize is the number of elements in p
    CMtx4x4<ET> & operator= (const ET &f);
    CMtx4x4<ET> & operator+= (const CMtx4x4<ET> &m);
    CMtx4x4<ET> & operator-= (const CMtx4x4<ET> &m);
    CMtx4x4<ET> & operator*= (const CMtx4x4<ET> &m);
    CMtx4x4<ET> & operator*= (const ET &f);
    CMtx4x4<ET> & operator/= (const ET &f);



protected:
    ET m_m[16];
};


typedef CVec4<int>  CVec4i;
typedef CVec4<float>  CVec4f;
typedef CVec4<double> CVec4d;
typedef CMtx4x4<float>  CMtx4x4f;
typedef CMtx4x4<double> CMtx4x4d;









#endif