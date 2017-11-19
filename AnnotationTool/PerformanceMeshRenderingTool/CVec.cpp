#include "CVec.h"



#include <cstdio>
#include <iostream>
#include <string.h>

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::~CVec()
//
//  Synopsis:  destructor
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET>::~CVec()
{
    if(m_p && !m_bWrap)
        delete m_p;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::Wrap()
//
//  Synopsis:  wrap around some memory
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
void CVec<ET>::Wrap(ET *p, int iSize)
{
    //ClearError();
    Free();
    m_p = p;
    m_iSize = iSize;
    m_bWrap = true;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::Create()
//
//  Synopsis:  create a vector of the specified size
//
//  Returns:   S_OK if no problems
//
//------------------------------------------------------------------------

template <class ET>
HRESULT CVec<ET>::Create(int iSize)
{
   // ClearError();

    if(m_bWrap)
        m_p = NULL;

    m_bWrap = false;

    if(iSize<=0)
    {
        Free();
        return 0;
    }

    if(m_p && iSize==m_iSize)
        return 0;

    if(m_p)
        delete m_p;

    m_p = new ET [iSize];
    if(m_p)
    {
        m_iSize = iSize;
        return 0;
    }

    m_iSize = 0;
    
    return 1;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::Free
//
//  Synopsis:  deallocates the memory and sets size to 0x0
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
void CVec<ET>::Free()
{
    m_iSize = 0;
    if(m_p && !m_bWrap)
        delete m_p;
    m_p = NULL;
    m_bWrap = false;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator=
//
//  Synopsis:  assigment operator from arbitrary vector types
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator= (const CVec<ET> &v)
{
    // if v has out of memory error then its size will be zero
    if(m_bWrap && Size()==v.Size() )
    {
        //SetError(v.GetError()); // pass on the error status

        //if(!IsError() && v.m_p && m_p)
            memcpy(m_p, v.m_p, v.Size() * sizeof(ET));
    }

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator=
//
//  Synopsis:  assigment operator from CVec2
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator= (const CVec2<ET> &v)
{
    if(m_bWrap && Size()==2 || SUCCEEDED(Create(2)))
        memcpy(m_p, v.Ptr(), 2 * sizeof(ET));

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator=
//
//  Synopsis:  assigment operator from CVec3
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator= (const CVec3<ET> &v)
{
    if(m_bWrap && Size()==3 || SUCCEEDED(Create(3)))
        memcpy(m_p, v.Ptr(), 3 * sizeof(ET));

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator=
//
//  Synopsis:  assigment operator from CVec4
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator= (const CVec4<ET> &v)
{
    if(m_bWrap && Size()==4 || SUCCEEDED(Create(4)))
        memcpy(m_p, v.Ptr(), 4 * sizeof(ET));

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator-
//
//  Synopsis:  unary minus operator
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator- () const
{
    CVec<ET> v(Size());
    //v.SetError(GetError());
    //if(!v.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            v(i) = -El(i);
    }
    return v;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator+
//
//  Synopsis:  vector addition
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator+ (const CVec<ET> &v) const
{
    CVec<ET> q(Size());
    //q.SetError(GetError());
    //q.SetError(v.GetError());
    //if(!q.IsError())
    {
        if(v.Size()!=Size())
        {
            q = (ET)0;
            //q.SetError(E_INVALIDARG);
            return q;
        }
        int i;
        for(i=0; i<Size(); i++)
            q(i) = El(i) + v[i];
    }
    return q;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator-
//
//  Synopsis:  vector substraction
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator- (const CVec<ET> &v) const
{
    CVec<ET> q(Size());
    //q.SetError(GetError());
    //q.SetError(v.GetError());
    //if(!q.IsError())
    {
        if(v.Size()!=Size())
        {
            q = (ET)0;
            //q.SetError(E_INVALIDARG);
            return q;
        }
        int i;
        for(i=0; i<Size(); i++)
            q(i) = El(i) - v[i];
    }
    return q;
}


template <class ET>
CMtx<ET> CVec<ET>::operator^ (const CVec<ET> &v) const
{
    CMtx<ET> m(Size(), v.Size());
    //m.SetError(GetError());
    //m.SetError(v.GetError());
    //if(m.IsError())
    //    return m;
    int i, j;
    for(i=0; i<Size(); i++)
        for(j=0; j<v.Size(); j++)
            m(i,j) = El(i) * v(j);
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator*
//
//  Synopsis:  dot product
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
ET CVec<ET>::operator* (const CVec<ET> &v) const
{
    if(v.Size()!=Size())
    {
        //SetError(E_INVALIDARG);
        return 0;
    }
    ET f = 0;
    int i;
    for(i=0; i<Size(); i++)
        f += El(i) *  v[i];

    return f;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator*
//
//  Synopsis:  rescale the vector
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator* (const ET &f) const
{
    CVec<ET> v(Size());
    //v.SetError(GetError());
    //if(!v.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            v(i) = El(i)*f;
    }
    return v;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator*
//
//  Synopsis:  right sided matrix multiply
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator* (const CMtx<ET> &m) const
{
    CVec<ET> v(m.Cols());
    //v.SetError(m.GetError());
    //v.SetError(GetError());

    //if(!v.IsError())
    {
        if(m.Rows()!=Size())
        {
            //v.SetError(E_INVALIDARG);
            v = (ET)0;
            return v;
        }

        int i, j;
        for(j=0; j<m.Cols(); j++)
        {
            ET f = 0;
            for(i=0; i<m.Rows(); i++)
                f += El(i) * m(i, j);
            v[j] = f;
        }
    }
    return v;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator/
//
//  Synopsis:  divide by a scalar
//
//  Returns:   the new CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CVec<ET>::operator/ (const ET &f) const
{
    CVec<ET> v(Size());
    //v.SetError(GetError());

    //if(!v.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            v(i) = El(i)/f;
    }
    return v;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::MagnitudeSq
//
//  Synopsis:  calculate the squared magnitude of the vector
//
//  Returns:   magnitude
//
//------------------------------------------------------------------------

template <class ET>
ET CVec<ET>::MagnitudeSq() const
{
    ET fMagSq = 0;
    int i;
    for(i=0; i<Size(); i++)
    {
        fMagSq += armsMagnitudeSq(El(i));
    }

    return fMagSq;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::Init
//
//  Synopsis:  initialize from an array pointer
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::Init(const ET *p, int iSize)
{
    if(p && Ptr() && iSize>0)
    {
        iSize = min(Size(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator=
//
//  Synopsis:  assignment from a scalar
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        El(i) = f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator+=
//
//  Synopsis:  add two vectors
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator+= (const CVec<ET> &v)
{ 
    //SetError(v.GetError());

    //if(!IsError())
    {
        int i;
        int iSize = min(Size(), v.Size());
        for(i=0; i<iSize; i++)
            El(i) += v[i];
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator-=
//
//  Synopsis:  subtract a vector
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator-= (const CVec<ET> &v)
{ 
    //SetError(v.GetError());

    //if(!IsError())
    {
        int i;
        int iSize = min(Size(), v.Size());
        for(i=0; i<iSize; i++)
            El(i) -= v[i];
    }

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator*=
//
//  Synopsis:  multiply by a scalar
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator*= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        El(i) *= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec<ET>::operator/=
//
//  Synopsis:  divide by a scalar
//
//  Returns:   the current CVec
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> & CVec<ET>::operator/= (const ET &f)
{ 
    int i;
    for(i=0; i<Size(); i++)
        El(i) /= f;
    return *this;
}




//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::Create
//
//  Synopsis:  create a matrix of the specified size
//
//  Returns:   S_OK if nothing goes wrong
//
//------------------------------------------------------------------------

template <class ET>
HRESULT CMtx<ET>::Create(int iRows, int iCols)
{
    //ClearError();

    if(m_bWrap)
        m_p = NULL;
    m_bWrap = false;

    if(iRows<=0 || iCols<=0)
    {
        Free();
        return 0;
    }

    if(m_p && iRows * iCols==Size())
    {
        m_iRows = iRows;
        m_iCols = iCols;
        return 0;
    }

    if(m_p)
        delete m_p;

    m_p = new ET [iRows * iCols];
    if(m_p)
    {
        m_iRows = iRows;
        m_iCols = iCols;
        return 0;
    }

    m_iRows = m_iCols = 0;
    
    return 1;//SetError(E_OUTOFMEMORY);
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::Free
//
//  Synopsis:  deallocates the memory and sets size to 0x0
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
void CMtx<ET>::Free()
{
    m_iRows = m_iCols = 0;
    if(m_p && !m_bWrap)
        delete m_p;
    m_bWrap = false;
    m_p = NULL;
}


//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::~CMtx()
//
//  Synopsis:  destructor
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET>::~CMtx()
{
    if(m_p && !m_bWrap)
        delete m_p;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::Wrap
//
//  Synopsis:  wraps existing memory
//
//  Returns:   none
//
//------------------------------------------------------------------------

template <class ET>
void CMtx<ET>::Wrap(ET *p, int iRows, int iCols)
{
    //ClearError();
    Free();
    m_p = p;
    m_iRows = iRows;
    m_iCols = iCols;
    m_bWrap = true;
}


//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator=
//
//  Synopsis:  assignment of matrix, creates new one with correct size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::operator= (const CMtx<ET> &m)
{
    // if m has out of memory then its size will be zero
    if(m_bWrap && Rows()==m.Rows() && Cols()==m.Cols() || SUCCEEDED(Create(m.Rows(), m.Cols())))
    {
        //SetError(m.GetError()); // pass on the error status

        //if(!IsError() && m.m_p && m_p)
            memcpy(m_p, m.m_p, m.Size() * sizeof(ET));
    }

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator=
//
//  Synopsis:  assignment of matrix, creates new one with correct size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::operator= (const CMtx2x2<ET> &m)
{
    if(m_bWrap && Rows()==2 && Cols()==2 || SUCCEEDED(Create(2,2)))
        memcpy(m_p, m.Ptr(), 4 * sizeof(ET));

    return *this;   
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator=
//
//  Synopsis:  assignment of matrix, creates new one with correct size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::operator= (const CMtx3x3<ET> &m)
{
    if(m_bWrap && Rows()==3 && Cols()==3 || SUCCEEDED(Create(3,3)))
        memcpy(m_p, m.Ptr(), 9 * sizeof(ET));

    return *this;   
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator=
//
//  Synopsis:  assignment of matrix, creates new one with correct size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::operator= (const CMtx4x4<ET> &m)
{
    if(m_bWrap && Rows()==4 && Cols()==4 || SUCCEEDED(Create(4,4)))
        memcpy(m_p, m.Ptr(), 16 * sizeof(ET));

    return *this;   
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::GetCols
//
//  Synopsis:  get a specified set of columns
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> CMtx<ET>::GetCols(const vector<int> &vCols) const
{
    CMtx<ET> m(Rows(), (int)vCols.size());
    //m.SetError(GetError());
    //if(!m.IsError())
    {
        int i,j;
        for(j=0; j<(int)vCols.size(); j++)
            for(i=0; i<Rows(); i++)
                m(i,j) = El(i, vCols[j]);
    }

    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::GetRows
//
//  Synopsis:  get a specified set of rows
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> CMtx<ET>::GetRows(const vector<int> &vRows) const
{
    CMtx<ET> m((int)vRows.size(), Cols());
    //m.SetError(GetError());
    //if(!m.IsError())
    {
        int i,j;
        for(j=0; j<Cols(); j++)
            for(i=0; i<(int)vRows.size(); i++)
                m(i,j) = El(vRows[i], j);
    }

    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::SetCols
//
//  Synopsis:  set a specified set of columns
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::SetCols(const vector<int> &vCols, const CMtx<ET> &m)
{
    //SetError(m.GetError());
    //if(!m.IsError())
    {
        int i,j;
        int iCols = min((int)vCols.size(), m.Cols());
        int iRows = min(Rows(), m.Rows());
        for(j=0; j<iCols; j++)
        {
            int k = vCols[j];
            if(k>=0 && k<Cols())
                for(i=0; i<iRows; i++)
                    El(i, k) = m(i, j);
        }
    }

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::SetRows
//
//  Synopsis:  set a specified set of rows
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::SetRows(const vector<int> &vRows, const CMtx<ET> &m)
{
    //SetError(m.GetError());
    //if(!m.IsError())
    {
        int i,j;
        int iRows = min((int)vRows.size(), m.Rows());
        int iCols = min(Cols(), m.Cols());
        for(i=0; i<iRows; i++)
        {
            int k = vRows[i];
            if(k>=0 && k<Rows())
                for(j=0; j<iCols; j++)
                    El(k, j) = m(i, j);
        }
    }

    return *this;
}



template <class ET>
CVec3<ET> CVec3<ET>::hsv ()
{
	float h = this->x;
	float s = this->y;
	float v = this->z;
	// From FvD
	if (s <= 0.0f)
		return CVec3<ET>(v,v,v);
	h = fmod(h, float(2.0f * M_PI));
	if (h < 0.0)
		h += (float)(2.0 * M_PI);
	h /= (float)(M_PI / 3.0);
	int i = int(floor(h));
	float f = h - i;
	float p = v * (1.0f - s);
	float q = v * (1.0f - (s*f));
	float t = v * (1.0f - (s*(1.0f-f)));
	switch(i) {
	case 0: return CVec3<ET>(v, t, p);
	case 1: return CVec3<ET>(q, v, p);
	case 2: return CVec3<ET>(p, v, t);
	case 3: return CVec3<ET>(p, q, v);
	case 4: return CVec3<ET>(t, p, v);
	default: return CVec3<ET>(v, p, q);
	}
}

template <class ET>
CVec3<ET> CVec3<ET>::srgb2hsv ()
{
	CVec3<ET> &v =*this;
	float V = max(max(v[0], v[1]), v[2]);
	if(V==0)
		return  CVec3<ET>(0, 0, 0);
	float diff = V - min(min(v[0], v[1]), v[2]);
	float S = diff / V;
	float H = 0.0f;
	if (S == 0.0f)
		return CVec3<ET>(H, S, V);
	if (V == v[0])
		H = (v[1] - v[2]) / diff;
	else if (V == v[1])
		H = (v[2] - v[0]) / diff + 2.0f;
	else
		H = (v[0] - v[1]) / diff + 4.0f;
	H *= float(M_PI / 3.0);
	if (H < 0.0f)
		H += float(2.0 * M_PI);

	H /= float(2.0 * M_PI);
	return CVec3<ET>(H, S, V);
}

template <class ET>
CMtx3x3<ET> CVec3<ET>::CrossMatrix () 
{

	CMtx3x3<ET> mat;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			mat(i,j) = 0;
		}
	}

	mat(0,1) = - z;
	mat(1,0) =  z;

	mat(0,2) =  y;
	mat(2,0) = - y;

	mat(1,2) = - x;
	mat(2,1) =  x;

	return mat;


}

template <class ET>
CMtx4x4<ET> CVec3<ET>::CrossMatrix (const CVec3<ET> v) 
{

	CMtx4x4<ET> mat;
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			mat(i,j) = 0;
		}
	}

	mat(0,1) = - z;
	mat(1,0) =  z;

	mat(0,2) =  y;
	mat(2,0) = - y;

	mat(1,2) = - x;
	mat(2,1) =  x;

	mat(0,3) = v[0];
	mat(1,3) = v[1];
	mat(2,3) = v[2];

	return mat;

}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator=
//
//  Synopsis:  assignment from a scalar
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::operator= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        Ptr()[i] = f;
    return *this;
}


//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator-
//
//  Synopsis:  unary minus operator
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator- () const
{
    CMtx<ET> m(Rows(),Cols());

    //m.SetError(GetError());
    //if(!m.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            m.Ptr()[i] = -Ptr()[i];
    }
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator+
//
//  Synopsis:  addition of two matrices
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator+ (const CMtx<ET> &m) const
{
    CMtx<ET> mnew(Rows(),Cols());
    //mnew.SetError(GetError());
    //mnew.SetError(m.GetError());

    //if(!mnew.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            mnew.Ptr()[i] = Ptr()[i];
        int iMinR = min(m.Rows(), Rows());
        int iMinC = min(m.Cols(), Cols());
        int j;
        for(i=0; i<iMinR; i++)
            for(j=0; j<iMinC; j++)
                mnew.El(i,j) += m.El(i,j);
    }
    return mnew;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator-
//
//  Synopsis:  matrix subtraction
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator- (const CMtx<ET> &m) const
{
    CMtx<ET> mnew(Rows(),Cols());

    //mnew.SetError(GetError());
    //mnew.SetError(m.GetError());
    //if(!mnew.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            mnew.Ptr()[i] = Ptr()[i];
        int iMinR = min(m.Rows(), Rows());
        int iMinC = min(m.Cols(), Cols());
        int j;
        for(i=0; i<iMinR; i++)
            for(j=0; j<iMinC; j++)
                mnew.El(i,j) -= m.El(i,j);
    }
    return mnew;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator*
//
//  Synopsis:  matrix multiplication
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator* (const CMtx<ET> &m) const
{
    CMtx<ET> mnew(Rows(),m.Cols());
    //mnew.SetError(GetError());
    //mnew.SetError(m.GetError());
    //if(!mnew.IsError())
    {
        if(Cols()!=m.Rows())
        {
            mnew = (ET)0;
            //mnew.SetError(E_INVALIDARG);
            return mnew;
        }

        int i,j,k;
        for(i=0; i<mnew.Rows(); i++)
            for(j=0; j<mnew.Cols(); j++)
            {
                ET f = (ET)0;
                for(k=0; k<Cols(); k++)
                    f += El(i,k) * m.El(k,j);
                mnew.El(i,j) = f;
            }
    }
    return mnew;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator*
//
//  Synopsis:  mutliply by a scalar
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator* (const ET &f) const
{
    CMtx<ET> m(Rows(),Cols());
    //m.SetError(GetError());
    //if(!m.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            m.Ptr()[i] = Ptr()[i] * f;
    }
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator*
//
//  Synopsis:  multiply by a vector (v2 = M * v1)
//
//  Returns:   the new vector
//
//------------------------------------------------------------------------

template <class ET> 
CVec<ET> CMtx<ET>::operator* (const CVec<ET> &v) const
{
    CVec<ET> vnew(Rows());

    //vnew.SetError(GetError());
    //vnew.SetError(v.GetError());
    //if(!vnew.IsError())
    {
        if(v.Size()!=Cols())
        {
            vnew = (ET)0;
            //vnew.SetError(E_INVALIDARG);
            return vnew;
        }
        int i, j;
        for(i=0; i<Rows(); i++)
        {
            ET f = (ET)0;
            for(j=0; j<Cols(); j++)
                f += El(i,j) * v[j];
            vnew[i] = f;
        }
    }
    return vnew;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator/
//
//  Synopsis:  divide by a scalar
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> CMtx<ET>::operator/ (const ET &f) const
{
    CMtx<ET> m(Rows(),Cols());
    //m.SetError(GetError());
    //if(!m.IsError())
    {
        int i;
        for(i=0; i<Size(); i++)
            m.Ptr()[i] = Ptr()[i] / f;
    }
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator==
//
//  Synopsis:  matrix comparison
//
//  Returns:   true or false
//
//------------------------------------------------------------------------

template <class ET> 
bool CMtx<ET>::operator== (const CMtx<ET> &m) const
{
    //if(IsError() || m.IsError() || m.Rows()!=Rows() || m.Cols()!=Cols())
    //    return false;
    int i;
    for(i=0; i<Size(); i++)
        if(m.Ptr()[i] != Ptr()[i])
            return false;
    return true;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator+=
//
//  Synopsis:  matrix addition
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> & CMtx<ET>::operator+= (const CMtx<ET> &m)
{
    //SetError(m.GetError());

    //if(!IsError())
    {
        int iMinR = min(m.Rows(), Rows());
        int iMinC = min(m.Cols(), Cols());
        int i,j;
        for(i=0; i<iMinR; i++)
            for(j=0; j<iMinC; j++)
                El(i,j) += m.El(i,j);
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator-=
//
//  Synopsis:  matrix subtraction
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> & CMtx<ET>::operator-= (const CMtx<ET> &m)
{
    //SetError(m.GetError());

    //if(!IsError())
    {  
        int iMinR = min(m.Rows(), Rows());
        int iMinC = min(m.Cols(), Cols());
        int i,j;
        for(i=0; i<iMinR; i++)
            for(j=0; j<iMinC; j++)
                El(i,j) -= m.El(i,j);
    }

    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator*=
//
//  Synopsis:  multiplication by a scalar
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> & CMtx<ET>::operator*= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        Ptr()[i] *= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::operator/=
//
//  Synopsis:  divide by a scalar
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET> 
CMtx<ET> & CMtx<ET>::operator/= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        Ptr()[i] /= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::Init
//
//  Synopsis:  initialize from an array pointer
//
//  Returns:   the current CMtx
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> & CMtx<ET>::Init(const ET *p, int iSize)
{
    if(p && Ptr() && iSize>0)
    {
        iSize = min(Rows() * Cols(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec2<ET>::operator^
//
//  Synopsis:  outer product operator
//
//  Returns:   the outer product matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx<ET> CMtx<ET>::T() const
{
    CMtx<ET> m(Cols(), Rows());
    //m.SetError(GetError());

    //if(!m.IsError())
    {
        int i, j;
        for(i=0; i<Rows(); i++)
            for(j=0; j<Cols(); j++)
                m.El(j,i) = El(i,j);
    }
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::GetRow
//
//  Synopsis:  makes a vector out of the specified row
//
//  Returns:   the new vector
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CMtx<ET>::GetRow(int i) const
{
   CVec<ET> v(Cols());
   //v.SetError(GetError());

   // if(!v.IsError())
    {
		 
        if(i<0 || i>=Rows())
		{
            v = (ET)0;
		}
        else
		{
			memcpy(v.Ptr(), Ptr() + i * Cols(), Cols() * sizeof(ET));
		}
	}
	
    return v;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx<ET>::GetCol
//
//  Synopsis:  makes a vector out of the specified column
//
//  Returns:   the new vector
//
//------------------------------------------------------------------------

template <class ET>
CVec<ET> CMtx<ET>::GetCol(int j) const
{
    CVec<ET> v(Rows());
    //v.SetError(GetError());

    //if(!v.IsError())
    {
        if(j<0 || j>=Cols())
            v = (ET)0;
        else
        {
            int i;
            for(i=0; i<Rows(); i++)
                v[i] = El(i, j);
        }
    }

    return v;
}


template <class ET>
inline CMtx2x2<ET> CVec2<ET>::operator^ (const CVec2<ET> &v) const
{
    CMtx2x2<ET> m;
    m(0,0) = El(0) * v[0];
    m(0,1) = El(0) * v[1];
    m(1,0) = El(1) * v[0];
    m(1,1) = El(1) * v[1];
    return m;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec2<ET>::operator*
//
//  Synopsis:  right hand size 2x2 matrix multiply
//
//  Returns:   the new CVec2
//
//------------------------------------------------------------------------

template <class ET>
inline CVec2<ET> CVec2<ET>::operator* (const CMtx2x2<ET> &m) const
{
    return CVec2<ET>(El(0) * m(0,0) + El(1) * m(1,0), 
                     El(0) * m(0,1) + El(1) * m(1,1));
}

template <class ET>
CVec2<ET> & CVec2<ET>::operator= (const CVec<ET> &v)
{
    if(v.Size()>0)
        El(0) = v[0];
    if(v.Size()>1)
        El(1) = v[1];
    return *this;
}


//+-----------------------------------------------------------------------
//
//  Member:    CVec2<ET>::operator*=
//
//  Synopsis:  right hand size 2x2 matrix multiply
//
//  Returns:   the current CVec2
//
//------------------------------------------------------------------------

template <class ET>
inline CVec2<ET> & CVec2<ET>::operator*= (const CMtx2x2<ET> &m)
{
    *this = *this * m;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator-
//
//  Synopsis:  unary minus
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> CMtx2x2<ET>::operator- () const
{
    CMtx2x2 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = -m_m[i];
    return mr;
}



//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator+
//
//  Synopsis:  addition
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> CMtx2x2<ET>::operator+ (const CMtx2x2<ET> &m) const
{
    CMtx2x2 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] + m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator-
//
//  Synopsis:  subtraction
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> CMtx2x2<ET>::operator- (const CMtx2x2<ET> &m) const
{
    CMtx2x2 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] - m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator*
//
//  Synopsis:  scaling by a constant
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> CMtx2x2<ET>::operator*(const ET &f) const
{
    CMtx2x2<ET> mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] * f;
    return mr;
}



//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator=
//
//  Synopsis:  set all elements to a constant
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator=(const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] = f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::Init
//
//  Synopsis:  initialize from an array pointer
//
//  Returns:   the current CMtx2x2
//
//------------------------------------------------------------------------

template <class ET>
CMtx2x2<ET> & CMtx2x2<ET>::Init(const ET *p, int iSize)
{
    if(p && iSize>0)
    {
        iSize = min(Rows() * Cols(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator+=
//
//  Synopsis:  add a matrix to the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator+= (const CMtx2x2<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] += m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator-=
//
//  Synopsis:  subtract a matrix to the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator-= (const CMtx2x2<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] -= m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator+=
//
//  Synopsis:  multiply two matrices RHS multiply by RHS matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator*= (const CMtx2x2<ET> &m)
{
    *this = *this * m;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator*=
//
//  Synopsis:  scale the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator*= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] *= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator/=
//
//  Synopsis:  divide the current matrix by a scale
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> & CMtx2x2<ET>::operator/= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] /= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator*
//
//  Synopsis:  matrix multiply
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx2x2<ET> CMtx2x2<ET>::operator*(const CMtx2x2<ET> &m) const
{
    CMtx2x2<ET> mr;
    mr(0,0) = El(0,0) * m(0,0) + El(0,1) * m(1,0);
    mr(0,1) = El(0,0) * m(0,1) + El(0,1) * m(1,1);
    mr(1,0) = El(1,0) * m(0,0) + El(1,1) * m(1,0);
    mr(1,1) = El(1,0) * m(0,1) + El(1,1) * m(1,1);
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx2x2<ET>::operator*
//
//  Synopsis:  multiply matrix by a 2x1 vector
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CVec2<ET> CMtx2x2<ET>::operator*(const CVec2<ET> &v) const
{
    CVec2<ET> vr;
    vr(0) = El(0,0) * v(0) + El(0,1) * v(1);
    vr(1) = El(1,0) * v(0) + El(1,1) * v(1);
    return vr;
}



template <class ET>
inline CMtx3x3<ET> CVec3<ET>::operator^ (const CVec3<ET> &v) const
{
    CMtx3x3<ET> m;
    m(0,0) = El(0) * v[0];
    m(0,1) = El(0) * v[1];
    m(0,2) = El(0) * v[2];
    m(1,0) = El(1) * v[0];
    m(1,1) = El(1) * v[1];
    m(1,2) = El(1) * v[2];
    m(2,0) = El(2) * v[0];
    m(2,1) = El(2) * v[1];
    m(2,2) = El(2) * v[2];
    return m;
}

template <class ET>
CVec3<ET> & CVec3<ET>::operator= (const CVec<ET> &v)
{
    if(v.Size()>0)
        El(0) = v[0];
    if(v.Size()>1)
        El(1) = v[1];
    if(v.Size()>2)
        El(2) = v[2];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec3<ET>::operator*
//
//  Synopsis:  right hand size 3x3 matrix multiply
//
//  Returns:   the new CVec3
//
//------------------------------------------------------------------------

template <class ET>
inline CVec3<ET> CVec3<ET>::operator* (const CMtx3x3<ET> &m) const
{
    return CVec3<ET>(El(0) * m(0,0) + El(1) * m(1,0) + El(2) * m(2,0),
        El(0) * m(0,1) + El(1) * m(1,1) + El(2) * m(2,1),
        El(0) * m(0,2) + El(1) * m(1,2) + El(2) * m(2,2));
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec3<ET>::operator*=
//
//  Synopsis:  right hand size 3x3 matrix multiply
//
//  Returns:   the current CVec3
//
//------------------------------------------------------------------------

template <class ET>
inline CVec3<ET> & CVec3<ET>::operator*= (const CMtx3x3<ET> &m)
{
    *this = *this * m;
    return *this;
}



template <class ET>
CMtx2x2<ET> CMtx3x3<ET>::DeleteRowCol(int iRow, int iCol) const
{
    CMtx2x2<ET> mr;
    if(iRow<0 || iRow>=Rows() || iCol<0 || iCol>=Cols())
        return mr=(ET)0;

    int i, j;
    int r, c;
    for(r=i=0; i<Rows(); i++)
        if(i!=iRow)
        {
            for(c=j=0; j<Cols(); j++)
                if(j!=iCol)
                    mr(r,c++) = El(i,j);
            r++;
        }
    return mr;
}




//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator==
//
//  Synopsis:  see if two matrices are equal
//
//  Returns:   true or false
//
//------------------------------------------------------------------------

template <class ET>
inline bool CMtx3x3<ET>::operator== (const CMtx3x3<ET> &m) const
{
    int i;
    for(i=0; i<9; i++)
        if(m_m[i]!=m.m_m[i])
            break;
    return i==9;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator=
//
//  Synopsis:  general assigment from another matrix of arbitrary size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx3x3<ET> & CMtx3x3<ET>::operator= (const CMtx<ET> &m)
{
    int iMinR = min(m.Rows(), Rows());
    int iMinC = min(m.Cols(), Cols());
    int i, j;
    for(i=0; i<iMinR; i++)
        for(j=0; j<iMinC; j++)
            El(i,j) = m.El(i,j);
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator/
//
//  Synopsis:  division by a constant
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::operator/ (const ET &f)  const
{
    CMtx3x3 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] / f;
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator-
//
//  Synopsis:  unary minus
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::operator- () const
{
    CMtx3x3 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = -m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator+
//
//  Synopsis:  addition
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::operator+ (const CMtx3x3<ET> &m) const
{
    CMtx3x3 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] + m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator-
//
//  Synopsis:  subtraction
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::operator- (const CMtx3x3<ET> &m) const
{
    CMtx3x3 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] - m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator*
//
//  Synopsis:  scaling by a constant
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::operator*(const ET &f) const
{
    CMtx3x3<ET> mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] * f;
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator=
//
//  Synopsis:  set all elements to a constant
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator=(const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] = f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::Init
//
//  Synopsis:  initialize from an array pointer
//
//  Returns:   the current CMtx3x3
//
//------------------------------------------------------------------------

template <class ET>
CMtx3x3<ET> & CMtx3x3<ET>::Init(const ET *p, int iSize)
{
    if(p && iSize>0)
    {
        iSize = min(Rows() * Cols(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator+=
//
//  Synopsis:  add a matrix to the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator+= (const CMtx3x3<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] += m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator-=
//
//  Synopsis:  subtract a matrix from the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator-= (const CMtx3x3<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] -= m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator*=
//
//  Synopsis:  multiply two matrices RHS multiply by RHS matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator*= (const CMtx3x3<ET> &m)
{
    *this = *this * m;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator*=
//
//  Synopsis:  scale the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator*= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] *= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator/=
//
//  Synopsis:  divide the current matrix by a scale
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::operator/= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] /= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator*
//
//  Synopsis:  matrix multiply
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx3x3<ET> CMtx3x3<ET>::operator*(const CMtx3x3<ET> &m) const
{
    CMtx3x3<ET> mr;
    mr(0,0) = El(0,0) * m(0,0) + El(0,1) * m(1,0) + El(0,2) * m(2,0);
    mr(0,1) = El(0,0) * m(0,1) + El(0,1) * m(1,1) + El(0,2) * m(2,1);
    mr(0,2) = El(0,0) * m(0,2) + El(0,1) * m(1,2) + El(0,2) * m(2,2);

    mr(1,0) = El(1,0) * m(0,0) + El(1,1) * m(1,0) + El(1,2) * m(2,0);
    mr(1,1) = El(1,0) * m(0,1) + El(1,1) * m(1,1) + El(1,2) * m(2,1);
    mr(1,2) = El(1,0) * m(0,2) + El(1,1) * m(1,2) + El(1,2) * m(2,2);

    mr(2,0) = El(2,0) * m(0,0) + El(2,1) * m(1,0) + El(2,2) * m(2,0);
    mr(2,1) = El(2,0) * m(0,1) + El(2,1) * m(1,1) + El(2,2) * m(2,1);
    mr(2,2) = El(2,0) * m(0,2) + El(2,1) * m(1,2) + El(2,2) * m(2,2);
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::operator*
//
//  Synopsis:  multiply matrix by a 3x1 vector
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CVec3<ET> CMtx3x3<ET>::operator*(const CVec3<ET> &v) const
{
    CVec3<ET> vr;
    vr(0) = El(0,0) * v(0) + El(0,1) * v(1) + El(0,2) * v(2);
    vr(1) = El(1,0) * v(0) + El(1,1) * v(1) + El(1,2) * v(2);
    vr(2) = El(2,0) * v(0) + El(2,1) * v(1) + El(2,2) * v(2);
    return vr;
}


template <class ET>
inline CMtx3x3<ET> CMtx3x3<ET>::T() const
{
    CMtx3x3<ET> mr;
    mr(0,0) = El(0,0);
    mr(0,1) = El(1,0);
    mr(0,2) = El(2,0);
    mr(1,0) = El(0,1);
    mr(1,1) = El(1,1);
    mr(1,2) = El(2,1);
    mr(2,0) = El(0,2);
    mr(2,1) = El(1,2);
    mr(2,2) = El(2,2);
    return mr;
}

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::T() const
{
    CMtx4x4<ET> mr;
    mr(0,0) = El(0,0);
    mr(0,1) = El(1,0);
    mr(0,2) = El(2,0);
	mr(0,3) = El(3,0);
    mr(1,0) = El(0,1);
    mr(1,1) = El(1,1);
    mr(1,2) = El(2,1);
	mr(1,3) = El(3,1);
    mr(2,0) = El(0,2);
    mr(2,1) = El(1,2);
    mr(2,2) = El(2,2);
	mr(2,3) = El(3,2);
    mr(3,0) = El(0,3);
    mr(3,1) = El(1,3);
    mr(3,2) = El(2,3);
	mr(3,3) = El(3,3);
    return mr;
}

template <class ET>
inline CMtx3x3<ET> & CMtx3x3<ET>::MakeI()
{
    m_m[0] = m_m[4] = m_m[8] = (ET)1.0;
    m_m[1] = m_m[2] = m_m[3] = m_m[5] = m_m[6] = m_m[7] = 0;
    return *this;
}

//
//  invert matrix (Gauss Jordan with partial pivoting)
//

template <class ET>
inline static void
ExchangeRows3(int iRow1, int iRow2, CMtx3x3<ET> &mat)
{
    ET f;

    f = mat[iRow1][0];  mat[iRow1][0] = mat[iRow2][0];  mat[iRow2][0] = f;
    f = mat[iRow1][1];  mat[iRow1][1] = mat[iRow2][1];  mat[iRow2][1] = f;
    f = mat[iRow1][2];  mat[iRow1][2] = mat[iRow2][2];  mat[iRow2][2] = f;
}

template <class ET>
inline static void
ScaleRow3(int iRow, ET f, CMtx3x3<ET> &mat)
{
    mat[iRow][0] *= f;
    mat[iRow][1] *= f;
    mat[iRow][2] *= f;
}

template <class ET>
inline static void
AddScaled3(int iRow1, ET f, int iRow2, CMtx3x3<ET> &mat)
{
    mat[iRow2][0] += f*mat[iRow1][0];
    mat[iRow2][1] += f*mat[iRow1][1];
    mat[iRow2][2] += f*mat[iRow1][2];
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx3x3<ET>::Inv
//
//  Synopsis:  invert the 3x3 matrix
//
//  Returns:   the new matrix. returns {0} if exactly singular
//
//------------------------------------------------------------------------

template <class ET>
CMtx3x3<ET> CMtx3x3<ET>::Inv() const
{
    CMtx3x3<ET> mr;
    mr.MakeI();
    CMtx3x3<ET> m = *this;

    int iRow, iCol, iPivotRow;
    ET fPivot, f;

    for (iCol = 0; iCol < Cols(); iCol++)
    {
        double fMax = 0;
        for (iRow = iCol; iRow < Rows(); iRow++)
        {
            // find maximum row element to pivot on
            double fm = armsReScalar(armsMagnitude(m(iRow,iCol)));
            if (fm > fMax)
            {
                iPivotRow = iRow;
                fMax = fm;
            }
        }
        if (fMax == 0)
        {
            // singular
            mr = (ET)0;
            return mr;
        }

        if (iPivotRow != iCol)
        {
            ExchangeRows3(iCol, iPivotRow, m);
            ExchangeRows3(iCol, iPivotRow, mr);
        }
        fPivot = m(iCol,iCol);
        ScaleRow3(iCol, ((ET)1)/fPivot, m);
        ScaleRow3(iCol, ((ET)1)/fPivot, mr);
        for (iRow = 0; iRow < Rows(); iRow++)
            if (iRow != iCol) 
            {
                f = -m(iRow,iCol);
                AddScaled3(iCol, f, iRow, m);
                AddScaled3(iCol, f, iRow, mr);
            }
    }

    return mr;
}


//template <class ET>
//CMtx3x3<ET> CMtx4x4<ET>::Submat3x3() 
//{
//	CMtx3x3<ET> retmat;
//	for (int i = 0; i < 3; i++)
//	{
//		for (int j = 0; j < 3; j++)
//			retmat(i, j) = (*this)(i, j);
//	}
//	return retmat;
//
//}




//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::Inv
//
//  Synopsis:  invert the 4x4 matrix using opencv function
//
//  Returns:   the new matrix. returns {0} if exactly singular
//
//------------------------------------------------------------------------

template <class ET>
CMtx4x4<ET> CMtx4x4<ET>::Inv() const
{
	cv::Mat mat(4, 4, CV_32F);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat.at<float>(i, j) = (*this)(i, j);
		}
	}

	cv::Mat invmat(4,4,CV_32F);
	cv::invert(mat, invmat);

	CMtx4x4<ET> mr;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mr(i, j) = invmat.at<float>(i, j);
		}
	}

	return mr;
}



template <class ET>
CVec4<ET> & CVec4<ET>::operator= (const CVec<ET> &v)
{
    if(v.Size()>0)
        El(0) = v[0];
    if(v.Size()>1)
        El(1) = v[1];
    if(v.Size()>2)
        El(2) = v[2];
    if(v.Size()>3)
        El(3) = v[3];
    return *this;
}

template <class ET>
CVec4<ET> & CVec4<ET>::Init(const ET *p, int iSize)
{
    if(p && iSize>0)
    {
        iSize = min(Size(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec4<ET>::operator^
//
//  Synopsis:  outer product operator
//
//  Returns:   the outer product matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CVec4<ET>::operator^ (const CVec4<ET> &v) const
{
    CMtx4x4<ET> m;
    m(0,0) = El(0) * v[0];
    m(0,1) = El(0) * v[1];
    m(0,2) = El(0) * v[2];
    m(0,3) = El(0) * v[3];
    m(1,0) = El(1) * v[0];
    m(1,1) = El(1) * v[1];
    m(1,2) = El(1) * v[2];
    m(1,3) = El(1) * v[3];
    m(2,0) = El(2) * v[0];
    m(2,1) = El(2) * v[1];
    m(2,2) = El(2) * v[2];
    m(2,3) = El(2) * v[3];
    m(3,0) = El(3) * v[0];
    m(3,1) = El(3) * v[1];
    m(3,2) = El(3) * v[2];
    m(3,3) = El(3) * v[3];
    return m;
}


//+-----------------------------------------------------------------------
//
//  Member:    CVec4<ET>::operator*
//
//  Synopsis:  right hand size 4x4 matrix multiply
//
//  Returns:   the new CVec4
//
//------------------------------------------------------------------------

template <class ET>
inline CVec4<ET> CVec4<ET>::operator* (const CMtx4x4<ET> &m) const
{
    return CVec4<ET>(El(0) * m(0,0) + El(1) * m(1,0) + El(2) * m(2,0) + El(3) * m(3,0),
        El(0) * m(0,1) + El(1) * m(1,1) + El(2) * m(2,1) + El(3) * m(3,1),
        El(0) * m(0,2) + El(1) * m(1,2) + El(2) * m(2,2) + El(3) * m(3,2),
        El(0) * m(0,3) + El(1) * m(1,3) + El(2) * m(2,3) + El(3) * m(3,3) );
}

//+-----------------------------------------------------------------------
//
//  Member:    CVec4<ET>::operator*=
//
//  Synopsis:  right hand size 4x4 matrix multiply
//
//  Returns:   the current CVec4
//
//------------------------------------------------------------------------

template <class ET>
inline CVec4<ET> & CVec4<ET>::operator*= (const CMtx4x4<ET> &m)
{
    *this = *this * m;
    return *this;
}



//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::MakeRotationX
//
//  Synopsis:  make a 4x4 rotation matrix from the angle given
//             rotation is about the X axis
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

//template <class ET>
//CMtx4x4<ET> & CMtx4x4<ET>::MakeRotationX(ET fTheta)
//{
//    ET fSin, fCos;
//    VtSinCos( fTheta, &fSin, &fCos );
//
//    El(0,0) = (ET)1; El(0,1) = (ET)0; El(0,2) = (ET)0;  El(0,3) = (ET)0; 
//    El(1,0) = (ET)0; El(1,1) = fCos; El(1,2) = -fSin; El(1,3) = (ET)0; 
//    El(2,0) = (ET)0; El(2,1) = fSin; El(2,2) = fCos;  El(2,3) = (ET)0; 
//    El(3,0) = (ET)0; El(3,1) = (ET)0; El(3,2) = (ET)0;  El(3,3) = (ET)1;
//    return *this;
//}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::MakeRotationY
//
//  Synopsis:  make a 4x4 rotation matrix from the angle given
//             rotation is about the Y axis
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

//template <class ET>
//CMtx4x4<ET> & CMtx4x4<ET>::MakeRotationY(ET fTheta)
//{
//    ET fSin, fCos;
//    VtSinCos( fTheta, &fSin, &fCos );
//
//    El(0,0) = fCos;  El(0,1) = (ET)0; El(0,2) = fSin;  El(0,3) = (ET)0; 
//    El(1,0) = (ET)0;  El(1,1) = (ET)1; El(1,2) = (ET)0;  El(1,3) = (ET)0; 
//    El(2,0) = -fSin; El(2,1) = (ET)0; El(2,2) = fCos;  El(2,3) = (ET)0; 
//    El(3,0) = (ET)0;  El(3,1) = (ET)0; El(3,2) = (ET)0;  El(3,3) = (ET)1;
//    return *this;
//}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::MakeRotationZ
//
//  Synopsis:  make a 4x4 rotation matrix from the angle given
//             rotation is about the Z axis
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

//template <class ET>
//CMtx4x4<ET> & CMtx4x4<ET>::MakeRotationZ(ET fTheta)
//{
//    ET fSin, fCos;
//    VtSinCos( fTheta, &fSin, &fCos );
//
//    El(0,0) = fCos; El(0,1) = -fSin; El(0,2) = (ET)0; El(0,3) = (ET)0; 
//    El(1,0) = fSin; El(1,1) = fCos;  El(1,2) = (ET)0; El(1,3) = (ET)0; 
//    El(2,0) = (ET)0; El(2,1) = (ET)0;  El(2,2) = (ET)1; El(2,3) = (ET)0; 
//    El(3,0) = (ET)0; El(3,1) = (ET)0;  El(3,2) = (ET)0; El(3,3) = (ET)1;
//    return *this;
//}


//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator=
//
//  Synopsis:  general assigment from another matrix of arbitrary size
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx4x4<ET> & CMtx4x4<ET>::operator= (const CMtx<ET> &m)
{
    int iMinR = min(m.Rows(), Rows());
    int iMinC = min(m.Cols(), Cols());
    int i, j;
    for(i=0; i<iMinR; i++)
        for(j=0; j<iMinC; j++)
            El(i,j) = m.El(i,j);
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator/
//
//  Synopsis:  division by a constant
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::operator/ (const ET &f)  const
{
    CMtx4x4 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] / f;
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator-
//
//  Synopsis:  unary minus
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::operator- () const
{
    CMtx4x4 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = -m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator+
//
//  Synopsis:  addition
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::operator+ (const CMtx4x4<ET> &m) const
{
    CMtx4x4 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] + m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator-
//
//  Synopsis:  subtraction
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::operator- (const CMtx4x4<ET> &m) const
{
    CMtx4x4 mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] - m.m_m[i];
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator*
//
//  Synopsis:  scaling by a constant
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> CMtx4x4<ET>::operator*(const ET &f) const
{
    CMtx4x4<ET> mr;
    int i;
    for(i=0; i<Size(); i++)
        mr.m_m[i] = m_m[i] * f;
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator=
//
//  Synopsis:  set all elements to a constant
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator=(const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] = f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::Init
//
//  Synopsis:  initialize from an array pointer
//
//  Returns:   the current CMtx4x4
//
//------------------------------------------------------------------------

template <class ET>
CMtx4x4<ET> & CMtx4x4<ET>::Init(const ET *p, int iSize)
{
    if(p && iSize>0)
    {
        iSize = min(Rows() * Cols(), iSize);
        memcpy(Ptr(), p, sizeof(ET) * iSize); // memcpy handles size = 0 case
    }
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator+=
//
//  Synopsis:  add a matrix to the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator+= (const CMtx4x4<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] += m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator-=
//
//  Synopsis:  subtract a matrix from the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator-= (const CMtx4x4<ET> &m)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] -= m.m_m[i];
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator*=
//
//  Synopsis:  multiply two matrices RHS multiply by RHS matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator*= (const CMtx4x4<ET> &m)
{
    *this = *this * m;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator*=
//
//  Synopsis:  scale the current matrix
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator*= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] *= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator/=
//
//  Synopsis:  divide the current matrix by a scale
//
//  Returns:   the current matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::operator/= (const ET &f)
{
    int i;
    for(i=0; i<Size(); i++)
        m_m[i] /= f;
    return *this;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator*
//
//  Synopsis:  matrix multiply
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
CMtx4x4<ET> CMtx4x4<ET>::operator*(const CMtx4x4<ET> &m) const
{
    CMtx4x4<ET> mr;
    mr(0,0) = El(0,0) * m(0,0) + El(0,1) * m(1,0) + El(0,2) * m(2,0) + El(0,3) * m(3,0);
    mr(0,1) = El(0,0) * m(0,1) + El(0,1) * m(1,1) + El(0,2) * m(2,1) + El(0,3) * m(3,1);
    mr(0,2) = El(0,0) * m(0,2) + El(0,1) * m(1,2) + El(0,2) * m(2,2) + El(0,3) * m(3,2);
    mr(0,3) = El(0,0) * m(0,3) + El(0,1) * m(1,3) + El(0,2) * m(2,3) + El(0,3) * m(3,3);

    mr(1,0) = El(1,0) * m(0,0) + El(1,1) * m(1,0) + El(1,2) * m(2,0) + El(1,3) * m(3,0);
    mr(1,1) = El(1,0) * m(0,1) + El(1,1) * m(1,1) + El(1,2) * m(2,1) + El(1,3) * m(3,1);
    mr(1,2) = El(1,0) * m(0,2) + El(1,1) * m(1,2) + El(1,2) * m(2,2) + El(1,3) * m(3,2);
    mr(1,3) = El(1,0) * m(0,3) + El(1,1) * m(1,3) + El(1,2) * m(2,3) + El(1,3) * m(3,3);

    mr(2,0) = El(2,0) * m(0,0) + El(2,1) * m(1,0) + El(2,2) * m(2,0) + El(2,3) * m(3,0);
    mr(2,1) = El(2,0) * m(0,1) + El(2,1) * m(1,1) + El(2,2) * m(2,1) + El(2,3) * m(3,1);
    mr(2,2) = El(2,0) * m(0,2) + El(2,1) * m(1,2) + El(2,2) * m(2,2) + El(2,3) * m(3,2);
    mr(2,3) = El(2,0) * m(0,3) + El(2,1) * m(1,3) + El(2,2) * m(2,3) + El(2,3) * m(3,3);

    mr(3,0) = El(3,0) * m(0,0) + El(3,1) * m(1,0) + El(3,2) * m(2,0) + El(3,3) * m(3,0);
    mr(3,1) = El(3,0) * m(0,1) + El(3,1) * m(1,1) + El(3,2) * m(2,1) + El(3,3) * m(3,1);
    mr(3,2) = El(3,0) * m(0,2) + El(3,1) * m(1,2) + El(3,2) * m(2,2) + El(3,3) * m(3,2);
    mr(3,3) = El(3,0) * m(0,3) + El(3,1) * m(1,3) + El(3,2) * m(2,3) + El(3,3) * m(3,3);
    return mr;
}

//+-----------------------------------------------------------------------
//
//  Member:    CMtx4x4<ET>::operator*
//
//  Synopsis:  multiply matrix by a 3x1 vector
//
//  Returns:   the new matrix
//
//------------------------------------------------------------------------

template <class ET>
inline CVec4<ET> CMtx4x4<ET>::operator*(const CVec4<ET> &v) const
{
    CVec4<ET> vr;
    vr(0) = El(0,0) * v(0) + El(0,1) * v(1) + El(0,2) * v(2) + El(0,3) * v(3);
    vr(1) = El(1,0) * v(0) + El(1,1) * v(1) + El(1,2) * v(2) + El(1,3) * v(3);
    vr(2) = El(2,0) * v(0) + El(2,1) * v(1) + El(2,2) * v(2) + El(2,3) * v(3);
    vr(3) = El(3,0) * v(0) + El(3,1) * v(1) + El(3,2) * v(2) + El(3,3) * v(3);
    return vr;
}




template <class ET>
inline CMtx4x4<ET> & CMtx4x4<ET>::MakeI()
{
    memset(m_m, 0, Size() * sizeof(ET));
    m_m[0] = m_m[5] = m_m[10] = m_m[15] = (ET)1.0;
    return *this;
}


//template <class ET>
//inline CVec3<ET> CMtx4x4<ET>::transform_point(const CVec3<ET> &v) const
//{
//	CVec3<ET> vr;
//
//	vr(0) = El(0, 0) * v(0) + El(0, 1) * v(1) + El(0, 2) * v(2) + El(0, 3);
//	vr(1) = El(1, 0) * v(0) + El(1, 1) * v(1) + El(1, 2) * v(2) + El(1, 3);
//	vr(2) = El(2, 0) * v(0) + El(2, 1) * v(1) + El(2, 2) * v(2) + El(2, 3);
//
//	return vr;
//}
//
//template <class ET>
//inline CVec3<ET> CMtx4x4<ET>::transform_vector(const CVec3<ET> &v) const
//{
//	CVec3<ET> vr;
//
//	vr(0) = El(0, 0) * v(0) + El(0, 1) * v(1) + El(0, 2) * v(2);
//	vr(1) = El(1, 0) * v(0) + El(1, 1) * v(1) + El(1, 2) * v(2);
//	vr(2) = El(2, 0) * v(0) + El(2, 1) * v(1) + El(2, 2) * v(2);
//
//	return vr;
//}






