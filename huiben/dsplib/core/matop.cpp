/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009-2010, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

/* ////////////////////////////////////////////////////////////////////
//
//  Mat basic operations: Copy, Set
//
// */

#include "precomp.hpp"

namespace libakaze
{

class MatOp_Identity : public MatOp
{
public:
    MatOp_Identity() {}
    virtual ~MatOp_Identity() {}

    bool elementWise(const MatExpr& /*expr*/) const { return true; }
    void assign(const MatExpr& expr, Mat& m, int type=-1) const;

    static void makeExpr(MatExpr& res, const Mat& m);
};

static MatOp_Identity g_MatOp_Identity;


class MatOp_AddEx : public MatOp
{
public:
    MatOp_AddEx() {}
    virtual ~MatOp_AddEx() {}

    bool elementWise(const MatExpr& /*expr*/) const { return true; }
    void assign(const MatExpr& expr, Mat& m, int type=-1) const;

    void add(const MatExpr& e1, const Scalar& s, MatExpr& res) const;
    void subtract(const Scalar& s, const MatExpr& expr, MatExpr& res) const;
    void multiply(const MatExpr& e1, double s, MatExpr& res) const;
    void divide(double s, const MatExpr& e, MatExpr& res) const;

    void transpose(const MatExpr& e1, MatExpr& res) const;
    void abs(const MatExpr& expr, MatExpr& res) const;

    static void makeExpr(MatExpr& res, const Mat& a, const Mat& b, double alpha, double beta, const Scalar& s=Scalar());
};

static MatOp_AddEx g_MatOp_AddEx;

class MatOp_Initializer : public MatOp
{
public:
    MatOp_Initializer() {}
    virtual ~MatOp_Initializer() {}

    bool elementWise(const MatExpr& /*expr*/) const { return false; }
    void assign(const MatExpr& expr, Mat& m, int type=-1) const;

    void multiply(const MatExpr& e, double s, MatExpr& res) const;

    static void makeExpr(MatExpr& res, int method, Size sz, int type, double alpha=1);
};

static MatOp_Initializer g_MatOp_Initializer;
static inline bool isIdentity(const MatExpr& e) { return e.op == &g_MatOp_Identity; }
static inline bool isAddEx(const MatExpr& e) { return e.op == &g_MatOp_AddEx; }
static inline bool isScaled(const MatExpr& e) { return isAddEx(e) && (!e.b.data || e.beta == 0) && e.s == Scalar(); }
static inline bool isInitializer(const MatExpr& e) { return e.op == &g_MatOp_Initializer; }

/////////////////////////////////////////////////////////////////////////////////////////////////////

bool MatOp::elementWise(const MatExpr& /*expr*/) const
{
    return false;
}



void MatOp::roi(const MatExpr& expr, const Range& rowRange, const Range& colRange, MatExpr& e) const
{
    if( elementWise(expr) )
    {
        e = MatExpr(expr.op, expr.flags, Mat(), Mat(), Mat(),
                    expr.alpha, expr.beta, expr.s);
        if(expr.a.data)
            e.a = expr.a(rowRange, colRange);
        if(expr.b.data)
            e.b = expr.b(rowRange, colRange);
        if(expr.c.data)
            e.c = expr.c(rowRange, colRange);
    }
    else
    {
        Mat m;
        expr.op->assign(expr, m);
        e = MatExpr(&g_MatOp_Identity, 0, m(rowRange, colRange), Mat(), Mat());
    }
}

void MatOp::diag(const MatExpr& expr, int d, MatExpr& e) const
{
    if( elementWise(expr) )
    {
        e = MatExpr(expr.op, expr.flags, Mat(), Mat(), Mat(),
                    expr.alpha, expr.beta, expr.s);
        if(expr.a.data)
            e.a = expr.a.diag(d);
        if(expr.b.data)
            e.b = expr.b.diag(d);
        if(expr.c.data)
            e.c = expr.c.diag(d);
    }
    else
    {
        Mat m;
        expr.op->assign(expr, m);
        e = MatExpr(&g_MatOp_Identity, 0, m.diag(d), Mat(), Mat());
    }
}


void MatOp::augAssignAdd(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m += temp;
}


void MatOp::augAssignSubtract(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m -= temp;
}


void MatOp::augAssignMultiply(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m *= temp;
}


void MatOp::augAssignDivide(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m /= temp;
}


void MatOp::augAssignAnd(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m &= temp;
}


void MatOp::augAssignOr(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m |= temp;
}


void MatOp::augAssignXor(const MatExpr& expr, Mat& m) const
{
    Mat temp;
    expr.op->assign(expr, temp);
    m ^= temp;
}


void MatOp::add(const MatExpr& e1, const MatExpr& e2, MatExpr& res) const
{
    if( this == e2.op )
    {
        double alpha = 1, beta = 1;
        Scalar s;
        Mat m1, m2;
        if( isAddEx(e1) && (!e1.b.data || e1.beta == 0) )
        {
            m1 = e1.a;
            alpha = e1.alpha;
            s = e1.s;
        }
        else
            e1.op->assign(e1, m1);

        if( isAddEx(e2) && (!e2.b.data || e2.beta == 0) )
        {
            m2 = e2.a;
            beta = e2.alpha;
            s += e2.s;
        }
        else
            e2.op->assign(e2, m2);
        MatOp_AddEx::makeExpr(res, m1, m2, alpha, beta, s);
    }
    else
        e2.op->add(e1, e2, res);
}


void MatOp::add(const MatExpr& expr1, const Scalar& s, MatExpr& res) const
{
    Mat m1;
    expr1.op->assign(expr1, m1);
    MatOp_AddEx::makeExpr(res, m1, Mat(), 1, 0, s);
}


void MatOp::subtract(const MatExpr& e1, const MatExpr& e2, MatExpr& res) const
{
    if( this == e2.op )
    {
        double alpha = 1, beta = -1;
        Scalar s;
        Mat m1, m2;
        if( isAddEx(e1) && (!e1.b.data || e1.beta == 0) )
        {
            m1 = e1.a;
            alpha = e1.alpha;
            s = e1.s;
        }
        else
            e1.op->assign(e1, m1);

        if( isAddEx(e2) && (!e2.b.data || e2.beta == 0) )
        {
            m2 = e2.a;
            beta = -e2.alpha;
            s -= e2.s;
        }
        else
            e2.op->assign(e2, m2);
        MatOp_AddEx::makeExpr(res, m1, m2, alpha, beta, s);
    }
    else
        e2.op->subtract(e1, e2, res);
}


void MatOp::subtract(const Scalar& s, const MatExpr& expr, MatExpr& res) const
{
    Mat m;
    expr.op->assign(expr, m);
    MatOp_AddEx::makeExpr(res, m, Mat(), -1, 0, s);
}


void MatOp::multiply(const MatExpr& e1, const MatExpr& e2, MatExpr& res, double scale) const
{
   1;
}


void MatOp::multiply(const MatExpr& expr, double s, MatExpr& res) const
{
    Mat m;
    expr.op->assign(expr, m);
    MatOp_AddEx::makeExpr(res, m, Mat(), s, 0);
}


void MatOp::divide(const MatExpr& e1, const MatExpr& e2, MatExpr& res, double scale) const
{
   1;
}


void MatOp::divide(double s, const MatExpr& expr, MatExpr& res) const
{
 1;
}


void MatOp::abs(const MatExpr& expr, MatExpr& res) const
{
   1;
}


void MatOp::transpose(const MatExpr& expr, MatExpr& res) const
{
   1;
}


void MatOp::matmul(const MatExpr& e1, const MatExpr& e2, MatExpr& res) const
{
  1;
}


void MatOp::invert(const MatExpr& expr, int method, MatExpr& res) const
{
   1;
}


Size MatOp::size(const MatExpr& expr) const
{
    return !expr.a.empty() ? expr.a.size() : expr.b.empty() ? expr.b.size() : expr.c.size();
}

int MatOp::type(const MatExpr& expr) const
{
    return !expr.a.empty() ? expr.a.type() : expr.b.empty() ? expr.b.type() : expr.c.type();
}
//////////////////////////////////////////////////////////////////////////////////////////////////

MatExpr::MatExpr(const Mat& m) : op(&g_MatOp_Identity), flags(0), a(m), b(Mat()), c(Mat()), alpha(1), beta(0), s(Scalar())
{
	1;
}

MatExpr operator - (const Mat& a, const Mat& b)
{
    MatExpr e;
    MatOp_AddEx::makeExpr(e, a, b, 1, -1);
    return e;
}

MatExpr operator - (const Mat& a, const Scalar& s)
{
    MatExpr e;
    MatOp_AddEx::makeExpr(e, a, Mat(), 1, 0, -s);
    return e;
}

MatExpr operator - (const Scalar& s, const Mat& a)
{
    MatExpr e;
    MatOp_AddEx::makeExpr(e, a, Mat(), -1, 0, s);
    return e;
}

MatExpr operator - (const MatExpr& e, const Mat& m)
{
    MatExpr en;
    e.op->subtract(e, MatExpr(m), en);
    return en;
}

MatExpr operator - (const Mat& m, const MatExpr& e)
{
    MatExpr en;
    e.op->subtract(MatExpr(m), e, en);
    return en;
}

MatExpr operator - (const MatExpr& e, const Scalar& s)
{
    MatExpr en;
    e.op->add(e, -s, en);
    return en;
}

MatExpr operator - (const Scalar& s, const MatExpr& e)
{
    MatExpr en;
    e.op->subtract(s, e, en);
    return en;
}

MatExpr operator - (const MatExpr& e1, const MatExpr& e2)
{
    MatExpr en;
    e1.op->subtract(e1, e2, en);
    return en;
}

MatExpr operator - (const Mat& m)
{
    MatExpr e;
    MatOp_AddEx::makeExpr(e, m, Mat(), -1, 0);
    return e;
}

MatExpr operator - (const MatExpr& e)
{
    MatExpr en;
    e.op->subtract(Scalar(0), e, en);
    return en;
}
Size MatExpr::size() const
{
   // if( isT(*this) || isInv(*this) )
   //     return Size(a.rows, a.cols);
//    if( isGEMM(*this) )
        return Size(b.cols, a.rows);
    //if( isSolve(*this) )
        return Size(b.cols, a.cols);
   // if( isInitializer(*this) )
    //    return a.size();
    return op ? op->size(*this) : Size();
}


int MatExpr::type() const
{
 //   if( isInitializer(*this) )
  //      return a.type();
 //   if( isCmp(*this) )
 //       return CV_8U;
    return op ? op->type(*this) : -1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

void MatOp_Identity::assign(const MatExpr& e, Mat& m, int _type) const
{
    if( _type == -1 || _type == e.a.type() )
        m = e.a;
    else
    {
        CV_Assert( CV_MAT_CN(_type) == e.a.channels() );
        e.a.convertTo(m, _type);
    }
}

inline void MatOp_Identity::makeExpr(MatExpr& res, const Mat& m)
{
    res = MatExpr(&g_MatOp_Identity, 0, m, Mat(), Mat(), 1, 0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void MatOp_AddEx::assign(const MatExpr& e, Mat& m, int _type) const
{
    Mat temp, &dst = _type == -1 || e.a.type() == _type ? m : temp;
    if( e.b.data )
    {
        if( e.s == Scalar() || !e.s.isReal() )
        {
            if( e.alpha == 1 )
            {
                if( e.beta == 1 )
                    libakaze::add(e.a, e.b, dst);
                else if( e.beta == -1 )
                    libakaze::subtract(e.a, e.b, dst);
                else
					1;
             //       libakaze::scaleAdd(e.b, e.beta, e.a, dst);
            }
            else if( e.beta == 1 )
            {
                if( e.alpha == -1 )
                    libakaze::subtract(e.b, e.a, dst);
                else
					1;
                //    libakaze::scaleAdd(e.a, e.alpha, e.b, dst);
            }
            else
               1;// libakaze::addWeighted(e.a, e.alpha, e.b, e.beta, 0, dst);

            if( !e.s.isReal() )
                libakaze::add(dst, e.s, dst);
        }
        else
           1;// libakaze::addWeighted(e.a, e.alpha, e.b, e.beta, e.s[0], dst);
    }
    else if( e.s.isReal() && (dst.data != m.data || fabs(e.alpha) != 1))
    {
        e.a.convertTo(m, _type, e.alpha, e.s[0]);
        return;
    }
    else if( e.alpha == 1 )
        libakaze::add(e.a, e.s, dst);
    else if( e.alpha == -1 )
        libakaze::subtract(e.s, e.a, dst);
    else
    {
        e.a.convertTo(dst, e.a.type(), e.alpha);
        libakaze::add(dst, e.s, dst);
    }

    if( dst.data != m.data )
        dst.convertTo(m, m.type());
}


void MatOp_AddEx::add(const MatExpr& e, const Scalar& s, MatExpr& res) const
{
    res = e;
    res.s += s;
}


void MatOp_AddEx::subtract(const Scalar& s, const MatExpr& e, MatExpr& res) const
{
    res = e;
    res.alpha = -res.alpha;
    res.beta = -res.beta;
    res.s = s - res.s;
}

void MatOp_AddEx::multiply(const MatExpr& e, double s, MatExpr& res) const
{
    res = e;
    res.alpha *= s;
    res.beta *= s;
    res.s *= s;
}

void MatOp_AddEx::divide(double s, const MatExpr& e, MatExpr& res) const
{
   // if( isScaled(e) )
     //   MatOp_Bin::makeExpr(res, '/', e.a, Mat(), s/e.alpha);
   // else
        MatOp::divide(s, e, res);
}


void MatOp_AddEx::transpose(const MatExpr& e, MatExpr& res) const
{
    if( isScaled(e) )
     1;//  MatOp_T::makeExpr(res, e.a, e.alpha);
    else
        MatOp::transpose(e, res);
}

void MatOp_AddEx::abs(const MatExpr& e, MatExpr& res) const
{
	
	1;
	/*
    if( (!e.b.data || e.beta == 0) && fabs(e.alpha) == 1 )
        MatOp_Bin::makeExpr(res, 'a', e.a, -e.s*e.alpha);
    else if( e.b.data && e.alpha + e.beta == 0 && e.alpha*e.beta == -1 )
        MatOp_Bin::makeExpr(res, 'a', e.a, e.b);
    else
        MatOp::abs(e, res);
		*/
}

inline void MatOp_AddEx::makeExpr(MatExpr& res, const Mat& a, const Mat& b, double alpha, double beta, const Scalar& s)
{
    res = MatExpr(&g_MatOp_AddEx, 0, a, b, Mat(), alpha, beta, s);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void MatOp_GEMM::assign(const MatExpr& e, Mat& m, int _type) const
{
    Mat temp, &dst = _type == -1 || _type == e.a.type() ? m : temp;

    libakaze::gemm(e.a, e.b, e.alpha, e.c, e.beta, dst, e.flags);
    if( dst.data != m.data )
        dst.convertTo(m, _type);
}

void MatOp_GEMM::add(const MatExpr& e1, const MatExpr& e2, MatExpr& res) const
{
    bool i1 = isIdentity(e1), i2 = isIdentity(e2);
    double alpha1 = i1 ? 1 : e1.alpha, alpha2 = i2 ? 1 : e2.alpha;

    if( isMatProd(e1) && (i2 || isScaled(e2) || isT(e2)) )
        MatOp_GEMM::makeExpr(res, (e1.flags & ~CV_GEMM_C_T)|(isT(e2) ? CV_GEMM_C_T : 0),
                             e1.a, e1.b, alpha1, e2.a, alpha2);
    else if( isMatProd(e2) && (i1 || isScaled(e1) || isT(e1)) )
        MatOp_GEMM::makeExpr(res, (e2.flags & ~CV_GEMM_C_T)|(isT(e1) ? CV_GEMM_C_T : 0),
                             e2.a, e2.b, alpha2, e1.a, alpha1);
    else if( this == e2.op )
        MatOp::add(e1, e2, res);
    else
        e2.op->add(e1, e2, res);
}

void MatOp_GEMM::subtract(const MatExpr& e1, const MatExpr& e2, MatExpr& res) const
{
    bool i1 = isIdentity(e1), i2 = isIdentity(e2);
    double alpha1 = i1 ? 1 : e1.alpha, alpha2 = i2 ? 1 : e2.alpha;

    if( isMatProd(e1) && (i2 || isScaled(e2) || isT(e2)) )
        MatOp_GEMM::makeExpr(res, (e1.flags & ~CV_GEMM_C_T)|(isT(e2) ? CV_GEMM_C_T : 0),
                             e1.a, e1.b, alpha1, e2.a, -alpha2);
    else if( isMatProd(e2) && (i1 || isScaled(e1) || isT(e1)) )
        MatOp_GEMM::makeExpr(res, (e2.flags & ~CV_GEMM_C_T)|(isT(e1) ? CV_GEMM_C_T : 0),
                            e2.a, e2.b, -alpha2, e1.a, alpha1);
    else if( this == e2.op )
        MatOp::subtract(e1, e2, res);
    else
        e2.op->subtract(e1, e2, res);
}

void MatOp_GEMM::multiply(const MatExpr& e, double s, MatExpr& res) const
{
    res = e;
    res.alpha *= s;
    res.beta *= s;
}

void MatOp_GEMM::transpose(const MatExpr& e, MatExpr& res) const
{
    res = e;
    res.flags = (!(e.flags & CV_GEMM_A_T) ? CV_GEMM_B_T : 0) |
                (!(e.flags & CV_GEMM_B_T) ? CV_GEMM_A_T : 0) |
                (!(e.flags & CV_GEMM_C_T) ? CV_GEMM_C_T : 0);
    swap(res.a, res.b);
}

inline void MatOp_GEMM::makeExpr(MatExpr& res, int flags, const Mat& a, const Mat& b,
                                 double alpha, const Mat& c, double beta)
{
    res = MatExpr(&g_MatOp_GEMM, flags, a, b, c, alpha, beta);
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////

void MatOp_Initializer::assign(const MatExpr& e, Mat& m, int _type) const
{
    if( _type == -1 )
        _type = e.a.type();
    m.create(e.a.size(), _type);
    if( e.flags == 'I' )
        setIdentity(m, Scalar(e.alpha));
    else if( e.flags == '0' )
        m = Scalar();
    else if( e.flags == '1' )
        m = Scalar(e.alpha);
    else
        CV_Error(CV_StsError, "Invalid matrix initializer type");
}

void MatOp_Initializer::multiply(const MatExpr& e, double s, MatExpr& res) const
{
    res = e;
    res.alpha *= s;
}

inline void MatOp_Initializer::makeExpr(MatExpr& res, int method, Size sz, int type, double alpha)
{
    res = MatExpr(&g_MatOp_Initializer, method, Mat(sz, type, (void*)0), Mat(), Mat(), alpha, 0);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////


MatExpr Mat::zeros(int rows, int cols, int type)
{
    MatExpr e;
    MatOp_Initializer::makeExpr(e, '0', Size(cols, rows), type);
    return e;
}

MatExpr Mat::zeros(Size size, int type)
{
    MatExpr e;
    MatOp_Initializer::makeExpr(e, '0', size, type);
    return e;
}


}

/* End of file. */
