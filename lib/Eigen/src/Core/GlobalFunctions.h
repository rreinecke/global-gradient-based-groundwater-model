// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010-2016 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_GLOBAL_FUNCTIONS_H
#define EIGEN_GLOBAL_FUNCTIONS_H

#ifdef EIGEN_PARSED_BY_DOXYGEN

#define EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(NAME,FUNCTOR,DOC_OP,DOC_DETAILS) \
  /** \returns an expression of the coefficient-wise DOC_OP of \a x

    DOC_DETAILS

    \sa <a href="group__CoeffwiseMathFunctions.html#cwisetable_##NAME">Math functions</a>, class CwiseUnaryOp
    */ \
  template<typename Derived> \
  inline const Eigen::CwiseUnaryOp<Eigen::internal::FUNCTOR<typename Derived::Scalar>, const Derived> \
  NAME(const Eigen::ArrayBase<Derived>& x);

#else

#define EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(NAME,FUNCTOR,DOC_OP,DOC_DETAILS) \
  template<typename Derived> \
  inline const Eigen::CwiseUnaryOp<Eigen::internal::FUNCTOR<typename Derived::Scalar>, const Derived> \
  (NAME)(const Eigen::ArrayBase<Derived>& x) { \
    return Eigen::CwiseUnaryOp<Eigen::internal::FUNCTOR<typename Derived::Scalar>, const Derived>(x.derived()); \
  }

#endif // EIGEN_PARSED_BY_DOXYGEN

#define EIGEN_ARRAY_DECLARE_GLOBAL_EIGEN_UNARY(NAME,FUNCTOR) \
  \
  template<typename Derived> \
  struct NAME##_retval<ArrayBase<Derived> > \
  { \
    typedef const Eigen::CwiseUnaryOp<Eigen::internal::FUNCTOR<typename Derived::Scalar>, const Derived> type; \
  }; \
  template<typename Derived> \
  struct NAME##_impl<ArrayBase<Derived> > \
  { \
    static inline typename NAME##_retval<ArrayBase<Derived> >::type run(const Eigen::ArrayBase<Derived>& x) \
    { \
      return typename NAME##_retval<ArrayBase<Derived> >::type(x.derived()); \
    } \
  };

namespace Eigen
{
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(real,scalar_real_op,real part,\sa ArrayBase::real)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(imag,scalar_imag_op,imaginary part,\sa ArrayBase::imag)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(conj,scalar_conjugate_op,complex conjugate,\sa ArrayBase::conjugate)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(inverse,scalar_inverse_op,inverse,\sa ArrayBase::inverse)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(sin,scalar_sin_op,sine,\sa ArrayBase::sin)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(cos,scalar_cos_op,cosine,\sa ArrayBase::cos)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(tan,scalar_tan_op,tangent,\sa ArrayBase::tan)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(atan,scalar_atan_op,arc-tangent,\sa ArrayBase::atan)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(asin,scalar_asin_op,arc-sine,\sa ArrayBase::asin)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(acos,scalar_acos_op,arc-consine,\sa ArrayBase::acos)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(sinh,scalar_sinh_op,hyperbolic sine,\sa ArrayBase::sinh)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(cosh,scalar_cosh_op,hyperbolic cosine,\sa ArrayBase::cosh)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(tanh,scalar_tanh_op,hyperbolic tangent,\sa ArrayBase::tanh)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(lgamma,scalar_lgamma_op,natural logarithm of the gamma function,\sa ArrayBase::lgamma)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(digamma,scalar_digamma_op,derivative of lgamma,\sa ArrayBase::digamma)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(erf,scalar_erf_op,error function,\sa ArrayBase::erf)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(erfc,scalar_erfc_op,complement error function,\sa ArrayBase::erfc)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(exp,scalar_exp_op,exponential,\sa ArrayBase::exp)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(log,scalar_log_op,natural logarithm,\sa Eigen::log10 DOXCOMMA ArrayBase::log)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(log1p,scalar_log1p_op,natural logarithm of 1 plus the value,\sa ArrayBase::log1p)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(log10,scalar_log10_op,base 10 logarithm,\sa Eigen::log DOXCOMMA ArrayBase::log)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(abs,scalar_abs_op,absolute value,\sa ArrayBase::abs DOXCOMMA MatrixBase::cwiseAbs)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(abs2,scalar_abs2_op,squared absolute value,\sa ArrayBase::abs2 DOXCOMMA MatrixBase::cwiseAbs2)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(arg,scalar_arg_op,complex argument,\sa ArrayBase::arg)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(sqrt,scalar_sqrt_op,square root,\sa ArrayBase::sqrt DOXCOMMA MatrixBase::cwiseSqrt)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(square,scalar_square_op,square (power 2),\sa Eigen::abs2 DOXCOMMA Eigen::pow DOXCOMMA ArrayBase::square)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(cube,scalar_cube_op,cube (power 3),\sa Eigen::pow DOXCOMMA ArrayBase::cube)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(round,scalar_round_op,nearest integer,\sa Eigen::floor DOXCOMMA Eigen::ceil DOXCOMMA ArrayBase::round)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(floor,scalar_floor_op,nearest integer not greater than the giben value,\sa Eigen::ceil DOXCOMMA ArrayBase::floor)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(ceil,scalar_ceil_op,nearest integer not less than the giben value,\sa Eigen::floor DOXCOMMA ArrayBase::ceil)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(isnan,scalar_isnan_op,not-a-number test,\sa Eigen::isinf DOXCOMMA Eigen::isfinite DOXCOMMA ArrayBase::isnan)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(isinf,scalar_isinf_op,infinite value test,\sa Eigen::isnan DOXCOMMA Eigen::isfinite DOXCOMMA ArrayBase::isinf)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(isfinite,scalar_isfinite_op,finite value test,\sa Eigen::isinf DOXCOMMA Eigen::isnan DOXCOMMA ArrayBase::isfinite)
  EIGEN_ARRAY_DECLARE_GLOBAL_UNARY(sign,scalar_sign_op,sign (or 0),\sa ArrayBase::sign)
  
  /** \returns an expression of the coefficient-wise power of \a x to the given constant \a exponent.
    *
    * \sa ArrayBase::pow()
    */
  template<typename Derived>
  inline const Eigen::CwiseUnaryOp<Eigen::internal::scalar_pow_op<typename Derived::Scalar>, const Derived>
  pow(const Eigen::ArrayBase<Derived>& x, const typename Derived::Scalar& exponent) {
    return x.derived().pow(exponent);
  }

  /** \returns an expression of the coefficient-wise power of \a x to the given array of \a exponents.
    *
    * This function computes the coefficient-wise power.
    *
    * Example: \include Cwise_array_power_array.cpp
    * Output: \verbinclude Cwise_array_power_array.out
    * 
    * \sa ArrayBase::pow()
    */
  template<typename Derived,typename ExponentDerived>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_binary_pow_op<typename Derived::Scalar, typename ExponentDerived::Scalar>, const Derived, const ExponentDerived>
  pow(const Eigen::ArrayBase<Derived>& x, const Eigen::ArrayBase<ExponentDerived>& exponents) 
  {
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_binary_pow_op<typename Derived::Scalar, typename ExponentDerived::Scalar>, const Derived, const ExponentDerived>(
      x.derived(),
      exponents.derived()
    );
  }
  
  /** \returns an expression of the coefficient-wise power of the scalar \a x to the given array of \a exponents.
    *
    * This function computes the coefficient-wise power between a scalar and an array of exponents.
    * Beaware that the scalar type of the input scalar \a x and the exponents \a exponents must be the same.
    *
    * Example: \include Cwise_scalar_power_array.cpp
    * Output: \verbinclude Cwise_scalar_power_array.out
    * 
    * \sa ArrayBase::pow()
    */
  template<typename Derived>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_binary_pow_op<typename Derived::Scalar, typename Derived::Scalar>, const typename Derived::ConstantReturnType, const Derived>
  pow(const typename Derived::Scalar& x, const Eigen::ArrayBase<Derived>& exponents) 
  {
    typename Derived::ConstantReturnType constant_x(exponents.rows(), exponents.cols(), x);
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_binary_pow_op<typename Derived::Scalar, typename Derived::Scalar>, const typename Derived::ConstantReturnType, const Derived>(
      constant_x,
      exponents.derived()
    );
  }
  
  /**
    * \brief Component-wise division of a scalar by array elements.
    **/
  template <typename Derived>
  inline const Eigen::CwiseUnaryOp<Eigen::internal::scalar_inverse_mult_op<typename Derived::Scalar>, const Derived>
    operator/(const typename Derived::Scalar& s, const Eigen::ArrayBase<Derived>& a)
  {
    return Eigen::CwiseUnaryOp<Eigen::internal::scalar_inverse_mult_op<typename Derived::Scalar>, const Derived>(
      a.derived(),
      Eigen::internal::scalar_inverse_mult_op<typename Derived::Scalar>(s)  
    );
  }

  /** \cpp11 \returns an expression of the coefficient-wise igamma(\a a, \a x) to the given arrays.
    *
    * This function computes the coefficient-wise incomplete gamma function.
    *
    * \note This function supports only float and double scalar types in c++11 mode. To support other scalar types,
    * or float/double in non c++11 mode, the user has to provide implementations of igammac(T,T) for any scalar
    * type T to be supported.
    *
    * \sa Eigen::igammac(), Eigen::lgamma()
    */
  template<typename Derived,typename ExponentDerived>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_igamma_op<typename Derived::Scalar>, const Derived, const ExponentDerived>
  igamma(const Eigen::ArrayBase<Derived>& a, const Eigen::ArrayBase<ExponentDerived>& x) 
  {
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_igamma_op<typename Derived::Scalar>, const Derived, const ExponentDerived>(
      a.derived(),
      x.derived()
    );
  }

  /** \cpp11 \returns an expression of the coefficient-wise igammac(\a a, \a x) to the given arrays.
    *
    * This function computes the coefficient-wise complementary incomplete gamma function.
    *
    * \note This function supports only float and double scalar types in c++11 mode. To support other scalar types,
    * or float/double in non c++11 mode, the user has to provide implementations of igammac(T,T) for any scalar
    * type T to be supported.
    *
    * \sa Eigen::igamma(), Eigen::lgamma()
    */
  template<typename Derived,typename ExponentDerived>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_igammac_op<typename Derived::Scalar>, const Derived, const ExponentDerived>
  igammac(const Eigen::ArrayBase<Derived>& a, const Eigen::ArrayBase<ExponentDerived>& x) 
  {
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_igammac_op<typename Derived::Scalar>, const Derived, const ExponentDerived>(
      a.derived(),
      x.derived()
    );
  }

  /** \cpp11 \returns an expression of the coefficient-wise polygamma(\a n, \a x) to the given arrays.
    *
    * It returns the \a n -th derivative of the digamma(psi) evaluated at \c x.
    *
    * \note This function supports only float and double scalar types in c++11 mode. To support other scalar types,
    * or float/double in non c++11 mode, the user has to provide implementations of polygamma(T,T) for any scalar
    * type T to be supported.
    *
    * \sa Eigen::digamma()
    */
  // * \warning Be careful with the order of the parameters: x.polygamma(n) is equivalent to polygamma(n,x)
  // * \sa ArrayBase::polygamma()
  template<typename DerivedN,typename DerivedX>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_polygamma_op<typename DerivedX::Scalar>, const DerivedN, const DerivedX>
  polygamma(const Eigen::ArrayBase<DerivedN>& n, const Eigen::ArrayBase<DerivedX>& x)
  {
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_polygamma_op<typename DerivedX::Scalar>, const DerivedN, const DerivedX>(
      n.derived(),
      x.derived()
    );
  }

  /** \cpp11 \returns an expression of the coefficient-wise betainc(\a x, \a a, \a b) to the given arrays.
    *
    * This function computes the regularized incomplete beta function (integral).
    *
    * \note This function supports only float and double scalar types in c++11 mode. To support other scalar types,
    * or float/double in non c++11 mode, the user has to provide implementations of betainc(T,T,T) for any scalar
    * type T to be supported.
    *
    * \sa Eigen::betainc(), Eigen::lgamma()
    */
  template<typename ArgADerived, typename ArgBDerived, typename ArgXDerived>
  inline const Eigen::CwiseTernaryOp<Eigen::internal::scalar_betainc_op<typename ArgXDerived::Scalar>, const ArgADerived, const ArgBDerived, const ArgXDerived>
  betainc(const Eigen::ArrayBase<ArgADerived>& a, const Eigen::ArrayBase<ArgBDerived>& b, const Eigen::ArrayBase<ArgXDerived>& x)
  {
    return Eigen::CwiseTernaryOp<Eigen::internal::scalar_betainc_op<typename ArgXDerived::Scalar>, const ArgADerived, const ArgBDerived, const ArgXDerived>(
      a.derived(),
      b.derived(),
      x.derived()
    );
  }


  /** \returns an expression of the coefficient-wise zeta(\a x, \a q) to the given arrays.
    *
    * It returns the Riemann zeta function of two arguments \a x and \a q:
    *
    * \param x is the exposent, it must be > 1
    * \param q is the shift, it must be > 0
    *
    * \note This function supports only float and double scalar types. To support other scalar types, the user has
    * to provide implementations of zeta(T,T) for any scalar type T to be supported.
    *
    * \sa ArrayBase::zeta()
    */
  template<typename DerivedX,typename DerivedQ>
  inline const Eigen::CwiseBinaryOp<Eigen::internal::scalar_zeta_op<typename DerivedX::Scalar>, const DerivedX, const DerivedQ>
  zeta(const Eigen::ArrayBase<DerivedX>& x, const Eigen::ArrayBase<DerivedQ>& q)
  {
    return Eigen::CwiseBinaryOp<Eigen::internal::scalar_zeta_op<typename DerivedX::Scalar>, const DerivedX, const DerivedQ>(
      x.derived(),
      q.derived()
    );
  }

  namespace internal
  {
    EIGEN_ARRAY_DECLARE_GLOBAL_EIGEN_UNARY(real,scalar_real_op)
    EIGEN_ARRAY_DECLARE_GLOBAL_EIGEN_UNARY(imag,scalar_imag_op)
    EIGEN_ARRAY_DECLARE_GLOBAL_EIGEN_UNARY(abs2,scalar_abs2_op)
  }
}

// TODO: cleanly disable those functions that are not supported on Array (numext::real_ref, internal::random, internal::isApprox...)

#endif // EIGEN_GLOBAL_FUNCTIONS_H
