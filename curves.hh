#include "vector.hh"

#include <cstddef>

using Point        = Vector;
using PointVector  = std::vector<Point>;
using VectorVector = std::vector<Vector>;
using DoubleVector = std::vector<double>;
using DoubleMatrix = std::vector<DoubleVector>;
using PointMatrix  = std::vector<PointVector>;

struct BezierCurve
{
  size_t n;                     // degree => n + 1 control points
  PointVector cp;               // control points

  static double bernstein(size_t i, size_t n, double u);
  Point evaluateOneByOne(double u) const;
  static void bernsteinAll(size_t n, double u, DoubleVector &coeff);
  Point evaluate(double u) const;
  Point evaluateWithCachedCofficients(const DoubleVector &coeff) const;
  Point evaluateByDeCasteljau(double u) const;
  void derivativeControlPoints(size_t d, PointMatrix &dcp) const;
  static void bernsteinAll(size_t n, double u, DoubleMatrix &coeff);
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
};

struct BSplineCurve
{
  size_t p;                     // degree
  size_t n;                     // n + 1 = cp.size()
  DoubleVector knots;           // first and last p+1 values are the same ("clamped")
  PointVector cp;               // knots.size() = cp.size() + p + 1

  size_t findSpan(double u) const;
  void basisFunctions(size_t i, double u, DoubleVector &coeff) const;
  Point evaluate(double u) const;
  void basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const;
  Point derivatives(double u, size_t d, VectorVector &der) const;
  void derivativeControlPoints(size_t d, size_t r1, size_t r2, PointMatrix &dcp) const;
  void basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const;
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
  size_t findSpanWithMultiplicity(double u, size_t &multi) const;
  Point evaluateByKnotInsertion(double u) const;
  Point evaluate2DRational(double u) const;
  static size_t binomial(size_t n, size_t k);
  Point derivatives2DRational(double u, size_t d, VectorVector &der) const;
  BSplineCurve insertKnot(double u, size_t k, size_t s, size_t r) const;
  BSplineCurve refineKnots(DoubleVector new_knots) const;
  Point projectPoint(const Point &p, double &u, double &distance,
                     size_t resolution, double distance_tol, double cosine_tol) const;
  static BSplineCurve interpolate(PointVector points, size_t degree);
  static BSplineCurve approximate(PointVector points, size_t degree, size_t n);
};
