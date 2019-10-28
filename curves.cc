#include "curves.hh"

#include <algorithm>
#include <cassert>
#include <limits>

#include <Eigen/LU>

double BezierCurve::bernstein(size_t i, size_t n, double u)
{
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t j = n; j >= k; --j)
      tmp[j] = tmp[j] * u1 + tmp[j-1] * u;
  return tmp[n];
}

Point BezierCurve::evaluateOneByOne(double u) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * bernstein(k, n, u);
  return p;
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector &coeff)
{
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

Point BezierCurve::evaluate(double u) const
{
  DoubleVector coeff; bernsteinAll(n, u, coeff);
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateWithCachedCofficients(const DoubleVector &coeff) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateByDeCasteljau(double u) const
{
  PointVector tmp = cp;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t i = 0; i <= n - k; ++i)
      tmp[i] = tmp[i] * u1 + tmp[i+1] * u;
  return tmp[0];
}

void BezierCurve::derivativeControlPoints(size_t d, PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  dcp[0] = cp;
  for (size_t k = 1; k <= d; ++k) {
    size_t tmp = n - k + 1;
    dcp[k].reserve(tmp);
    for (size_t i = 0; i <= n - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp);
  }
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleMatrix &coeff)
{
  coeff.clear(); coeff.resize(n + 1);
  coeff[0].push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    coeff[j].reserve(j + 1);
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[j-1][k];
      coeff[j].push_back(saved + tmp * u1);
      saved = tmp * u;
    }
    coeff[j].push_back(saved);
  }
}

Point BezierCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
  size_t du = std::min(d, n);
  der.clear(); der.reserve(d + 1);
  DoubleMatrix coeff; bernsteinAll(n, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, dcp);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= n - k; ++j)
      der[k] += dcp[k][j] * coeff[n-k][j];
  }
  for (size_t k = n + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

size_t BSplineCurve::findSpan(double u) const
{
  if (u == knots[n+1])
    return n;
  return (std::upper_bound(knots.begin() + p + 1, knots.end(), u) - knots.begin()) - 1;
}

void BSplineCurve::basisFunctions(size_t i, double u, DoubleVector &coeff) const
{
  coeff.clear(); coeff.reserve(p + 1);
  coeff.push_back(1.0);
  DoubleVector left(p + 1), right(p + 1);
  for (size_t j = 1; j <= p; ++j) {
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      double tmp = coeff[r] / (right[r+1] + left[j-r]);
      coeff[r] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    coeff.push_back(saved);
  }
}

Point BSplineCurve::evaluate(double u) const
{
  double span = findSpan(u);
  DoubleVector coeff; basisFunctions(span, u, coeff);
  Point point(0.0, 0.0, 0.0);
  for (size_t i = 0; i <= p; ++i)
    point += cp[span - p + i] * coeff[i];
  return point;
}

void BSplineCurve::basisFunctionDerivatives(size_t i, double u, size_t d, DoubleMatrix &der) const
{
  der.clear(); der.resize(d + 1);
  DoubleVector left(p + 1), right(p + 1), a[2];
  a[0].resize(p + 1); a[1].resize(p + 1);
  DoubleMatrix ndu(p + 1);
  ndu[0].resize(p + 1); ndu[0][0] = 1.0;
  for (size_t j = 1; j <= p; ++j) {
    ndu[j].resize(p + 1);
    left[j] = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      // lower triangle
      ndu[j][r] = right[r+1] + left[j-r];
      double tmp = ndu[r][j-1] / ndu[j][r];
      // upper triangle
      ndu[r][j] = saved + tmp * right[r+1];
      saved = tmp * left[j-r];
    }
    ndu[j][j] = saved;
  }
  for (size_t j = 0; j <= p; ++j)
    der[0].push_back(ndu[j][p]);
  for (size_t r = 0; r <= p; ++r) {
    size_t s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for (size_t k = 1; k <= d; ++k) {
      double dd = 0.0;
      int rk = r - k;
      int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk+1][rk];
        dd = a[s2][0] * ndu[rk][pk];
      }
      size_t j1 = rk >= -1 ? 1 : -rk;
      size_t j2 = (int)r - 1 <= pk ? k - 1 : p - r;
      for (size_t j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
        dd += a[s2][j] * ndu[rk+j][pk];
      }
      if (r <= (size_t)pk) {
        a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
        dd += a[s2][k] * ndu[r][pk];
      }
      der[k].push_back(dd);
      std::swap(s1, s2);
    }
  }
  size_t r = p;
  for (size_t k = 1; k <= d; ++k) {
    for (size_t j = 0; j <= p; ++j)
      der[k][j] *= r;
    r *= p - k;
  }
}

Point BSplineCurve::derivatives(double u, size_t d, VectorVector &der) const
{
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix nder; basisFunctionDerivatives(span, u, du, nder);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p; ++j)
      der[k] += cp[span-p+j] * nder[k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

void BSplineCurve::derivativeControlPoints(size_t d, size_t r1, size_t r2, PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  size_t r = r2 - r1;
  dcp[0].reserve(r + 1);
  for (size_t i = 0; i <= r; ++i)
    dcp[0].push_back(cp[r1+i]);
  for (size_t k = 1; k <= d; ++k) {
    dcp[k].reserve(r + 1 - k);
    size_t tmp = p - k + 1;
    for (size_t i = 0; i <= r - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp / (knots[r1+i+p+1] - knots[r1+i+k]));
  }
}

void BSplineCurve::basisFunctionsAll(size_t i, double u, DoubleMatrix &coeff) const
{
  coeff.clear(); coeff.resize(p + 1);
  coeff[0].push_back(1.0);
  DoubleVector left(p + 1), right(p + 1);
  for (size_t j = 1; j <= p; ++j) {
    coeff[j].reserve(j + 1);
    left[j]  = u - knots[i+1-j];
    right[j] = knots[i+j] - u;
    double saved = 0.0;
    for (size_t r = 0; r < j; ++r) {
      double tmp = coeff[j-1][r] / (right[r+1] + left[j-r]);
      coeff[j].push_back(saved + tmp * right[r+1]);
      saved = tmp * left[j-r];
    }
    coeff[j].push_back(saved);
  }
}

Point BSplineCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
  size_t du = std::min(d, p);
  der.clear();
  size_t span = findSpan(u);
  DoubleMatrix coeff; basisFunctionsAll(span, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, span - p, span, dcp);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= p - k; ++j)
      der[k] += dcp[k][j] * coeff[p-k][j];
  }
  for (size_t k = p + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}

size_t BSplineCurve::findSpanWithMultiplicity(double u, size_t &multi) const
{
  auto range = std::equal_range(knots.begin(), knots.end(), u);
  multi = range.second - range.first;

  if (u == knots[n+1])
    return n;
  return (range.second - knots.begin()) - 1;
}

Point BSplineCurve::evaluateByKnotInsertion(double u) const
{
  if (u == knots[0])
    return cp[0];
  if (u == knots[n+p+1])
    return cp[n];
  size_t s, k = findSpanWithMultiplicity(u, s), r = p - s;
  PointVector tmp; tmp.reserve(r + 1);
  std::copy_n(cp.begin() + k - p, r + 1, std::back_inserter(tmp));
  for (size_t j = 1; j <= r; ++j)
    for (size_t i = 0; i <= r - j; ++i) {
      double alpha = (u - knots[k-p+j+i]) / (knots[i+k+1] - knots[k-p+j+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
  return tmp[0];
}

Point BSplineCurve::evaluate2DRational(double u) const
{
  Point p = evaluate(u);
  return Point(p.x / p.z, p.y / p.z, 1.0);
}

size_t BSplineCurve::binomial(size_t n, size_t k)
{
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

Point BSplineCurve::derivatives2DRational(double u, size_t d, VectorVector &der) const
{
  der.clear(); der.reserve(d + 1);
  VectorVector der3d; derivativesByControlPoints(u, d, der3d);
  for (size_t k = 0; k <= d; ++k) {
    Vector v = der3d[k];
    for (size_t i = 1; i <= k; ++i)
      v = v - der[k-i] * der3d[i].z * binomial(k, i);
    der.push_back(v / der3d[0].z);
  }
  return der[0];
}

BSplineCurve BSplineCurve::insertKnot(double u, size_t k, size_t s, size_t r) const
{
  PointVector tmp; tmp.reserve(p - s + 1);

  BSplineCurve result;
  result.p = p; result.n = n + r;

  result.knots.reserve(knots.size() + r);
  std::copy_n(knots.begin(), k + 1, std::back_inserter(result.knots));
  std::fill_n(std::back_inserter(result.knots), r, u);
  std::copy(knots.begin() + k + 1, knots.end(), std::back_inserter(result.knots));

  result.cp.resize(cp.size() + r);
  std::copy_n(cp.begin(), k - p + 1, result.cp.begin());
  std::copy(cp.begin() + k - s, cp.end(), result.cp.begin() + r + k - s);

  std::copy_n(cp.begin() + k - p, p - s + 1, std::back_inserter(tmp));
  size_t L;
  for (size_t j = 1; j <= r; ++j) {
    L = k - p + j;
    for (size_t i = 0; i <= p - j - s; ++i) {
      double alpha = (u - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
      tmp[i] = tmp[i+1] * alpha + tmp[i] * (1.0 - alpha);
    }
    result.cp[L] = tmp[0];
    result.cp[k+r-j-s] = tmp[p-j-s];
  }
  std::copy_n(tmp.begin() + 1, p - s - 1 - r, result.cp.begin() + L + 1);

  return result;
}

BSplineCurve BSplineCurve::refineKnots(DoubleVector new_knots) const
{
  size_t r = new_knots.size();
  size_t a = findSpan(new_knots[0]);
  size_t b = findSpan(new_knots[r-1]) + 1;

  BSplineCurve result;
  result.p = p; result.n = n + r;
  result.knots.resize(knots.size() + r);
  result.cp.resize(cp.size() + r);

  std::copy_n(knots.begin(), a + 1, result.knots.begin());
  std::copy(knots.begin() + b + p, knots.end(), result.knots.begin() + b + p + r);
  std::copy_n(cp.begin(), a - p + 1, result.cp.begin());
  std::copy(cp.begin() + b - 1, cp.end(), result.cp.begin() + b - 1 + r);

  size_t i = b + p - 1,         // next position in knots
         k = b + p + r,         // next position in result.knots
         j = r;                 // next position in new_knots
  do {
    --j; --k;
    for (; new_knots[j] <= knots[i] && i > a; --i, --k) {
      result.knots[k] = knots[i];
      result.cp[k-p-1] = cp[i-p-1];
    }
    result.cp[k-p-1] = result.cp[k-p];
    for (size_t l = 1; l <= p; ++l) {
      size_t index = k - p + l;
      double alpha = result.knots[k+l] - new_knots[j];
      if (fabs(alpha) == 0.0)
        result.cp[index-1] = result.cp[index];
      else {
        alpha /= result.knots[k+l] - knots[i-p+l];
        result.cp[index-1] = result.cp[index-1] * alpha +
                             result.cp[index] * (1.0 - alpha);
      }
    }
    result.knots[k] = new_knots[j];
  } while (j > 0);
  return result;
}

Point BSplineCurve::projectPoint(const Point &point, double &u, double &distance,
                                 size_t resolution, double distance_tol, double cosine_tol) const
{
  // If we know that point is on the curve,
  // we could check only the spans of the convex hull point is in
  double span_min = knots[p], span_max = knots[n+1];
  distance = std::numeric_limits<double>::max();
  for (size_t i = 0; i < resolution; ++i) {
    double v = span_min + (span_max - span_min) * (double)i / (double)(resolution - 1);
    double d = (evaluate(v) - point).norm();
    if (d < distance) {
      distance = d;
      u = v;
    }
  }

  VectorVector der;
  derivativesByControlPoints(u, 2, der);
  Vector deviation = der[0] - point;

  while (distance > distance_tol) {
    double scaled_error = der[1] * deviation;
    double cosine_err = fabs(scaled_error) / (der[1].norm() * distance);
    if (cosine_err < cosine_tol)
      break;
    
    double old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, span_min), span_max);

    if ((der[1] * (u - old)).norm() < distance_tol)
      break;

    derivativesByControlPoints(u, 2, der);
    deviation = der[0] - point;
    distance = deviation.norm();
  }

  return der[0];
}

BSplineCurve BSplineCurve::interpolate(PointVector points, size_t p)
{
  size_t n = points.size() - 1;
  size_t m = n + p + 1;
  DoubleVector u; u.reserve(n + 1);
  u.push_back(0.0);
  double total = 0.0;
  for (size_t i = 1; i <= n; ++i) {
    u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    total += u.back();
  }
  for (size_t i = 1; i < n; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(m + 1);
  result.cp.reserve(n + 1);

  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for (size_t j = 1; j <= n - p; ++j) {
    double t = 0.0;
    for (size_t i = j; i <= j + p - 1; ++i)
      t += u[i];
    result.knots.push_back(t / p);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);

  Eigen::MatrixXd A(n + 1, n + 1); A.setZero();
  for (size_t i = 0; i <= n; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for (size_t j = 0; j <= p; ++j)
      A(i, span - p + j) = coeff[j];
  }
  Eigen::MatrixXd b(n + 1, 3);
  for (size_t i = 0; i <= n; ++i) {
    b(i, 0) = points[i].x;
    b(i, 1) = points[i].y;
    b(i, 2) = points[i].z;
  }
  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t i = 0; i <= n; ++i)
    result.cp.emplace_back(x(i, 0), x(i, 1), x(i, 2));

  return result;
}

// Simple version - no constraints, weights or derivatives as input.
// Approximates `points' with a B-spline of degree `p' with `n'+1 control points.
// This is essentially the same algorithm as `interpolate', just the knotvector generation differs.
BSplineCurve BSplineCurve::approximate(PointVector points, size_t p, size_t n)
{
  size_t m = points.size() - 1;
  assert(n <= m);

  BSplineCurve result;
  result.p = p; result.n = n;
  result.knots.reserve(n + p + 2);
  result.cp.reserve(n + 1);

  // Generate parameters for the data points by the centripetal method
  DoubleVector u; u.reserve(m + 1);
  u.push_back(0.0);
  double total = 0.0;
  for (size_t i = 1; i <= m; ++i) {
    u.push_back(std::sqrt((points[i] - points[i-1]).norm()));
    total += u.back();
  }
  for (size_t i = 1; i < m; ++i)
    u[i] = u[i-1] + u[i] / total;
  u.back() = 1.0;

  // Generate knotvector
  std::fill_n(std::back_inserter(result.knots), p + 1, 0.0);
  for (size_t j = 1; j <= n - p; ++j) {
    double d = (double)(m + 1) / (n - p + 1);
    size_t i = d * j;
    double alpha = d * j - i;
    double knot = (1.0 - alpha) * u[i-1] + alpha * u[i];
    result.knots.push_back(knot);
  }
  std::fill_n(std::back_inserter(result.knots), p + 1, 1.0);
  
  // Set up the equation
  Eigen::MatrixXd N(m + 1, n + 1); N.setZero();
  Eigen::MatrixXd S(m + 1, 3);
  for (size_t i = 0; i <= m; ++i) {
    size_t span = result.findSpan(u[i]);
    DoubleVector coeff; result.basisFunctions(span, u[i], coeff);
    for (size_t j = 0; j <= p; ++j)
      N(i, span - p + j) = coeff[j];
    S(i, 0) = points[i].x;
    S(i, 1) = points[i].y;
    S(i, 2) = points[i].z;
  }

  // Fill the control points
  Eigen::MatrixXd x = N.fullPivLu().solve(S);
  for (size_t i = 0; i <= n; ++i)
    result.cp.emplace_back(x(i, 0), x(i, 1), x(i, 2));

  return result;
}
