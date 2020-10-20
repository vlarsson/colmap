// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "base/triangulation.h"

#include "base/essential_matrix.h"
#include "base/pose.h"

namespace colmap {

Eigen::Vector3d TriangulatePoint(const Eigen::Matrix3x4d& proj_matrix1,
                                 const Eigen::Matrix3x4d& proj_matrix2,
                                 const Eigen::Vector2d& point1,
                                 const Eigen::Vector2d& point2) {
  Eigen::Matrix4d A;

  A.row(0) = point1(0) * proj_matrix1.row(2) - proj_matrix1.row(0);
  A.row(1) = point1(1) * proj_matrix1.row(2) - proj_matrix1.row(1);
  A.row(2) = point2(0) * proj_matrix2.row(2) - proj_matrix2.row(0);
  A.row(3) = point2(1) * proj_matrix2.row(2) - proj_matrix2.row(1);

  Eigen::JacobiSVD<Eigen::Matrix4d> svd(A, Eigen::ComputeFullV);

  return svd.matrixV().col(3).hnormalized();
}

std::vector<Eigen::Vector3d> TriangulatePoints(
    const Eigen::Matrix3x4d& proj_matrix1,
    const Eigen::Matrix3x4d& proj_matrix2,
    const std::vector<Eigen::Vector2d>& points1,
    const std::vector<Eigen::Vector2d>& points2) {
  CHECK_EQ(points1.size(), points2.size());

  std::vector<Eigen::Vector3d> points3D(points1.size());

  for (size_t i = 0; i < points3D.size(); ++i) {
    points3D[i] =
        TriangulatePoint(proj_matrix1, proj_matrix2, points1[i], points2[i]);
  }

  return points3D;
}

Eigen::Vector3d TriangulateMultiViewPoint(
    const std::vector<Eigen::Matrix3x4d>& proj_matrices,
    const std::vector<Eigen::Vector2d>& points) {
  CHECK_EQ(proj_matrices.size(), points.size());

  Eigen::Matrix4d A = Eigen::Matrix4d::Zero();

  for (size_t i = 0; i < points.size(); i++) {
    const Eigen::Vector3d point = points[i].homogeneous().normalized();
    const Eigen::Matrix3x4d term =
        proj_matrices[i] - point * point.transpose() * proj_matrices[i];
    A += term.transpose() * term;
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigen_solver(A);

  return eigen_solver.eigenvectors().col(0).hnormalized();
}

Eigen::Vector3d TriangulateOptimalPoint(const Eigen::Matrix3x4d& proj_matrix1,
                                        const Eigen::Matrix3x4d& proj_matrix2,
                                        const Eigen::Vector2d& point1,
                                        const Eigen::Vector2d& point2) {
  const Eigen::Matrix3d E =
      EssentialMatrixFromAbsolutePoses(proj_matrix1, proj_matrix2);

  Eigen::Vector2d optimal_point1;
  Eigen::Vector2d optimal_point2;
  FindOptimalImageObservations(E, point1, point2, &optimal_point1,
                               &optimal_point2);

  return TriangulatePoint(proj_matrix1, proj_matrix2, optimal_point1,
                          optimal_point2);
}

std::vector<Eigen::Vector3d> TriangulateOptimalPoints(
    const Eigen::Matrix3x4d& proj_matrix1,
    const Eigen::Matrix3x4d& proj_matrix2,
    const std::vector<Eigen::Vector2d>& points1,
    const std::vector<Eigen::Vector2d>& points2) {
  std::vector<Eigen::Vector3d> points3D(points1.size());

  for (size_t i = 0; i < points3D.size(); ++i) {
    points3D[i] =
        TriangulatePoint(proj_matrix1, proj_matrix2, points1[i], points2[i]);
  }

  return points3D;
}

double CalculateTriangulationAngle(const Eigen::Vector3d& proj_center1,
                                   const Eigen::Vector3d& proj_center2,
                                   const Eigen::Vector3d& point3D) {
  const double baseline_length_squared =
      (proj_center1 - proj_center2).squaredNorm();

  const double ray_length_squared1 = (point3D - proj_center1).squaredNorm();
  const double ray_length_squared2 = (point3D - proj_center2).squaredNorm();

  // Using "law of cosines" to compute the enclosing angle between rays.
  const double denominator =
      2.0 * std::sqrt(ray_length_squared1 * ray_length_squared2);
  if (denominator == 0.0) {
    return 0.0;
  }
  const double nominator =
      ray_length_squared1 + ray_length_squared2 - baseline_length_squared;
  const double angle = std::abs(std::acos(nominator / denominator));

  // Triangulation is unstable for acute angles (far away points) and
  // obtuse angles (close points), so always compute the minimum angle
  // between the two intersecting rays.
  return std::min(angle, M_PI - angle);
}

std::vector<double> CalculateTriangulationAngles(
    const Eigen::Vector3d& proj_center1, const Eigen::Vector3d& proj_center2,
    const std::vector<Eigen::Vector3d>& points3D) {
  // Baseline length between camera centers.
  const double baseline_length_squared =
      (proj_center1 - proj_center2).squaredNorm();

  std::vector<double> angles(points3D.size());

  for (size_t i = 0; i < points3D.size(); ++i) {
    // Ray lengths from cameras to point.
    const double ray_length_squared1 =
        (points3D[i] - proj_center1).squaredNorm();
    const double ray_length_squared2 =
        (points3D[i] - proj_center2).squaredNorm();

    // Using "law of cosines" to compute the enclosing angle between rays.
    const double denominator =
        2.0 * std::sqrt(ray_length_squared1 * ray_length_squared2);
    if (denominator == 0.0) {
      angles[i] = 0.0;
      continue;
    }
    const double nominator =
        ray_length_squared1 + ray_length_squared2 - baseline_length_squared;
    const double angle = std::abs(std::acos(nominator / denominator));

    // Triangulation is unstable for acute angles (far away points) and
    // obtuse angles (close points), so always compute the minimum angle
    // between the two intersecting rays.
    angles[i] = std::min(angle, M_PI - angle);
  }

  return angles;
}



// Calculate angle in radians between the two rays of a triangulated line.
double CalculateLineTriangulationAngle(const Eigen::Vector3d& proj_center1,
                                   const Eigen::Vector3d& proj_center2,
                                   const std::pair<Eigen::Vector3d, Eigen::Vector3d>& line3D) {

  Eigen::Vector3d ray1 = line3D.first - proj_center1;
  Eigen::Vector3d ray2 = line3D.first - proj_center2;

  // Remove component parallel to the line
  Eigen::Vector3d v = (line3D.first - line3D.second).normalized();
  ray1 = ray1 - v.dot(ray1) * v;
  ray2 = ray2 - v.dot(ray2) * v;

  // these are now the vectors that are orthogonal to the line and go to the projection centers
  // compute the angle between these
  const double cos_angle = std::max(-1.0,std::min(1.0, ray1.normalized().dot(ray2.normalized())));
  const double angle = std::abs(std::acos(cos_angle));
  // Triangulation is unstable for acute angles (far away points) and
  // obtuse angles (close points), so always compute the minimum angle
  // between the two intersecting rays.
  return std::min(angle, M_PI - angle);
}


std::vector<double> CalculateLineTriangulationAngles(
    const Eigen::Vector3d& proj_center1, const Eigen::Vector3d& proj_center2,
    const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& lines3D) {
      
  std::vector<double> angles(lines3D.size());
  for(size_t i = 0; i < lines3D.size(); ++i) {
    angles[i] = CalculateLineTriangulationAngle(proj_center1, proj_center2, lines3D[i]);
  }

  return angles;
}


// Triangulate 3D line from multiple line segments.
// Calls AdjustLineMultiView on the result
// 
//
// @param proj_matrices       Projection matrices of multi-view observations.
// @param line_segment        Image observations from multi-views.
//
// @return                    Estimated 3D line.
std::pair<Eigen::Vector3d, Eigen::Vector3d> TriangulateMultiViewLine(
    const std::vector<Eigen::Matrix3x4d>& proj_matrices,
    const std::vector<Eigen::Vector2d>& points1,
    const std::vector<Eigen::Vector2d>& points2) {


  CHECK_EQ(proj_matrices.size(), points1.size());
  CHECK_EQ(points1.size(), points2.size());

  Eigen::Matrix4d A;
  A.setZero();

  for(size_t i = 0; i < proj_matrices.size(); ++i) {

    Eigen::Vector3d line2d = points1[i].homogeneous().cross(points2[i].homogeneous()).normalized();

    // backprojected plane
    Eigen::Vector4d p = proj_matrices[i].transpose() * line2d;
    p = p / p.topRows<3>().norm();
    A +=  p * p.transpose();
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigen_solver(A);
  Eigen::Matrix<double, 4, 2> basis = eigen_solver.eigenvectors().leftCols<2>();
  
  // transform basis such that the line is represented
  //    X(lambda) = X0 + lambda * X1

  Eigen::Matrix2d H;
  H << basis(3,0), basis(3,1),
       basis(3,1),-basis(3,0);
  basis = basis * H;

  Eigen::Vector3d X0, X1;
  X0 = basis.block<3,1>(0,0) / basis(3,0);
  X1 = X0 + basis.block<3,1>(0,1).normalized();

  std::pair<Eigen::Vector3d, Eigen::Vector3d> xyz = 
    AdjustLineMultiView(std::make_pair(X0,X1), proj_matrices, points1, points2);

  return xyz;
}

// Ensures that the line segment in 3D corresponds to the hull of the union of 
// the 2D line segments.
std::pair<Eigen::Vector3d, Eigen::Vector3d> AdjustLineMultiView(
    const std::pair<Eigen::Vector3d, Eigen::Vector3d> &line3d, 
    const std::vector<Eigen::Matrix3x4d>& proj_matrices,
    const std::vector<Eigen::Vector2d>& points1,
    const std::vector<Eigen::Vector2d>& points2) {

  CHECK_EQ(proj_matrices.size(), points1.size());
  CHECK_EQ(points1.size(), points2.size());

  // Write line as X0 + t * X1
  Eigen::Vector3d X0 = line3d.first;
  Eigen::Vector3d X1 = line3d.second - X0;
    
  // compute closest point on the line for each line segment endpoint
  std::vector<double> parameter_values;
  parameter_values.reserve(2*proj_matrices.size());
  for(size_t i = 0; i < proj_matrices.size(); ++i) {


    // Transform X0 and X1 into the target frame
    // Note that X1 is a direction here and is just rotated
    Eigen::Vector3d Z0 = proj_matrices[i] * X0.homogeneous();
    Eigen::Vector3d Z1 = proj_matrices[i].leftCols<3>() * X1;

    // Z0 + t * Z1 = lambda * x
    // [Z1 -x] * [t; -lambda] = -Z0

    Eigen::Matrix<double,3,2> M;
    M << Z1, -points1[i].homogeneous();

    Eigen::Vector2d z = M.colPivHouseholderQr().solve(-Z0);
    parameter_values.push_back(z(0));

    M.col(1) = -points2[i].homogeneous();
    z = M.colPivHouseholderQr().solve(-Z0);
    parameter_values.push_back(z(0));
  }

  auto mm = std::minmax_element(parameter_values.begin(), parameter_values.end());
  
  // We take the maximum and minimum parameter values as the line end points
  Eigen::Vector3d line_point1 = X0 + X1 * (*mm.first);
  Eigen::Vector3d line_point2 = X0 + X1 * (*mm.second);

  return std::make_pair(line_point1, line_point2);
}

}  // namespace colmap
