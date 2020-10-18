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

#include "estimators/triangulation.h"

#include <Eigen/Geometry>

#include "base/projection.h"
#include "base/triangulation.h"
#include "estimators/essential_matrix.h"
#include "optim/combination_sampler.h"
#include "optim/loransac.h"
#include "util/logging.h"
#include "util/math.h"

namespace colmap {

void TriangulationEstimator::SetMinTriAngle(const double min_tri_angle) {
  CHECK_GE(min_tri_angle, 0);
  min_tri_angle_ = min_tri_angle;
}

void TriangulationEstimator::SetResidualType(const ResidualType residual_type) {
  residual_type_ = residual_type;
}

std::vector<TriangulationEstimator::M_t> TriangulationEstimator::Estimate(
    const std::vector<X_t>& point_data,
    const std::vector<Y_t>& pose_data) const {
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());

  if (point_data.size() == 2) {
    // Two-view triangulation.

    const M_t xyz = TriangulatePoint(
        pose_data[0].proj_matrix, pose_data[1].proj_matrix,
        point_data[0].point_normalized, point_data[1].point_normalized);

    if (HasPointPositiveDepth(pose_data[0].proj_matrix, xyz) &&
        HasPointPositiveDepth(pose_data[1].proj_matrix, xyz) &&
        CalculateTriangulationAngle(pose_data[0].proj_center,
                                    pose_data[1].proj_center,
                                    xyz) >= min_tri_angle_) {
      return std::vector<M_t>{xyz};
    }
  } else {
    // Multi-view triangulation.

    std::vector<Eigen::Matrix3x4d> proj_matrices;
    proj_matrices.reserve(point_data.size());
    std::vector<Eigen::Vector2d> points;
    points.reserve(point_data.size());
    for (size_t i = 0; i < point_data.size(); ++i) {
      proj_matrices.push_back(pose_data[i].proj_matrix);
      points.push_back(point_data[i].point_normalized);
    }

    const M_t xyz = TriangulateMultiViewPoint(proj_matrices, points);

    // Check for cheirality constraint.
    for (const auto& pose : pose_data) {
      if (!HasPointPositiveDepth(pose.proj_matrix, xyz)) {
        return std::vector<M_t>();
      }
    }

    // Check for sufficient triangulation angle.
    for (size_t i = 0; i < pose_data.size(); ++i) {
      for (size_t j = 0; j < i; ++j) {
        const double tri_angle = CalculateTriangulationAngle(
            pose_data[i].proj_center, pose_data[j].proj_center, xyz);
        if (tri_angle >= min_tri_angle_) {
          return std::vector<M_t>{xyz};
        }
      }
    }
  }

  return std::vector<M_t>();
}

void TriangulationEstimator::Residuals(const std::vector<X_t>& point_data,
                                       const std::vector<Y_t>& pose_data,
                                       const M_t& xyz,
                                       std::vector<double>* residuals) const {
  CHECK_EQ(point_data.size(), pose_data.size());

  residuals->resize(point_data.size());

  for (size_t i = 0; i < point_data.size(); ++i) {
    if (residual_type_ == ResidualType::REPROJECTION_ERROR) {
      (*residuals)[i] = CalculateSquaredReprojectionError(
          point_data[i].point, xyz, pose_data[i].proj_matrix,
          *pose_data[i].camera);
    } else if (residual_type_ == ResidualType::ANGULAR_ERROR) {
      const double angular_error = CalculateNormalizedAngularError(
          point_data[i].point_normalized, xyz, pose_data[i].proj_matrix);
      (*residuals)[i] = angular_error * angular_error;
    }
  }
}

bool EstimateTriangulation(
    const EstimateTriangulationOptions& options,
    const std::vector<TriangulationEstimator::PointData>& point_data,
    const std::vector<TriangulationEstimator::PoseData>& pose_data,
    std::vector<char>* inlier_mask, Eigen::Vector3d* xyz) {
  CHECK_NOTNULL(inlier_mask);
  CHECK_NOTNULL(xyz);
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());
  options.Check();

  // Robustly estimate track using LORANSAC.
  LORANSAC<TriangulationEstimator, TriangulationEstimator,
           InlierSupportMeasurer, CombinationSampler>
      ransac(options.ransac_options);
  ransac.estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.estimator.SetResidualType(options.residual_type);
  ransac.local_estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.local_estimator.SetResidualType(options.residual_type);
  const auto report = ransac.Estimate(point_data, pose_data);
  if (!report.success) {
    return false;
  }

  *inlier_mask = report.inlier_mask;
  *xyz = report.model;

  return report.success;
}

/////////////////////////////////////////////////////////////
// Line triangulation below

void LineTriangulationEstimator::SetMinTriAngle(const double min_tri_angle) {
  CHECK_GE(min_tri_angle, 0);
  min_tri_angle_ = min_tri_angle;
}


void LineTriangulationEstimator::SetResidualType(const TriangulationEstimator::ResidualType residual_type) {
  residual_type_ = residual_type;
}


std::vector<LineTriangulationEstimator::M_t> LineTriangulationEstimator::Estimate(
    const std::vector<X_t>& line_data,
    const std::vector<Y_t>& pose_data) const {
  CHECK_GE(line_data.size(), 2);
  CHECK_EQ(line_data.size(), pose_data.size());
  
  Eigen::Matrix4d A;
  A.setZero();

  for(int i = 0; i < line_data.size(); ++i) {

    Eigen::Vector3d line2d = line_data[i].point1_normalized.homogeneous().cross(
            line_data[i].point2_normalized.homogeneous()).normalized();

    // backprojected plane
    Eigen::Vector4d p = pose_data[i].proj_matrix.transpose() * line2d;
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
  X1 = basis.block<3,1>(0,1).normalized();
  
  // compute closest point on the line for each line segment endpoint
  std::vector<double> parameter_values;
  for(int i = 0; i < line_data.size(); ++i) {
    Eigen::Vector3d Z0 = pose_data[i].proj_matrix * X0.homogeneous();
    Eigen::Vector3d Z1 = pose_data[i].proj_matrix.leftCols<3>() * X1;

    Eigen::Matrix<double,3,2> M;
    M << Z1, line_data[i].point1_normalized.homogeneous();

    Eigen::Vector2d z = M.colPivHouseholderQr().solve(Z0);
    parameter_values.push_back(z(0));

    M.col(1) = line_data[i].point1_normalized.homogeneous();
    z = M.colPivHouseholderQr().solve(Z0);
    parameter_values.push_back(z(0));
  }

  auto mm = std::minmax_element(parameter_values.begin(), parameter_values.end());
  
  Eigen::Vector3d line_point1 = X0 + X1 * (*mm.first);
  Eigen::Vector3d line_point2 = X0 + X1 * (*mm.second);

  return {std::make_pair(line_point1, line_point2)};
}

void LineTriangulationEstimator::Residuals(const std::vector<X_t>& line_data,
                                       const std::vector<Y_t>& pose_data,
                                       const M_t& xyz,
                                       std::vector<double>* residuals) const {
  CHECK_EQ(line_data.size(), pose_data.size());

  residuals->resize(line_data.size());

  for (size_t i = 0; i < line_data.size(); ++i) {
    if (residual_type_ == TriangulationEstimator::ResidualType::REPROJECTION_ERROR) {
      (*residuals)[i] = CalculateSquaredLineReprojectionError(
        line_data[i].point1, line_data[i].point2, xyz.first, xyz.second, 
        pose_data[i].proj_matrix, *pose_data[i].camera);
    } else if (residual_type_ == TriangulationEstimator::ResidualType::ANGULAR_ERROR) {
      const double angular_error = CalculateNormalizedLineAngularError(
           line_data[i].point1_normalized, line_data[i].point2_normalized,
           xyz.first, xyz.second, pose_data[i].proj_matrix);
  
      (*residuals)[i] = angular_error * angular_error;
    }
  }
}

bool EstimateLineTriangulation(
    const EstimateTriangulationOptions& options,
    const std::vector<LineTriangulationEstimator::LineData>& line_data,
    const std::vector<LineTriangulationEstimator::PoseData>& pose_data,
    std::vector<char>* inlier_mask, std::pair<Eigen::Vector3d, Eigen::Vector3d>* xyz) {
  CHECK_NOTNULL(inlier_mask);
  CHECK_NOTNULL(xyz);
  CHECK_GE(line_data.size(), 2);
  CHECK_EQ(line_data.size(), pose_data.size());
  options.Check();

  // Robustly estimate track using LORANSAC.
  LORANSAC<LineTriangulationEstimator, LineTriangulationEstimator,
           InlierSupportMeasurer, CombinationSampler>
      ransac(options.ransac_options);
  ransac.estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.estimator.SetResidualType(options.residual_type);
  ransac.local_estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.local_estimator.SetResidualType(options.residual_type);
  const auto report = ransac.Estimate(line_data, pose_data);
  if (!report.success) {
    return false;
  }

  *inlier_mask = report.inlier_mask;
  *xyz = report.model;

  return report.success;
}

}  // namespace colmap
