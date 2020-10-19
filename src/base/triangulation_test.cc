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

#define TEST_NAME "base/triangulation"
#include "util/testing.h"

#include <Eigen/Core>

#include "base/similarity_transform.h"
#include "base/triangulation.h"

using namespace colmap;

BOOST_AUTO_TEST_CASE(TestTriangulatePoint) {
  std::vector<Eigen::Vector3d> points3D(6);
  points3D[0] = Eigen::Vector3d(0, 0.1, 0.1);
  points3D[1] = Eigen::Vector3d(0, 1, 3);
  points3D[2] = Eigen::Vector3d(0, 1, 2);
  points3D[3] = Eigen::Vector3d(0.01, 0.2, 3);
  points3D[4] = Eigen::Vector3d(-1, 0.1, 1);
  points3D[5] = Eigen::Vector3d(0.1, 0.1, 0.2);

  Eigen::Matrix3x4d proj_matrix1 = Eigen::MatrixXd::Identity(3, 4);

  for (double qz = 0; qz < 1; qz += 0.2) {
    for (double tx = 0; tx < 10; tx += 2) {
      SimilarityTransform3 tform(1, Eigen::Vector4d(0.2, 0.3, 0.4, qz),
                                 Eigen::Vector3d(tx, 2, 3));

      Eigen::Matrix3x4d proj_matrix2 = tform.Matrix().topLeftCorner<3, 4>();

      for (size_t i = 0; i < points3D.size(); ++i) {
        const Eigen::Vector3d& point3D = points3D[i];
        const Eigen::Vector4d point3D1(point3D(0), point3D(1), point3D(2), 1);
        Eigen::Vector3d point2D1 = proj_matrix1 * point3D1;
        Eigen::Vector3d point2D2 = proj_matrix2 * point3D1;
        point2D1 /= point2D1(2);
        point2D2 /= point2D2(2);

        const Eigen::Vector2d point2D1_N(point2D1(0), point2D1(1));
        const Eigen::Vector2d point2D2_N(point2D2(0), point2D2(1));

        const Eigen::Vector3d tri_point3D = TriangulatePoint(
            proj_matrix1, proj_matrix2, point2D1_N, point2D2_N);

        BOOST_CHECK((point3D - tri_point3D).norm() < 1e-10);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestCalculateTriangulationAngle) {
  const Eigen::Vector3d tvec1(0, 0, 0);
  const Eigen::Vector3d tvec2(0, 1, 0);

  BOOST_CHECK_CLOSE(
      CalculateTriangulationAngle(tvec1, tvec2, Eigen::Vector3d(0, 0, 100)),
      0.009999666687, 1e-8);
  BOOST_CHECK_CLOSE(
      CalculateTriangulationAngle(tvec1, tvec2, Eigen::Vector3d(0, 0, 50)),
      0.019997333973, 1e-8);
  BOOST_CHECK_CLOSE(CalculateTriangulationAngles(
                        tvec1, tvec2, {Eigen::Vector3d(0, 0, 50)})[0],
                    0.019997333973, 1e-8);
}



BOOST_AUTO_TEST_CASE(TestTriangulateLine) {
  std::vector<Eigen::Vector3d> points3D(6);
  points3D[0] = Eigen::Vector3d(0, 0.1, 0.1);
  points3D[1] = Eigen::Vector3d(0.2, 1, 3);
  points3D[2] = Eigen::Vector3d(0, 1, 2);
  points3D[3] = Eigen::Vector3d(0.01, 0.2, 3);
  points3D[4] = Eigen::Vector3d(-1, 0.1, 1);
  points3D[5] = Eigen::Vector3d(0.1, 0.1, 0.2);

  Eigen::Matrix3x4d proj_matrix1 = Eigen::MatrixXd::Identity(3, 4);

  for (double qz = 0.2; qz < 1; qz += 0.2) {
    for (double tx = 2; tx < 10; tx += 2) {
      SimilarityTransform3 tform(1, Eigen::Vector4d(0.2, 0.3, 0.4, qz),
                                 Eigen::Vector3d(tx, 2, 3));

      Eigen::Matrix3x4d proj_matrix2 = tform.Matrix().topLeftCorner<3, 4>();

      for (size_t i = 0; i < points3D.size()-1; ++i) {
        const Eigen::Vector3d& point3D_1 = points3D[i];
        const Eigen::Vector3d& point3D_2 = points3D[i+1];
        
        std::vector<Eigen::Vector2d> points2D_1(2);
        std::vector<Eigen::Vector2d> points2D_2(2);

        points2D_1[0] = (proj_matrix1 * point3D_1.homogeneous()).hnormalized();
        points2D_1[1] = (proj_matrix2 * point3D_1.homogeneous()).hnormalized();

        points2D_2[0] = (proj_matrix1 * point3D_2.homogeneous()).hnormalized();
        points2D_2[1] = (proj_matrix2 * point3D_2.homogeneous()).hnormalized();
        
        const std::pair<Eigen::Vector3d, Eigen::Vector3d> tri = 
            TriangulateMultiViewLine({proj_matrix1, proj_matrix2},points2D_1, points2D_2);

        // We cant be sure which point ends up as the first or second here
        const double err1 = (tri.first - point3D_1).norm() + (tri.second - point3D_2).norm();
        const double err2 = (tri.second - point3D_1).norm() + (tri.first - point3D_2).norm();
              

        BOOST_CHECK(std::min(err1, err2) < 1e-10);
      }
    }
  }
}


BOOST_AUTO_TEST_CASE(TestTriangulateLineCorrectExtent) {
  std::vector<Eigen::Vector3d> points3D(6);  
  points3D[0] = Eigen::Vector3d(0, 0.1, 0.1);
  points3D[1] = Eigen::Vector3d(0.2, 1, 3);
  points3D[2] = Eigen::Vector3d(0, 1, 2);
  points3D[3] = Eigen::Vector3d(0.01, 0.2, 3);
  points3D[4] = Eigen::Vector3d(-1, 0.1, 1);
  points3D[5] = Eigen::Vector3d(0.1, 0.1, 0.2);

  // interpolation is
  //    X = (1-t) * X0 + t * X1
  // For each test the maximal extent should be [X0, X1]
  // so one segment needs to have 0 and one 1

  std::vector<std::pair<double,double>> ivals1 = {
    {0.5, 1.0},  {0.25, 0.5}, {1.0, 0.5}, {1.0, 0.9}
  };
  std::vector<std::pair<double,double>> ivals2 = {
    {0.0, 0.5},  {0.0, 1.0},  {0.25, 0.0}, {0.0, 0.15}
  };

  Eigen::Matrix3x4d proj_matrix1 = Eigen::MatrixXd::Identity(3, 4);

  SimilarityTransform3 tform0(1, Eigen::Vector4d(0.1, 0.01, 0.01, 0.02),
                                 Eigen::Vector3d(0.1, 0.001, -0.1));

  proj_matrix1 = tform0.Matrix().topLeftCorner<3,4>();

  for (double qz = 0.2; qz < 1; qz += 0.2) {
    for (double tx = 2; tx < 10; tx += 2) {
      SimilarityTransform3 tform(1, Eigen::Vector4d(0.2, 0.3, 0.4, qz),
                                 Eigen::Vector3d(tx, 2, 3));

      Eigen::Matrix3x4d proj_matrix2 = tform.Matrix().topLeftCorner<3, 4>();

      for (size_t i = 0; i < points3D.size()-1; ++i) {
        const Eigen::Vector3d& X0 = points3D[i];
        const Eigen::Vector3d& X1 = points3D[i+1];
        for(size_t k = 0; k < ivals1.size(); ++k) {

          // End points in the first image
          Eigen::Vector3d X0_1 = (1 - ivals1[k].first) * X0 + ivals1[k].first * X1;
          Eigen::Vector3d X1_1 = (1 - ivals1[k].second) * X0 + ivals1[k].second * X1;

          // End points in the second image
          Eigen::Vector3d X0_2 = (1 - ivals2[k].first) * X0 + ivals2[k].first * X1;          
          Eigen::Vector3d X1_2 = (1 - ivals2[k].second) * X0 + ivals2[k].second * X1;
          
          std::vector<Eigen::Vector2d> points2D_1(2);
          std::vector<Eigen::Vector2d> points2D_2(2);

          points2D_1[0] = (proj_matrix1 * X0_1.homogeneous()).hnormalized();
          points2D_2[0] = (proj_matrix1 * X1_1.homogeneous()).hnormalized();

          points2D_1[1] = (proj_matrix2 * X0_2.homogeneous()).hnormalized();          
          points2D_2[1] = (proj_matrix2 * X1_2.homogeneous()).hnormalized();
          
          const std::pair<Eigen::Vector3d, Eigen::Vector3d> tri = 
              TriangulateMultiViewLine({proj_matrix1, proj_matrix2},points2D_1, points2D_2);

          // We cant be sure which point ends up as the first or second here
          const double err1 = (tri.first - X0).norm() + (tri.second - X1).norm();
          const double err2 = (tri.second - X0).norm() + (tri.first - X1).norm();
          
          BOOST_CHECK(std::min(err1, err2) < 1e-8);
        }
      }
    }
  }
}