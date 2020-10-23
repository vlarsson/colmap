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

#ifndef COLMAP_SRC_BASE_LINE3D_H_
#define COLMAP_SRC_BASE_LINE3D_H_

#include <vector>
#include <utility>
#include <Eigen/Core>

#include "base/track.h"
#include "util/logging.h"
#include "util/types.h"

namespace colmap {

// The 3D line class represents a triangulated line
// Note that while this represents a potentially infinite line in 3D
// we represent it in terms of 2 3D points. These points are coincide with
// the maximally observed extent of the line (i.e. union of all 2d line segments)
class Line3D {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Line3D();

  // The point coordinate in world space.  
  inline const Eigen::Vector3d& XYZ1() const;
  inline const Eigen::Vector3d& XYZ2() const;
  inline Eigen::Vector3d& XYZ1();
  inline Eigen::Vector3d& XYZ2();
  inline double XYZ1(const size_t idx) const;
  inline double XYZ2(const size_t idx) const;
  inline double& XYZ1(const size_t idx);
  inline double& XYZ2(const size_t idx);
  inline double X1() const;
  inline double Y1() const;
  inline double Z1() const;
  inline double X2() const;
  inline double Y2() const;
  inline double Z2() const;
  inline void SetXYZ(const Eigen::Vector3d& xyz1, const Eigen::Vector3d& xyz2);
  
  // The RGB color of the point.
  inline const Eigen::Vector3ub& Color() const;
  inline Eigen::Vector3ub& Color();
  inline uint8_t Color(const size_t idx) const;
  inline uint8_t& Color(const size_t idx);
  inline void SetColor(const Eigen::Vector3ub& color);

  // The mean reprojection error in image space.
  inline double Error() const;
  inline bool HasError() const;
  inline void SetError(const double error);

  inline const class Track& Track() const;
  inline class Track& Track();
  inline void SetTrack(const class Track& track);

 private:
  // Two points representing the line
  Eigen::Vector3d xyz1_;
  Eigen::Vector3d xyz2_;

  // The color of the point in the range [0, 255].
  Eigen::Vector3ub color_;

  // The mean reprojection error in pixels.
  double error_;

  // The track of the point as a list of image observations.
  class Track track_;
};

// Fits a line to a vector of points. The returned endpoints matches the extent of the original points
std::pair<Eigen::Vector3d, Eigen::Vector3d> FitLineToPoints(const std::vector<Eigen::Vector3d> &points);

// Projects a point to a line
Eigen::Vector3d ProjectPointToLine(const Eigen::Vector3d &point, const std::pair<Eigen::Vector3d, Eigen::Vector3d> &line);

// Finds the parameter for the convex combination of the end points, i.e.
//        t * X1 + (1-t) * X2 = projection
// thus if t \in [0,1], we have a point inside the current line segment
double ProjectPointToLineParameter(const Eigen::Vector3d &point, const std::pair<Eigen::Vector3d, Eigen::Vector3d> &line);


////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

const Eigen::Vector3d& Line3D::XYZ1() const { return xyz1_; }
const Eigen::Vector3d& Line3D::XYZ2() const { return xyz2_; }

Eigen::Vector3d& Line3D::XYZ1() { return xyz1_; }
Eigen::Vector3d& Line3D::XYZ2() { return xyz2_; }

double Line3D::XYZ1(const size_t idx) const { return xyz1_(idx); }
double Line3D::XYZ2(const size_t idx) const { return xyz2_(idx); }

double& Line3D::XYZ1(const size_t idx) { return xyz1_(idx); }
double& Line3D::XYZ2(const size_t idx) { return xyz2_(idx); }

double Line3D::X1() const { return xyz1_.x(); }
double Line3D::Y1() const { return xyz1_.y(); }
double Line3D::Z1() const { return xyz1_.z(); }

double Line3D::X2() const { return xyz2_.x(); }
double Line3D::Y2() const { return xyz2_.y(); }
double Line3D::Z2() const { return xyz2_.z(); }

void Line3D::SetXYZ(const Eigen::Vector3d& xyz1, const Eigen::Vector3d& xyz2) {
   xyz1_ = xyz1;
   xyz2_ = xyz2;
}


const Eigen::Vector3ub& Line3D::Color() const { return color_; }

Eigen::Vector3ub& Line3D::Color() { return color_; }

uint8_t Line3D::Color(const size_t idx) const { return color_(idx); }

uint8_t& Line3D::Color(const size_t idx) { return color_(idx); }

void Line3D::SetColor(const Eigen::Vector3ub& color) { color_ = color; }

double Line3D::Error() const { return error_; }

bool Line3D::HasError() const { return error_ != -1.0; }

void Line3D::SetError(const double error) { error_ = error; }

const class Track& Line3D::Track() const { return track_; }

class Track& Line3D::Track() {
  return track_;
}

void Line3D::SetTrack(const class Track& track) { track_ = track; }

}  // namespace colmap

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(colmap::Line3D)

#endif  // COLMAP_SRC_BASE_LINE3D_H_
