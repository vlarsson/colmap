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

#ifndef COLMAP_SRC_BASE_LINE2D_H_
#define COLMAP_SRC_BASE_LINE2D_H_

#include <Eigen/Core>

#include "util/alignment.h"
#include "util/types.h"

namespace colmap {

// 2D Line class corresponds to a line segment in an image. It may or may not have a
// corresponding 3D line if it is part of a triangulated track.
class Line2D {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Line2D();

  // The coordinate in image space in pixels.
  inline const std::pair<Eigen::Vector2d, Eigen::Vector2d> XY() const;
  inline const Eigen::Vector2d& XY1() const;
  inline const Eigen::Vector2d& XY2() const;
  inline double X1() const;
  inline double Y1() const;
  inline double X2() const;
  inline double Y2() const;
  
  inline void SetXY(const Eigen::Vector2d& xy1, const Eigen::Vector2d& xy2);

  // The identifier of the observed 3D point. If the image point does not
  // observe a 3D point, the identifier is `kInvalidPoint3Did`.
  inline point3D_t Line3DId() const;
  inline bool HasLine3D() const;
  inline void SetLine3DId(const point3D_t point3D_id);

 private:
  // The image coordinates in pixels, starting at upper left corner with 0.
  Eigen::Vector2d xy1_, xy2_;

  // The identifier of the 3D line. If the 2D line is not part of a 3D line
  // track the identifier is `kInvalidPoint3DId` and `HasLine3D() = false`.
  point3D_t line3D_id_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

const std::pair<Eigen::Vector2d, Eigen::Vector2d> Line2D::XY() const { return std::make_pair(xy1_,xy2_); }

const Eigen::Vector2d& Line2D::XY1() const { return xy1_; }
const Eigen::Vector2d& Line2D::XY2() const { return xy2_; }

double Line2D::X1() const { return xy1_.x(); }
double Line2D::Y1() const { return xy1_.y(); }
double Line2D::X2() const { return xy2_.x(); }
double Line2D::Y2() const { return xy2_.y(); }

void Line2D::SetXY(const Eigen::Vector2d& xy1, const Eigen::Vector2d& xy2) { xy1_ = xy1; xy2_ = xy2; }

point3D_t Line2D::Line3DId() const { return line3D_id_; }

bool Line2D::HasLine3D() const { return line3D_id_ != kInvalidPoint3DId; }

void Line2D::SetLine3DId(const point3D_t line3D_id) {
  line3D_id_ = line3D_id;
}

}  // namespace colmap

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(colmap::Line2D)

#endif  // COLMAP_SRC_BASE_LINE2D_H_
