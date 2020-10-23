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

#include "base/line3d.h"
#include <Eigen/Dense>

namespace colmap {

Line3D::Line3D() : xyz1_(0.0, 0.0, 0.0), xyz2_(0.0, 0.0, 0.0), color_(0, 0, 0), error_(-1.0) {}

std::pair<Eigen::Vector3d, Eigen::Vector3d> FitLineToPoints(const std::vector<Eigen::Vector3d> &points) {

    Eigen::Vector3d m;
    Eigen::Matrix3d tensor;

    m.setZero();
    tensor.setZero();
    for(const Eigen::Vector3d &pt : points) {
        m += pt;
    }
    m /= points.size();

    // Compute centered structure tensor
    for(const Eigen::Vector3d &pt : points) {
        tensor += (pt - m) * (pt - m).transpose();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(tensor);    
    // the eigenvalues are sorted in increasing order
    Eigen::Vector3d v = eig.eigenvectors().col(2);

    Eigen::Vector3d X0 = m;
    Eigen::Vector3d X1 = m + v;


    double min_p = std::numeric_limits<double>::max();
    double max_p = std::numeric_limits<double>::min();
    for(const Eigen::Vector3d &pt : points) {
        double p = ProjectPointToLineParameter(pt, {X0, X1});
        min_p = std::min(min_p, p);
        max_p = std::max(max_p, p);
    }

    return {min_p * X0 + (1-min_p) * X1, min_p * X0 + (1-max_p) * X1};
}

Eigen::Vector3d ProjectPointToLine(const Eigen::Vector3d &point, const std::pair<Eigen::Vector3d, Eigen::Vector3d> &line) {    
    double t = ProjectPointToLineParameter(point, line);
    return t * line.first + (1-t) * line.second;
}

double ProjectPointToLineParameter(const Eigen::Vector3d &point, const std::pair<Eigen::Vector3d, Eigen::Vector3d> &line) {
    // t * X1 + (1-t) * X2 = point
    // (X1-X2) * t = point - X2
    Eigen::Vector3d V = line.first - line.second;
    return V.dot(point - line.second) / V.squaredNorm();    
}

}  // namespace colmap
