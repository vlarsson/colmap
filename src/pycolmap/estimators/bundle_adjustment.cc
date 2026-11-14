#include "colmap/estimators/bundle_adjustment.h"

#include "pycolmap/helpers.h"
#include "pycolmap/pybind11_extension.h"
#include "pycolmap/utils.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace colmap;
using namespace pybind11::literals;
namespace py = pybind11;

void BindBundleAdjuster(py::module& m) {
  IsPyceresAvailable();  // Try to import pyceres to populate the docstrings.

  using BACfg = BundleAdjustmentConfig;
  py::class_<BACfg> PyBundleAdjustmentConfig(m, "BundleAdjustmentConfig");
  PyBundleAdjustmentConfig.def(py::init<>())
      .def("num_images", &BACfg::NumImages)
      .def("num_points", &BACfg::NumPoints)
      .def("num_constant_cam_intrinsics", &BACfg::NumConstantCamIntrinsics)
      .def("num_constant_cam_poses", &BACfg::NumConstantCamPoses)
      .def("num_constant_cam_positions", &BACfg::NumConstantCamPositions)
      .def("num_variable_points", &BACfg::NumVariablePoints)
      .def("num_constant_points", &BACfg::NumConstantPoints)
      .def("num_residuals", &BACfg::NumResiduals, "reconstruction"_a)
      .def("add_image", &BACfg::AddImage, "image_id"_a)
      .def("has_image", &BACfg::HasImage, "image_id"_a)
      .def("remove_image", &BACfg::RemoveImage, "image_id"_a)
      .def("set_constant_cam_intrinsics",
           &BACfg::SetConstantCamIntrinsics,
           "camera_id"_a)
      .def("set_variable_cam_intrinsics",
           &BACfg::SetVariableCamIntrinsics,
           "camera_id"_a)
      .def("has_constant_cam_intrinsics",
           &BACfg::HasConstantCamIntrinsics,
           "camera_id"_a)
      .def("set_constant_cam_pose", &BACfg::SetConstantCamPose, "image_id"_a)
      .def("set_variable_cam_pose", &BACfg::SetVariableCamPose, "image_id"_a)
      .def("has_constant_cam_pose", &BACfg::HasConstantCamPose, "image_id"_a)
      .def("set_constant_cam_positions",
           &BACfg::SetConstantCamPositions,
           "image_id"_a,
           "idxs"_a)
      .def("remove_variable_cam_positions",
           &BACfg::RemoveConstantCamPositions,
           "image_id"_a)
      .def("has_constant_cam_positions",
           &BACfg::HasConstantCamPositions,
           "image_id"_a)
      .def("add_variable_point", &BACfg::AddVariablePoint, "point3D_id"_a)
      .def("add_constant_point", &BACfg::AddConstantPoint, "point3D_id"_a)
      .def("has_point", &BACfg::HasPoint, "point3D_id"_a)
      .def("has_variable_point", &BACfg::HasVariablePoint, "point3D_id"_a)
      .def("has_constant_point", &BACfg::HasConstantPoint, "point3D_id"_a)
      .def("remove_variable_point", &BACfg::RemoveVariablePoint, "point3D_id"_a)
      .def("remove_constant_point", &BACfg::RemoveConstantPoint, "point3D_id"_a)
      .def_property_readonly("constant_intrinsics", &BACfg::ConstantIntrinsics)
      .def_property_readonly("image_ids", &BACfg::Images)
      .def_property_readonly("variable_point3D_ids", &BACfg::VariablePoints)
      .def_property_readonly("constant_point3D_ids", &BACfg::ConstantPoints)
      .def_property_readonly("constant_cam_poses", &BACfg::ConstantCamPoses)
      .def(
          "constant_cam_positions", &BACfg::ConstantCamPositions, "image_id"_a);
  MakeDataclass(PyBundleAdjustmentConfig);

  using BAOpts = BundleAdjustmentOptions;
  auto PyBALossFunctionType =
      py::enum_<BAOpts::LossFunctionType>(m, "LossFunctionType")
          .value("TRIVIAL", BAOpts::LossFunctionType::TRIVIAL)
          .value("SOFT_L1", BAOpts::LossFunctionType::SOFT_L1)
          .value("CAUCHY", BAOpts::LossFunctionType::CAUCHY);
  AddStringToEnumConstructor(PyBALossFunctionType);

  auto PyBundleAdjustmentOptions =
      py::class_<BAOpts>(m, "BundleAdjustmentOptions")
          .def(py::init<>())
          .def("create_loss_function", &BAOpts::CreateLossFunction)
          .def("create_solver_options",
               &BAOpts::CreateSolverOptions,
               "config"_a,
               "problem"_a)
          .def_readwrite("loss_function_type",
                         &BAOpts::loss_function_type,
                         "Loss function types: Trivial (non-robust) and Cauchy "
                         "(robust) loss.")
          .def_readwrite("loss_function_scale",
                         &BAOpts::loss_function_scale,
                         "Scaling factor determines residual at which "
                         "robustification takes place.")
          .def_readwrite("loss_function_magnitude",
                         &BAOpts::loss_function_magnitude,
                         "Magnitude determines the scaling of the residual.")
          .def_readwrite("refine_focal_length",
                         &BAOpts::refine_focal_length,
                         "Whether to refine the focal length parameter group.")
          .def_readwrite(
              "refine_principal_point",
              &BAOpts::refine_principal_point,
              "Whether to refine the principal point parameter group.")
          .def_readwrite("refine_extra_params",
                         &BAOpts::refine_extra_params,
                         "Whether to refine the extra parameter group.")
          .def_readwrite("refine_extrinsics",
                         &BAOpts::refine_extrinsics,
                         "Whether to refine the extrinsic parameter group.")
          .def_readwrite("print_summary",
                         &BAOpts::print_summary,
                         "Whether to print a final summary.")
          .def_readwrite("use_gpu",
                         &BAOpts::use_gpu,
                         "Whether to use Ceres' CUDA linear algebra library, "
                         "if available.")
          .def_readwrite("gpu_index",
                         &BAOpts::gpu_index,
                         "Which GPU to use for solving the problem.")
          .def_readwrite(
              "min_num_residuals_for_cpu_multi_threading",
              &BAOpts::min_num_residuals_for_cpu_multi_threading,
              "Minimum number of residuals to enable multi-threading. Note "
              "that single-threaded is typically better for small bundle "
              "adjustment problems due to the overhead of threading.")
          .def_readwrite("min_num_images_gpu_solver",
                         &BAOpts::min_num_images_gpu_solver,
                         "Minimum number of images to use the GPU solver.")
          .def_readwrite("max_num_images_direct_dense_cpu_solver",
                         &BAOpts::max_num_images_direct_dense_cpu_solver,
                         "Threshold to switch between direct, sparse, and "
                         "iterative solvers.")
          .def_readwrite("max_num_images_direct_sparse_cpu_solver",
                         &BAOpts::max_num_images_direct_sparse_cpu_solver,
                         "Threshold to switch between direct, sparse, and "
                         "iterative solvers.")
          .def_readwrite("max_num_images_direct_dense_gpu_solver",
                         &BAOpts::max_num_images_direct_dense_gpu_solver,
                         "Threshold to switch between direct, sparse, and "
                         "iterative solvers.")
          .def_readwrite("max_num_images_direct_sparse_gpu_solver",
                         &BAOpts::max_num_images_direct_sparse_gpu_solver,
                         "Threshold to switch between direct, sparse, and "
                         "iterative solvers.")
          .def_readwrite("solver_options",
                         &BAOpts::solver_options,
                         "Options for the Ceres solver. Using this member "
                         "requires having PyCeres installed.");
  MakeDataclass(PyBundleAdjustmentOptions);

  using PosePriorBAOpts = PosePriorBundleAdjustmentOptions;
  auto PyPosePriorBundleAdjustmentOptions =
      py::class_<PosePriorBAOpts>(m, "PosePriorBundleAdjustmentOptions")
          .def(py::init<>())
          .def_readwrite("use_robust_loss_on_prior_position",
                         &PosePriorBAOpts::use_robust_loss_on_prior_position,
                         "Whether to use a robust loss on prior locations.")
          .def_readwrite("prior_position_loss_scale",
                         &PosePriorBAOpts::prior_position_loss_scale,
                         "Threshold on the residual for the robust loss (chi2 "
                         "for 3DOF at 95% = 7.815).")
          .def_readwrite("ransac_max_error",
                         &PosePriorBAOpts::ransac_max_error,
                         "Maximum RANSAC error for Sim3 alignment.");
  MakeDataclass(PyPosePriorBundleAdjustmentOptions);

  class PyBundleAdjuster : public BundleAdjuster {
   public:
    PyBundleAdjuster(BundleAdjustmentOptions options,
                     BundleAdjustmentConfig config)
        : BundleAdjuster(std::move(options), std::move(config)) {}

    ceres::Solver::Summary Solve() override {
      PYBIND11_OVERRIDE_PURE(ceres::Solver::Summary, BundleAdjuster, Solve);
    }

    std::shared_ptr<ceres::Problem>& Problem() override {
      PYBIND11_OVERRIDE_PURE(
          std::shared_ptr<ceres::Problem>&, BundleAdjuster, Problem);
    }
  };

  py::class_<BundleAdjuster, PyBundleAdjuster>(m, "BundleAdjuster")
      .def(py::init<BundleAdjustmentOptions, BundleAdjustmentConfig>(),
           "options"_a,
           "config"_a)
      .def("solve", &BundleAdjuster::Solve)
      .def_property_readonly("problem", &BundleAdjuster::Problem)
      .def_property_readonly("options", &BundleAdjuster::Options)
      .def_property_readonly("config", &BundleAdjuster::Config);

  m.def("create_default_bundle_adjuster",
        CreateDefaultBundleAdjuster,
        "options"_a,
        "config"_a,
        "reconstruction"_a);

  m.def("create_pose_prior_bundle_adjuster",
        CreatePosePriorBundleAdjuster,
        "options"_a,
        "prior_options"_a,
        "config"_a,
        "pose_priors"_a,
        "reconstruction"_a);
        
  m.def("create_depth_bundle_adjuster",
        [](ceres::Problem* problem,
            image_t image_id,
            const std::vector<point3D_t>& point3D_ids,
            const std::vector<double>& depths,
            const std::vector<double>& loss_magnitudes,
            const std::vector<double>& loss_params,
            BAOpts::LossFunctionType loss_type,
            py::array_t<double> shift_scale,
            Reconstruction& reconstruction,
            bool logloss,
            bool fix_shift,
            bool fix_scale) {

            auto buf = shift_scale.request();
            if (buf.ndim != 1 || buf.shape[0] != 2)
                throw std::runtime_error("shift_scale must have exactly 2 elements.");

            double* shift_scale_ptr = static_cast<double*>(buf.ptr);

            DepthPriorBundleAdjuster(problem, image_id, point3D_ids, depths,
                                            loss_magnitudes, loss_params, loss_type,
                                            shift_scale_ptr, reconstruction,
                                            logloss, fix_shift, fix_scale);
        },
        "problem"_a,
        "image_id"_a,
        "point3D_ids"_a,
        "depths"_a,
        "loss_magnitudes"_a,
        "loss_params"_a,
        "loss_type"_a,
        "shift_scale"_a,
        "reconstruction"_a,
        "logloss"_a = false,
        "fix_shift"_a = false,
        "fix_scale"_a = false);
  m.def("create_depth_basis_bundle_adjuster",
        [](ceres::Problem* problem,
           image_t image_id,
           const std::vector<point3D_t>& point3D_ids,
           const std::vector<double>& depths,
           py::array_t<double, py::array::c_style | py::array::forcecast> basis,
           const std::vector<double>& loss_magnitudes,
           const std::vector<double>& loss_params,
           BAOpts::LossFunctionType loss_type,
           py::array_t<double> theta,
           Reconstruction& reconstruction) {
          auto basis_buf = basis.request();
          if (basis_buf.ndim != 2 ||
              basis_buf.shape[0] != static_cast<ssize_t>(point3D_ids.size()))
            throw std::runtime_error("basis must have shape (num_points, rank).");
          const int rank = static_cast<int>(basis_buf.shape[1]);
          auto theta_buf = theta.request();
          if (theta_buf.ndim != 1 || theta_buf.shape[0] != rank + 1)
            throw std::runtime_error("theta must have rank + 1 elements.");
          if (!(theta.flags() & py::array::c_style) || !theta.writeable())
            throw std::runtime_error("theta must be a writeable contiguous float64 array.");
          const double* basis_ptr = static_cast<const double*>(basis_buf.ptr);
          std::vector<double> basis_vec(basis_ptr,
                                        basis_ptr + basis_buf.shape[0] * rank);
          DepthBasisBundleAdjuster(problem, image_id, point3D_ids, depths,
                                   basis_vec, rank, loss_magnitudes,
                                   loss_params, loss_type,
                                   static_cast<double*>(theta_buf.ptr),
                                   reconstruction);
        },
        "problem"_a,
        "image_id"_a,
        "point3D_ids"_a,
        "depths"_a,
        "basis"_a,
        "loss_magnitudes"_a,
        "loss_params"_a,
        "loss_type"_a,
        "theta"_a,
        "reconstruction"_a);
}
