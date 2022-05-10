/* This header file contains the utilities functions and identities
** for the 3D Lie Group. The utilities mainly include
** identities and operations of SO(3) and SE(3) vector and matrices.
**
**
** Author: MT
** Creation Date: 2022-May-05
** Previous Edit: 2022-May-08
*/


#pragma once
#ifndef UTILS_INCLUDE_H
#define UTILS_INCLUDE_H

// include 
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include <math.h>
#include <string.h>
#include <algorithm>

using namespace std;

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;


////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Below are SO(3) operations ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T>
inline void unitQuaternionToAngleAxis(const T *unit_quaternion, T* angle_axis) {
    // Note the input unit quaternion should have a format of:
    // [qw, qx, qy, qz], that is the scalar part should be the first element
    // The output angle_axis will be a 1D array with 3 elements with the format of
    // [p1, p2, p3], the three values for the angle-axis rotation, where, norm is the angle
    

    // // Naive implementation - expensive and unstable
    // const T &half_theta = arccos(unit_quaternion[0]);
    // const T &sin_half_theta = sin(half_theta);
    // const T &theta = 2.0 * half_theta;
    // angle_axis[0] = unit_quaternion[1] / sin_half_theta * theta; // v1
    // angle_axis[1] = unit_quaternion[2] / sin_half_theta * theta; // v2
    // angle_axis[2] = unit_quaternion[3] / sin_half_theta * theta; // v3


    // Implementation reference: Google's Ceres-solver - stable
    const T &q1 = unit_quaternion[1];
    const T &q2 = unit_quaternion[2];
    const T &q3 = unit_quaternion[3];
    const T sin_squared_half_theta = q1 * q1 + q2 * q2 + q3 * q3;

    // For quaternions representing non-zero rotation, the conversion is numerically stable
    if (sin_squared_half_theta > T(std::numeric_limits<double>::epsilon())) {
        const T sin_half_theta = sqrt(sin_squared_half_theta);
        const T &cos_half_theta = unit_quaternion[0];
        
        // If cos_theta is negative, theta is greater than pi/2, which
        // means that angle for the angle-axis vector (i.e., 2 * theta)
        // would be greater than pi

        // syntax: (<condition> ? <value if true> : <value if false>)
        const T theta = T(2.0) * ((cos_half_theta < 0.0)
                                      ? atan2(-sin_half_theta, -cos_half_theta)
                                      : atan2(sin_half_theta, cos_half_theta));

        const T theta_by_norm = theta / sin_half_theta;

        angle_axis[0] = q1 * theta_by_norm;
        angle_axis[1] = q2 * theta_by_norm;
        angle_axis[2] = q3 * theta_by_norm;
    } else {
        // For zero rotation, sqrt() will produce NaN in derivative since the argument is zero
        // By approximating the derivative with a Taylor serise and truncating at the first term
        // (first order approximation), the value of the first derivative will be computed
        // correctly.
        const T k(2.0);
        angle_axis[0] = q1 * k;
        angle_axis[1] = q2 * k;
        angle_axis[2] = q3 * k;
    }
} // end of unitQuaternionToAngleAxis()


template<typename T>
inline void angleAxisToUnitQuaternion(const T* angle_axis, T* unit_quaternion){
    // The input angle_axis should be a 1D array with three elements with the format of
    // [p1, p2, p3], the three values for the angle-axis rotation, where, norm is the angle    
    // Note the output unit quaternion will have a format of:
    // [qw, qx, qy, qz], that is the scalar part should be the first element
    const T &p1 = angle_axis[0];
    const T &p2 = angle_axis[1];
    const T &p3 = angle_axis[2];
    const T rot_angle_squared = p1*p1 + p2*p2 + p3*p3;
    const T rot_angle = sqrt(rot_angle_squared); // norm of the angle-axis vector, note: sqrt(0.0)=0.0, sqrt(-0.0)=-0.0
    
    if (rot_angle_squared > T(0.0) && rot_angle > T(std::numeric_limits<double>::epsilon())) 
    { // when rotation angle is not zero
        const T half_angle = rot_angle * T(0.5);
        const T sin_half_angle_over_angle = sin(half_angle) / rot_angle; // sin(phi) or sin(theta)
        unit_quaternion[0] = cos(half_angle);                // qw
        unit_quaternion[1] = p1 * sin_half_angle_over_angle; // qx
        unit_quaternion[2] = p2 * sin_half_angle_over_angle; // qy
        unit_quaternion[3] = p3 * sin_half_angle_over_angle; // qz
    } else { // reference: "rotation.h" from Google's Ceres solver
        // By approximating with a Taylor series and truncating at one term
        const T sin_half_angle_over_angle = T(0.5); // note: lim(sin(x/2) / x) at x->0 is 0.5
        unit_quaternion[0] = T(1.0);                         // qw
        unit_quaternion[1] = p1 * sin_half_angle_over_angle; // qx
        unit_quaternion[2] = p2 * sin_half_angle_over_angle; // qy
        unit_quaternion[3] = p3 * sin_half_angle_over_angle; // qz
    }
} // end of angleAxisToUnitQuaternion()



// The trace < 0 part refers to Google's Ceres-solver's "rotation.h" file.
// According to Ceres, the algorithm further referenced "Quaternion Calculus and Fast Animation",
// Ken Shoemake, 1987, SIGGRAPH course notes
template<typename T>
inline void rotationMatrixArrayToUnitQuaternion(const T rotation_matrix_vec[9], T unit_quaternion[4]) {
    // The input rotation matrix should be a 9x1 vector representing the raveled form of a 
    // row major 3x3 rotation matrix on SO(3). Note: Row major order!
    // Note the output unit quaternion should have a format of:
    // [qw, qx, qy, qz], that is the scalar part should be the first element
    const T &r11 = rotation_matrix_vec[0]; 
    const T &r12 = rotation_matrix_vec[1]; 
    const T &r13 = rotation_matrix_vec[2]; // end of first row
    const T &r21 = rotation_matrix_vec[3];
    const T &r22 = rotation_matrix_vec[4];
    const T &r23 = rotation_matrix_vec[5]; // end of second row
    const T &r31 = rotation_matrix_vec[6];
    const T &r32 = rotation_matrix_vec[7];
    const T &r33 = rotation_matrix_vec[8]; // end of third row

    T rot_trace = r11 + r22 + r33; // the trace of the rotation matrix
    if (rot_trace >= 0.0) {
        T sqrt_trace_plus_1 = sqrt(rot_trace + 1); // an intermediate term to be used later
        T twice_sqrt_trace_plus_1 = T(2.0) * sqrt_trace_plus_1;

        unit_quaternion[0] = sqrt_trace_plus_1 * 0.5; // the real part of the unit quaternion, qw
        unit_quaternion[1] = (r32 - r23) / twice_sqrt_trace_plus_1; // qx
        unit_quaternion[2] = (r13 - r31) / twice_sqrt_trace_plus_1; // qy
        unit_quaternion[3] = (r21 - r12) / twice_sqrt_trace_plus_1; // qz

    } else { // trace is smaller than 0.0
        const T R_arr[3][3] = {{r11,r12,r13}, {r21,r22,r23}, {r31,r32,r33}}; // Rotation mat as 2D array, row major
        int i = 0;
        if (r22 > r11) {
            i = 1;
        }

        if (r33 > R_arr[i][i]) {
            i = 2;
        }

        const int j = (i + 1) % 3;
        const int k = (j + 1) % 3;
        T t = sqrt(R_arr[i][i] - R_arr[j][j] - R_arr[k][k] + T(1.0));
        unit_quaternion[i + 1] = T(0.5) * t;
        t = T(0.5) / t;
        unit_quaternion[0] = (R_arr[k][j] - R_arr[j][k]) * t; // real part of unit quaternion
        unit_quaternion[j + 1] = (R_arr[j][i] + R_arr[i][j]) * t;
        unit_quaternion[k + 1] = (R_arr[k][i] + R_arr[i][k]) * t;
    }
} // end of rotationMatrixArrayToUnitQuaternion()


// Not using a MatrixAdapter to accept matrix input
// To use MatrixAdapter to accept matrix input, see:
// https://github.com/ceres-solver/ceres-solver/blob/master/include/ceres/rotation.h
// Here uses a 9D dimensional array to store the output hat matrix of input 3D vector
template<typename T>
inline void vecHat(const T vec[3], T vec_hat[9]) {
    // Note the input should be a 3D vector
    // The output matrix will be stored in a 1D array with 9 elements in ROW MAJOR order.
    // Again, the output is in ROW MAJOR order.
    vec_hat[1] = -vec[2];
    vec_hat[2] =  vec[1];
    vec_hat[3] =  vec[2];
    vec_hat[5] = -vec[0];
    vec_hat[6] = -vec[1];
    vec_hat[7] =  vec[0];
    vec_hat[0] = 0.0;
    vec_hat[4] = 0.0;
    vec_hat[8] = 0.0;
} // end of vecHat()


// the expPhiHat() funciton performs the exponential map over the Phi_hat matrix
// to obtain the corresponding rotation matrix for phi
void expPhiHat(const Eigen::Matrix3d &phi_hat_mat, Eigen::Matrix3d &R_mat) {
    // The input to the function should be 3x3 Eigen matrix with double precision values
    // The output is the corresponding 3x3 Eigen matrix to the input phi hat matrix
    Eigen::Vector3d rot_vec(phi_hat_mat.coeff(2,1), phi_hat_mat.coeff(0,2), phi_hat_mat.coeff(1,0)); // extract the rotation vector
    double rot_angle = sqrt(rot_vec.transpose() * rot_vec); // this gives the rotation angle in [radian]
    if (rot_angle > std::numeric_limits<double>::epsilon()) {
        // This part only performs if the rotation angle is not (approaching) zero
        Eigen::Vector3d rot_axis = rot_vec.normalized();
        Eigen::Matrix3d axis_hat;
        axis_hat <<          0.0, -rot_axis[2],  rot_axis[1], 
                     rot_axis[2],          0.0, -rot_axis[0],
                    -rot_axis[1],  rot_axis[0],          0.0;
        double cos_angle = cos(rot_angle);
        R_mat = cos_angle * Eigen::Matrix3d::Identity() 
                + (1 - cos_angle) * rot_axis * rot_axis.transpose()
                + sin(rot_angle) * axis_hat;
    } else {
        // If the theta is smaller than the numeric limit, using an approximation as in Barfoot's SO(3) table
        R_mat = Eigen::Matrix3d::Identity() + phi_hat_mat;
    }
} // end of expPhiHat()


// the phiToSO3() function first converts the angle-axis vector, phi, to the phi hat matrix,
// then perform the exponential map over the phi_hat matrix to obtain the rotation matrix
void phiToSO3(const Eigen::Vector3d &phi_vec, Eigen::Matrix3d &R_mat) {
    // Note the input phi_vec vector should be an Eigen 3D vector as angle-axis rotation vector
    // The output is the corresponding 3x3 Eigen matrix to the input phi hat matrix.
    double phi_hat[9]; // 1D array with 9 elements in ROW MAJOR order of phi_hat matrix
    vecHat(phi_vec.data(), phi_hat); // obtain the phi_hat matrix stored in row-major order
    expPhiHat(Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(phi_hat), R_mat); // get the rotation mat
} // end of phiToSO3


// the lnVeeToPhi() function performs the logarithm map for a Rotation matrix on SO(3),
// computing a 3 x 1 tangent vector, phi (i.e., an angle-axis representation), from a rotation
void lnVeeToPhi(const Eigen::Matrix3d &R_mat, Eigen::Vector3d &phi) {
    // The input to the function should be a 3x3 Eigen matrix with double precision values.
    // The output is the corresponding 3x1 tangential lie algebra vector, phi.

    // The workflow follow the "rotation.h" file used by Google's Ceres solver: Rot_mat -> unit_quaternion -> angle_axis
    // the rotation matrix will be first raveled into an 1-D array, 
    // then converted to unit quaternion, finally from quaternion to rotation vector
    double R_vec[9], unit_qua[4], phi_vec[3]; // define the necessary
    // to store data row by row, by default Eigen stores things in col major order, so needs to transpose the Rotation mat
    Eigen::Map<Eigen::MatrixXd>(R_vec, R_mat.cols(), R_mat.rows()) = R_mat.transpose();
    rotationMatrixArrayToUnitQuaternion(R_vec, unit_qua); // convert the rotation matrix vector to corresponding unit quaternion
    unitQuaternionToAngleAxis(unit_qua, phi_vec); // convert the unit quaternion into a 3x1 rotation vector
    phi << phi_vec[0], phi_vec[1], phi_vec[2]; // assign values to the Eigen Vector
} // end of lnVeeToPhi


// the leftJacobianSO3() function computes the left Jacobian w.r.t. a SO(3) rotation vector.
void leftJacobianSO3(const Eigen::Vector3d &phi_vec, Eigen::Matrix3d &left_Jac_mat) {
    // The input to the function should be a 3D vector of the angle-axis rotation vector with double precision values
    // The output matrix will be the corresponding Jacobian of the rotation vector in a 3x3 matrix with double precision values
    double squared_norm = phi_vec.transpose() * phi_vec;
    double rot_angle = sqrt(squared_norm); // this is the angle of rotation, scalar phi
    double sin_phi = sin(rot_angle);

    if (sin_phi > std::numeric_limits<double>::epsilon()) {
        // This part is only performed when the rotation angle is greater than the numeric limit to avoid division by zero
        Eigen::Vector3d axis_vec = phi_vec.normalized(); // obtain the rotation axis
        double axis_hat[9];
        vecHat(axis_vec.data(), axis_hat); // obtain the skew-symmetric matrix of the rotation axis

        // calculate several intermediate terms to be used
        double one_over_phi = 1.0 / rot_angle; // 1 / phi
        double sin_phi_by_phi = sin_phi * one_over_phi; // sin(phi) / phi

        // calculate the left Jacobian - following formula in Barfoot's SO(3) table
        left_Jac_mat =  sin_phi_by_phi * Eigen::Matrix3d::Identity()
                      + (1.0 - sin_phi_by_phi) * axis_vec * axis_vec.transpose()
                      + ((1.0-cos(rot_angle)) * one_over_phi) * Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(axis_hat);

    } else { // when the rotation angle is smaller than the numeric limit, e.g., near the identity
        double phi_hat[9];
        vecHat(phi_vec.data(), phi_hat); // obtain the skew-symmetric matrix of the rotation vector
        left_Jac_mat =   Eigen::Matrix3d::Identity() 
                       + 0.5 * Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(phi_hat);
    }
} // end of leftJacobianSO3()


// the invLeftJacobianSO3() function computes the inversion of left Jacobian w.r.t. a SO(3) rotation vector.
void invLeftJacobianSO3(const Eigen::Vector3d &phi_vec, Eigen::Matrix3d &inv_Jac_mat) {
    // The input to the function should be a 3D vector of the angle-axis rotation vector with double precision values
    // The output matrix will be the inversion of a Jacobian of a SO(3) rotation vector in a 3x3 matrix with double precision values
    double squared_norm = phi_vec.transpose() * phi_vec;
    double rot_angle = sqrt(squared_norm); // this is the angle of rotation, scalar phi
    double sin_phi = sin(rot_angle);

    if (sin_phi > std::numeric_limits<double>::epsilon()) {
        // This part is only performed when the rotation angle is greater than the numeric limit to avoid division by zero
        Eigen::Vector3d axis_vec = phi_vec.normalized(); // obtain the rotation axis
        double axis_hat[9];
        vecHat(axis_vec.data(), axis_hat); // obtain the skew-symmetric matrix of the rotation axis

        // calculate several intermediate terms to be used
        double half_phi = rot_angle * 0.5; // half angle
        double cot_half_phi = (1.0 + cos(rot_angle)) / sin_phi; // the half angle formula for cotangent
        double half_phi_times_cot = half_phi * cot_half_phi; // (phi/2) * cot(phi/2)

        // calculate the inversion of the left Jacobian
        inv_Jac_mat =   half_phi_times_cot * Eigen::Matrix3d::Identity()
                      + (1.0 - half_phi_times_cot) * axis_vec * axis_vec.transpose()
                      - half_phi * Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(axis_hat);

    } else { // when the rotation angle is smaller than the numeric limit, e.g., near the identity
        double phi_hat[9];
        vecHat(phi_vec.data(), phi_hat); // obtain the skew-symmetric matrix of the rotation vector
        inv_Jac_mat = Eigen::Matrix3d::Identity() 
                     - 0.5 * Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(phi_hat);
    }
} // end of invLeftJacobianSO3()


// the eigenRotationMatrixToUnitQuaternion() function wraps the two templated function
// to take 3x3 Eigen rotation matrix as input and outputs a Eigen 4D vector of unit quaternion
void eigenRotationMatrixToUnitQuaternion(const Eigen::Matrix3d &R_mat, Vector4d &unit_quaternion) {
    // The input of the function should be a Rotation matrix as Eigen::Matrix3d
    // The output is the corresponding unit quaternion in the order of [qw, qx, qy, qz],
    // where qw is the scalar part
    double R_vec[9], unit_qua[4]; // store the rotation matrix array and unit quaternion vector array
    // to store data row by row, by default Eigen stores things in col major order, so needs to transpose the Rotation mat
    Eigen::Map<Eigen::MatrixXd>(R_vec, R_mat.cols(), R_mat.rows()) = R_mat.transpose();
    rotationMatrixArrayToUnitQuaternion(R_vec, unit_qua); // convert the rotation matrix vector to corresponding unit quaternion
    unit_quaternion << unit_qua[0], unit_qua[1], unit_qua[2], unit_qua[3];
} // end of eigenRotationMatrixToUnitQuaternion()



////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Below are SE(3) operations ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


// Not using a MatrixAdapter to accept matrix input
// To use MatrixAdapter to accept matrix input, see:
// https://github.com/ceres-solver/ceres-solver/blob/master/include/ceres/rotation.h
// Here uses a 1D array with 16 elements to store the output hat matrix of zeta
template<typename T>
inline void zetaHat(const T zeta_vec[6], T zeta_hat[16]) {
    // Note the input zeta_vec vector should have a format of:
    // [r1, r2, r3, p1, p2, p3], that is the translation-related part should be the first in the front
    // angle-axis rotation vector comes in the back
    // The output matrix will be stored in a 1D array with 16 elements in ROW MAJOR order.
    // Again, the output is in ROW MAJOR order.
    memset(zeta_hat, T(0.0), 16*sizeof(*zeta_hat)); // restet the 16D array with zeros

    T tmp_vec[3] = {zeta_vec[3], zeta_vec[4], zeta_vec[5]};
    T tmp_hat[9] = {T(0.0)};
    vecHat(tmp_vec, tmp_hat);
    
    // allocate the return matrix, row major
    zeta_hat[0] = tmp_hat[0];  zeta_hat[1] = tmp_hat[1];  zeta_hat[2] = tmp_hat[2];   zeta_hat[3] = zeta_vec[0];
    zeta_hat[4] = tmp_hat[3];  zeta_hat[5] = tmp_hat[4];  zeta_hat[6] = tmp_hat[5];   zeta_hat[7] = zeta_vec[1];
    zeta_hat[8] = tmp_hat[6];  zeta_hat[9] = tmp_hat[7];  zeta_hat[10] = tmp_hat[8];  zeta_hat[11] = zeta_vec[2];
} // end of zetaHat


// the expZetaHat() funciton performs the exponential map over the zeta_hat matrix
// to obtain the corresponding transformation matrix for zeta
void expZetaHat(const Matrix4d &zeta_hat_mat, Matrix4d &T_mat) {
    // The input to the function should be 4x4 Eigen matrix with double precision values.
    // The output is the corresponding 4x4 Eigen matrix to input zeta hat matrix.
    Eigen::Vector3d rot_vec(zeta_hat_mat.coeff(2,1), zeta_hat_mat.coeff(0,2), zeta_hat_mat.coeff(1,0)); // extract the rotation vector
    double rot_angle_squared = rot_vec.transpose() * rot_vec;
    double rot_angle = sqrt(rot_angle_squared);
    if (rot_angle * rot_angle_squared > std::numeric_limits<double>::epsilon()) { // smallest denominator should be greater than numeric limit
        T_mat = Matrix4d::Identity() + zeta_hat_mat 
                + (1 - cos(rot_angle))/rot_angle_squared * zeta_hat_mat * zeta_hat_mat
                + (rot_angle - sin(rot_angle))/(rot_angle_squared * rot_angle) * zeta_hat_mat*zeta_hat_mat*zeta_hat_mat;
    } else {
        T_mat = Matrix4d::Identity() + zeta_hat_mat; 
    }
} // end of expZetaHat


// the zetaToSE3() function first converts the lie algebra tangent vector, zeta to the zeta hat matrix,
// then perform the exponential map over the zeta_hat matrix to obtain an transformation matrix
void zetaToSE3(const Vector6d &zeta_vec, Matrix4d &T_mat) {
    // Note the input zeta_vec vector should have a format of:
    // [tx, ty, tz, r1, r2, r3], that is the translation part should be the first in the front
    // angle-axis rotation vector comes in the back
    // The output is the corresponding 4x4 Eigen matrix to input zeta hat matrix.
    double zeta_hat[16]; // 1D array with 16 elements in ROW MAJOR order of zeta_hat matrix
    zetaHat(zeta_vec.data(), zeta_hat); // obtain the zeta_hat matrix stored in row-major order
    expZetaHat(Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::RowMajor>>(zeta_hat), T_mat); // get the transformation mat
} // zetaToSE3


// the lnVeeToZeta() function performs the logarithm map for a transformation matrix on SE(3),
// computing a 6 x 1 tangent vector, zeta, from a transformation
void lnVeeToZeta(const Matrix4d &T_mat, Vector6d &zeta) {
    // The input to the function should be a 4x4 Eigen matrix with double precision values.
    // The output is the corresponding 6x1 lie algebra vector, zeta.

    // 1. extract the rotation matrix and rho vector from the input transformation matrix
    Eigen::Matrix3d R_mat = T_mat.block<3, 3>(0, 0);
    Eigen::Vector3d Jac_rho_vec = T_mat.block<3, 1>(0, 3);

    // 2. recover 3x1 axis-angle vector from the extracted rotation matrix
    // the rotation matrix will be first raveled into an 1-D array, 
    // then converted to unit quaternion, finally from quaternion to rotation vector
    // (similar to the approach used by Google's Ceres-solver)
    double R_vec[9], unit_qua[4], phi_vec[3]; // define the necessary
    // to store data row by row, by default Eigen stores things in col major order, so needs to transpose the Rotation mat
    Eigen::Map<Eigen::MatrixXd>(R_vec, R_mat.cols(), R_mat.rows()) = R_mat.transpose();
    rotationMatrixArrayToUnitQuaternion(R_vec, unit_qua); // convert the rotation matrix vector to corresponding unit quaternion
    unitQuaternionToAngleAxis(unit_qua, phi_vec); // convert the unit quaternion into a 3x1 rotation vector

    // 3. compute the inversion of left jacobian of phi_vec on SO(3)
    Eigen::Matrix3d inv_left_Jac_mat; // define the inversion of left Jacobian w.r.t. the rotation vector
    invLeftJacobianSO3(Eigen::Vector3d(phi_vec[0], phi_vec[1], phi_vec[2]), inv_left_Jac_mat); // get inverted left Jac mat

    // 4. left multiply Jac_rho_vec by the inverse of left Jacobian to obtain the rho vector
    Eigen::Vector3d rho_vec = inv_left_Jac_mat * Jac_rho_vec; // J^{-1} * (J*rho) = rho

    // 5. construct the zeta from phi_vec and rho_vec vectors
    zeta.block<3, 1>(0, 0) = rho_vec; // translation vector at the top
    zeta.block<3, 1>(3, 0) = Eigen::Vector3d(phi_vec[0], phi_vec[1], phi_vec[2]);
} // end of lnVeeToZeta


// Not using a MatrixAdapter to accept matrix input
// Here uses a 1D array with 36 elements to store the output curly hat matrix of zeta
template<typename T>
inline void zetaCurlyHat(const T zeta_vec[6], T zeta_curly_hat[36]) {
    // Note the input zeta_vec vector should have a format of:
    // [tx, ty, tz, r1, r2, r3], that is the translation part should be the first in the front
    // angle-axis rotation vector comes in the back
    // The output matrix will be stored in a 1D array with 36 elements in ROW MAJOR order.
    // Again, the output is in ROW MAJOR order.
    memset(zeta_curly_hat, T(0.0), 36*sizeof(*zeta_curly_hat)); // restet the 16D array with zeros

    T phi_vec[3] = {zeta_vec[3], zeta_vec[4], zeta_vec[5]}; // extract the angle-axis rotation vector, phi_vec
    T phi_hat[9] = {T(0.0)}; // initialize the skew-symmetric matrix with zeros
    vecHat(phi_vec, phi_hat); // get the skew-symmetric matrix for phi_vec
    
    T rho_vec[3] = {zeta_vec[0], zeta_vec[1], zeta_vec[2]}; // extract the rho vector
    T rho_hat[9] = {T(0.0)}; // initialize the skew-symmetric matrix with zeros
    vecHat(rho_vec, rho_hat); // get the skew-symmetric matrix for rho_vec

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            zeta_curly_hat[i*6+j]   = phi_hat[i*3+j];
            zeta_curly_hat[i*6+3+j] = rho_hat[i*3+j];
            zeta_curly_hat[21+i*6+j]= phi_hat[i*3+j];
        }
    }
} // end of zetaCurlyHat


// the AdjointSE3() function convert a 4x4 SE(3) transformation matrix into a 6x6 adjointed transformation matrix
void AdjointSE3(const Matrix4d &T_mat, Matrix6d &Adj_T_mat) {
    // The input to the function should be a 4x4 Eigen transformation matrix with double precision values
    // The output is the corresponding 6x6 adjointed form of the transformation matrix
    Eigen::Matrix3d R_mat = T_mat.block<3, 3>(0, 0); // extract the rotation matrix, [3 x 3]
    // extract the translation vector, which is equivalent to the left Jacobian of rotation times the rho vector
    Eigen::Vector3d t_vec = T_mat.block<3, 1>(0, 3); // [3 x 1]
    double t_vec_hat[9]; // an array to store the skew-symmetric matrix for t_vec
    vecHat(t_vec.data(), t_vec_hat); // get the skew-symmetric form of t_vec, data stored in a ROW-MAJOR 1D array with 9 elements

    // Start allocating the matrix blocks
    Adj_T_mat.setZero(); // reset the Adjoint matrix to zeros before value allocation

    Adj_T_mat.block<3,3>(0,0) = R_mat; // top-left block of the Adjoint matrix is the rotation matrix
    Adj_T_mat.block<3,3>(3,3) = R_mat; // bottom-right block of the Adjoint matrix is the rotation matrix
    // top-right 3x3 block of the Adjoint is (J*rho)^{hat} * R_mat - Barfoot's SE(3) Table
    Adj_T_mat.block<3,3>(0,3) = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(t_vec_hat) * R_mat;
} // end of AdjointSE3


// the invAdjointSE3() function convert a 4x4 SE(3) transformation matrix 
// into a 6x6 inverted adjointed transformation matrix associated with the initial transformation matrix
void invAdjointSE3(const Matrix4d &T_mat, Matrix6d &inv_Adj_T_mat) {
    // The input to the function should be a 4x4 Eigen transformation matrix with double precision values
    // The output is the corresponding inversion of the 6x6 adjointed form of the transformation matrix
    Eigen::Matrix3d R_mat = T_mat.block<3, 3>(0, 0); // extract the rotation matrix, [3 x 3]
    Eigen::Matrix3d R_mat_transpose = R_mat.transpose(); // the transpose of the rotation matrix
    // extract the translation vector, which is equivalent to the left Jacobian of rotation times the rho vector
    Eigen::Vector3d t_vec = T_mat.block<3, 1>(0, 3); // [3 x 1]
    double t_vec_hat[9]; // an array to store the skew-symmetric matrix for t_vec
    vecHat(t_vec.data(), t_vec_hat); // get the skew-symmetric form of t_vec, data stored in a ROW-MAJOR 1D array with 9 elements

    // Start allocating the matrix blocks
    inv_Adj_T_mat.setZero(); // reset the inverted Adjoint matrix to zeros before value allocation

    inv_Adj_T_mat.block<3,3>(0,0) = R_mat_transpose; // top-left block is the tranposed rotation matrix
    inv_Adj_T_mat.block<3,3>(3,3) = R_mat_transpose; // bottom-right block is the tranposed rotation matrix
    // top-right 3x3 block of the inverted Adjoint is - R_mat^T * (J*rho)^{hat} - Barfoot's SE(3) Table
    inv_Adj_T_mat.block<3,3>(0,3) = - R_mat_transpose * Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(t_vec_hat);
} // end of invAdjointSE3


// the leftJacobianSE3() function computes the left Jacobian w.r.t. a SE(3) lie algebra vector, zeta.
void leftJacobianSE3(const Vector6d &zeta_vec, Matrix6d &left_Jac_mat) {
    // The input to the function should be a 6D vector of the lie algebra tangent vector with double precision values
    // The output matrix will be the corresponding Jacobian of the tangent vector in a 6x6 matrix with double precisions
    // Extract the rotation vector, phi_vec
    Eigen::Vector3d phi_vec(zeta_vec[3], zeta_vec[4], zeta_vec[5]); // rotation vector is behind the rho vector, [3 x 1]
    double rot_angle_squared = phi_vec.transpose() * phi_vec;
    double rot_angle = sqrt(rot_angle_squared); // rotaion angle, phi - scalar
    double rot_angle_pow5 = rot_angle_squared * rot_angle_squared * rot_angle; // rot_angle^5

    // the zeta_curly_hat matrix will be used in both cases
    double zeta_curly_hat[36]; // the curly hat form of the zeta vector - a raveled 6x6 matrix
    zetaCurlyHat(zeta_vec.data(), zeta_curly_hat); // the zeta_curly_hat variable stores the ROW-MAJOR matrix elements
    Matrix6d zeta_curly_hat_mat = Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(zeta_curly_hat);

    if (rot_angle_pow5 > std::numeric_limits<double>::epsilon()) {
        // This part is run only if (rot_angle)^5 is greater than numeric limits to avoid dividing by a very small value
        
        // define some intermediate terms
        double sin_phi = sin(rot_angle); // sin(phi)
        double cos_phi = cos(rot_angle); // cos(phi)
        double coeff_1 = (4.0 - rot_angle*sin_phi - 4.0*cos_phi) / (2.0 * rot_angle_squared);
        double coeff_2 = (4.0*rot_angle - 5.0*sin_phi + rot_angle*cos_phi) / (2.0 * rot_angle_squared * rot_angle);
        double coeff_3 = (2.0 - rot_angle*sin_phi - 2.0*cos_phi) / (2.0 * rot_angle_squared * rot_angle_squared);
        double coeff_4 = (2.0*rot_angle - 3.0*sin_phi + rot_angle*cos_phi) / (2.0 * rot_angle_pow5);
        Matrix6d hat_mat_pow2 = zeta_curly_hat_mat * zeta_curly_hat_mat; // power 2 of zeta_curly_hat_mat

        // calculate the output left jacobian matrix
        left_Jac_mat =   Matrix6d::Identity()
                       + coeff_1 * zeta_curly_hat_mat
                       + coeff_2 * hat_mat_pow2
                       + coeff_3 * hat_mat_pow2 * zeta_curly_hat_mat
                       + coeff_4 * hat_mat_pow2 * hat_mat_pow2;

    } else { // if the rot_angle is small
        left_Jac_mat =   Matrix6d::Identity() 
                       + 0.5 * zeta_curly_hat_mat;
    }
} // end of leftJacobianSE3()


// the invLeftJacobianSE3() function computes the inversion of left Jacobian w.r.t. a SE(3) lie algebra vector, zeta.
void invLeftJacobianSE3(const Vector6d &zeta_vec, Matrix6d &inv_Jac_mat) {
    // The input to the function should be a 6D vector of the lie algebra tangent vector with double precision values
    // The output matrix will be the corresponding inversion of Jacobian of the tangent vector in a 6x6 matrix
    // Extract the rotation vector, phi_vec
    Eigen::Vector3d phi_vec(zeta_vec[3], zeta_vec[4], zeta_vec[5]); // rotation vector is behind the rho vector, [3 x 1]

    // 1. Call the invLeftJacobianSO3() function to obtain the inverted left Jacobian matrix in SO(3) based on the phi_vec
    Eigen::Matrix3d inv_left_Jac_so3; // [3 x 3]
    invLeftJacobianSO3(phi_vec, inv_left_Jac_so3); // J_so3^(-1)

    // 2. Call the leftJacobianSE3() function to obtain the left Jacobian matrix in SE(3) based on the 6D zeta
    Matrix6d left_Jac_se3;
    leftJacobianSE3(zeta_vec, left_Jac_se3); // this make easy to get the Q matrix
    Eigen::Matrix3d Q_mat = left_Jac_se3.block<3,3>(0,3); // the upper-right block is the Q matrix

    // 3. Construct the inv_left_jacobian matirx
    inv_Jac_mat.setZero(); // reset the Jacobian matrix
    inv_Jac_mat.block<3, 3>(0, 0) = inv_left_Jac_so3; // top-left block as J_so3^(-1)
    inv_Jac_mat.block<3, 3>(3, 3) = inv_left_Jac_so3; // bottom-right block as J_so3^(-1)
    inv_Jac_mat.block<3, 3>(0, 3) = - inv_left_Jac_so3 * Q_mat * inv_left_Jac_so3; // top-right block -J^(-1) * Q * J^(-1)
} // end of invLeftJacobianSE3



#endif // UTILS_INCLUDE_H
