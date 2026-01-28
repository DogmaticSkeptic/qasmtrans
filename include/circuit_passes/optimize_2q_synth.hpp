#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "../IR/circuit.hpp"
#include "../IR/gate.hpp"
#include "../QASMTransPrimitives.hpp"
#include "optimize_1q.hpp"

#ifdef QASMTRANS_USE_EIGEN
#include <Eigen/Dense>
#endif

namespace QASMTrans
{
    namespace optimize
    {
#ifdef QASMTRANS_USE_EIGEN
        namespace detail_2q_synth
        {
            using Complex = std::complex<double>;
            using Mat2 = Eigen::Matrix2cd;
            using Mat4 = Eigen::Matrix4cd;

            constexpr double kTol = 1e-12;

            inline Mat2 to_eigen(const optimize::detail::Mat2 &in)
            {
                Mat2 out;
                out(0, 0) = in.m[0][0];
                out(0, 1) = in.m[0][1];
                out(1, 0) = in.m[1][0];
                out(1, 1) = in.m[1][1];
                return out;
            }

            inline optimize::detail::Mat2 to_mat2(const Mat2 &in)
            {
                optimize::detail::Mat2 out{};
                out.m[0][0] = in(0, 0);
                out.m[0][1] = in(0, 1);
                out.m[1][0] = in(1, 0);
                out.m[1][1] = in(1, 1);
                return out;
            }

            inline Mat4 kron(const Mat2 &a, const Mat2 &b)
            {
                Mat4 out;
                out.block<2, 2>(0, 0) = a(0, 0) * b;
                out.block<2, 2>(0, 2) = a(0, 1) * b;
                out.block<2, 2>(2, 0) = a(1, 0) * b;
                out.block<2, 2>(2, 2) = a(1, 1) * b;
                return out;
            }

            inline Mat4 cx_matrix(bool control_first)
            {
                Mat4 out = Mat4::Zero();
                if (control_first)
                {
                    out(0, 0) = 1.0;
                    out(1, 1) = 1.0;
                    out(2, 3) = 1.0;
                    out(3, 2) = 1.0;
                }
                else
                {
                    out(0, 0) = 1.0;
                    out(2, 2) = 1.0;
                    out(1, 3) = 1.0;
                    out(3, 1) = 1.0;
                }
                return out;
            }

            inline Mat4 cz_matrix()
            {
                Mat4 out = Mat4::Identity();
                out(3, 3) = -1.0;
                return out;
            }

            inline Mat4 swap_matrix()
            {
                Mat4 out = Mat4::Zero();
                out(0, 0) = 1.0;
                out(1, 2) = 1.0;
                out(2, 1) = 1.0;
                out(3, 3) = 1.0;
                return out;
            }

            inline Mat4 iswap_matrix()
            {
                Mat4 out = Mat4::Zero();
                const Complex I(0.0, 1.0);
                out(0, 0) = 1.0;
                out(1, 2) = I;
                out(2, 1) = I;
                out(3, 3) = 1.0;
                return out;
            }

            inline Mat4 ecr_matrix()
            {
                const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
                const Complex I(0.0, 1.0);
                Mat4 out;
                out << Complex(0.0, 0.0), inv_sqrt2, Complex(0.0, 0.0), I * inv_sqrt2,
                    inv_sqrt2, Complex(0.0, 0.0), -I * inv_sqrt2, Complex(0.0, 0.0),
                    Complex(0.0, 0.0), I * inv_sqrt2, Complex(0.0, 0.0), inv_sqrt2,
                    -I * inv_sqrt2, Complex(0.0, 0.0), inv_sqrt2, Complex(0.0, 0.0);
                return out;
            }

            inline Mat2 x_matrix()
            {
                Mat2 out;
                out << 0.0, 1.0,
                    1.0, 0.0;
                return out;
            }

            inline Mat2 y_matrix()
            {
                const Complex I(0.0, 1.0);
                Mat2 out;
                out << 0.0, -I,
                    I, 0.0;
                return out;
            }

            inline Mat2 z_matrix()
            {
                Mat2 out;
                out << 1.0, 0.0,
                    0.0, -1.0;
                return out;
            }

            inline Mat4 rxx_matrix(double theta)
            {
                const Complex I(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat4 xx = kron(x_matrix(), x_matrix());
                return c * Mat4::Identity() - I * s * xx;
            }

            inline Mat4 ryy_matrix(double theta)
            {
                const Complex I(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat4 yy = kron(y_matrix(), y_matrix());
                return c * Mat4::Identity() - I * s * yy;
            }

            inline Mat4 rzz_matrix(double theta)
            {
                const Complex I(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat4 zz = kron(z_matrix(), z_matrix());
                return c * Mat4::Identity() - I * s * zz;
            }

            inline Mat4 rzx_matrix(double theta, bool control_first)
            {
                const Complex I(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat4 zx = control_first ? kron(z_matrix(), x_matrix()) : kron(x_matrix(), z_matrix());
                return c * Mat4::Identity() - I * s * zx;
            }

            inline Mat4 controlled_from_u(const Mat2 &u, bool control_first)
            {
                Mat4 out = Mat4::Identity();
                if (control_first)
                {
                    out.block<2, 2>(2, 2) = u;
                    return out;
                }
                Mat4 swap = swap_matrix();
                Mat4 base = Mat4::Identity();
                base.block<2, 2>(2, 2) = u;
                return swap * base * swap;
            }

            inline bool two_qubit_gate_matrix(const Gate &gate, IdxType q0, IdxType q1, Mat4 &out)
            {
                const bool control_first = gate.ctrl == q0 && gate.qubit == q1;
                const bool control_second = gate.ctrl == q1 && gate.qubit == q0;
                switch (gate.op_name)
                {
                case OP::CX:
                    if (!control_first && !control_second)
                        return false;
                    out = cx_matrix(control_first);
                    return true;
                case OP::CZ:
                    out = cz_matrix();
                    return true;
                case OP::SWAP:
                    out = swap_matrix();
                    return true;
                case OP::ISWAP:
                    out = iswap_matrix();
                    return true;
                case OP::ECR:
                    out = control_first ? ecr_matrix() : (swap_matrix() * ecr_matrix() * swap_matrix());
                    return true;
                case OP::RXX:
                    out = rxx_matrix(gate.theta);
                    return true;
                case OP::RYY:
                    out = ryy_matrix(gate.theta);
                    return true;
                case OP::RZZ:
                    out = rzz_matrix(gate.theta);
                    return true;
                case OP::RZX:
                    out = rzx_matrix(gate.theta, control_first);
                    return true;
                case OP::ZZ:
                    out = rzz_matrix(PI / 2.0);
                    return true;
                case OP::CRX:
                {
                    if (!control_first && !control_second)
                        return false;
                    optimize::detail::Mat2 rx = optimize::detail::rx_matrix(gate.theta);
                    out = controlled_from_u(to_eigen(rx), control_first);
                    return true;
                }
                case OP::CRY:
                {
                    if (!control_first && !control_second)
                        return false;
                    optimize::detail::Mat2 ry = optimize::detail::ry_matrix(gate.theta);
                    out = controlled_from_u(to_eigen(ry), control_first);
                    return true;
                }
                case OP::CRZ:
                {
                    if (!control_first && !control_second)
                        return false;
                    optimize::detail::Mat2 rz = optimize::detail::rz_matrix(gate.theta);
                    out = controlled_from_u(to_eigen(rz), control_first);
                    return true;
                }
                case OP::CP:
                {
                    if (!control_first && !control_second)
                        return false;
                    Mat2 phase;
                    phase << 1.0, 0.0, 0.0, std::exp(Complex(0.0, gate.theta));
                    out = controlled_from_u(phase, control_first);
                    return true;
                }
                default:
                    break;
                }
                return false;
            }

            inline bool gate_matrix_4x4(const Gate &gate, IdxType q0, IdxType q1, Mat4 &out)
            {
                if (gate.has_custom_name())
                {
                    return false;
                }
                if (gate.ctrl >= 0 && gate.extra < 0)
                {
                    return two_qubit_gate_matrix(gate, q0, q1, out);
                }
                if (gate.ctrl < 0 && gate.extra < 0)
                {
                    optimize::detail::Mat2 mat{};
                    if (!optimize::detail::gate_matrix(gate, mat))
                    {
                        return false;
                    }
                    Mat2 oneq = to_eigen(mat);
                    if (gate.qubit == q0)
                    {
                        out = kron(oneq, Mat2::Identity());
                        return true;
                    }
                    if (gate.qubit == q1)
                    {
                        out = kron(Mat2::Identity(), oneq);
                        return true;
                    }
                }
                return false;
            }

            inline bool unitary_from_gates(const std::vector<Gate> &gates, IdxType q0, IdxType q1, Mat4 &out)
            {
                out = Mat4::Identity();
                for (const auto &gate : gates)
                {
                    Mat4 mat;
                    if (!gate_matrix_4x4(gate, q0, q1, mat))
                    {
                        return false;
                    }
                    out = mat * out;
                }
                return true;
            }

            inline double unitary_overlap(const Mat4 &a, const Mat4 &b)
            {
                Complex trace = (a.adjoint() * b).trace();
                return std::abs(trace) / 4.0;
            }

            inline bool is_two_qubit_gate(const Gate &gate)
            {
                return gate.ctrl >= 0 && gate.qubit >= 0 && gate.extra < 0;
            }

            inline Mat4 magic_basis()
            {
                Mat4 out;
                const Complex I(0.0, 1.0);
                out << 1.0, I, 0.0, 0.0,
                    0.0, 0.0, I, 1.0,
                    0.0, 0.0, I, -1.0,
                    1.0, -I, 0.0, 0.0;
                return out;
            }

            inline Mat4 magic_basis_dagger()
            {
                Mat4 out;
                out << 0.5, 0.0, 0.0, 0.5,
                    Complex(0.0, -0.5), 0.0, 0.0, Complex(0.0, 0.5),
                    0.0, Complex(0.0, -0.5), Complex(0.0, -0.5), 0.0,
                    0.0, 0.5, -0.5, 0.0;
                return out;
            }

            enum class MagicBasisTransform
            {
                Into,
                OutOf
            };

            inline void dump_unitary(std::ostream &os, const Mat4 &mat, const std::string &label)
            {
                os << label << ":\n";
                for (int r = 0; r < 4; ++r)
                {
                    for (int c = 0; c < 4; ++c)
                    {
                        const Complex &val = mat(r, c);
                        os << "(" << val.real() << "," << val.imag() << ")";
                        if (c < 3)
                        {
                            os << " ";
                        }
                    }
                    os << "\n";
                }
            }

            inline void dump_block_gates(std::ostream &os,
                                         const std::vector<Gate> &gates,
                                         const std::vector<std::size_t> &positions)
            {
                os << "block_gates:\n";
                for (std::size_t pos : positions)
                {
                    if (pos >= gates.size())
                    {
                        continue;
                    }
                    const Gate &g = gates[pos];
                    os << "  idx=" << pos
                       << " op=" << g.lower_name()
                       << " ctrl=" << g.ctrl
                       << " qubit=" << g.qubit
                       << " theta=" << g.theta
                       << " phi=" << g.phi
                       << " lam=" << g.lam
                       << "\n";
                }
            }

            inline void dump_synth_failure(const Mat4 &target,
                                           bool synth_ok,
                                           const Mat4 *synth,
                                           double overlap,
                                           bool reverse_attempted,
                                           bool reverse_ok,
                                           double overlap_rev,
                                           IdxType q0,
                                           IdxType q1,
                                           const std::vector<Gate> &gates,
                                           const std::vector<std::size_t> &positions)
            {
                const char *dump_dir = std::getenv("QASMTRANS_2Q_SYNTH_DUMP");
                if (!dump_dir || dump_dir[0] == '\0')
                {
                    return;
                }
                static int dump_count = 0;
                static constexpr int kMaxDumps = 3;
                if (dump_count >= kMaxDumps)
                {
                    return;
                }
                std::filesystem::path dir(dump_dir);
                std::error_code ec;
                std::filesystem::create_directories(dir, ec);
                std::filesystem::path out_path = dir / ("synth_fail_" + std::to_string(dump_count) + ".txt");
                std::ofstream out(out_path);
                if (!out)
                {
                    return;
                }
                out << "q0=" << q0 << " q1=" << q1 << "\n";
                out << "synth_ok=" << (synth_ok ? "true" : "false") << " overlap=" << overlap << "\n";
                out << "reverse_attempted=" << (reverse_attempted ? "true" : "false")
                    << " reverse_ok=" << (reverse_ok ? "true" : "false")
                    << " overlap_rev=" << overlap_rev << "\n";
                dump_block_gates(out, gates, positions);
                dump_unitary(out, target, "target_unitary");
                if (synth_ok && synth)
                {
                    dump_unitary(out, *synth, "synth_unitary");
                }
                dump_count += 1;
            }

            inline Mat4 magic_basis_transform(const Mat4 &unitary, MagicBasisTransform direction)
            {
                static const Mat4 b = magic_basis();
                static const Mat4 b_dag = magic_basis_dagger();
                if (direction == MagicBasisTransform::OutOf)
                {
                    return b_dag * unitary * b;
                }
                return b * unitary * b_dag;
            }

            inline Mat2 ipx_matrix()
            {
                const Complex I(0.0, 1.0);
                Mat2 out;
                out << 0.0, I,
                    I, 0.0;
                return out;
            }

            inline Mat2 ipy_matrix()
            {
                Mat2 out;
                out << 0.0, 1.0,
                    -1.0, 0.0;
                return out;
            }

            inline Mat2 ipz_matrix()
            {
                const Complex I(0.0, 1.0);
                Mat2 out;
                out << I, 0.0,
                    0.0, -I;
                return out;
            }

            inline double mod_positive(double value, double period)
            {
                double out = std::fmod(value, period);
                if (out < 0.0)
                {
                    out += period;
                }
                return out;
            }

            inline bool approx_equal(const Mat4 &a, const Mat4 &b, double tol)
            {
                return (a - b).cwiseAbs().maxCoeff() < tol;
            }

            inline bool decompose_two_qubit_product_gate(const Mat4 &su4,
                                                         Mat2 &left,
                                                         Mat2 &right,
                                                         double &phase)
            {
                Mat2 r = su4.block<2, 2>(0, 0);
                Complex det_r = r.determinant();
                if (std::abs(det_r) < 0.1)
                {
                    r = su4.block<2, 2>(2, 0);
                    det_r = r.determinant();
                }
                if (std::abs(det_r) < 0.1)
                {
                    return false;
                }
                r /= std::sqrt(det_r);
                Mat2 r_t_conj = r.conjugate().transpose();
                Mat4 temp = kron(Mat2::Identity(), r_t_conj);
                Mat4 tmp = su4 * temp;
                Mat2 l;
                l(0, 0) = tmp(0, 0);
                l(0, 1) = tmp(0, 2);
                l(1, 0) = tmp(2, 0);
                l(1, 1) = tmp(2, 2);
                Complex det_l = l.determinant();
                if (std::abs(det_l) < 0.9)
                {
                    return false;
                }
                l /= std::sqrt(det_l);
                left = l;
                right = r;
                phase = std::arg(det_l) / 2.0;
                return true;
            }

            struct WeylDecomposition
            {
                double a = 0.0;
                double b = 0.0;
                double c = 0.0;
                double global_phase = 0.0;
                Mat2 K1l;
                Mat2 K1r;
                Mat2 K2l;
                Mat2 K2r;
                Mat4 unitary;
            };

            inline WeylDecomposition weyl_decompose(const Mat4 &unitary)
            {
                WeylDecomposition out;
                out.unitary = unitary;
                Complex det_u = unitary.determinant();
                Complex det_pow = std::pow(det_u, 0.25);
                Mat4 u = unitary / det_pow;
                double global_phase = std::arg(det_u) / 4.0;

                Mat4 u_p = magic_basis_transform(u, MagicBasisTransform::OutOf);
                Mat4 m2 = u_p.transpose() * u_p;

                std::mt19937_64 rng(2023);
                std::normal_distribution<double> dist(0.0, 1.0);
                bool diag_found = false;
                bool decomp_found = false;
                std::array<double, 3> cs_final{};
                Mat2 K1l, K1r, K2l, K2r;
                double global_phase_final = 0.0;
                for (int i = 0; i < 100; ++i)
                {
                    double rand_a = (i == 0) ? 1.2602066112249388 : dist(rng);
                    double rand_b = (i == 0) ? 0.22317849046722027 : dist(rng);
                    Eigen::Matrix4d m2_real = rand_a * m2.real() + rand_b * m2.imag();
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> es(m2_real);
                    if (es.info() != Eigen::Success)
                    {
                        continue;
                    }
                    Eigen::Matrix4d p_real = es.eigenvectors();
                    Mat4 p_inner = p_real.cast<Complex>();
                    Eigen::Matrix<Complex, 4, 1> d_inner = (p_inner.transpose() * m2 * p_inner).diagonal();
                    Mat4 diag = Mat4::Zero();
                    diag.diagonal() = d_inner;
                    Mat4 compare = p_inner * diag * p_inner.transpose();
                    if (!approx_equal(compare, m2, 1.0e-12))
                    {
                        continue;
                    }
                    diag_found = true;

                    std::array<double, 4> dargs{};
                    for (int idx = 0; idx < 4; ++idx)
                    {
                        dargs[idx] = -std::arg(d_inner(idx)) / 2.0;
                    }
                    dargs[3] = -dargs[0] - dargs[1] - dargs[2];
                    std::array<double, 3> cs{};
                    for (int idx = 0; idx < 3; ++idx)
                    {
                        cs[idx] = mod_positive((dargs[idx] + dargs[3]) / 2.0, 2.0 * PI);
                    }
                    std::array<double, 3> cstemp{};
                    for (int idx = 0; idx < 3; ++idx)
                    {
                        double v = mod_positive(cs[idx], PI / 2.0);
                        cstemp[idx] = std::min(v, PI / 2.0 - v);
                    }
                    std::array<int, 3> order{0, 1, 2};
                    std::sort(order.begin(), order.end(), [&](int lhs, int rhs)
                              { return cstemp[lhs] < cstemp[rhs]; });
                    std::array<int, 3> perm{order[1], order[2], order[0]};
                    std::array<double, 3> cs_new{};
                    std::array<double, 4> d_new{};
                    for (int idx = 0; idx < 3; ++idx)
                    {
                        cs_new[idx] = cs[perm[idx]];
                        d_new[idx] = dargs[perm[idx]];
                    }
                    d_new[3] = dargs[3];
                    cs = cs_new;
                    dargs = d_new;

                    Mat4 p_try = p_inner;
                    Mat4 p_orig = p_try;
                    for (int idx = 0; idx < 3; ++idx)
                    {
                        p_try.col(idx) = p_orig.col(perm[idx]);
                    }
                    if (p_try.determinant().real() < 0.0)
                    {
                        p_try.col(3) = -p_try.col(3);
                    }

                    Mat4 temp = Mat4::Zero();
                    for (int idx = 0; idx < 4; ++idx)
                    {
                        temp(idx, idx) = std::exp(Complex(0.0, dargs[idx]));
                    }
                    Mat4 k1 = magic_basis_transform(u_p * p_try * temp, MagicBasisTransform::Into);
                    Mat4 k2 = magic_basis_transform(p_try.transpose(), MagicBasisTransform::Into);

                    Mat2 cand_K1l, cand_K1r, cand_K2l, cand_K2r;
                    double phase_l = 0.0;
                    double phase_r = 0.0;
                    if (!decompose_two_qubit_product_gate(k1, cand_K1l, cand_K1r, phase_l) ||
                        !decompose_two_qubit_product_gate(k2, cand_K2l, cand_K2r, phase_r))
                    {
                        continue;
                    }

                    double cand_global_phase = global_phase + phase_l + phase_r;

                    Mat2 ipx = ipx_matrix();
                    Mat2 ipy = ipy_matrix();
                    Mat2 ipz = ipz_matrix();

                    if (cs[0] > PI / 2.0)
                    {
                        cs[0] -= 3.0 * PI / 2.0;
                        cand_K1l = cand_K1l * ipy;
                        cand_K1r = cand_K1r * ipy;
                        cand_global_phase += PI / 2.0;
                    }
                    if (cs[1] > PI / 2.0)
                    {
                        cs[1] -= 3.0 * PI / 2.0;
                        cand_K1l = cand_K1l * ipx;
                        cand_K1r = cand_K1r * ipx;
                        cand_global_phase += PI / 2.0;
                    }
                    int conjs = 0;
                    if (cs[0] > PI / 4.0)
                    {
                        cs[0] = PI / 2.0 - cs[0];
                        cand_K1l = cand_K1l * ipy;
                        cand_K2r = ipy * cand_K2r;
                        conjs += 1;
                        cand_global_phase -= PI / 2.0;
                    }
                    if (cs[1] > PI / 4.0)
                    {
                        cs[1] = PI / 2.0 - cs[1];
                        cand_K1l = cand_K1l * ipx;
                        cand_K2r = ipx * cand_K2r;
                        conjs += 1;
                        cand_global_phase += PI / 2.0;
                        if (conjs == 1)
                        {
                            cand_global_phase -= PI;
                        }
                    }
                    if (cs[2] > PI / 2.0)
                    {
                        cs[2] -= 3.0 * PI / 2.0;
                        cand_K1l = cand_K1l * ipz;
                        cand_K1r = cand_K1r * ipz;
                        cand_global_phase += PI / 2.0;
                        if (conjs == 1)
                        {
                            cand_global_phase -= PI;
                        }
                    }
                    if (conjs == 1)
                    {
                        cs[2] = PI / 2.0 - cs[2];
                        cand_K1l = cand_K1l * ipz;
                        cand_K2r = ipz * cand_K2r;
                        cand_global_phase += PI / 2.0;
                    }
                    if (cs[2] > PI / 4.0)
                    {
                        cs[2] -= PI / 2.0;
                        cand_K1l = cand_K1l * ipz;
                        cand_K1r = cand_K1r * ipz;
                        cand_global_phase -= PI / 2.0;
                    }

                    cs_final = cs;
                    K1l = cand_K1l;
                    K1r = cand_K1r;
                    K2l = cand_K2l;
                    K2r = cand_K2r;
                    global_phase_final = cand_global_phase;
                    decomp_found = true;
                    break;
                }
                if (!diag_found)
                {
                    throw std::runtime_error("TwoQubitWeylDecomposition failed to diagonalize M2");
                }
                if (!decomp_found)
                {
                    throw std::runtime_error("TwoQubitWeylDecomposition failed to decompose k1/k2");
                }

                out.a = cs_final[1];
                out.b = cs_final[0];
                out.c = cs_final[2];
                out.global_phase = global_phase_final;
                out.K1l = K1l;
                out.K1r = K1r;
                out.K2l = K2l;
                out.K2r = K2r;
                return out;
            }

            inline double trace_to_fidelity(const Complex &trace)
            {
                return (4.0 + std::norm(trace)) / 20.0;
            }

            inline int best_basis_count(double a, double b, double c, double basis_b, double basis_fidelity)
            {
                Complex trace0(4.0 * (std::cos(a) * std::cos(b) * std::cos(c)),
                               4.0 * (std::sin(a) * std::sin(b) * std::sin(c)));
                Complex trace1(4.0 * (std::cos(PI / 4.0 - a) * std::cos(basis_b - b) * std::cos(c)),
                               4.0 * (std::sin(PI / 4.0 - a) * std::sin(basis_b - b) * std::sin(c)));
                Complex trace2(4.0 * std::cos(c), 0.0);
                Complex trace3(4.0, 0.0);
                std::array<double, 4> fid{
                    trace_to_fidelity(trace0),
                    trace_to_fidelity(trace1) * std::pow(basis_fidelity, 1),
                    trace_to_fidelity(trace2) * std::pow(basis_fidelity, 2),
                    trace_to_fidelity(trace3) * std::pow(basis_fidelity, 3)};

                int best = 0;
                double best_val = fid[0];
                for (int i = 1; i < 4; ++i)
                {
                    if (fid[i] > best_val + 1e-15)
                    {
                        best_val = fid[i];
                        best = i;
                    }
                }
                return best;
            }

            inline std::vector<Gate> synthesize_one_qubit(const Mat2 &unitary,
                                                          IdxType qubit,
                                                          const std::shared_ptr<Chip> &chip,
                                                          const std::unordered_set<std::string> *basis_gates,
                                                          bool pulse_only,
                                                          optimize::detail::EulerBasis pulse_basis)
            {
                optimize::detail::Mat2 mat = to_mat2(unitary);
                bool use_errors = optimize::detail::has_error_data(chip, qubit);
                double best_error = std::numeric_limits<double>::infinity();
                std::vector<Gate> best_seq;

                if (pulse_only)
                {
                    double theta = 0.0, phi = 0.0, lam = 0.0, phase = 0.0;
                    optimize::detail::angles_from_unitary(mat, pulse_basis, theta, phi, lam, phase);
                    best_seq = optimize::detail::generate_circuit(pulse_basis, theta, phi, lam, phase, qubit, basis_gates);
                    return best_seq;
                }

                const std::vector<optimize::detail::EulerBasis> bases = optimize::detail::possible_bases(basis_gates);
                for (auto basis : bases)
                {
                    double theta = 0.0, phi = 0.0, lam = 0.0, phase = 0.0;
                    optimize::detail::angles_from_unitary(mat, basis, theta, phi, lam, phase);
                    std::vector<Gate> candidate = optimize::detail::generate_circuit(basis, theta, phi, lam, phase, qubit, basis_gates);
                    double err = optimize::detail::compute_sequence_error(candidate, chip, qubit, use_errors);
                    if (err < best_error || (std::abs(err - best_error) < 1e-12 && candidate.size() < best_seq.size()))
                    {
                        best_error = err;
                        best_seq = std::move(candidate);
                    }
                }
                return best_seq;
            }

            struct BasisContext
            {
                WeylDecomposition basis;
                double basis_b = 0.0;
                Mat2 u0l, u0r, u1l, u1ra, u1rb, u2la, u2lb, u2ra, u2rb, u3l, u3r;
                Mat2 q0l, q0r, q1la, q1lb, q1ra, q1rb, q2l, q2r;
            };

            enum class SynthFailReason : int
            {
                None = 0,
                BasisContext = 1,
                WeylDiag = 2,
                WeylK = 3,
                WeylOther = 4,
            };

            inline void set_fail_reason(int *out, SynthFailReason reason)
            {
                if (out)
                {
                    *out = static_cast<int>(reason);
                }
            }

            inline BasisContext build_basis_context(const Mat4 &basis_matrix)
            {
                BasisContext ctx;
                ctx.basis = weyl_decompose(basis_matrix);
                ctx.basis_b = ctx.basis.b;

                const Complex I(0.0, 1.0);
                const Complex M_I(0.0, -1.0);
                const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
                const double b = ctx.basis.b;

                const Complex temp1(0.5, -0.5);
                Mat2 k11l;
                k11l << temp1 * (M_I * std::exp(Complex(0.0, -b))),
                    temp1 * std::exp(Complex(0.0, -b)),
                    temp1 * (M_I * std::exp(Complex(0.0, b))),
                    temp1 * (-std::exp(Complex(0.0, b)));
                Mat2 k11r;
                k11r << inv_sqrt2 * (I * std::exp(Complex(0.0, -b))),
                    inv_sqrt2 * (-std::exp(Complex(0.0, -b))),
                    inv_sqrt2 * std::exp(Complex(0.0, b)),
                    inv_sqrt2 * (M_I * std::exp(Complex(0.0, b)));

                Mat2 k12l;
                k12l << Complex(0.5, 0.5), Complex(0.5, 0.5),
                    Complex(-0.5, 0.5), Complex(0.5, -0.5);
                Mat2 k12r;
                k12r << Complex(0.0, inv_sqrt2), Complex(inv_sqrt2, 0.0),
                    Complex(-inv_sqrt2, 0.0), Complex(0.0, -inv_sqrt2);

                Mat2 k32l_k21l;
                k32l_k21l << inv_sqrt2 * Complex(1.0, std::cos(2.0 * b)),
                    inv_sqrt2 * (I * std::sin(2.0 * b)),
                    inv_sqrt2 * (I * std::sin(2.0 * b)),
                    inv_sqrt2 * Complex(1.0, -std::cos(2.0 * b));

                const Complex temp2(0.5, 0.5);
                Mat2 k21r;
                k21r << temp2 * (M_I * std::exp(Complex(0.0, -2.0 * b))),
                    temp2 * std::exp(Complex(0.0, -2.0 * b)),
                    temp2 * (I * std::exp(Complex(0.0, 2.0 * b))),
                    temp2 * std::exp(Complex(0.0, 2.0 * b));

                Mat2 k22l;
                k22l << inv_sqrt2, -inv_sqrt2,
                    inv_sqrt2, inv_sqrt2;

                Mat2 k22r;
                k22r << 0.0, 1.0,
                    -1.0, 0.0;

                Mat2 k31l;
                k31l << inv_sqrt2 * std::exp(Complex(0.0, -b)),
                    inv_sqrt2 * std::exp(Complex(0.0, -b)),
                    inv_sqrt2 * -std::exp(Complex(0.0, b)),
                    inv_sqrt2 * std::exp(Complex(0.0, b));

                Mat2 k31r;
                k31r << I * std::exp(Complex(0.0, b)), 0.0,
                    0.0, M_I * std::exp(Complex(0.0, -b));

                Mat2 k32r;
                k32r << temp2 * std::exp(Complex(0.0, b)),
                    temp2 * -std::exp(Complex(0.0, -b)),
                    temp2 * (M_I * std::exp(Complex(0.0, b))),
                    temp2 * (M_I * std::exp(Complex(0.0, -b)));

                Mat2 k1ld = ctx.basis.K1l.conjugate().transpose();
                Mat2 k1rd = ctx.basis.K1r.conjugate().transpose();
                Mat2 k2ld = ctx.basis.K2l.conjugate().transpose();
                Mat2 k2rd = ctx.basis.K2r.conjugate().transpose();

                ctx.u0l = k31l * k1ld;
                ctx.u0r = k31r * k1rd;
                ctx.u1l = k2ld * k32l_k21l * k1ld;
                ctx.u1ra = k2rd * k32r;
                ctx.u1rb = k21r * k1rd;
                ctx.u2la = k2ld * k22l;
                ctx.u2lb = k11l * k1ld;
                ctx.u2ra = k2rd * k22r;
                ctx.u2rb = k11r * k1rd;
                ctx.u3l = k2ld * k12l;
                ctx.u3r = k2rd * k12r;

                Mat2 ipz = ipz_matrix();
                ctx.q0l = k12l.conjugate().transpose() * k1ld;
                ctx.q0r = k12r.conjugate().transpose() * ipz * k1rd;
                ctx.q1la = k2ld * k11l.conjugate().transpose();
                ctx.q1lb = k11l * k1ld;
                ctx.q1ra = k2rd * ipz * k11r.conjugate().transpose();
                ctx.q1rb = k11r * k1rd;
                ctx.q2l = k2ld * k12l;
                ctx.q2r = k2rd * k12r;
                return ctx;
            }

            struct BasisSpec
            {
                std::string name;
                OP op;
                Mat4 matrix;
                bool directional = false;
                bool is_param = false;
                double fixed_theta = 0.0;
            };

            inline BasisSpec choose_basis_spec(const std::unordered_set<std::string> *basis_gates)
            {
                auto has_gate = [&](const char *name) {
                    return basis_gates && basis_gates->find(name) != basis_gates->end();
                };

                if (has_gate("ecr"))
                {
                    return {"ecr", OP::ECR, ecr_matrix(), true, false, 0.0};
                }
                if (has_gate("cx"))
                {
                    return {"cx", OP::CX, cx_matrix(true), true, false, 0.0};
                }
                if (has_gate("cz"))
                {
                    return {"cz", OP::CZ, cz_matrix(), false, false, 0.0};
                }
                if (has_gate("iswap"))
                {
                    return {"iswap", OP::ISWAP, iswap_matrix(), false, false, 0.0};
                }
                if (has_gate("rxx"))
                {
                    return {"rxx", OP::RXX, rxx_matrix(PI / 2.0), false, true, PI / 2.0};
                }
                if (has_gate("ryy"))
                {
                    return {"ryy", OP::RYY, ryy_matrix(PI / 2.0), false, true, PI / 2.0};
                }
                if (has_gate("rzz"))
                {
                    return {"rzz", OP::RZZ, rzz_matrix(PI / 2.0), false, true, PI / 2.0};
                }
                if (has_gate("rzx"))
                {
                    return {"rzx", OP::RZX, rzx_matrix(PI / 2.0, true), true, true, PI / 2.0};
                }

                return {"cx", OP::CX, cx_matrix(true), true, false, 0.0};
            }

            inline std::vector<Mat2> decomp0(const WeylDecomposition &target)
            {
                return {target.K1r * target.K2r, target.K1l * target.K2l};
            }

            inline std::vector<Mat2> decomp1(const BasisContext &basis, const WeylDecomposition &target)
            {
                Mat2 k2rd = basis.basis.K2r.conjugate().transpose();
                Mat2 k2ld = basis.basis.K2l.conjugate().transpose();
                Mat2 k1rd = basis.basis.K1r.conjugate().transpose();
                Mat2 k1ld = basis.basis.K1l.conjugate().transpose();
                return {
                    k2rd * target.K2r,
                    k2ld * target.K2l,
                    target.K1r * k1rd,
                    target.K1l * k1ld};
            }

            inline std::vector<Mat2> decomp2(const BasisContext &basis, const WeylDecomposition &target)
            {
                Mat2 rz_b = Mat2::Zero();
                rz_b << std::exp(Complex(0.0, -target.b)),
                    0.0,
                    0.0,
                    std::exp(Complex(0.0, target.b));
                Mat2 rz_a = Mat2::Zero();
                rz_a << std::exp(Complex(0.0, -target.a)),
                    0.0,
                    0.0,
                    std::exp(Complex(0.0, target.a));

                return {
                    basis.q2r * target.K2r,
                    basis.q2l * target.K2l,
                    basis.q1ra * rz_b * basis.q1rb,
                    basis.q1la * rz_a.conjugate() * basis.q1lb,
                    target.K1r * basis.q0r,
                    target.K1l * basis.q0l};
            }

            inline std::vector<Mat2> decomp3(const BasisContext &basis, const WeylDecomposition &target)
            {
                Mat2 rz_b = Mat2::Zero();
                rz_b << std::exp(Complex(0.0, -target.b)),
                    0.0,
                    0.0,
                    std::exp(Complex(0.0, target.b));
                Mat2 rz_a = Mat2::Zero();
                rz_a << std::exp(Complex(0.0, -target.a)),
                    0.0,
                    0.0,
                    std::exp(Complex(0.0, target.a));
                Mat2 rz_c = Mat2::Zero();
                rz_c << std::exp(Complex(0.0, -target.c)),
                    0.0,
                    0.0,
                    std::exp(Complex(0.0, target.c));

                return {
                    basis.u3r * target.K2r,
                    basis.u3l * target.K2l,
                    basis.u2ra * rz_b * basis.u2rb,
                    basis.u2la * rz_a.conjugate() * basis.u2lb,
                    basis.u1ra * rz_c.conjugate() * basis.u1rb,
                    basis.u1l,
                    target.K1r * basis.u0r,
                    target.K1l * basis.u0l};
            }

            inline bool synthesize_basis_block(const Mat4 &unitary,
                                               IdxType q0,
                                               IdxType q1,
                                               const std::shared_ptr<Chip> &chip,
                                               const std::unordered_set<std::string> *basis_gates,
                                               std::vector<Gate> &out,
                                               int &num_2q,
                                               int *fail_reason)
            {
                set_fail_reason(fail_reason, SynthFailReason::None);
                const BasisSpec basis_spec = choose_basis_spec(basis_gates);

                static bool cx_ready = false;
                static bool cx_failed = false;
                static BasisContext cx_ctx;
                static bool cz_ready = false;
                static bool cz_failed = false;
                static BasisContext cz_ctx;
                static bool iswap_ready = false;
                static bool iswap_failed = false;
                static BasisContext iswap_ctx;
                static bool ecr_ready = false;
                static bool ecr_failed = false;
                static BasisContext ecr_ctx;
                static bool rxx_ready = false;
                static bool rxx_failed = false;
                static BasisContext rxx_ctx;
                static bool ryy_ready = false;
                static bool ryy_failed = false;
                static BasisContext ryy_ctx;
                static bool rzz_ready = false;
                static bool rzz_failed = false;
                static BasisContext rzz_ctx;
                static bool rzx_ready = false;
                static bool rzx_failed = false;
                static BasisContext rzx_ctx;

                auto build_once = [&](bool &ready, bool &failed, BasisContext &ctx) -> const BasisContext * {
                    if (!ready && !failed)
                    {
                        try
                        {
                            ctx = build_basis_context(basis_spec.matrix);
                            ready = true;
                        }
                        catch (const std::exception &e)
                        {
                            failed = true;
                            std::cerr << "[2q_synth] basis_context_failed basis=" << basis_spec.name
                                      << " error=" << e.what() << std::endl;
                        }
                    }
                    if (!ready || failed)
                    {
                        return nullptr;
                    }
                    return &ctx;
                };

                const BasisContext *basis = nullptr;
                if (basis_spec.name == "cx")
                    basis = build_once(cx_ready, cx_failed, cx_ctx);
                else if (basis_spec.name == "cz")
                    basis = build_once(cz_ready, cz_failed, cz_ctx);
                else if (basis_spec.name == "iswap")
                    basis = build_once(iswap_ready, iswap_failed, iswap_ctx);
                else if (basis_spec.name == "ecr")
                    basis = build_once(ecr_ready, ecr_failed, ecr_ctx);
                else if (basis_spec.name == "rxx")
                    basis = build_once(rxx_ready, rxx_failed, rxx_ctx);
                else if (basis_spec.name == "ryy")
                    basis = build_once(ryy_ready, ryy_failed, ryy_ctx);
                else if (basis_spec.name == "rzz")
                    basis = build_once(rzz_ready, rzz_failed, rzz_ctx);
                else if (basis_spec.name == "rzx")
                    basis = build_once(rzx_ready, rzx_failed, rzx_ctx);

                if (!basis)
                {
                    set_fail_reason(fail_reason, SynthFailReason::BasisContext);
                    return false;
                }

                WeylDecomposition target;
                try
                {
                    target = weyl_decompose(unitary);
                }
                catch (const std::exception &e)
                {
                    const std::string msg = e.what();
                    if (msg.find("diagonalize") != std::string::npos)
                    {
                        set_fail_reason(fail_reason, SynthFailReason::WeylDiag);
                    }
                    else if (msg.find("k1/k2") != std::string::npos)
                    {
                        set_fail_reason(fail_reason, SynthFailReason::WeylK);
                    }
                    else
                    {
                        set_fail_reason(fail_reason, SynthFailReason::WeylOther);
                    }
                    return false;
                }
                const int best_nbasis = best_basis_count(target.a, target.b, target.c, basis->basis_b, 1.0);
                num_2q = best_nbasis;
                std::vector<Mat2> decomposition;
                switch (best_nbasis)
                {
                case 0:
                    decomposition = decomp0(target);
                    break;
                case 1:
                    decomposition = decomp1(*basis, target);
                    break;
                case 2:
                    decomposition = decomp2(*basis, target);
                    break;
                case 3:
                default:
                    decomposition = decomp3(*basis, target);
                    break;
                }

                bool pulse_only = false;
                optimize::detail::EulerBasis pulse_basis = optimize::detail::EulerBasis::ZSX;
                if (basis_spec.name == "cx")
                {
                    bool has_sx = basis_gates && basis_gates->find("sx") != basis_gates->end();
                    bool has_rz = basis_gates && basis_gates->find("rz") != basis_gates->end();
                    bool has_x = basis_gates && basis_gates->find("x") != basis_gates->end();
                    if (has_sx && has_rz)
                    {
                        pulse_only = true;
                        pulse_basis = has_x ? optimize::detail::EulerBasis::ZSXX : optimize::detail::EulerBasis::ZSX;
                    }
                }

                std::vector<Gate> synthesized;
                synthesized.reserve(32);
                for (int i = 0; i < best_nbasis; ++i)
                {
                    // Decomposition order is (right, left); apply left to q0 and right to q1.
                    std::vector<Gate> left = synthesize_one_qubit(decomposition[2 * i + 1], q0, chip, basis_gates, pulse_only, pulse_basis);
                    std::vector<Gate> right = synthesize_one_qubit(decomposition[2 * i], q1, chip, basis_gates, pulse_only, pulse_basis);
                    synthesized.insert(synthesized.end(), left.begin(), left.end());
                    synthesized.insert(synthesized.end(), right.begin(), right.end());
                    Gate basis_gate(basis_spec.op, q1, q0, -1, 2);
                    if (basis_spec.is_param)
                    {
                        basis_gate.theta = basis_spec.fixed_theta;
                    }
                    synthesized.push_back(basis_gate);
                }
                std::vector<Gate> tail_left = synthesize_one_qubit(decomposition[2 * best_nbasis + 1], q0, chip, basis_gates, pulse_only, pulse_basis);
                std::vector<Gate> tail_right = synthesize_one_qubit(decomposition[2 * best_nbasis], q1, chip, basis_gates, pulse_only, pulse_basis);
                synthesized.insert(synthesized.end(), tail_left.begin(), tail_left.end());
                synthesized.insert(synthesized.end(), tail_right.begin(), tail_right.end());

                out = std::move(synthesized);
                return true;
            }

            inline bool gate_uses_only(const Gate &gate, IdxType q0, IdxType q1)
            {
                if (gate.qubit >= 0 && gate.qubit != q0 && gate.qubit != q1)
                {
                    return false;
                }
                if (gate.ctrl >= 0 && gate.ctrl != q0 && gate.ctrl != q1)
                {
                    return false;
                }
                if (gate.extra >= 0)
                {
                    return false;
                }
                return true;
            }

            inline bool gate_overlaps(const Gate &gate, IdxType q0, IdxType q1)
            {
                if (gate.qubit == q0 || gate.qubit == q1)
                {
                    return true;
                }
                if (gate.ctrl == q0 || gate.ctrl == q1)
                {
                    return true;
                }
                if (gate.extra == q0 || gate.extra == q1)
                {
                    return true;
                }
                return false;
            }

            inline bool is_nonunitary_gate(const Gate &gate)
            {
                return gate.op_name == OP::M || gate.op_name == OP::MA || gate.op_name == OP::RESET;
            }
        } // namespace detail_2q_synth
#endif

        inline void synthesize_2q_blocks(std::shared_ptr<Circuit> circuit,
                                         const std::shared_ptr<Chip> &chip,
                                         const std::unordered_set<std::string> *basis_gates)
        {
#ifndef QASMTRANS_USE_EIGEN
            (void)circuit;
            (void)chip;
            (void)basis_gates;
#else
            if (!circuit)
            {
                return;
            }
            std::vector<Gate> gates = circuit->get_gates();
            if (gates.size() < 2)
            {
                return;
            }

            struct Counters
            {
                std::size_t blocks_seen = 0;
                std::size_t blocks_supported = 0;
                std::size_t blocks_synthesized = 0;
                std::size_t blocks_no_improve = 0;
                std::size_t blocks_unsupported = 0;
                std::size_t blocks_fail_verify = 0;
                std::size_t blocks_rev_used = 0;
                std::size_t total_twoq_before = 0;
                std::size_t total_twoq_after = 0;
                std::size_t total_oneq_before = 0;
                std::size_t total_oneq_after = 0;
                std::size_t blocks_fail_basis = 0;
                std::size_t blocks_fail_weyl_diag = 0;
                std::size_t blocks_fail_weyl_k = 0;
                std::size_t blocks_fail_weyl_other = 0;
                std::size_t blocks_best_lt = 0;
                std::size_t blocks_best_eq = 0;
                std::size_t blocks_best_gt = 0;
                std::size_t blocks_twoq_1 = 0;
                std::size_t blocks_twoq_2 = 0;
                std::size_t blocks_twoq_3 = 0;
                std::size_t blocks_twoq_4p = 0;
            };
            Counters stats;

            std::vector<char> remove_gate(gates.size(), 0);
            std::vector<std::vector<Gate>> replacements(gates.size());
            std::vector<char> block_assigned(gates.size(), 0);

            std::size_t idx = 0;
            while (idx < gates.size())
            {
                if (remove_gate[idx] || block_assigned[idx])
                {
                    ++idx;
                    continue;
                }
                const Gate &gate = gates[idx];
                if (!detail_2q_synth::is_two_qubit_gate(gate))
                {
                    ++idx;
                    continue;
                }
                const IdxType q0 = gate.ctrl;
                const IdxType q1 = gate.qubit;
                std::vector<std::size_t> block_positions;
                block_positions.reserve(16);
                std::size_t twoq_count = 0;
                std::size_t oneq_count = 0;
                bool supported = true;
                // Extend backwards to include leading 1q gates on q0/q1.
                std::size_t back = idx;
                while (back > 0)
                {
                    std::size_t pos = back - 1;
                    if (remove_gate[pos] || block_assigned[pos])
                    {
                        break;
                    }
                    const Gate &cand = gates[pos];
                    if (detail_2q_synth::is_nonunitary_gate(cand))
                    {
                        break;
                    }
                    if (detail_2q_synth::gate_uses_only(cand, q0, q1))
                    {
                        detail_2q_synth::Mat4 mat;
                        if (!detail_2q_synth::gate_matrix_4x4(cand, q0, q1, mat))
                        {
                            supported = false;
                            break;
                        }
                        back = pos;
                        continue;
                    }
                    if (detail_2q_synth::gate_overlaps(cand, q0, q1))
                    {
                        break;
                    }
                    back = pos;
                }

                if (!supported)
                {
                    ++stats.blocks_unsupported;
                    ++idx;
                    continue;
                }

                std::size_t pos = back;
                while (pos < gates.size())
                {
                    if (remove_gate[pos] || block_assigned[pos])
                    {
                        ++pos;
                        continue;
                    }
                    const Gate &cand = gates[pos];
                    if (detail_2q_synth::gate_uses_only(cand, q0, q1))
                    {
                        detail_2q_synth::Mat4 mat;
                        if (!detail_2q_synth::gate_matrix_4x4(cand, q0, q1, mat))
                        {
                            supported = false;
                            break;
                        }
                        block_positions.push_back(pos);
                        if (detail_2q_synth::is_two_qubit_gate(cand))
                        {
                            ++twoq_count;
                        }
                        else if (cand.ctrl < 0 && cand.extra < 0)
                        {
                            ++oneq_count;
                        }
                        ++pos;
                        continue;
                    }
                    if (detail_2q_synth::gate_overlaps(cand, q0, q1))
                    {
                        break;
                    }
                    if (detail_2q_synth::is_nonunitary_gate(cand))
                    {
                        break;
                    }
                    ++pos;
                }
                ++stats.blocks_seen;
                if (!supported || twoq_count < 1 || block_positions.size() < 2)
                {
                    if (!supported)
                    {
                        ++stats.blocks_unsupported;
                    }
                    ++idx;
                    continue;
                }
                ++stats.blocks_supported;
                stats.total_twoq_before += twoq_count;
                stats.total_oneq_before += oneq_count;
                if (twoq_count == 1)
                {
                    ++stats.blocks_twoq_1;
                }
                else if (twoq_count == 2)
                {
                    ++stats.blocks_twoq_2;
                }
                else if (twoq_count == 3)
                {
                    ++stats.blocks_twoq_3;
                }
                else
                {
                    ++stats.blocks_twoq_4p;
                }

                detail_2q_synth::Mat4 unitary = detail_2q_synth::Mat4::Identity();
                for (std::size_t pos_idx = 0; pos_idx < block_positions.size(); ++pos_idx)
                {
                    std::size_t pos = block_positions[pos_idx];
                    detail_2q_synth::Mat4 mat;
                    if (!detail_2q_synth::gate_matrix_4x4(gates[pos], q0, q1, mat))
                    {
                        supported = false;
                        break;
                    }
                    unitary = mat * unitary;
                }
                if (!supported)
                {
                    ++idx;
                    continue;
                }

                std::vector<Gate> synthesized;
                int new_twoq = 0;
                int fail_reason = 0;
                bool ok = detail_2q_synth::synthesize_basis_block(unitary, q0, q1, chip, basis_gates, synthesized, new_twoq, &fail_reason);
                if (!ok)
                {
                    ++stats.blocks_unsupported;
                    switch (static_cast<detail_2q_synth::SynthFailReason>(fail_reason))
                    {
                    case detail_2q_synth::SynthFailReason::BasisContext:
                        ++stats.blocks_fail_basis;
                        break;
                    case detail_2q_synth::SynthFailReason::WeylDiag:
                        ++stats.blocks_fail_weyl_diag;
                        break;
                    case detail_2q_synth::SynthFailReason::WeylK:
                        ++stats.blocks_fail_weyl_k;
                        break;
                    case detail_2q_synth::SynthFailReason::WeylOther:
                        ++stats.blocks_fail_weyl_other;
                        break;
                    case detail_2q_synth::SynthFailReason::None:
                    default:
                        break;
                    }
                    ++idx;
                    continue;
                }

                const double overlap_tol = 1.0 - 1e-6;
                detail_2q_synth::Mat4 unitary_synth;
                bool verify_ok = detail_2q_synth::unitary_from_gates(synthesized, q0, q1, unitary_synth);
                double overlap = verify_ok ? detail_2q_synth::unitary_overlap(unitary, unitary_synth) : 0.0;
                bool reverse_attempted = false;
                bool reverse_ok = false;
                double overlap_rev = 0.0;
                if (!verify_ok || overlap < overlap_tol)
                {
                    detail_2q_synth::Mat4 unitary_rev = detail_2q_synth::Mat4::Identity();
                    bool reverse_supported = true;
                    for (auto it = block_positions.rbegin(); it != block_positions.rend(); ++it)
                    {
                        detail_2q_synth::Mat4 mat;
                        if (!detail_2q_synth::gate_matrix_4x4(gates[*it], q0, q1, mat))
                        {
                            reverse_supported = false;
                            break;
                        }
                        unitary_rev = mat * unitary_rev;
                    }
                    reverse_attempted = true;
                    if (reverse_supported)
                    {
                        std::vector<Gate> synthesized_rev;
                        int new_twoq_rev = 0;
                        int fail_reason_rev = 0;
                        bool ok_rev = detail_2q_synth::synthesize_basis_block(
                            unitary_rev, q0, q1, chip, basis_gates, synthesized_rev, new_twoq_rev, &fail_reason_rev);
                        if (ok_rev)
                        {
                            detail_2q_synth::Mat4 unitary_synth_rev;
                            bool verify_rev = detail_2q_synth::unitary_from_gates(synthesized_rev, q0, q1, unitary_synth_rev);
                            overlap_rev = verify_rev ? detail_2q_synth::unitary_overlap(unitary_rev, unitary_synth_rev) : 0.0;
                            reverse_ok = verify_rev;
                            if (overlap_rev > overlap)
                            {
                                synthesized = std::move(synthesized_rev);
                                new_twoq = new_twoq_rev;
                                overlap = overlap_rev;
                                verify_ok = verify_rev;
                                ++stats.blocks_rev_used;
                            }
                        }
                    }
                }
                if (!verify_ok || overlap < overlap_tol)
                {
                    detail_2q_synth::dump_synth_failure(
                        unitary,
                        verify_ok,
                        verify_ok ? &unitary_synth : nullptr,
                        overlap,
                        reverse_attempted,
                        reverse_ok,
                        overlap_rev,
                        q0,
                        q1,
                        gates,
                        block_positions);
                    ++stats.blocks_fail_verify;
                    ++idx;
                    continue;
                }
                std::size_t new_oneq = 0;
                for (const auto &g : synthesized)
                {
                    if (g.ctrl < 0 && g.extra < 0)
                    {
                        ++new_oneq;
                    }
                }

                if (static_cast<std::size_t>(new_twoq) < twoq_count)
                {
                    ++stats.blocks_best_lt;
                }
                else if (static_cast<std::size_t>(new_twoq) == twoq_count)
                {
                    ++stats.blocks_best_eq;
                }
                else
                {
                    ++stats.blocks_best_gt;
                }

                bool improve = (static_cast<std::size_t>(new_twoq) < twoq_count) ||
                               (static_cast<std::size_t>(new_twoq) == twoq_count && new_oneq < oneq_count);
                if (!improve)
                {
                    ++stats.blocks_no_improve;
                    stats.total_twoq_after += static_cast<std::size_t>(new_twoq);
                    stats.total_oneq_after += new_oneq;
                    ++idx;
                    continue;
                }
                stats.total_twoq_after += static_cast<std::size_t>(new_twoq);
                stats.total_oneq_after += new_oneq;
                ++stats.blocks_synthesized;

                for (std::size_t pos_idx = 0; pos_idx < block_positions.size(); ++pos_idx)
                {
                    remove_gate[block_positions[pos_idx]] = 1;
                    block_assigned[block_positions[pos_idx]] = 1;
                }
                for (auto &new_gate : synthesized)
                {
                    new_gate.inherit_logical_metadata(gate);
                }
                replacements[block_positions.front()] = std::move(synthesized);
                ++idx;
            }

            std::vector<Gate> out;
            out.reserve(gates.size());
            for (std::size_t i = 0; i < gates.size(); ++i)
            {
                if (!replacements[i].empty())
                {
                    for (auto &g : replacements[i])
                    {
                        out.push_back(g);
                    }
                }
                if (remove_gate[i])
                {
                    continue;
                }
                out.push_back(gates[i]);
            }
            circuit->set_gates(out);

            std::cout << "[2q_synth] blocks_seen=" << stats.blocks_seen
                      << " blocks_supported=" << stats.blocks_supported
                      << " blocks_synthesized=" << stats.blocks_synthesized
                      << " blocks_no_improve=" << stats.blocks_no_improve
                      << " blocks_unsupported=" << stats.blocks_unsupported
                      << " blocks_fail_verify=" << stats.blocks_fail_verify
                      << " blocks_rev_used=" << stats.blocks_rev_used
                      << " fail_basis=" << stats.blocks_fail_basis
                      << " fail_weyl_diag=" << stats.blocks_fail_weyl_diag
                      << " fail_weyl_k=" << stats.blocks_fail_weyl_k
                      << " fail_weyl_other=" << stats.blocks_fail_weyl_other
                      << " best_lt=" << stats.blocks_best_lt
                      << " best_eq=" << stats.blocks_best_eq
                      << " best_gt=" << stats.blocks_best_gt
                      << " twoq_blocks_1=" << stats.blocks_twoq_1
                      << " twoq_blocks_2=" << stats.blocks_twoq_2
                      << " twoq_blocks_3=" << stats.blocks_twoq_3
                      << " twoq_blocks_4p=" << stats.blocks_twoq_4p
                      << " twoq_before=" << stats.total_twoq_before
                      << " twoq_after=" << stats.total_twoq_after
                      << " oneq_before=" << stats.total_oneq_before
                      << " oneq_after=" << stats.total_oneq_after
                      << std::endl;
#endif
        }
    } // namespace optimize
} // namespace QASMTrans
