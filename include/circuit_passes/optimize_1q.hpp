#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

#include "../IR/chip.hpp"
#include "../IR/circuit.hpp"
#include "../IR/gate.hpp"
#include "../QASMTransPrimitives.hpp"

namespace QASMTrans
{
    namespace optimize
    {
        namespace detail
        {
            constexpr double kAngleEps = 1e-12;

            struct Mat2
            {
                std::complex<double> m[2][2];
            };

            inline Mat2 identity()
            {
                Mat2 out{};
                out.m[0][0] = {1.0, 0.0};
                out.m[0][1] = {0.0, 0.0};
                out.m[1][0] = {0.0, 0.0};
                out.m[1][1] = {1.0, 0.0};
                return out;
            }

            inline Mat2 matmul(const Mat2 &a, const Mat2 &b)
            {
                Mat2 out{};
                out.m[0][0] = a.m[0][0] * b.m[0][0] + a.m[0][1] * b.m[1][0];
                out.m[0][1] = a.m[0][0] * b.m[0][1] + a.m[0][1] * b.m[1][1];
                out.m[1][0] = a.m[1][0] * b.m[0][0] + a.m[1][1] * b.m[1][0];
                out.m[1][1] = a.m[1][0] * b.m[0][1] + a.m[1][1] * b.m[1][1];
                return out;
            }

            inline Mat2 rz_matrix(double theta)
            {
                const std::complex<double> i(0.0, 1.0);
                const std::complex<double> phase = std::exp(-i * theta / 2.0);
                Mat2 out{};
                out.m[0][0] = phase;
                out.m[0][1] = {0.0, 0.0};
                out.m[1][0] = {0.0, 0.0};
                out.m[1][1] = std::conj(phase);
                return out;
            }

            inline Mat2 rx_matrix(double theta)
            {
                const std::complex<double> i(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat2 out{};
                out.m[0][0] = {c, 0.0};
                out.m[0][1] = -i * s;
                out.m[1][0] = -i * s;
                out.m[1][1] = {c, 0.0};
                return out;
            }

            inline Mat2 ry_matrix(double theta)
            {
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                Mat2 out{};
                out.m[0][0] = {c, 0.0};
                out.m[0][1] = {-s, 0.0};
                out.m[1][0] = {s, 0.0};
                out.m[1][1] = {c, 0.0};
                return out;
            }

            inline Mat2 phase_matrix(double theta)
            {
                const std::complex<double> i(0.0, 1.0);
                Mat2 out{};
                out.m[0][0] = {1.0, 0.0};
                out.m[0][1] = {0.0, 0.0};
                out.m[1][0] = {0.0, 0.0};
                out.m[1][1] = std::exp(i * theta);
                return out;
            }

            inline Mat2 u_matrix(double theta, double phi, double lam)
            {
                const std::complex<double> i(0.0, 1.0);
                const double c = std::cos(theta / 2.0);
                const double s = std::sin(theta / 2.0);
                const std::complex<double> e_phi = std::exp(i * phi);
                const std::complex<double> e_lam = std::exp(i * lam);
                const std::complex<double> e_phi_lam = std::exp(i * (phi + lam));
                Mat2 out{};
                out.m[0][0] = {c, 0.0};
                out.m[0][1] = -e_lam * s;
                out.m[1][0] = e_phi * s;
                out.m[1][1] = e_phi_lam * c;
                return out;
            }

            inline Mat2 prx_matrix(double theta, double phi)
            {
                return matmul(rz_matrix(phi), matmul(rx_matrix(theta), rz_matrix(-phi)));
            }

            inline bool gate_matrix(const Gate &gate, Mat2 &out)
            {
                switch (gate.op_name)
                {
                case OP::ID:
                    out = identity();
                    return true;
                case OP::X:
                    out = rx_matrix(PI);
                    return true;
                case OP::Y:
                    out = matmul(rz_matrix(PI / 2.0), matmul(rx_matrix(PI), rz_matrix(-PI / 2.0)));
                    return true;
                case OP::Z:
                    out = rz_matrix(PI);
                    return true;
                case OP::H:
                {
                    const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
                    Mat2 h{};
                    h.m[0][0] = {inv_sqrt2, 0.0};
                    h.m[0][1] = {inv_sqrt2, 0.0};
                    h.m[1][0] = {inv_sqrt2, 0.0};
                    h.m[1][1] = {-inv_sqrt2, 0.0};
                    out = h;
                    return true;
                }
                case OP::S:
                    out = phase_matrix(PI / 2.0);
                    return true;
                case OP::SDG:
                    out = phase_matrix(-PI / 2.0);
                    return true;
                case OP::T:
                    out = phase_matrix(PI / 4.0);
                    return true;
                case OP::TDG:
                    out = phase_matrix(-PI / 4.0);
                    return true;
                case OP::RI:
                    out = phase_matrix(gate.theta);
                    return true;
                case OP::RX:
                    out = rx_matrix(gate.theta);
                    return true;
                case OP::RY:
                    out = ry_matrix(gate.theta);
                    return true;
                case OP::RZ:
                    out = rz_matrix(gate.theta);
                    return true;
                case OP::SX:
                    out = rx_matrix(PI / 2.0);
                    return true;
                case OP::P:
                    out = phase_matrix(gate.theta);
                    return true;
                case OP::U:
                    out = u_matrix(gate.theta, gate.phi, gate.lam);
                    return true;
                case OP::PRX:
                    out = prx_matrix(gate.theta, gate.phi);
                    return true;
                default:
                    break;
                }
                return false;
            }

            inline double mod_2pi(double angle, double atol)
            {
                double wrapped = std::fmod(angle + PI, 2.0 * PI);
                if (wrapped < 0.0)
                {
                    wrapped += 2.0 * PI;
                }
                wrapped -= PI;
                if (std::abs(wrapped - PI) < atol)
                {
                    return -PI;
                }
                return wrapped;
            }

            inline bool near_zero(double value, double atol)
            {
                return std::abs(value) < atol;
            }

            inline void params_zyz(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                const std::complex<double> det = mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0];
                const double det_arg = std::arg(det);
                phase = 0.5 * det_arg;
                theta = 2.0 * std::atan2(std::abs(mat.m[1][0]), std::abs(mat.m[0][0]));
                const double ang1 = std::arg(mat.m[1][1]);
                const double ang2 = std::arg(mat.m[1][0]);
                phi = ang1 + ang2 - det_arg;
                lam = ang1 - ang2;
            }

            inline void params_zxz(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                params_zyz(mat, theta, phi, lam, phase);
                phi += PI / 2.0;
                lam -= PI / 2.0;
            }

            inline void params_u3(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                params_zyz(mat, theta, phi, lam, phase);
                phase -= 0.5 * (phi + lam);
            }

            inline void params_u1x(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                params_zyz(mat, theta, phi, lam, phase);
                phase -= 0.5 * (theta + phi + lam);
            }

            inline void params_xyx(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                Mat2 mat_zyz{};
                mat_zyz.m[0][0] = 0.5 * (mat.m[0][0] + mat.m[0][1] + mat.m[1][0] + mat.m[1][1]);
                mat_zyz.m[0][1] = 0.5 * (mat.m[0][0] - mat.m[0][1] + mat.m[1][0] - mat.m[1][1]);
                mat_zyz.m[1][0] = 0.5 * (mat.m[0][0] + mat.m[0][1] - mat.m[1][0] - mat.m[1][1]);
                mat_zyz.m[1][1] = 0.5 * (mat.m[0][0] - mat.m[0][1] - mat.m[1][0] + mat.m[1][1]);

                params_zyz(mat_zyz, theta, phi, lam, phase);
                const double new_phi = mod_2pi(phi + PI, 0.0);
                const double new_lam = mod_2pi(lam + PI, 0.0);
                phase += (new_phi + new_lam - phi - lam) / 2.0;
                phi = new_phi;
                lam = new_lam;
            }

            inline void params_xzx(const Mat2 &mat, double &theta, double &phi, double &lam, double &phase)
            {
                const std::complex<double> det = mat.m[0][0] * mat.m[1][1] - mat.m[0][1] * mat.m[1][0];
                phase = std::log(det).imag() / 2.0;
                const std::complex<double> sqrt_det = std::sqrt(det);

                Mat2 mat_zyz{};
                mat_zyz.m[0][0] = {((mat.m[0][0] / sqrt_det).real()), ((mat.m[1][0] / sqrt_det).imag())};
                mat_zyz.m[0][1] = {((mat.m[1][0] / sqrt_det).real()), ((mat.m[0][0] / sqrt_det).imag())};
                mat_zyz.m[1][0] = {-(mat.m[1][0] / sqrt_det).real(), (mat.m[0][0] / sqrt_det).imag()};
                mat_zyz.m[1][1] = {((mat.m[0][0] / sqrt_det).real()), -(mat.m[1][0] / sqrt_det).imag()};

                double phase_zxz = 0.0;
                params_zxz(mat_zyz, theta, phi, lam, phase_zxz);
                phase += phase_zxz;
            }

            enum class AxisGate
            {
                RZ,
                RX,
                RY
            };

            inline void emit_axis_gate(std::vector<Gate> &out, AxisGate axis, double angle, IdxType qubit, double atol)
            {
                if (near_zero(angle, atol))
                {
                    return;
                }
                switch (axis)
                {
                case AxisGate::RZ:
                    out.emplace_back(OP::RZ, qubit, -1, -1, 1, angle);
                    break;
                case AxisGate::RX:
                    out.emplace_back(OP::RX, qubit, -1, -1, 1, angle);
                    break;
                case AxisGate::RY:
                    out.emplace_back(OP::RY, qubit, -1, -1, 1, angle);
                    break;
                }
            }

            inline std::vector<Gate> circuit_kak(double theta,
                                                 double phi,
                                                 double lam,
                                                 double phase,
                                                 AxisGate k_gate,
                                                 AxisGate a_gate,
                                                 IdxType qubit,
                                                 bool simplify,
                                                 double atol)
            {
                std::vector<Gate> circuit;
                double local_atol = simplify ? atol : -1.0;

                if (std::abs(theta) < local_atol)
                {
                    lam += phi;
                    lam = mod_2pi(lam, local_atol);
                    emit_axis_gate(circuit, k_gate, lam, qubit, local_atol);
                    return circuit;
                }
                if (std::abs(theta - PI) < local_atol)
                {
                    lam -= phi;
                    phi = 0.0;
                }
                if (std::abs(mod_2pi(lam + PI, local_atol)) < local_atol ||
                    std::abs(mod_2pi(phi + PI, local_atol)) < local_atol)
                {
                    lam += PI;
                    theta = -theta;
                    phi += PI;
                }
                lam = mod_2pi(lam, local_atol);
                emit_axis_gate(circuit, k_gate, lam, qubit, local_atol);
                emit_axis_gate(circuit, a_gate, theta, qubit, local_atol);
                phi = mod_2pi(phi, local_atol);
                emit_axis_gate(circuit, k_gate, phi, qubit, local_atol);
                (void)phase;
                return circuit;
            }

            inline std::vector<Gate> circuit_u(double theta,
                                               double phi,
                                               double lam,
                                               double phase,
                                               IdxType qubit,
                                               bool simplify,
                                               double atol,
                                               const std::string &custom_name)
            {
                std::vector<Gate> circuit;
                const double local_atol = simplify ? atol : -1.0;
                const double phi_mod = mod_2pi(phi, local_atol);
                const double lam_mod = mod_2pi(lam, local_atol);
                if (!simplify || std::abs(theta) > local_atol || std::abs(phi_mod) > local_atol || std::abs(lam_mod) > local_atol)
                {
                    Gate gate(OP::U, qubit, -1, -1, 1, theta, phi_mod, lam_mod);
                    if (!custom_name.empty())
                    {
                        gate.set_custom_name(custom_name);
                    }
                    circuit.push_back(gate);
                }
                (void)phase;
                return circuit;
            }

            inline std::vector<Gate> circuit_u321(double theta,
                                                  double phi,
                                                  double lam,
                                                  double phase,
                                                  IdxType qubit,
                                                  bool simplify,
                                                  double atol,
                                                  const std::unordered_set<std::string> *basis)
            {
                std::vector<Gate> circuit;
                double local_atol = simplify ? atol : -1.0;
                if (std::abs(theta) < local_atol)
                {
                    const double tot = mod_2pi(phi + lam, local_atol);
                    if (std::abs(tot) > local_atol)
                    {
                        std::string custom = "u1";
                        if (basis && basis->find("p") != basis->end())
                        {
                            circuit.emplace_back(OP::P, qubit, -1, -1, 1, tot);
                        }
                        else
                        {
                            Gate gate(OP::U, qubit, -1, -1, 1, 0.0, 0.0, tot);
                            if (!custom.empty())
                            {
                                gate.set_custom_name(custom);
                            }
                            circuit.push_back(gate);
                        }
                    }
                }
                else if (std::abs(theta - PI / 2.0) < local_atol)
                {
                    Gate gate(OP::U, qubit, -1, -1, 1, PI / 2.0, mod_2pi(phi, local_atol), mod_2pi(lam, local_atol));
                    gate.set_custom_name("u2");
                    circuit.push_back(gate);
                }
                else
                {
                    Gate gate(OP::U, qubit, -1, -1, 1, theta, mod_2pi(phi, local_atol), mod_2pi(lam, local_atol));
                    gate.set_custom_name("u3");
                    circuit.push_back(gate);
                }
                (void)phase;
                return circuit;
            }

            enum class PhaseKind
            {
                RZ,
                P
            };

            enum class XKind
            {
                SX,
                RX
            };

            inline void emit_phase(std::vector<Gate> &out, PhaseKind kind, double angle, IdxType qubit, double atol)
            {
                const double wrapped = mod_2pi(angle, atol);
                if (near_zero(wrapped, atol))
                {
                    return;
                }
                if (kind == PhaseKind::RZ)
                {
                    out.emplace_back(OP::RZ, qubit, -1, -1, 1, wrapped);
                }
                else
                {
                    out.emplace_back(OP::P, qubit, -1, -1, 1, wrapped);
                }
            }

            inline void emit_x(std::vector<Gate> &out, XKind kind, IdxType qubit)
            {
                if (kind == XKind::SX)
                {
                    out.emplace_back(OP::SX, qubit);
                }
                else
                {
                    out.emplace_back(OP::RX, qubit, -1, -1, 1, PI / 2.0);
                }
            }

            inline std::vector<Gate> circuit_psx_like(double theta,
                                                      double phi,
                                                      double lam,
                                                      double phase,
                                                      IdxType qubit,
                                                      bool simplify,
                                                      double atol,
                                                      PhaseKind p_kind,
                                                      XKind x_kind,
                                                      bool allow_x_pi)
            {
                std::vector<Gate> circuit;
                double local_atol = simplify ? atol : -1.0;

                if (std::abs(theta) < local_atol)
                {
                    emit_phase(circuit, p_kind, lam + phi, qubit, local_atol);
                    return circuit;
                }
                if (std::abs(theta - PI / 2.0) < local_atol)
                {
                    emit_phase(circuit, p_kind, lam - PI / 2.0, qubit, local_atol);
                    emit_x(circuit, x_kind, qubit);
                    emit_phase(circuit, p_kind, phi + PI / 2.0, qubit, local_atol);
                    return circuit;
                }
                if (std::abs(theta - PI) < local_atol)
                {
                    phi -= lam;
                    lam = 0.0;
                }
                if (std::abs(mod_2pi(lam + PI, local_atol)) < local_atol ||
                    std::abs(mod_2pi(phi, local_atol)) < local_atol)
                {
                    lam += PI;
                    theta = -theta;
                    phi += PI;
                }
                theta += PI;
                phi += PI;

                emit_phase(circuit, p_kind, lam, qubit, local_atol);
                if (allow_x_pi && std::abs(mod_2pi(theta, local_atol)) < local_atol)
                {
                    circuit.emplace_back(OP::X, qubit);
                }
                else
                {
                    emit_x(circuit, x_kind, qubit);
                    emit_phase(circuit, p_kind, theta, qubit, local_atol);
                    emit_x(circuit, x_kind, qubit);
                }
                emit_phase(circuit, p_kind, phi, qubit, local_atol);
                (void)phase;
                return circuit;
            }

            inline std::vector<Gate> circuit_rr(double theta,
                                                double phi,
                                                double lam,
                                                double phase,
                                                IdxType qubit,
                                                bool simplify,
                                                double atol)
            {
                std::vector<Gate> circuit;
                double local_atol = simplify ? atol : -1.0;
                if (std::abs(mod_2pi((phi + lam) / 2.0, local_atol)) < local_atol)
                {
                    if (std::abs(theta) > local_atol)
                    {
                        circuit.emplace_back(OP::PRX, qubit, -1, -1, 1, theta, mod_2pi(PI / 2.0 + phi, local_atol));
                    }
                }
                else
                {
                    if (std::abs(theta - PI) > local_atol)
                    {
                        circuit.emplace_back(OP::PRX, qubit, -1, -1, 1, theta - PI, mod_2pi(PI / 2.0 - lam, local_atol));
                    }
                    circuit.emplace_back(OP::PRX, qubit, -1, -1, 1, PI, mod_2pi(0.5 * (phi - lam + PI), local_atol));
                }
                (void)phase;
                return circuit;
            }

            inline bool basis_contains(const std::unordered_set<std::string> &basis, const std::string &gate)
            {
                if (basis.find(gate) != basis.end())
                {
                    return true;
                }
                if (gate == "u1")
                {
                    return basis.find("p") != basis.end() || basis.find("u") != basis.end();
                }
                if (gate == "u2" || gate == "u3")
                {
                    return basis.find("u") != basis.end();
                }
                if (gate == "u")
                {
                    return basis.find("u") != basis.end() || basis.find("u3") != basis.end() ||
                           basis.find("u2") != basis.end() || basis.find("u1") != basis.end();
                }
                if (gate == "p")
                {
                    return basis.find("u1") != basis.end();
                }
                if (gate == "r")
                {
                    return basis.find("prx") != basis.end();
                }
                return false;
            }

            inline bool basis_supports(const std::unordered_set<std::string> *basis,
                                       const std::vector<std::string> &required)
            {
                if (!basis || basis->empty())
                {
                    return true;
                }
                for (const auto &gate : required)
                {
                    if (!basis_contains(*basis, gate))
                    {
                        return false;
                    }
                }
                return true;
            }

            enum class EulerBasis
            {
                U3,
                U321,
                U,
                PSX,
                U1X,
                RR,
                ZYZ,
                ZXZ,
                XZX,
                XYX,
                ZSXX,
                ZSX
            };

            struct BasisSpec
            {
                EulerBasis basis;
                const char *name;
                std::vector<std::string> required;
            };

            inline std::vector<EulerBasis> possible_bases(const std::unordered_set<std::string> *basis)
            {
                static const std::vector<BasisSpec> kSpecs = {
                    {EulerBasis::U3, "U3", {"u3"}},
                    {EulerBasis::U321, "U321", {"u3", "u2", "u1"}},
                    {EulerBasis::U, "U", {"u"}},
                    {EulerBasis::PSX, "PSX", {"p", "sx"}},
                    {EulerBasis::U1X, "U1X", {"u1", "rx"}},
                    {EulerBasis::RR, "RR", {"r"}},
                    {EulerBasis::ZYZ, "ZYZ", {"rz", "ry"}},
                    {EulerBasis::ZXZ, "ZXZ", {"rz", "rx"}},
                    {EulerBasis::XZX, "XZX", {"rx", "rz"}},
                    {EulerBasis::XYX, "XYX", {"rx", "ry"}},
                    {EulerBasis::ZSXX, "ZSXX", {"rz", "sx", "x"}},
                    {EulerBasis::ZSX, "ZSX", {"rz", "sx"}},
                };

                std::vector<EulerBasis> out;
                out.reserve(kSpecs.size());
                for (const auto &spec : kSpecs)
                {
                    if (basis_supports(basis, spec.required))
                    {
                        out.push_back(spec.basis);
                    }
                }
                const bool has_zsxx = std::find(out.begin(), out.end(), EulerBasis::ZSXX) != out.end();
                if (has_zsxx)
                {
                    out.erase(std::remove(out.begin(), out.end(), EulerBasis::ZSX), out.end());
                }
                const bool has_u321 = std::find(out.begin(), out.end(), EulerBasis::U321) != out.end();
                if (has_u321)
                {
                    out.erase(std::remove(out.begin(), out.end(), EulerBasis::U3), out.end());
                }
                return out;
            }

            inline void angles_from_unitary(const Mat2 &unitary,
                                            EulerBasis basis,
                                            double &theta,
                                            double &phi,
                                            double &lam,
                                            double &phase)
            {
                switch (basis)
                {
                case EulerBasis::U3:
                case EulerBasis::U321:
                case EulerBasis::U:
                    params_u3(unitary, theta, phi, lam, phase);
                    break;
                case EulerBasis::PSX:
                case EulerBasis::ZSX:
                case EulerBasis::ZSXX:
                case EulerBasis::U1X:
                    params_u1x(unitary, theta, phi, lam, phase);
                    break;
                case EulerBasis::RR:
                case EulerBasis::ZYZ:
                    params_zyz(unitary, theta, phi, lam, phase);
                    break;
                case EulerBasis::ZXZ:
                    params_zxz(unitary, theta, phi, lam, phase);
                    break;
                case EulerBasis::XYX:
                    params_xyx(unitary, theta, phi, lam, phase);
                    break;
                case EulerBasis::XZX:
                    params_xzx(unitary, theta, phi, lam, phase);
                    break;
                }
            }

            inline std::vector<Gate> generate_circuit(EulerBasis basis,
                                                      double theta,
                                                      double phi,
                                                      double lam,
                                                      double phase,
                                                      IdxType qubit,
                                                      const std::unordered_set<std::string> *basis_gates)
            {
                switch (basis)
                {
                case EulerBasis::ZYZ:
                    return circuit_kak(theta, phi, lam, phase, AxisGate::RZ, AxisGate::RY, qubit, true, kAngleEps);
                case EulerBasis::ZXZ:
                    return circuit_kak(theta, phi, lam, phase, AxisGate::RZ, AxisGate::RX, qubit, true, kAngleEps);
                case EulerBasis::XZX:
                    return circuit_kak(theta, phi, lam, phase, AxisGate::RX, AxisGate::RZ, qubit, true, kAngleEps);
                case EulerBasis::XYX:
                    return circuit_kak(theta, phi, lam, phase, AxisGate::RX, AxisGate::RY, qubit, true, kAngleEps);
                case EulerBasis::U3:
                    return circuit_u(theta, phi, lam, phase, qubit, true, kAngleEps,
                                     basis_gates && basis_gates->find("u") == basis_gates->end() ? "u3" : "");
                case EulerBasis::U:
                    return circuit_u(theta, phi, lam, phase, qubit, true, kAngleEps, "");
                case EulerBasis::U321:
                    return circuit_u321(theta, phi, lam, phase, qubit, true, kAngleEps, basis_gates);
                case EulerBasis::PSX:
                    return circuit_psx_like(theta, phi, lam, phase, qubit, true, kAngleEps, PhaseKind::P, XKind::SX, false);
                case EulerBasis::U1X:
                    return circuit_psx_like(theta, phi, lam, phase, qubit, true, kAngleEps, PhaseKind::P, XKind::RX, false);
                case EulerBasis::ZSX:
                    return circuit_psx_like(theta, phi, lam, phase, qubit, true, kAngleEps, PhaseKind::RZ, XKind::SX, false);
                case EulerBasis::ZSXX:
                    return circuit_psx_like(theta, phi, lam, phase, qubit, true, kAngleEps, PhaseKind::RZ, XKind::SX, true);
                case EulerBasis::RR:
                    return circuit_rr(theta, phi, lam, phase, qubit, true, kAngleEps);
                }
                return {};
            }

            inline bool is_single_qubit_gate(const Gate &gate)
            {
                return gate.ctrl < 0 && gate.extra < 0 && gate.qubit >= 0;
            }

            inline bool has_error_data(const std::shared_ptr<Chip> &chip, IdxType qubit)
            {
                if (!chip)
                {
                    return false;
                }
                if (qubit < 0 || qubit >= static_cast<IdxType>(chip->single_qubit_errors.size()))
                {
                    return false;
                }
                return !chip->single_qubit_errors[static_cast<std::size_t>(qubit)].empty();
            }

            inline double lookup_error(const std::shared_ptr<Chip> &chip, IdxType qubit, const std::string &gate_name)
            {
                if (!chip)
                {
                    return 0.0;
                }
                if (qubit < 0 || qubit >= static_cast<IdxType>(chip->single_qubit_errors.size()))
                {
                    return 0.0;
                }
                const auto &err_map = chip->single_qubit_errors[static_cast<std::size_t>(qubit)];
                auto it = err_map.find(gate_name);
                if (it == err_map.end())
                {
                    return 0.0;
                }
                return it->second;
            }

            inline double compute_sequence_error(const std::vector<Gate> &sequence,
                                                 const std::shared_ptr<Chip> &chip,
                                                 IdxType qubit,
                                                 bool use_errors)
            {
                if (!use_errors)
                {
                    return static_cast<double>(sequence.size());
                }
                double fidelity = 1.0;
                for (const auto &gate : sequence)
                {
                    const std::string name = gate.lower_name();
                    fidelity *= (1.0 - lookup_error(chip, qubit, name));
                }
                return 1.0 - fidelity;
            }

            inline bool gate_in_basis(const Gate &gate, const std::unordered_set<std::string> *basis)
            {
                if (!basis || basis->empty())
                {
                    return true;
                }
                const std::string name = gate.lower_name();
                return basis->find(name) != basis->end();
            }

            inline IdxType compute_qubit_capacity(const std::vector<Gate> &gates)
            {
                IdxType max_index = -1;
                for (const auto &gate : gates)
                {
                    if (gate.qubit > max_index)
                        max_index = gate.qubit;
                    if (gate.ctrl > max_index)
                        max_index = gate.ctrl;
                    if (gate.extra > max_index)
                        max_index = gate.extra;
                }
                return max_index + 1;
            }
        } // namespace detail

        inline void optimize_1q_gates_decomposition(std::shared_ptr<Circuit> circuit,
                                                    const std::shared_ptr<Chip> &chip,
                                                    const std::unordered_set<std::string> *basis_gates)
        {
            if (!circuit)
            {
                return;
            }
            const std::vector<Gate> gates = circuit->get_gates();
            if (gates.empty())
            {
                return;
            }

            const IdxType capacity = detail::compute_qubit_capacity(gates);
            if (capacity <= 0)
            {
                return;
            }

            std::vector<std::vector<IdxType>> current_runs(static_cast<std::size_t>(capacity));
            std::vector<std::vector<IdxType>> runs;
            runs.reserve(gates.size());

            auto flush_run = [&](IdxType qubit)
            {
                if (qubit < 0 || qubit >= capacity)
                {
                    return;
                }
                auto &run = current_runs[static_cast<std::size_t>(qubit)];
                if (!run.empty())
                {
                    runs.push_back(run);
                    run.clear();
                }
            };

            for (IdxType idx = 0; idx < static_cast<IdxType>(gates.size()); ++idx)
            {
                const Gate &gate = gates[static_cast<std::size_t>(idx)];
                std::vector<IdxType> touched;
                if (gate.qubit >= 0)
                    touched.push_back(gate.qubit);
                if (gate.ctrl >= 0)
                    touched.push_back(gate.ctrl);
                if (gate.extra >= 0)
                    touched.push_back(gate.extra);
                std::sort(touched.begin(), touched.end());
                touched.erase(std::unique(touched.begin(), touched.end()), touched.end());

                for (IdxType qubit : touched)
                {
                    const bool eligible = detail::is_single_qubit_gate(gate) && gate.qubit == qubit && !gate.has_custom_name();
                    if (!eligible)
                    {
                        flush_run(qubit);
                        continue;
                    }
                    detail::Mat2 mat;
                    if (!detail::gate_matrix(gate, mat))
                    {
                        flush_run(qubit);
                        continue;
                    }
                    current_runs[static_cast<std::size_t>(qubit)].push_back(idx);
                }
            }

            for (IdxType qubit = 0; qubit < capacity; ++qubit)
            {
                flush_run(qubit);
            }

            if (runs.empty())
            {
                return;
            }

            const bool use_errors = detail::has_error_data(chip, 0) ||
                                    std::any_of(runs.begin(), runs.end(), [&](const auto &run) {
                                        if (run.empty())
                                            return false;
                                        IdxType q = gates[static_cast<std::size_t>(run.front())].qubit;
                                        return detail::has_error_data(chip, q);
                                    });

            std::vector<char> remove_gate(gates.size(), 0);
            std::vector<std::vector<Gate>> replacements(gates.size());

            const std::vector<detail::EulerBasis> bases = detail::possible_bases(basis_gates);

            for (const auto &run : runs)
            {
                if (run.empty())
                {
                    continue;
                }
                const Gate &first_gate = gates[static_cast<std::size_t>(run.front())];
                const IdxType qubit = first_gate.qubit;

                detail::Mat2 unitary = detail::identity();
                bool supported = true;
                for (IdxType gate_idx : run)
                {
                    detail::Mat2 mat;
                    if (!detail::gate_matrix(gates[static_cast<std::size_t>(gate_idx)], mat))
                    {
                        supported = false;
                        break;
                    }
                    unitary = detail::matmul(mat, unitary);
                }
                if (!supported)
                {
                    continue;
                }

                double best_error = std::numeric_limits<double>::infinity();
                std::vector<Gate> best_seq;
                for (detail::EulerBasis basis : bases)
                {
                    double theta = 0.0, phi = 0.0, lam = 0.0, phase = 0.0;
                    detail::angles_from_unitary(unitary, basis, theta, phi, lam, phase);
                    std::vector<Gate> candidate = detail::generate_circuit(basis, theta, phi, lam, phase, qubit, basis_gates);
                    const double err = detail::compute_sequence_error(candidate, chip, qubit, use_errors);
                    if (err < best_error || (std::abs(err - best_error) < 1e-12 && candidate.size() < best_seq.size()))
                    {
                        best_error = err;
                        best_seq = std::move(candidate);
                    }
                }

                if (best_seq.empty() && bases.empty())
                {
                    continue;
                }

                const double old_error = detail::compute_sequence_error([
                    &]() -> std::vector<Gate> {
                        std::vector<Gate> seq;
                        seq.reserve(run.size());
                        for (IdxType gate_idx : run)
                        {
                            seq.push_back(gates[static_cast<std::size_t>(gate_idx)]);
                        }
                        return seq;
                    }(),
                    chip,
                    qubit,
                    use_errors);

                bool outside_basis = false;
                if (basis_gates && !basis_gates->empty())
                {
                    for (IdxType gate_idx : run)
                    {
                        if (!detail::gate_in_basis(gates[static_cast<std::size_t>(gate_idx)], basis_gates))
                        {
                            outside_basis = true;
                            break;
                        }
                    }
                }

                const bool old_is_identity = std::abs(old_error) < 1e-12 && !run.empty();
                const bool new_is_identity = best_seq.empty();
                const bool replace = outside_basis ||
                                     (best_error < old_error) ||
                                     (new_is_identity && !old_is_identity);
                if (!replace)
                {
                    continue;
                }

                for (IdxType gate_idx : run)
                {
                    remove_gate[static_cast<std::size_t>(gate_idx)] = 1;
                }

                for (auto &gate : best_seq)
                {
                    gate.inherit_logical_metadata(first_gate);
                }
                replacements[static_cast<std::size_t>(run.front())] = std::move(best_seq);
            }

            std::vector<Gate> out;
            out.reserve(gates.size());
            for (std::size_t idx = 0; idx < gates.size(); ++idx)
            {
                if (!replacements[idx].empty())
                {
                    for (auto &gate : replacements[idx])
                    {
                        out.push_back(gate);
                    }
                }
                if (remove_gate[idx])
                {
                    continue;
                }
                out.push_back(gates[idx]);
            }

            circuit->set_gates(out);
        }

        inline void consolidate_single_qubit_chains(std::shared_ptr<Circuit> circuit)
        {
            optimize_1q_gates_decomposition(circuit, nullptr, nullptr);
        }
    } // namespace optimize
} // namespace QASMTrans
