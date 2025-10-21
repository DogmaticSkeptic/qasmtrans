#pragma once

#include "../QASMTransPrimitives.hpp"
#include "../IR/gate.hpp"

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>

using namespace QASMTrans;
using namespace std;

namespace QASMTrans
{
    extern std::unordered_set<std::string> g_device_basis_gates;
extern std::unordered_map<std::string, std::string> g_merged_gate_aliases;
}

inline std::string toLowerCase(const char *name)
{
    std::string lower = name ? std::string(name) : std::string{};
    std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });
    return lower;
}

inline std::string lookupMergedAlias(const Gate &gate)
{
    if (g_merged_gate_aliases.empty())
    {
        return {};
    }
    std::vector<std::string> logical_candidates;
    logical_candidates.push_back(gate.lower_base_name());
    logical_candidates.push_back(gate.lower_name());
    if (!gate.logical_label.empty())
    {
        std::string label = gate.logical_label;
        std::transform(label.begin(), label.end(), label.begin(), [](unsigned char c)
                       { return static_cast<char>(std::tolower(c)); });
        logical_candidates.push_back(std::move(label));
    }
    for (const auto &candidate : logical_candidates)
    {
        if (candidate.empty())
        {
            continue;
        }
        auto it = g_merged_gate_aliases.find(candidate);
        if (it != g_merged_gate_aliases.end())
        {
            return it->second;
        }
    }
    return {};
}

Gate BasicRZ(ValType theta, IdxType qubit)
{
    Gate G(OP::RZ, qubit, -1, -1, 1, theta);
    return G;
}
Gate BasicSX(IdxType qubit)
{
    Gate G(OP::SX, qubit);
    return G;
}
Gate BasicX(IdxType qubit)
{
    Gate G(OP::X, qubit);
    return G;
}
Gate BasicCX(IdxType ctrl, IdxType qubit)
{
    Gate G(OP::CX, qubit, ctrl, -1, 2);
    return G;
}
vector<Gate> decomposeHadamard(IdxType qubit)
{
    vector<Gate> decomposedGates;
    // Assuming the Gate constructor takes name, control qubit, target qubit, and angle (in that order)
    // Also assuming single qubit gates use -1 or similar for non-applicable target/control qubits
    // sx x rz(-pi/2) sx x
    Gate rzgate = BasicRZ(-PI / 2, qubit);
    Gate sxgate = BasicSX(qubit);
    Gate xgate = BasicX(qubit);
    decomposedGates.push_back(xgate);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(rzgate);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(xgate);
    return decomposedGates;
}
vector<Gate> decomposeT(IdxType qubit)
{
    vector<Gate> decomposedGates;
    // Applies the π/8 gate to a single qubit.
    Gate rzgate = BasicRZ(PI / 4, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}
vector<Gate> decomposeTdg(IdxType qubit)
{
    vector<Gate> decomposedGates;
    // Applies the π/8 gate to a single qubit.
    Gate rzgate = BasicRZ(-PI / 4, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}
vector<Gate> decomposeZ(IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate rzgate = BasicRZ(PI, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}
vector<Gate> decomposeY(IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate sxgate = BasicSX(qubit);
    decomposedGates.push_back(sxgate);
    vector<Gate> decomposez = decomposeZ(qubit);
    decomposedGates.insert(decomposedGates.end(), decomposez.begin(), decomposez.end());
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(sxgate);
    return decomposedGates;
}
vector<Gate> decomposeRx(ValType theta, IdxType qubit)
{
    vector<Gate> decomposedGates;
    vector<Gate> decomposeh = decomposeHadamard(qubit);
    decomposedGates.insert(decomposedGates.end(), decomposeh.begin(), decomposeh.end());
    Gate rzgate = BasicRZ(theta, qubit);
    decomposedGates.push_back(rzgate);
    decomposedGates.insert(decomposedGates.end(), decomposeh.begin(), decomposeh.end());
    return decomposedGates;
}
vector<Gate> decomposeP(ValType theta, IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate rzgate = BasicRZ(theta, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}

vector<Gate> decomposeRI(ValType theta, IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate rzgate = BasicRZ(2 * theta, qubit);
    decomposedGates.push_back(rzgate);
    vector<Gate> decomposeh = decomposeZ(qubit);
    decomposedGates.insert(decomposedGates.end(), decomposeh.begin(), decomposeh.end());
    return decomposedGates;
}
vector<Gate> decomposeRy(ValType theta, IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate sxgate = BasicSX(qubit);
    Gate rzgate = BasicRZ(theta, qubit);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(rzgate);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(sxgate);
    decomposedGates.push_back(sxgate);
    return decomposedGates;
}
vector<Gate> decomposeS(IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate rzgate = BasicRZ(PI / 2, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}
vector<Gate> decomposeSdg(IdxType qubit)
{
    vector<Gate> decomposedGates;
    Gate rzgate = BasicRZ(-PI / 2, qubit);
    decomposedGates.push_back(rzgate);
    return decomposedGates;
}
vector<Gate> decomposeU(ValType theta, ValType phi, ValType lam, IdxType qubit)
{
    vector<Gate> decomposedGates;
    if (lam != 0)
    {
        decomposedGates.push_back(BasicRZ(lam, qubit));
    }
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(theta + PI, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(3 * PI + phi, qubit));

    return decomposedGates;
}
vector<Gate> decomposeCZ(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    vector<Gate> decomposery = decomposeHadamard(qubit);
    decomposedGates.insert(decomposedGates.end(), decomposery.begin(), decomposery.end());
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.insert(decomposedGates.end(), decomposery.begin(), decomposery.end());
    return decomposedGates;
}
vector<Gate> decomposeCY(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(-PI / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCH(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(-PI, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI * 3 / 4, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 4, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    return decomposedGates;
}
vector<Gate> decomposeCS(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 4, ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI / 4, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 4, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCSDG(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicRZ(PI / 4, ctrl));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI / 4, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 4, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCT(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 8, ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI / 8, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 8, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCTDG(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(-PI / 8, ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 8, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI / 8, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCRX(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicRZ(theta / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-theta / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    return decomposedGates;
}
vector<Gate> decomposeRXX(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, ctrl));
    decomposedGates.push_back(BasicSX(ctrl));
    decomposedGates.push_back(BasicRZ(PI / 2, ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(theta, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, ctrl));
    decomposedGates.push_back(BasicSX(ctrl));
    decomposedGates.push_back(BasicRZ(PI / 2, ctrl));
    return decomposedGates;
}
vector<Gate> decomposeRYY(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicSX(ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(theta, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(-PI, qubit));
    decomposedGates.push_back(BasicRZ(-PI, ctrl));
    decomposedGates.push_back(BasicSX(ctrl));
    decomposedGates.push_back(BasicRZ(-PI, ctrl));
    return decomposedGates;
}
vector<Gate> decomposeRZZ(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(theta, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCRY(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI + theta / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(3 * PI, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI - theta / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(3 * PI, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCRZ(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    // cout<<"theta is"<<theta<<endl;
    // cout<<"phi is"<<phi<<endl;
    // cout<<"lam is"<<lam<<endl;
    // cout<<"gamma is"<<gamma<<endl;
    decomposedGates.push_back(BasicRZ(theta / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-theta / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCSX(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicRZ(PI / 4, ctrl));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-PI / 4, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(3 * PI / 4, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI / 2, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCP(ValType theta, IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicRZ(theta / 2, ctrl));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-theta / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(theta / 2, qubit));
    return decomposedGates;
}
vector<Gate> decomposeCU(ValType theta, ValType phi, ValType lam, ValType gamma,
                         IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    // cout<<"theta is"<<theta<<endl;
    // cout<<"phi is"<<phi<<endl;
    // cout<<"lam is"<<lam<<endl;
    // cout<<"gamma is"<<gamma<<endl;
    decomposedGates.push_back(BasicRZ(gamma, ctrl));
    decomposedGates.push_back(BasicRZ(lam / 2 + phi / 2, ctrl));
    decomposedGates.push_back(BasicRZ(lam / 2 - phi / 2, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicRZ(-lam / 2 - phi / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI - theta / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(3 * PI, qubit));
    decomposedGates.push_back(BasicCX(ctrl, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(PI + theta / 2, qubit));
    decomposedGates.push_back(BasicSX(qubit));
    decomposedGates.push_back(BasicRZ(3 * PI + phi, qubit));
    return decomposedGates;
}
vector<Gate> decomposeSWAP(IdxType qubit, IdxType ctrl)
{
    vector<Gate> decomposedGates;
    Gate cxgate = BasicCX(ctrl, qubit);
    Gate cxgate2 = BasicCX(qubit, ctrl);
    decomposedGates.push_back(cxgate);
    decomposedGates.push_back(cxgate2);
    decomposedGates.push_back(cxgate);
    return decomposedGates;
}
// gate ccx a,b,c
// {
//   h c;
//   cx b,c; tdg c;
//   cx a,c; t c;
//   cx b,c; tdg c;
//   cx a,c; t b; t c; h c;
//   cx a,b; t a; tdg b;
//   cx a,b;
// }
vector<Gate> decomposeCCX(IdxType a, IdxType b, IdxType c)
{
    vector<Gate> decomposedGates;

    decomposedGates.push_back(Gate(OP::H, c));
    decomposedGates.push_back(BasicCX(b, c));
    decomposedGates.push_back(Gate(OP::TDG, c));
    decomposedGates.push_back(BasicCX(a, c));
    decomposedGates.push_back(Gate(OP::T, c));
    decomposedGates.push_back(BasicCX(b, c));
    decomposedGates.push_back(Gate(OP::TDG, c));
    decomposedGates.push_back(BasicCX(a, c));
    decomposedGates.push_back(Gate(OP::T, b));
    decomposedGates.push_back(Gate(OP::T, c));
    decomposedGates.push_back(Gate(OP::H, c));
    decomposedGates.push_back(BasicCX(a, b));
    decomposedGates.push_back(Gate(OP::T, a));
    decomposedGates.push_back(Gate(OP::TDG, b));
    decomposedGates.push_back(BasicCX(a, b));
    return decomposedGates;
}

vector<Gate> decomposeRCCX(IdxType a, IdxType b, IdxType c)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, PI / 2, 0, PI));
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, 0, 0, PI / 4));
    decomposedGates.push_back(BasicCX(b, c));
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, 0, 0, -PI / 4));
    decomposedGates.push_back(BasicCX(a, c));
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, 0, 0, PI / 4));
    decomposedGates.push_back(BasicCX(b, c));
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, 0, 0, -PI / 4));
    decomposedGates.push_back(Gate(OP::U, c, -1, -1, 1, PI / 2, 0, PI));
    return decomposedGates;
}
vector<Gate> decomposeCSWAP(IdxType a, IdxType b, IdxType c)
{
    vector<Gate> decomposedGates;
    decomposedGates.push_back(BasicCX(c, b));
    vector<Gate> decomposeccx = decomposeCCX(a, b, c);
    decomposedGates.insert(decomposedGates.end(), decomposeccx.begin(), decomposeccx.end());
    decomposedGates.push_back(BasicCX(c, b));
    return decomposedGates;
}
void Decompose_three_to_two(shared_ptr<Circuit> circuit)
{
    vector<Gate> circuit_gates = circuit->get_gates();
    vector<Gate> decomposedGates;
    for (Gate g : circuit_gates)
    {
        std::string gate_name_lower = g.lower_name();
        if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
        {
            decomposedGates.push_back(g);
            continue;
        }
        if (g.n_qubits > 2)
        {
            // std::cout<<"find three-qubit gates"<<std::endl;
            // print gate and control, target , extra qubit
            // std::cout<<"gate name is"<<OP_NAMES[g.op_name];
            // std::cout<<"gate control is"<<g.ctrl;
            // std::cout<<"gate target is"<<g.qubit;
            // std::cout<<"gate extra is"<<g.extra<<std::endl;
            if (g.name_equals("CSWAP"))
            {
                vector<Gate> Decomposed_gates = decomposeCSWAP(g.qubit, g.ctrl, g.extra);
                decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
            }
            else if (g.name_equals("CCX"))
            {
                vector<Gate> Decomposed_gates = decomposeCCX(g.qubit, g.ctrl, g.extra);
                decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
            }
            else if (g.name_equals("RCCX"))
            {
                vector<Gate> Decomposed_gates = decomposeRCCX(g.qubit, g.ctrl, g.extra);
                decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
            }
        }
        else
        {
            decomposedGates.push_back(g);
        }
    }

    circuit->set_gates(decomposedGates);
    // print circuit gates
    //  if (debug_level > 1) {
    //      vector<Gate> circuit_gates2 = circuit->get_gates();
    //      for (Gate g : circuit_gates2)
    //      {
    //          std::cout<<"gate name is"<<OP_NAMES[g.op_name];
    //          std::cout<<"gate control is"<<g.ctrl;
    //          std::cout<<"gate target is"<<g.qubit;
    //          std::cout<<"gate extra is"<<g.extra;
    //          //print angle
    //          std::cout<<"gate angle is"<<g.theta<<","<<g.phi<<","<<g.lam<<std::endl;

    //     }
    // }

    return;
}
void Decompose(shared_ptr<Circuit> circuit, IdxType mode)
{
    vector<Gate> circuit_gates = circuit->get_gates();
    vector<Gate> decomposedGates;
    for (Gate g : circuit_gates)
    {
        auto append_with_metadata = [&](const std::vector<Gate> &src)
        {
            size_t local_before = decomposedGates.size();
            decomposedGates.insert(decomposedGates.end(), src.begin(), src.end());
            size_t after = decomposedGates.size();
            for (size_t idx = local_before; idx < after; ++idx)
            {
                decomposedGates[idx].inherit_logical_metadata(g);
            }
        };
        auto push_with_metadata = [&](Gate gate)
        {
            gate.inherit_logical_metadata(g);
            decomposedGates.push_back(gate);
        };
        std::string alias = lookupMergedAlias(g);
        if (!alias.empty())
        {
            g.set_custom_name(alias);
            push_with_metadata(g);
            continue;
        }
        std::string gate_name_lower = g.lower_name();
        if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
        {
            push_with_metadata(g);
            continue;
        }
        size_t before_size = decomposedGates.size();
        if (g.name_equals("H"))
        {
            vector<Gate> Decomposed_gates = decomposeHadamard(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("T"))
        {
            vector<Gate> Decomposed_gates = decomposeT(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("Z"))
        {
            vector<Gate> Decomposed_gates = decomposeZ(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("TDG"))
        {
            vector<Gate> Decomposed_gates = decomposeTdg(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("Y"))
        {
            vector<Gate> Decomposed_gates = decomposeY(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("S"))
        {
            vector<Gate> Decomposed_gates = decomposeS(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("SDG"))
        {
            vector<Gate> Decomposed_gates = decomposeSdg(g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("RX"))
        {
            vector<Gate> Decomposed_gates = decomposeRx(g.theta, g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("RY"))
        {
            vector<Gate> Decomposed_gates = decomposeRy(g.theta, g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());

        } //^ double check RI gate later
        else if (g.name_equals("RI"))
        {
            vector<Gate> Decomposed_gates = decomposeRI(g.theta, g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("P"))
        {
            vector<Gate> Decomposed_gates = decomposeP(g.theta, g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("U"))
        {
            vector<Gate> Decomposed_gates = decomposeU(g.theta, g.phi, g.lam, g.qubit);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CZ"))
        {
            vector<Gate> Decomposed_gates = decomposeCZ(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CY"))
        {
            vector<Gate> Decomposed_gates = decomposeCY(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CH"))
        {
            vector<Gate> Decomposed_gates = decomposeCH(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CS"))
        {
            vector<Gate> Decomposed_gates = decomposeCS(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CSDG"))
        {
            vector<Gate> Decomposed_gates = decomposeCSDG(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CT"))
        {
            vector<Gate> Decomposed_gates = decomposeCT(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CTDG"))
        {
            vector<Gate> Decomposed_gates = decomposeCTDG(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CRX"))
        {
            vector<Gate> Decomposed_gates = decomposeCRX(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CRY"))
        {
            vector<Gate> Decomposed_gates = decomposeCRY(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CRZ"))
        {
            vector<Gate> Decomposed_gates = decomposeCRZ(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CSX"))
        {
            vector<Gate> Decomposed_gates = decomposeCSX(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CP"))
        {
            vector<Gate> Decomposed_gates = decomposeCP(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CU"))
        {
            vector<Gate> Decomposed_gates = decomposeCU(g.theta, g.phi, g.lam, g.gamma, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("RXX"))
        {
            vector<Gate> Decomposed_gates = decomposeRXX(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("RYY"))
        {
            vector<Gate> Decomposed_gates = decomposeRYY(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("RZZ"))
        {
            vector<Gate> Decomposed_gates = decomposeRZZ(g.theta, g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("SWAP"))
        {
            vector<Gate> Decomposed_gates = decomposeSWAP(g.qubit, g.ctrl);
            decomposedGates.insert(decomposedGates.end(), Decomposed_gates.begin(), Decomposed_gates.end());
        }
        else if (g.name_equals("CX"))
        {
            decomposedGates.push_back(g);
        }
        else if (g.name_equals("RZ"))
        {
            // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
            decomposedGates.push_back(BasicRZ(g.theta, g.qubit));
        }
        else if (g.name_equals("SX"))
        {
            decomposedGates.push_back(g);
        }
        else if (g.name_equals("X"))
        {
            decomposedGates.push_back(g);
        }
        else if (g.name_equals("MA"))
        {
            decomposedGates.push_back(g);
        }
        else if (g.name_equals("ID"))
        {
            decomposedGates.push_back(g);
        }
        else if (g.name_equals("RESET"))
        {
            decomposedGates.push_back(g);
        }
        else
        {
            cout << "Error: cannot find this gate: " << endl;
            cout << "Gate " << g.name() << " not supported" << endl;
            decomposedGates.push_back(g);
        }
        size_t after_size = decomposedGates.size();
        for (size_t idx = before_size; idx < after_size; ++idx)
        {
            decomposedGates[idx].inherit_logical_metadata(g);
        }
    }
    if (mode == 0)
    {
        circuit->set_gates(decomposedGates);
        return;
    }
    else if (mode == 1)
    {
        vector<Gate> decomposedGates_IonQ;
        for (Gate g : decomposedGates)
        {
            size_t before_size = decomposedGates_IonQ.size();
            if (g.has_custom_name())
            {
                decomposedGates_IonQ.push_back(g);
                size_t after_size = decomposedGates_IonQ.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_IonQ[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string alias = lookupMergedAlias(g);
            if (!alias.empty())
            {
                g.set_custom_name(alias);
                decomposedGates_IonQ.push_back(g);
                size_t after_size = decomposedGates_IonQ.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_IonQ[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string gate_name_lower = g.lower_name();
            if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
            {
                decomposedGates_IonQ.push_back(g);
                size_t after_size = decomposedGates_IonQ.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_IonQ[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            if (g.name_equals("RZ"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_IonQ.push_back(BasicRZ(g.theta, g.qubit));
            }
            else if (g.name_equals("SX"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_IonQ.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI / 2));
            }
            else if (g.name_equals("X"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_IonQ.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI));
            }
            else if (g.name_equals("CX"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_IonQ.push_back(Gate(OP::RY, g.qubit, -1, -1, 1, PI / 2));
                decomposedGates_IonQ.push_back(Gate(OP::RXX, g.qubit, g.ctrl, -1, 2, PI / 2));
                decomposedGates_IonQ.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, -PI / 2));
                decomposedGates_IonQ.push_back(Gate(OP::RX, g.ctrl, -1, -1, 1, -PI / 2));
                decomposedGates_IonQ.push_back(Gate(OP::RY, g.qubit, -1, -1, 1, -PI / 2));
            }
            size_t after_size = decomposedGates_IonQ.size();
            for (size_t idx = before_size; idx < after_size; ++idx)
            {
                decomposedGates_IonQ[idx].inherit_logical_metadata(g);
            }
        }
        circuit->set_gates(decomposedGates_IonQ);
        return;
    }
    else if (mode == 2)
    {
        vector<Gate> decomposedGates_Quantinuum;
        for (Gate g : decomposedGates)
        {
            size_t before_size = decomposedGates_Quantinuum.size();
            if (g.has_custom_name())
            {
                decomposedGates_Quantinuum.push_back(g);
                size_t after_size = decomposedGates_Quantinuum.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quantinuum[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string alias = lookupMergedAlias(g);
            if (!alias.empty())
            {
                g.set_custom_name(alias);
                decomposedGates_Quantinuum.push_back(g);
                size_t after_size = decomposedGates_Quantinuum.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quantinuum[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string gate_name_lower = g.lower_name();
            if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
            {
                decomposedGates_Quantinuum.push_back(g);
                size_t after_size = decomposedGates_Quantinuum.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quantinuum[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            if (g.name_equals("RZ"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Quantinuum.push_back(BasicRZ(g.theta, g.qubit));
            }
            else if (g.name_equals("SX"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Quantinuum.push_back(Gate(OP::U, g.qubit, -1, -1, 1, PI / 2));
            }
            else if (g.name_equals("X"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Quantinuum.push_back(Gate(OP::U, g.qubit, -1, -1, 1, PI));
            }
            else if (g.name_equals("CX"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Quantinuum.push_back(Gate(OP::U, g.qubit, -1, -1, 1, -PI / 2, PI / 2));
                decomposedGates_Quantinuum.push_back(Gate(OP::ZZ, g.qubit, g.ctrl, -1, 2, PI / 2));
                decomposedGates_Quantinuum.push_back(Gate(OP::RZ, g.ctrl, -1, -1, 1, -PI / 2));
                decomposedGates_Quantinuum.push_back(Gate(OP::U, g.qubit, -1, -1, 1, PI / 2, PI));
                decomposedGates_Quantinuum.push_back(Gate(OP::RZ, g.ctrl, -1, -1, 1, -PI / 2));
            }
            size_t after_size = decomposedGates_Quantinuum.size();
            for (size_t idx = before_size; idx < after_size; ++idx)
            {
                decomposedGates_Quantinuum[idx].inherit_logical_metadata(g);
            }
        }
        circuit->set_gates(decomposedGates_Quantinuum);
        return;
    }
    else if (mode == 3)
    {
        vector<Gate> decomposedGates_Rigetti;
        for (Gate g : decomposedGates)
        {
            size_t before_size = decomposedGates_Rigetti.size();
            if (g.has_custom_name())
            {
                decomposedGates_Rigetti.push_back(g);
                size_t after_size = decomposedGates_Rigetti.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Rigetti[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string alias = lookupMergedAlias(g);
            if (!alias.empty())
            {
                g.set_custom_name(alias);
                decomposedGates_Rigetti.push_back(g);
                size_t after_size = decomposedGates_Rigetti.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Rigetti[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string gate_name_lower = g.lower_name();
            if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
            {
                decomposedGates_Rigetti.push_back(g);
                size_t after_size = decomposedGates_Rigetti.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Rigetti[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            if (g.name_equals("RZ"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Rigetti.push_back(BasicRZ(g.theta, g.qubit));
            }
            else if (g.name_equals("SX"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Rigetti.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI / 2));
            }
            else if (g.name_equals("X"))
            {
                // cout<<"gate name is"<<OP_NAMES[g.op_name]<<"angle is"<<g.theta<<endl;
                decomposedGates_Rigetti.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI));
            }
            else if (g.name_equals("CX"))
            {
                IdxType ctrl = g.ctrl;
                IdxType target = g.qubit;
                decomposedGates_Rigetti.push_back(Gate(OP::RZ, target, -1, -1, 1, PI / 2));
                decomposedGates_Rigetti.push_back(Gate(OP::ISWAP, target, ctrl, -1, 2));
                decomposedGates_Rigetti.push_back(Gate(OP::RX, ctrl, -1, -1, 1, -PI / 2));
                decomposedGates_Rigetti.push_back(Gate(OP::RZ, ctrl, -1, -1, 1, PI / 2));
                decomposedGates_Rigetti.push_back(Gate(OP::ISWAP, target, ctrl, -1, 2));
                decomposedGates_Rigetti.push_back(Gate(OP::RX, target, -1, -1, 1, -PI / 2));
                decomposedGates_Rigetti.push_back(Gate(OP::RZ, ctrl, -1, -1, 1, PI / 2));
            }
            size_t after_size = decomposedGates_Rigetti.size();
            for (size_t idx = before_size; idx < after_size; ++idx)
            {
                decomposedGates_Rigetti[idx].inherit_logical_metadata(g);
            }
        }
        circuit->set_gates(decomposedGates_Rigetti);
        return;
    }
    else if (mode == 4)
    {
        vector<Gate> decomposedGates_Quafu;
        for (Gate g : decomposedGates)
        {
            size_t before_size = decomposedGates_Quafu.size();
            if (g.has_custom_name())
            {
                decomposedGates_Quafu.push_back(g);
                size_t after_size = decomposedGates_Quafu.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quafu[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string alias = lookupMergedAlias(g);
            if (!alias.empty())
            {
                g.set_custom_name(alias);
                decomposedGates_Quafu.push_back(g);
                size_t after_size = decomposedGates_Quafu.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quafu[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            std::string gate_name_lower = g.lower_name();
            if (g_device_basis_gates.find(gate_name_lower) != g_device_basis_gates.end())
            {
                decomposedGates_Quafu.push_back(g);
                size_t after_size = decomposedGates_Quafu.size();
                for (size_t idx = before_size; idx < after_size; ++idx)
                {
                    decomposedGates_Quafu[idx].inherit_logical_metadata(g);
                }
                continue;
            }
            if (g.name_equals("RZ"))
            {
                decomposedGates_Quafu.push_back(BasicRZ(g.theta, g.qubit));
            }
            else if (g.name_equals("SX"))
            {
                decomposedGates_Quafu.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI / 2));
            }
            else if (g.name_equals("X"))
            {
                decomposedGates_Quafu.push_back(Gate(OP::RX, g.qubit, -1, -1, 1, PI));
            }
            else if (g.name_equals("CX"))
            {
                decomposedGates_Quafu.push_back(Gate(OP::H, g.qubit));
                decomposedGates_Quafu.push_back(Gate(OP::CZ, g.qubit, g.ctrl, 2));
                decomposedGates_Quafu.push_back(Gate(OP::H, g.qubit));
            }
            size_t after_size = decomposedGates_Quafu.size();
            for (size_t idx = before_size; idx < after_size; ++idx)
            {
                decomposedGates_Quafu[idx].inherit_logical_metadata(g);
            }
        }
        circuit->set_gates(decomposedGates_Quafu);
        return;
    }
}
