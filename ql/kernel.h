/**
 * @file   kernel.h
 * @date   04/2017
 * @author Nader Khammassi
 *         Imran Ashraf
 * @brief  openql kernel
 */

#ifndef QL_KERNEL_H
#define QL_KERNEL_H

#include <sstream>
#include <algorithm>
#include <iterator>

#include "ql/json.h"
#include "ql/utils.h"
#include "ql/options.h"
#include "ql/gate.h"
#include "ql/optimizer.h"
#include "ql/ir.h"

#define PI M_PI

#include "mapper.h"

#ifndef __disable_lemon__
#include "scheduler.h"
#endif // __disable_lemon__

namespace ql
{
// un-comment it to decompose
// #define DECOMPOSE

/**
 * quantum_kernel
 */
class quantum_kernel
{
public:

    quantum_kernel(std::string name) : name(name), iterations(1) {}

    quantum_kernel(std::string name, ql::quantum_platform& platform) : name(name), iterations(1)
    {
        gate_definition = platform.instruction_map;
        qubit_number = platform.qubit_number;
        cycle_time = platform.cycle_time;
    }

    void loop(size_t it)
    {
        iterations = it;
    }

    void identity(size_t qubit)
    {
        gate("identity", {qubit} );
    }

    void i(size_t qubit)
    {
        gate("identity", {qubit} );
    }

    void hadamard(size_t qubit)
    {
        gate("hadamard", {qubit} );
    }

    void h(size_t qubit)
    {
        hadamard(qubit);
    }

    void rx(size_t qubit, double angle)
    {
        std::string gname("rx");
        // to do : rotation decomposition
        c.push_back(new ql::rx(qubit,angle));
    }

    void ry(size_t qubit, double angle)
    {
        std::string gname("ry");
        // to do : rotation decomposition
        c.push_back(new ql::ry(qubit,angle));
    }

    void rz(size_t qubit, double angle)
    {
        std::string gname("rz");
        // to do : rotation decomposition
        c.push_back(new ql::rz(qubit,angle));
    }

    void s(size_t qubit)
    {
        gate("s", {qubit} );
    }

    void sdag(size_t qubit)
    {
        gate("sdag", {qubit} );
    }

    void t(size_t qubit)
    {
        gate("t", {qubit} );
    }

    void tdag(size_t qubit)
    {
        gate("tdag", {qubit} );
    }

    void x(size_t qubit)
    {
        gate("x", {qubit} );
    }

    void y(size_t qubit)
    {
        gate("y", {qubit} );
    }

    void z(size_t qubit)
    {
        gate("z", {qubit} );
    }

    void rx90(size_t qubit)
    {
        gate("rx90", {qubit} );
    }

    void mrx90(size_t qubit)
    {
        gate("mrx90", {qubit} );
    }

    void rx180(size_t qubit)
    {
        gate("rx180", {qubit} );
    }

    void ry90(size_t qubit)
    {
        gate("ry90", {qubit} );
    }

    void mry90(size_t qubit)
    {
        gate("mry90", {qubit} );
    }

    void ry180(size_t qubit)
    {
        gate("ry180", {qubit} );
    }

    void measure(size_t qubit)
    {
        gate("measure", {qubit} );
    }

    void prepz(size_t qubit)
    {
        gate("prepz", {qubit} );
    }

    void cnot(size_t qubit1, size_t qubit2)
    {
        gate("cnot", {qubit1, qubit2} );
    }

    void cz(size_t qubit1, size_t qubit2)
    {
        gate("cz", {qubit1, qubit2} );
    }

    void cphase(size_t qubit1, size_t qubit2)
    {
        gate("cphase", {qubit1, qubit2} );
    }

    void toffoli(size_t qubit1, size_t qubit2, size_t qubit3)
    {
        // TODO add custom gate check if needed
        c.push_back(new ql::toffoli(qubit1, qubit2, qubit3));
    }

    void wait(std::vector<size_t> qubits, size_t duration)
    {
        gate("wait", qubits, duration );
    }

    void display()
    {
        c.push_back(new ql::display());
    }

    /**
     * add clifford
     */
    void clifford(int id, size_t qubit=0)
    {
        switch (id)
        {
        case 0 :
            break;                                          //  ['I']
        case 1 :
            ry90(qubit);
            rx90(qubit);
            break;                //  ['Y90', 'X90']
        case 2 :
            mrx90(qubit);
            mry90(qubit);
            break;              //  ['mX90', 'mY90']
        case 3 :
            rx180(qubit);
            break;                            //  ['X180']
        case 4 :
            mry90(qubit);
            mrx90(qubit);
            break;              //  ['mY90', 'mX90']
        case 5 :
            rx90(qubit);
            mry90(qubit);
            break;               //  ['X90', 'mY90']
        case 6 :
            ry180(qubit);
            break;                            //  ['Y180']
        case 7 :
            mry90(qubit);
            rx90(qubit);
            break;               //  ['mY90', 'X90']
        case 8 :
            rx90(qubit);
            ry90(qubit);
            break;                //  ['X90', 'Y90']
        case 9 :
            rx180(qubit);
            ry180(qubit);
            break;              //  ['X180', 'Y180']
        case 10:
            ry90(qubit);
            mrx90(qubit);
            break;               //  ['Y90', 'mX90']
        case 11:
            mrx90(qubit);
            ry90(qubit);
            break;               //  ['mX90', 'Y90']
        case 12:
            ry90(qubit);
            rx180(qubit);
            break;               //  ['Y90', 'X180']
        case 13:
            mrx90(qubit);
            break;                            //  ['mX90']
        case 14:
            rx90(qubit);
            mry90(qubit);
            mrx90(qubit);
            break; //  ['X90', 'mY90', 'mX90']
        case 15:
            mry90(qubit);
            break;                            //  ['mY90']
        case 16:
            rx90(qubit);
            break;                             //  ['X90']
        case 17:
            rx90(qubit);
            ry90(qubit);
            rx90(qubit);
            break;   //  ['X90', 'Y90', 'X90']
        case 18:
            mry90(qubit);
            rx180(qubit);
            break;              //  ['mY90', 'X180']
        case 19:
            rx90(qubit);
            ry180(qubit);
            break;               //  ['X90', 'Y180']
        case 20:
            rx90(qubit);
            mry90(qubit);
            rx90(qubit);
            break;  //  ['X90', 'mY90', 'X90']
        case 21:
            ry90(qubit);
            break;                             //  ['Y90']
        case 22:
            mrx90(qubit);
            ry180(qubit);
            break;              //  ['mX90', 'Y180']
        case 23:
            rx90(qubit);
            ry90(qubit);
            mrx90(qubit);
            break;  //  ['X90', 'Y90', 'mX90']
        default:
            break;
        }
    }

    bool add_default_gate_if_available(std::string gname, std::vector<size_t> qubits,
                                       size_t duration=0, double angle=0.0)
    {
        bool result=false;

        bool is_one_qubit_gate = (gname == "identity") || (gname == "i")
                                 || (gname == "hadamard") || (gname == "h")
                                 || (gname == "pauli_x") || (gname == "pauli_y") || (gname == "pauli_z")
                                 || (gname == "x") || (gname == "y") || (gname == "z")
                                 || (gname == "s") || (gname == "sdag")
                                 || (gname == "t") || (gname == "tdag")
                                 || (gname == "rx") || (gname == "ry") || (gname == "rz")
                                 || (gname == "rx90") || (gname == "mrx90") || (gname == "rx180")
                                 || (gname == "ry90") || (gname == "mry90") || (gname == "ry180")
                                 || (gname == "measure") || (gname == "prepz");

        bool is_two_qubit_gate = (gname == "cnot")
                                 || (gname == "cz") || (gname == "cphase")
                                 || (gname == "swap");

        bool is_multi_qubit_gate = (gname == "toffoli")
                                   || (gname == "wait") || (gname == "barrier");

        if(is_one_qubit_gate)
        {
            if( qubits.size() != 1 )
                return false;
        }
        else if(is_two_qubit_gate)
        {
            if( qubits.size() != 2 )
                return false;
            if( qubits[0] == qubits[1] )
                return false;
        }
        else if(is_multi_qubit_gate)
        {
            // by default wait will be applied to all qubits
        }
        else
        {
            return false;
        }

        if( gname == "identity" || gname == "i" )
        {
            c.push_back(new ql::identity(qubits[0]) );
            result = true;
        }
        else if( gname == "hadamard" || gname == "h" )
        {
            c.push_back(new ql::hadamard(qubits[0]) );
            result = true;
        }
        else if( gname == "pauli_x" || gname == "x" )
        {
            c.push_back(new ql::pauli_x(qubits[0]) );
            result = true;
        }
        else if( gname == "pauli_y" || gname == "y" )
        {
            c.push_back(new ql::pauli_y(qubits[0]) );
            result = true;
        }
        else if( gname == "pauli_z" || gname == "z" )
        {
            c.push_back(new ql::pauli_z(qubits[0]) );
            result = true;
        }
        else if( gname == "s" || gname == "phase" )
        {
            c.push_back(new ql::phase(qubits[0]) );
            result = true;
        }
        else if( gname == "sdag" || gname == "phasedag" )
        {
            c.push_back(new ql::phasedag(qubits[0]) );
            result = true;
        }
        else if( gname == "t" )
        {
            c.push_back(new ql::t(qubits[0]) );
            result = true;
        }
        else if( gname == "tdag" )
        {
            c.push_back(new ql::tdag(qubits[0]) );
            result = true;
        }
        else if( gname == "rx" )
        {
            c.push_back(new ql::rx(qubits[0], angle));
            result = true;
        }
        else if( gname == "ry" )
        {
            c.push_back(new ql::ry(qubits[0], angle));
            result = true;
        }
        else if( gname == "rz" )
        {
            c.push_back(new ql::rz(qubits[0], angle));
            result = true;
        }
        else if( gname == "rx90" )
        {
            c.push_back(new ql::rx90(qubits[0]) );
            result = true;
        }
        else if( gname == "mrx90" )
        {
            c.push_back(new ql::mrx90(qubits[0]) );
            result = true;
        }
        else if( gname == "rx180" )
        {
            c.push_back(new ql::rx180(qubits[0]) );
            result = true;
        }
        else if( gname == "ry90" )
        {
            c.push_back(new ql::ry90(qubits[0]) );
            result = true;
        }
        else if( gname == "mry90" )
        {
            c.push_back(new ql::mry90(qubits[0]) );
            result = true;
        }
        else if( gname == "ry180" )
        {
            c.push_back(new ql::ry180(qubits[0]) );
            result = true;
        }
        else if( gname == "measure" )
        {
            c.push_back(new ql::measure(qubits[0]) );
            result = true;
        }
        else if( gname == "prepz" )
        {
            c.push_back(new ql::prepz(qubits[0]) );
            result = true;
        }
        else if( gname == "cnot" )
        {
            c.push_back(new ql::cnot(qubits[0], qubits[1]) );
            result = true;
        }
        else if( gname == "cz" || gname == "cphase" )
        {
            c.push_back(new ql::cphase(qubits[0], qubits[1]) );
            result = true;
        }
        else if( gname == "toffoli" )
        {
            c.push_back(new ql::toffoli(qubits[0], qubits[1], qubits[2]) );
            result = true;
        }
        else if( gname == "swap" )
        {
            c.push_back(new ql::swap(qubits[0], qubits[1]) );
            result = true;
        }
        else if( gname == "barrier")
        {
            c.push_back(new ql::wait(qubits, 0, 0));
            result = true;
        }
        else if( gname == "wait")
        {
            size_t duration_in_cycles = std::ceil(static_cast<float>(duration)/cycle_time);
            c.push_back(new ql::wait(qubits, duration, duration_in_cycles));
            result = true;
        }
        else result = false;

        return result;
    }

    bool add_custom_gate_if_available(std::string & gname, std::vector<size_t> qubits,
                                      size_t duration=0, double angle=0.0)
    {
        bool added = false;
        // first check if a specialized custom gate is available
        std::string instr = gname + " ";
        if(qubits.size() > 0)
        {
            for (size_t i=0; i<(qubits.size()-1); ++i)
                instr += "q" + std::to_string(qubits[i]) + ",";
            if(qubits.size() >= 1) // to make if work with gates without operands
                instr += "q" + std::to_string(qubits[qubits.size()-1]);
        }

        std::map<std::string,custom_gate*>::iterator it = gate_definition.find(instr);
        if (it != gate_definition.end())
        {
            custom_gate* g = new custom_gate(*(it->second));
            for(auto & qubit : qubits)
                g->operands.push_back(qubit);
            if(duration>0) g->duration = duration;
            g->angle = angle;
            added = true;
            c.push_back(g);
        }
        else
        {
            // otherwise, check if there is a parameterized custom gate (i.e. not specialized for arguments)
            std::map<std::string,custom_gate*>::iterator it = gate_definition.find(gname);
            if (it != gate_definition.end())
            {
                custom_gate* g = new custom_gate(*(it->second));
                for(auto & qubit : qubits)
                    g->operands.push_back(qubit);
                if(duration>0) g->duration = duration;
                g->angle = angle;
                added = true;
                c.push_back(g);
            }
        }

        if(added)
        {
            DOUT("custom gate added for " << gname);
        }
        else
        {
            DOUT("custom gate not added for " << gname);
        }

        return added;
    }

    void get_decomposed_ins( ql::composite_gate * gptr, std::vector<std::string> & sub_instructons )
    {
        auto & sub_gates = gptr->gs;
        DOUT("composite ins: " << gptr->name);
        for(auto & agate : sub_gates)
        {
            std::string & sub_ins = agate->name;
            DOUT("  sub ins: " << sub_ins);
            auto it = gate_definition.find(sub_ins);
            if( it != gate_definition.end() )
            {
                sub_instructons.push_back(sub_ins);
            }
            else
            {
                throw ql::exception("[x] error : ql::kernel::gate() : gate decomposition not available for '"+sub_ins+"'' in the target platform !",false);
            }
        }
    }

    bool add_spec_decomposed_gate_if_available(std::string gate_name, std::vector<size_t> all_qubits)
    {
        bool added = false;
        DOUT("Checking if specialized decomposition is available for " << gate_name);
        std::string instr_parameterized = gate_name + " ";
        size_t i;
        if(all_qubits.size() > 0)
        {
            for(i=0; i<all_qubits.size()-1; i++)
            {
                instr_parameterized += "q" + std::to_string(all_qubits[i]) + " ";
            }
            if(all_qubits.size() >= 1)
            {
                instr_parameterized += "q" + std::to_string(all_qubits[i]);
            }
        }
        DOUT("decomposed specialized instruction name: " << instr_parameterized);

        auto it = gate_definition.find(instr_parameterized);
        if( it != gate_definition.end() )
        {
            DOUT("specialized composite gate found for " << instr_parameterized);
            composite_gate * gptr = (composite_gate *)(it->second);
            if( __composite_gate__ == gptr->type() )
            {
                DOUT("composite gate type");
            }
            else
            {
                DOUT("Not a composite gate type");
                return false;
            }


            std::vector<std::string> sub_instructons;
            get_decomposed_ins( gptr, sub_instructons );
            for(auto & sub_ins : sub_instructons)
            {
                DOUT("Adding sub ins: " << sub_ins);
                std::replace( sub_ins.begin(), sub_ins.end(), ',', ' ');
                DOUT(" after comma removal, sub ins: " << sub_ins);
                std::istringstream iss(sub_ins);

                std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},
                                                 std::istream_iterator<std::string>{} };

                std::vector<size_t> this_gate_qubits;
                std::string & sub_ins_name = tokens[0];

                for(size_t i=1; i<tokens.size(); i++)
                {
                    DOUT("tokens[i] : " << tokens[i]);
                    auto sub_str_token = tokens[i].substr(1);
                    DOUT("sub_str_token[i] : " << sub_str_token);
                    this_gate_qubits.push_back( stoi( tokens[i].substr(1) ) );
                }

                DOUT( ql::utils::to_string<size_t>(this_gate_qubits, "actual qubits of this gate:") );

                // custom gate check
                bool custom_added = add_custom_gate_if_available(sub_ins_name, this_gate_qubits);
                if(!custom_added)
                {
                    if(ql::options::get("use_default_gates") == "yes")
                    {
                        // default gate check
                        DOUT("adding default gate for " << sub_ins_name);
                        bool default_available = add_default_gate_if_available(sub_ins_name, this_gate_qubits);
                        if( default_available )
                        {
                            WOUT("added default gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                        }
                        else
                        {
                            EOUT("unknown gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                            throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+sub_ins_name+"' with " +ql::utils::to_string(this_gate_qubits,"qubits")+" is not supported by the target platform !",false);
                        }
                    }
                    else
                    {
                        EOUT("unknown gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                        throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+sub_ins_name+"' with " +ql::utils::to_string(this_gate_qubits,"qubits")+" is not supported by the target platform !",false);
                    }
                }
            }
            added = true;
        }
        else
        {
            DOUT("composite gate not found for " << instr_parameterized);
        }

        return added;
    }


    bool add_param_decomposed_gate_if_available(std::string gate_name, std::vector<size_t> all_qubits)
    {
        bool added = false;
        DOUT("Checking if parameterized decomposition is available for " << gate_name);
        std::string instr_parameterized = gate_name + " ";
        size_t i;
        if(all_qubits.size() > 0)
        {
            for(i=0; i<all_qubits.size()-1; i++)
            {
                instr_parameterized += "%" + std::to_string(i) + " ";
            }
            if(all_qubits.size() >= 1)
            {
                instr_parameterized += "%" + std::to_string(i);
            }
        }
        DOUT("decomposed parameterized instruction name: " << instr_parameterized);

        // check for composite ins
        auto it = gate_definition.find(instr_parameterized);
        if( it != gate_definition.end() )
        {
            DOUT("parameterized composite gate found for " << instr_parameterized);
            composite_gate * gptr = (composite_gate *)(it->second);
            if( __composite_gate__ == gptr->type() )
            {
                DOUT("composite gate type");
            }
            else
            {
                DOUT("Not a composite gate type");
                return false;
            }

            std::vector<std::string> sub_instructons;
            get_decomposed_ins( gptr, sub_instructons );
            for(auto & sub_ins : sub_instructons)
            {
                DOUT("Adding sub ins: " << sub_ins);
                std::replace( sub_ins.begin(), sub_ins.end(), ',', ' ');
                DOUT(" after comma removal, sub ins: " << sub_ins);
                std::istringstream iss(sub_ins);

                std::vector<std::string> tokens{ std::istream_iterator<std::string>{iss},
                                                 std::istream_iterator<std::string>{} };

                std::vector<size_t> this_gate_qubits;
                std::string & sub_ins_name = tokens[0];

                for(size_t i=1; i<tokens.size(); i++)
                {
                    this_gate_qubits.push_back( all_qubits[ stoi( tokens[i].substr(1) ) ] );
                }

                DOUT( ql::utils::to_string<size_t>(this_gate_qubits, "actual qubits of this gate:") );

                // custom gate check
                bool custom_added = add_custom_gate_if_available(sub_ins_name, this_gate_qubits);
                if(!custom_added)
                {
                    if(ql::options::get("use_default_gates") == "yes")
                    {
                        // default gate check
                        DOUT("adding default gate for " << sub_ins_name);
                        bool default_available = add_default_gate_if_available(sub_ins_name, this_gate_qubits);
                        if( default_available )
                        {
                            WOUT("added default gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                        }
                        else
                        {
                            EOUT("unknown gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                            throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+sub_ins_name+"' with " +ql::utils::to_string(this_gate_qubits,"qubits")+" is not supported by the target platform !",false);
                        }
                    }
                    else
                    {
                        EOUT("unknown gate '" << sub_ins_name << "' with " << ql::utils::to_string(this_gate_qubits,"qubits") );
                        throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+sub_ins_name+"' with " +ql::utils::to_string(this_gate_qubits,"qubits")+" is not supported by the target platform !",false);
                    }
                }
            }
            added = true;
        }
        else
        {
            DOUT("composite gate not found for " << instr_parameterized);
        }
        return added;
    }


    /**
     * custom 1 qubit gate.
     */
    void gate(std::string gname, size_t q0)
    {
        gate(gname, std::vector<size_t> {q0});
    }

    /**
     * custom 2 qubits gate
     */
    void gate(std::string gname, size_t q0, size_t q1)
    {
        gate(gname, std::vector<size_t> {q0, q1});
    }


    /**
     * custom gate with arbitrary number of operands
     */
    void gate(std::string gname, std::vector<size_t> qubits = {}, size_t duration=0, double angle = 0.0)
    {
        for(auto & qno : qubits)
        {
            DOUT("qno : " << qno);
            if( qno >= qubit_number )
            {
                EOUT("Number of qubits in platform: " << std::to_string(qubit_number) << ", specified qubit numbers out of range for gate: '" << gname << "' with " << ql::utils::to_string(qubits,"qubits") );
                throw ql::exception("[x] error : ql::kernel::gate() : Number of qubits in platform: "+std::to_string(qubit_number)+", specified qubit numbers out of range for gate '"+gname+"' with " +ql::utils::to_string(qubits,"qubits")+" !",false);
            }
        }

        // check if specialized composite gate is available
        // if not, check if parameterized composite gate is available
        // if not, check if a specialized custom gate is available
        // if not, check if a parameterized custom gate is available
        // if not, check if a default gate is available
        // if not, then error

        str::lower_case(gname);
        DOUT("Adding gate : " << gname << " with " << ql::utils::to_string(qubits,"qubits"));

        // specialized composite gate check
        DOUT("trying to add specialized decomposed gate for: " << gname);
        bool spec_decom_added = add_spec_decomposed_gate_if_available(gname, qubits);
        if(spec_decom_added)
        {
            DOUT("specialized decomposed gates added for " << gname);
        }
        else
        {
            // parameterized composite gate check
            DOUT("trying to add parameterized decomposed gate for: " << gname);
            bool param_decom_added = add_param_decomposed_gate_if_available(gname, qubits);
            if(param_decom_added)
            {
                DOUT("decomposed gates added for " << gname);
            }
            else
            {
                // specialized/parameterized custom gate check
                DOUT("adding custom gate for " << gname);
                bool custom_added = add_custom_gate_if_available(gname, qubits, duration, angle);
                if(!custom_added)
                {
                    if(ql::options::get("use_default_gates") == "yes")
                    {
                        // default gate check (which is always parameterized)
                        DOUT("adding default gate for " << gname);

                        bool default_available = add_default_gate_if_available(gname, qubits, duration);
                        if( default_available )
                        {
                            WOUT("default gate added for " << gname);
                        }
                        else
                        {
                            EOUT("unknown gate '" << gname << "' with " << ql::utils::to_string(qubits,"qubits") );
                            throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+gname+"' with " +ql::utils::to_string(qubits,"qubits")+" is not supported by the target platform !",false);
                        }
                    }
                    else
                    {
                        EOUT("unknown gate '" << gname << "' with " << ql::utils::to_string(qubits,"qubits") );
                        throw ql::exception("[x] error : ql::kernel::gate() : the gate '"+gname+"' with " +ql::utils::to_string(qubits,"qubits")+" is not supported by the target platform !",false);
                    }
                }
                else
                {
                    DOUT("custom gate added for " << gname);
                }
            }
        }
        DOUT("");
    }


    /**
     * qasm
     */
    std::string qasm()
    {
        std::stringstream ss;
        ss << "." << name;
        if (iterations > 1)
            ss << "(" << iterations << ") \n";
        else
            ss << "\n";
        for (size_t i=0; i<c.size(); ++i)
        {
            ss << "   " << c[i]->qasm() << "\n";
            // std::cout << c[i]->qasm() << std::endl;
        }
        return ss.str();
    }

    /**
     * micro code
     */
    std::string micro_code()
    {
        std::stringstream ss;
        // ss << "." << name;
        // if (iterations > 1)
        // ss << "(" << iterations << ")\n";
        // else
        // ss << "\n";
        for (size_t i=0; i<c.size(); ++i)
        {
            ss << c[i]->micro_code() << "\n";
            // std::cout << c[i]->qasm() << std::endl;
        }
        return ss.str();
    }


    void optimize()
    {
        ql::rotations_merging rm;
        if (contains_measurements(c))
        {
            // decompose the circuit
            std::vector<circuit*> cs = split_circuit(c);
            std::vector<circuit > cs_opt;
            for (size_t i=0; i<cs.size(); ++i)
            {
                if (!contains_measurements(*cs[i]))
                {
                    circuit opt = rm.optimize(*cs[i]);
                    cs_opt.push_back(opt);
                }
                else
                    cs_opt.push_back(*cs[i]);
            }
            // for (int i=0; i<cs_opt.size(); ++i)
            // print(cs_opt[i]);
            c.clear( );
            for (size_t i=0; i<cs_opt.size(); ++i)
                for (size_t j=0; j<cs_opt[i].size(); j++)
                    c.push_back(cs_opt[i][j]);
        }
        else
        {
            c = rm.optimize(c);
        }

    }

    void decompose_toffoli()
    {
        DOUT("decompose_toffoli()");
        for( auto cit = c.begin(); cit != c.end(); ++cit )
        {
            auto g = *cit;
            ql::gate_type_t gtype = g->type();
            std::vector<size_t> goperands = g->operands;

            ql::quantum_kernel toff_kernel("toff_kernel");
            toff_kernel.gate_definition = gate_definition;
            toff_kernel.qubit_number = qubit_number;
            toff_kernel.cycle_time = cycle_time;

            if( __toffoli_gate__ == gtype )
            {
                size_t cq1 = goperands[0];
                size_t cq2 = goperands[1];
                size_t tq = goperands[2];
                auto opt = ql::options::get("decompose_toffoli");
                if ( opt == "AM" )
                {
                    toff_kernel.controlled_cnot_AM(tq, cq1, cq2);
                }
                else
                {
                    toff_kernel.controlled_cnot_NC(tq, cq1, cq2);
                }
                ql::circuit& toff_ckt = toff_kernel.get_circuit();
                cit = c.erase(cit);
                cit = c.insert(cit, toff_ckt.begin(), toff_ckt.end());
            }
        }
        DOUT("decompose_toffoli() [Done] ");
    }

    void map(Mapper& mapper, std::string& map_in_qasm, std::string& map_out_qasm)
    {
        DOUT("Mapping kernel: " << name);

        map_in_qasm += qasm();
        mapper.MapCircuit(c);                // maps circuit in current context
        map_out_qasm += qasm();
    }

    void schedule(size_t qubits, quantum_platform platform, std::string& sched_qasm, std::string& sched_dot)
    {
        std::string scheduler = ql::options::get("scheduler");

#ifndef __disable_lemon__
        IOUT( scheduler << " scheduling the quantum kernel '" << name << "'...");
        // DOUT("Qasm at start of kernel::schedule size=" << c.size() << ":");
        // DOUT(qasm());
        // DOUT("Qasm at start of kernel::schedule END");

        Scheduler sched;
        sched.Init(qubits, c, platform);
        // sched.Print();
        // sched.PrintMatrix();
        // sched.PrintDot();

        if("ASAP" == scheduler)
        {
            // sched.PrintScheduleASAP();
            // sched.PrintDotScheduleASAP();
            // sched_dot = sched.GetDotScheduleASAP();
            // sched.PrintQASMScheduledASAP();
            ql::ir::bundles_t bundles = sched.schedule_asap();
            sched_qasm = ql::ir::qasm(bundles);

        }
        else if("ALAP" == scheduler)
        {
            // sched.PrintScheduleALAP();
            // sched.PrintDotScheduleALAP();
            // sched_dot = sched.GetDotScheduleALAP();
            // sched.PrintQASMScheduledALAP();
            ql::ir::bundles_t bundles = sched.schedule_alap();
            sched_qasm = ql::ir::qasm(bundles);
        }
        else
        {
            EOUT("Unknown scheduler");
        }
#endif // __disable_lemon__
    }

    std::vector<circuit*> split_circuit(circuit x)
    {
        IOUT("circuit decomposition in basic blocks ... ");
        std::vector<circuit*> cs;
        cs.push_back(new circuit());
        for (size_t i=0; i<x.size(); i++)
        {
            if ((x[i]->type() == __prepz_gate__) || (x[i]->type() == __measure_gate__))
            {
                cs.push_back(new circuit());
                cs.back()->push_back(x[i]);
                cs.push_back(new circuit());
            }
            else
            {
                cs.back()->push_back(x[i]);
            }
        }
        IOUT("circuit decomposion done (" << cs.size() << ").");
        /*
           for (int i=0; i<cs.size(); ++i)
           {
           println(" |-- circuit " << i);
           print(*(cs[i]));
           }
         */
        return cs;
    }

    /**
     * detect measurements and qubit preparations
     */
    bool contains_measurements(circuit x)
    {
        for (size_t i=0; i<x.size(); i++)
        {
            if (x[i]->type() == __measure_gate__)
                return true;
            if (x[i]->type() == __prepz_gate__)
                return true;
        }
        return false;
    }

    /**
     * detect unoptimizable gates
     */
    bool contains_unoptimizable_gates(circuit x)
    {
        for (size_t i=0; i<x.size(); i++)
        {
            if (x[i]->type() == __measure_gate__)
                return true;
            if (x[i]->type() == __prepz_gate__)
                return true;
            if (!(x[i]->optimization_enabled))
                return true;
        }
        return false;
    }

    /**
     * load custom instructions from a json file
     */
    int load_custom_instructions(std::string file_name="instructions.json")
    {
        load_instructions(gate_definition,file_name);
        return 0;
    }

    /**
     * debug
     */
    void print_gates_definition()
    {
        for (std::map<std::string,custom_gate*>::iterator i=gate_definition.begin(); i!=gate_definition.end(); i++)
        {
            COUT("[-] gate '" << i->first << "'");
            COUT(" |- qumis : \n" << i->second->micro_code());
        }
    }

    std::string get_gates_definition()
    {
        std::stringstream ss;

        for (std::map<std::string,custom_gate*>::iterator i=gate_definition.begin(); i!=gate_definition.end(); i++)
        {
            ss << i->first << '\n';
        }
        return ss.str();
    }

    /**
     * name getter
     */
    std::string get_name()
    {
        return name;
    }

    /**
     * circuit getter
     */
    circuit& get_circuit()
    {
        return c;
    }

    void controlled_x(size_t tq, size_t cq)
    {
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        cnot(cq, tq);
    }
    void controlled_y(size_t tq, size_t cq)
    {
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        sdag(tq);
        cnot(cq, tq);
        s(tq);
    }
    void controlled_z(size_t tq, size_t cq)
    {
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        hadamard(tq);
        cnot(cq, tq);
        hadamard(tq);
    }
    void controlled_h(size_t tq, size_t cq)
    {
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        s(tq);
        hadamard(tq);
        t(tq);
        cnot(cq, tq);
        tdag(tq);
        hadamard(tq);
        sdag(tq);
    }
    void controlled_i(size_t tq, size_t cq)
    {
        // well, basically you dont need to do anything for it :‑)
    }

    void controlled_s(size_t tq, size_t cq)
    {
        // cphase(cq, tq);

        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits

        cnot(tq, cq);
        tdag(cq);
        cnot(tq, cq);
        t(cq);
        t(tq);
    }

    void controlled_sdag(size_t tq, size_t cq)
    {
        // based on: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits

        tdag(cq);
        tdag(tq);
        cnot(tq, cq);
        t(cq);
        cnot(tq, cq);
    }

    void controlled_t(size_t tq, size_t cq)
    {
        WOUT("Controlled-T implementation requires an ancilla");
        WOUT("At the moment, Qubit 0 is used as ancilla");
        WOUT("This will change when Qubit allocater is implemented");
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        size_t aq = 0; // TODO at the moment qubit 0 is used as ancilla

        cnot(cq, tq);
        hadamard(aq);
        sdag(cq);
        cnot(tq, aq);
        cnot(aq, cq);
        t(cq);
        tdag(aq);
        cnot(tq, cq);
        cnot(tq, aq);
        t(cq);
        tdag(aq);
        cnot(aq, cq);
        h(cq);
        t(cq);
        h(cq);
        cnot(aq, cq);
        tdag(cq);
        t(aq);
        cnot(tq, aq);
        cnot(tq, cq);
        t(aq);
        tdag(cq);
        cnot(aq, cq);
        s(cq);
        cnot(tq, aq);
        cnot(cq, tq);
        h(aq);
    }

    void controlled_tdag(size_t tq, size_t cq)
    {
        WOUT("Controlled-Tdag implementation requires an ancilla");
        WOUT("At the moment, Qubit 0 is used as ancilla");
        WOUT("This will change when Qubit allocater is implemented");
        // from: https://arxiv.org/pdf/1206.0758v3.pdf
        // A meet-in-the-middle algorithm for fast synthesis
        // of depth-optimal quantum circuits
        size_t aq = 0; // TODO at the moment qubit 0 is used as ancilla

        h(aq);
        cnot(cq, tq);
        sdag(cq);
        cnot(tq, aq);
        cnot(aq, cq);
        t(cq);
        cnot(tq, cq);
        tdag(aq);
        cnot(tq, aq);
        t(cq);
        tdag(aq);
        cnot(aq, cq);
        h(cq);
        tdag(cq);
        h(cq);
        cnot(aq, cq);
        tdag(cq);
        t(aq);
        cnot(tq, aq);
        cnot(tq, cq);
        tdag(cq);
        t(aq);
        cnot(aq, cq);
        s(cq);
        cnot(tq, aq);
        cnot(cq, tq);
        hadamard(aq);
    }

    void controlled_ix(size_t tq, size_t cq)
    {
        // from: https://arxiv.org/pdf/1210.0974.pdf
        // Quantum circuits of T-depth one
        cnot(cq, tq);
        s(cq);
    }

    // toffoli decomposition
    // from: https://arxiv.org/pdf/1210.0974.pdf
    // Quantum circuits of T-depth one
    void controlled_cnot_AM(size_t tq, size_t cq1, size_t cq2)
    {
        h(tq);
        t(cq1);
        t(cq2);
        t(tq);
        cnot(cq2, cq1);
        cnot(tq, cq2);
        cnot(cq1, tq);
        tdag(cq2);
        cnot(cq1, cq2);
        tdag(cq1);
        tdag(cq2);
        tdag(tq);
        cnot(tq, cq2);
        cnot(cq1, tq);
        cnot(cq2, cq1);
        h(tq);
    }

    // toffoli decomposition
    // Neilsen and Chuang
    void controlled_cnot_NC(size_t tq, size_t cq1, size_t cq2)
    {
        h(tq);
        cnot(cq2,tq);
        tdag(tq);
        cnot(cq1,tq);
        t(tq);
        cnot(cq2,tq);
        tdag(tq);
        cnot(cq1,tq);
        tdag(cq2);
        t(tq);
        cnot(cq1,cq2);
        h(tq);
        tdag(cq2);
        cnot(cq1,cq2);
        t(cq1);
        s(cq2);
    }

    void controlled_swap(size_t tq1, size_t tq2, size_t cq)
    {
        // from: https://arxiv.org/pdf/1210.0974.pdf
        // Quantum circuits of T-depth one
        cnot(tq2, tq1);
        cnot(cq, tq1);
        h(tq2);
        t(cq);
        tdag(tq1);
        t(tq2);
        cnot(tq2, tq1);
        cnot(cq, tq2);
        t(tq1);
        cnot(cq, tq1);
        tdag(tq2);
        tdag(tq1);
        cnot(cq, tq2);
        cnot(tq2, tq1);
        t(tq1);
        h(tq2);
        cnot(tq2, tq1);
    }
    void controlled_rx(size_t tq, size_t cq, double theta)
    {
        rx(tq, theta/2);
        cz(cq, tq);
        rx(tq, -theta/2);
        cz(cq, tq);
    }
    void controlled_ry(size_t tq, size_t cq, double theta)
    {
        ry(tq, theta/2);
        cnot(cq, tq);
        ry(tq, -theta/2);
        cnot(cq, tq);
    }
    void controlled_rz(size_t tq, size_t cq, double theta)
    {
        rz(tq, theta/2);
        cnot(cq, tq);
        rz(tq, -theta/2);
        cnot(cq, tq);
    }


    void controlled(ql::quantum_kernel *k, std::vector<size_t> control_qubits)
    {
        DOUT("Generating controlled kernel ... ");

        ql::circuit& ckt = k->get_circuit();
        for( auto & g : ckt )
        {
            std::string gname = g->name;
            ql::gate_type_t gtype = g->type();
            std::vector<size_t> goperands = g->operands;
            DOUT("Generating controlled gate for " << gname);
            DOUT("Type : " << gtype);
            if( __pauli_x_gate__ == gtype  || __rx180_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_x(tq, cq);
            }
            else if( __pauli_y_gate__ == gtype  || __ry180_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_y(tq, cq);
            }
            else if( __pauli_z_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_z(tq, cq);
            }
            else if( __hadamard_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_h(tq, cq);
            }
            else if( __identity_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_i(tq, cq);
            }
            else if( __t_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_t(tq, cq);
            }
            else if( __tdag_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_tdag(tq, cq);
            }
            else if( __phase_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_s(tq, cq);
            }
            else if( __phasedag_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_sdag(tq, cq);
            }
            else if( __cnot_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq1 = control_qubits[0];
                size_t cq2 = control_qubits[1];
                auto opt = ql::options::get("decompose_toffoli");
                if ( opt == "AM" )
                {
                    controlled_cnot_AM(tq, cq1, cq2);
                }
                else
                {
                    controlled_cnot_NC(tq, cq1, cq2);
                }
            }
            else if( __swap_gate__ == gtype )
            {
                size_t tq1 = goperands[0];
                size_t tq2 = goperands[1];
                size_t cq = control_qubits[0];
                controlled_swap(tq1, tq2, cq);
            }
            else if( __rx_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rz(tq, cq, g->angle);
            }
            else if( __ry_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rz(tq, cq, g->angle);
            }
            else if( __rz_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rz(tq, cq, g->angle);
            }
            else if( __rx90_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rx(tq, cq, PI/2);
            }
            else if( __mrx90_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rx(tq, cq, -1*PI/2);
            }
            else if( __rx180_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_rx(tq, cq, PI);
                // controlled_x(tq, cq);
            }
            else if( __ry90_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_ry(tq, cq, PI/4);
            }
            else if( __mry90_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_ry(tq, cq, -1*PI/4);
            }
            else if( __ry180_gate__ == gtype )
            {
                size_t tq = goperands[0];
                size_t cq = control_qubits[0];
                controlled_ry(tq, cq, PI);
                // controlled_y(tq, cq);
            }
            else
            {
                EOUT("Controlled version of gate '" << gname << "' not defined !");
                throw ql::exception("[x] error : ql::kernel::controlled : Controlled version of gate '"+gname+"' not defined ! ",false);
            }
        }
    }

public:
    std::string name;
    circuit     c;
    size_t      iterations;
    size_t      qubit_number;
    size_t      cycle_time;
    std::map<std::string,custom_gate*> gate_definition;
};




} // namespace ql

#endif // QL_KERNEL_H
