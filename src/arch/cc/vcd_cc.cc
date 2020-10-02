/**
 * @file    vcd_cc.cc
 * @date    20201001
 * @author  Wouter Vlothuizen (wouter.vlothuizen@tno.nl)
 * @brief   handle VCD file generation for the CC backend
 * @note
 */

#include "vcd_cc.h"
#include "options.h"

// NB: parameters qubit_number and cycle_time originate from OpenQL variable 'platform'
void vcd_cc::programStart(int qubit_number, int cycle_time, int maxGroups, const settings_cc &settings)
{
    this->cycle_time = cycle_time;
    kernelStartTime = 0;

    // define header
    vcd.start();

    // define kernel variable
    vcd.scope(vcd.ST_MODULE, "kernel");
    vcdVarKernel = vcd.registerVar("kernel", Vcd::VT_STRING);
    vcd.upscope();

    // define qubit variables
    vcd.scope(vcd.ST_MODULE, "qubits");
    vcdVarQubit.resize(qubit_number);
    for(size_t q=0; q<qubit_number; q++) {
        std::string name = "q"+std::to_string(q);
        vcdVarQubit[q] = vcd.registerVar(name, Vcd::VT_STRING);
    }
    vcd.upscope();

    // define signal variables
    size_t instrsUsed = settings.getInstrumentsSize();
    vcd.scope(vcd.ST_MODULE, "sd.signal");
    vcdVarSignal.assign(instrsUsed, std::vector<int>(maxGroups, {0}));
    for(size_t instrIdx=0; instrIdx<instrsUsed; instrIdx++) {
        const json &instrument = settings.getInstrumentAtIdx(instrIdx);                  // NB: always exists
        std::string instrumentPath = SS2S("instruments["<<instrIdx<<"]");       // for JSON error reporting
        std::string instrumentName = json_get<std::string>(instrument, "name", instrumentPath);
        const json qubits = json_get<const json>(instrument, "qubits", instrumentPath);
        for(size_t group=0; group<qubits.size(); group++) {
            std::string name = instrumentName+"-"+std::to_string(group);
            vcdVarSignal[instrIdx][group] = vcd.registerVar(name, Vcd::VT_STRING);
        }
    }
    vcd.upscope();

    // define codeword variables
    vcd.scope(vcd.ST_MODULE, "codewords");
    vcdVarCodeword.resize(qubit_number);
    for(size_t instrIdx=0; instrIdx<instrsUsed; instrIdx++) {
        const json &instrument = settings.getInstrumentAtIdx(instrIdx);         // NB: always exists
        std::string instrumentPath = SS2S("instruments["<<instrIdx<<"]");       // for JSON error reporting
        std::string instrumentName = json_get<std::string>(instrument, "name", instrumentPath);
        vcdVarCodeword[instrIdx] = vcd.registerVar(instrumentName, Vcd::VT_STRING);
    }
    vcd.upscope();
}


void vcd_cc::programFinish(const std::string &progName)
{
    // generate VCD
    vcd.finish();

    // write VCD to file
    std::string file_name(ql::options::get("output_dir") + "/" + progName + ".vcd");
    IOUT("Writing Value Change Dump to " << file_name);
    ql::utils::write_file(file_name, vcd.getVcd());
}


void vcd_cc::kernelFinish(const std::string &kernelName, size_t durationInCycles)
{
    // NB: timing starts anew for every kernel
    unsigned int durationInNs = durationInCycles*cycle_time;
    vcd.change(vcdVarKernel, kernelStartTime, kernelName);          // start of kernel
    vcd.change(vcdVarKernel, kernelStartTime + durationInNs, "");   // end of kernel
    kernelStartTime += durationInNs;
}


void vcd_cc::bundleFinishGroup(size_t startCycle, unsigned int durationInNs, uint32_t groupDigOut, const std::string &signalValue, int instrIdx, int group)
{
    // generate signal output for group
    unsigned int startTime = kernelStartTime + startCycle*cycle_time;
    int var = vcdVarSignal[instrIdx][group];
    std::string val = SS2S(groupDigOut) + "=" + signalValue;
    vcd.change(var, startTime, val);                    // start of signal
    vcd.change(var, startTime+durationInNs, "");        // end of signal
}


void vcd_cc::bundleFinish(size_t startCycle, uint32_t digOut, size_t maxDurationInCycles, int instrIdx)
{
    // generate codeword output for instrument
    unsigned int startTime = kernelStartTime + startCycle*cycle_time;
    unsigned int durationInNs = maxDurationInCycles*cycle_time;
    int var = vcdVarCodeword[instrIdx];
    std::string val = SS2S("0x" << std::hex << std::setfill('0') << std::setw(8) << digOut);
    vcd.change(var, startTime, val);                // start of signal
    vcd.change(var, startTime+durationInNs, "");    // end of signal
}


void vcd_cc::customgate(const std::string &iname, const std::vector<size_t> &qops, size_t startCycle, size_t durationInNs)
{
    // generate qubit VCD output
    unsigned int startTime = kernelStartTime + startCycle*cycle_time;
    for(size_t i=0; i<qops.size(); i++) {
        int var = vcdVarQubit[qops[i]];
        std::string name = iname;                       // FIXME: improve name for 2q gates
        vcd.change(var, startTime, name);               // start of instruction
        vcd.change(var, startTime+durationInNs, "");    // end of instruction
    }
}