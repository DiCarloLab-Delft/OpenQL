/**
 * @file    datapath_cc.h
 * @date    20201119
 * @author  Wouter Vlothuizen (wouter.vlothuizen@tno.nl)
 * @brief   handling of Central Controller datapath (input MUX, Distributed Shared Memory, output PL)
 * @note
 */

#pragma once

#include "options_cc.h"
#include "gate.h"
#include "utils/vec.h"
//#include "codegen_cc.h"

#include <string>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <utils/logger.h>

namespace ql {

namespace { using namespace utils; }

// NB: types shared with codegen_cc
typedef struct {
	cond_type_t condition;
	Vec<UInt> cond_operands;
	uint32_t groupDigOut;
} tCondGateInfo;											// information for conditional gate on single instrument group
using tCondGateMap = std::map<int, tCondGateInfo>;			// NB: key is instrument group



class datapath_cc
{
public: // types

public: // functions
    datapath_cc() = default;
    ~datapath_cc() = default;

	void programStart();
	void programFinish();

	int allocateSmBit(size_t breg_operand, size_t instrIdx);
	int getSmBit(size_t breg_operand, size_t instrIdx);
	int getOrAssignMux(size_t instrIdx);
	int getOrAssignPl(size_t instrIdx);
	static int getSizeTag(int numReadouts);
	void emitPl(int pl, int smAddr, const tCondGateMap &condGateMap, size_t instrIdx, int slot);

	std::string getDatapathSection() { return datapathSection.str(); }

	void comment(const std::string &cmnt) {
	    datapathSection << cmnt << std::endl;
	}

// FIXME private:
	void emit(const std::string &sel, const std::string &statement, const std::string &comment="") {
		datapathSection << std::setw(16) << sel << std::setw(16) << statement << std::setw(24) << comment << std::endl;
	}
	void emit(int sel, const std::string &statement, const std::string &comment="") {
		emit(QL_SS2S("[" << sel << "]"), statement, comment);
	}

private:    // vars
	std::stringstream datapathSection;                          // the data path configuration generated
}; // class

} // namespace ql
