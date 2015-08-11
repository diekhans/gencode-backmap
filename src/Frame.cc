#include "Frame.hh"
#include <stdexcept>
#include "typeOps.hh"


/*
 * Frame string names.
 */
const string Frame::INVALID_NAME = "invalid";
const string Frame::F0_NAME = "frame0";
const string Frame::F1_NAME = "frame1";
const string Frame::F2_NAME = "frame2";
const string Frame::NONE_NAME = "none";

/* create a frame from a GFF-style phase */
Frame Frame::fromPhase(int phase) {
    switch (phase) {
    case 0:
        return Frame(F0);
    case 1:
        return Frame(F2);
    case 2:
        return Frame(F1);
    default:
        throw logic_error("Invalid phase: " + ::toString(phase));
    }
}

/* create a frame from a GFF-style phase string */
Frame Frame::fromPhaseStr(const string& phase) {
    if (phase == ".") {
        return Frame(NONE);
    } else if (phase == "0") {
        return Frame(F0);
    } else if (phase == "1") {
        return Frame(F2);
    } else if (phase == "2") {
        return Frame(F1);
    } else {
        throw logic_error("Invalid phase: " + phase);
    }
}

/*
 * Convert a frame number to a GFF file frame number.
 */
int Frame::toPhase() const {
    switch (fVal) {
    case 0:
        return 0;
    case 1:
        return 2;
    case 2:
        return 1;
    default:
        return -1;
    }
}

/* get a GFF-static phase tring for a frame */
string Frame::toPhaseStr() const {
    switch (fVal) {
    case 0:
        return "0";
    case 1:
        return "2";
    case 2:
        return "1";
    default:
        return ".";
    }
}


/** 
 * Get a single character abbreviation for the frame
 * Returned as a string  rather than a char to avoid problems using it
 * in string expressions.
 */
string Frame::getAbbrev() const {
    switch (fVal) {
    case F0:
        return "0";
    case F1:
        return "1";
    case F2:
        return "2";
    case NONE:
        return "N";
    case INVALID:
    default:
        return "I";
        }
}

/*
 * Convert a frame to a symbolic name.
 */
const string& Frame::toString() const {
    assertLegal();
    switch (fVal) {
      case INVALID:
      default:
        return INVALID_NAME;
      case F0:
        return F0_NAME;
      case F1:
        return F1_NAME;
      case F2:
        return F2_NAME;
      case NONE:
        return NONE_NAME;
    }
}

/*
 * Determine distance between from this frame to another.
 */
int Frame::diff(Frame otherFrame) const {
    assertValid();
    if ((*this == INVALID) || (*this == NONE)) {
        throw logic_error("Feature::frameDifference: invalid *this");
    }
    if ((otherFrame == INVALID) || (otherFrame == NONE)) {
        throw logic_error("Feature::frameDifference: invalid otherFrame");
    }

    Frame frame = *this;
    int diff = 0;
    while (frame != otherFrame) {
        frame = frame.incr();
        diff++;
    }
    return diff;
}



/*
 * Local Variables:
 * mode: c++
 * c-basic-off: 4
 * End:
 */

