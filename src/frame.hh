/*
 * FILE: Frame.h
 * AUTHOR: Mark Diekhans
 * CREATE DATE: Frame of an arc
 * PROJECT: G-Known
 * DESCRIPTION: 28 Feb 2003
 * VERSION: $Revision$
 *
 * Copyright 1998-2001, The Regents of the University of California
 *
 * Departments of Computer Engineering and Computer Science
 * Jack Baskin School of Engineering
 * University of California, Santa Cruz, CA 95064
 */
#ifndef FRAME_H
#define FRAME_H
#include <typeinfo>
#include <string>
#include <assert.h>
using namespace std;

/*
 * Object representing a frame.Frame is maintained as 0, 1 or 2 throughout
 * each gene and is designated as <CODE>NO_FRAME</CODE> between genes.  
 * <br>
 * Note: FRAME_0, FRAME_1, and FRAME_2 are defined to have the respective numeric
 * values 0, 1, 2 and code may depend on this.
 */
class Frame {
    public:
    /* constants for value of frame */
    typedef enum {
        INVALID = -1,
        F0 = 0,
        F1 = 1,
        F2 = 2,
        NONE = 3
    } FrameEnum;

    /** Frame string names. */
    static const string INVALID_NAME;
    static const string F0_NAME;
    static const string F1_NAME;
    static const string F2_NAME;
    static const string NONE_NAME;

    private:
    /* value of this frame object */
    FrameEnum fVal;

    /* check for a legal value (which includes INVALID) */
    void assertLegal() const {
        assert((INVALID <= fVal) && (fVal <= NONE));
    }

    /* check for a valid value (which excluds INVALID) */
    void assertValid() const {
        assert((F0 <= fVal) && (fVal <= NONE));
    }

    public:
    /* Constructor, defaults to invalid value */
    Frame(FrameEnum val = INVALID):
        fVal(val) {
        assertLegal();
    }

    /* Constructor from int */
    Frame(int val):
        fVal(FrameEnum(val)) {
        assertLegal();
    }

    /* Copy constructor */
    Frame(const Frame& src):
        fVal(src.fVal) {
        assertLegal();
    }

    /* create a frame from a GFF-style phase */
    static Frame fromPhase(int phase);

    /* create a frame from a GFF-style phase string */
    static Frame fromPhaseStr(const string& phase);

    /* get a GFF-static phase for a frame */
    int toPhase() const;

    /* get a GFF-static phase tring for a frame */
    string toPhaseStr() const;

    /* conversion to int */
    operator int() const {
        return fVal;
    }

    /** Assignment */
    Frame& operator=(const Frame& frame) {
        fVal = frame.fVal;
        assertLegal();
        return *this;
    }

    /** Comparison */
    bool operator==(const Frame& frame) const {
        return (fVal == frame.fVal);
    }
    bool operator!=(const Frame& frame) const {
        return (fVal != frame.fVal);
    }
    bool operator==(FrameEnum frame) const {
        return (fVal == frame);
    }
    bool operator!=(FrameEnum frame) const {
        return (fVal != frame);
    }

    /** Determine if this is a real frame (0, 1, or 2) */
    bool isRealFrame() const {
        return (fVal == F0) ||  (fVal == F1) ||  (fVal == F2);
    }

    /** Determine distance between from this frame to another. */
    int diff(Frame otherFrame) const;

    /**
     * Increment the frame value, returning a new Frame.  That is F0->F1,
     * F1->F2, F2->F0.  NONE -> NONE.
     */
    Frame incr(int amount = 1) const {
        assertValid();
        switch (fVal) {
        case F0:
        case F1:
        case F2:
            if (amount >= 0) {
                return Frame(((fVal+amount) % 3));
            } else {
                int amount3 = (-amount) % 3;
                return Frame(((fVal - (amount - amount3)) % 3));
            }
        case INVALID:
        case NONE:
        default:
            return Frame(fVal);;
        }
    }

    /** 
     * Get the next frame, but don't cycle.  After frame 2, NONE is returned.
     * NONE always returns NONE. 
     */
    Frame& operator++(int) {
        assertValid();
        if ((F0 <= fVal) && (fVal <= F2)) {
            fVal = FrameEnum(((int)fVal)+1);
        }
        return *this;
    }

    /** Get a single character abbreviation for the frame. */
    string getAbbrev() const;

    /** Convert a frame to a symbolic name. */
    const string& toString() const;
};

#endif

/*
 * Local Variables:
 * mode: c++
 * c-basic-off: 4
 * End:
 */

