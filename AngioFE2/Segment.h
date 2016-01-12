#pragma once
#include <FECore/vec3d.h>

//-----------------------------------------------------------------------------
// A helper class for locating points on the grid using an element number and
// a natural coordinates
class GridPoint
{
public:
	int		nelem;		// element number
	vec3d	q;			// natural coordinates

public:
	GridPoint() { nelem = -1; }
	GridPoint(int ne, vec3d& r) { nelem = ne; q = r; }
};

//-----------------------------------------------------------------------------
// Microvessels are represented by a collection of line segments. 
// Growth is represented by the addition of new segments onto the 
// active tips of exisiting segments. Within the simulation, these 
// line segments are found in the Segment class.
class Segment  
{
public:
	// status flags
	enum {
		BC_DEAD     = 1,		// segment is "dead" because it reached a boundar condition (?)
		INIT_BRANCH = 2,		// this is an initial sprout that is allowed to branch
		ANAST       = 4,		// this segment underwent anastimoses
		INIT_SPROUT = 8			// this is an intitial segment
	};

	// struct defining a tip of the segment
	class TIP
	{
	public:
		vec3d		rt;			// current position of segment tip
		bool		bactive;	// flag if tip is active
		int			bdyf_id;	// ID of the body force
		int			BC;			// something to do with body forces?
		GridPoint	pt;			// point in grid where this tip lies

	public:
		TIP();
	};

public:
    Segment();
	virtual ~Segment();
    
	// update the segment data
	// Call this each time the position of the nodes has changed
    void Update();

	// get the segment length
	double length() const { return m_length; }

	// get the unit vector
	vec3d uvect() { return m_uvect; }

	// return one of the tip ends
	TIP& tip(int i) { return m_tip[i]; }

	// return one of the tip ends
	const TIP& tip(int i) const { return m_tip[i]; }

	// add a flag
	void SetFlagOn(unsigned int nflag) { m_nflag |= nflag; }

	// turn a flag off
	void SetFlagOff(unsigned int nflag) { m_nflag &= ~nflag; }

	// get the status of a flag
	bool GetFlag(unsigned int nflag) { return ((m_nflag & nflag) != 0); }

	// get all the flags
	unsigned int GetFlags() { return m_nflag; }

	// set the time of birth
	void SetTimeOfBirth(double t) { m_TofBirth = t; }

	// get the time of birth
	double GetTimeOfBirth() { return m_TofBirth; }

public:
	int m_nseed;		// Label that indicates which initial fragment the segment orginated from
	int m_nvessel;		// Label that indicates which vessel the segment belongs to
	int m_nid;			// segment id (unique zero-based ID)

	// TODO: I don't think there is any reason why segments should be killed and this
	//       is probably a way to clean up after bugs. Better solution: fix bugs and remove this.
	bool mark_of_death;			// this flag marks the segment as dead and should be removed.
	int death_label;			// cause of death

private:
	TIP				m_tip[2];		// the two end tips
    double			m_length;       // Length of the segment.
    vec3d			m_uvect;		// unit direction vector
	unsigned int	m_nflag;		// status flag
	double			m_TofBirth;		// Time point at which the segment was created
};
