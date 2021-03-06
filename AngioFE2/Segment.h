#pragma once
#include <FECore/vec3d.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
// A helper class for locating points on the grid using an element number and
// a natural coordinates
class GridPoint
{
public:
	int		nelem;		// element id
	int		elemindex = -1;
	FESolidDomain *		ndomain=nullptr;    //domain id
	vec3d	q;			// natural coordinates
	vec3d	r;			// spatial position

	GridPoint() : nelem(-1) { }
	GridPoint(int ne, vec3d& k) : q(k), nelem(ne) { }
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
		BC_DEAD     = 1,		// segment is "dead" because it reached a boundary condition (?)
		INIT_BRANCH = 2,		// this is an initial sprout that is allowed to branch
		ANAST       = 4,		// this segment underwent anastomosis
		INIT_SPROUT = 8			// this is an intitial segment
	};

	// struct defining a tip of the segment
	class TIP
	{
	public:
		bool		bactive;	// flag if tip is active
		int			BC;			// something to do with body forces?
		GridPoint	pt;			// point in grid where this tip lies
		vec3d		u;			// sprout force vector
		Segment * parent;//segment that contains this tip
		Segment * connected = nullptr;//the segment that is connected to this tip
		int		nseed;		// seed where this tip started from (TODO: I want to remove this)
		int		nvessel;	// vessel this seed is part of

		double length_to_branch = 0.0;//determines how much farther the tip must grow before it branches
		double wait_time_to_branch = 0.0;//determines how long a segment waits before the branch grows

		const vec3d& pos() const { return pt.r; }

		TIP();
	};

    Segment();
	Segment(const Segment &obj);
	void operator = (Segment& seg);
	~Segment();
    
	// update the segment data
	// Call this each time the position of the nodes has changed
    void Update();

	// get the segment length
	double length() const { return m_length; }

	// get the unit vector
	const vec3d& uvect() const { return m_uvect; }

	// return one of the tip ends
	const TIP& tip_c(int i) const { return m_tip[i]; }

	TIP& tip(int i) { return m_tip[i]; }

	// add a flag
	void SetFlagOn(unsigned int nflag) { m_nflag |= nflag; }

	// turn a flag off
	void SetFlagOff(unsigned int nflag) { m_nflag &= ~nflag; }

	// get the status of a flag
	bool GetFlag(unsigned int nflag) const { return ((m_nflag & nflag) != 0); }

	// get all the flags
	unsigned int GetFlags() const { return m_nflag; }

	// set the time of birth
	void SetTimeOfBirth(double t) { m_TofBirth = t; }

	// get the time of birth
	double GetTimeOfBirth() const { return m_TofBirth; }

	// set the seed value
	void seed(int n)
	{
		m_nseed = n;
		m_tip[0].nseed = m_tip[1].nseed = n;
	}

	// get the seed value
	int seed() const { return m_nseed; }

	// set the vessel value
	void vessel(int n)
	{
		m_nvessel = n;
		m_tip[0].nvessel = m_tip[1].nvessel = n;
	}

	// get the vessel value
	int vessel() const { return m_nvessel; }

	int m_nid;			// segment id (unique zero-based ID)
	int m_nvessel;		// Label that indicates which vessel the segment belongs to
	double expected_length = 0.0;//the total growth this segment is supposed to take in the current timestep
private:
	TIP				m_tip[2];		// the two end tips
    double			m_length;       // Length of the segment.
    vec3d			m_uvect;		// unit direction vector
	unsigned int	m_nflag;		// status flag
	double			m_TofBirth;		// Time point at which the segment was created

	int m_nseed;		// Label that indicates which initial fragment the segment orginated from
	
};
