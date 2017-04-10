#include "FiberManager.h"
#include <FEBioMech/FESPRProjection.h>
#include "FEAngioMaterial.h"

vec3d FiberManager::GetFiberDirection(GridPoint & pt, double& lambda)
{
	assert(pt.elemindex >= 0);
	assert(pt.ndomain != nullptr);

	// get the element
	FEElement * elem = &pt.ndomain->ElementRef(pt.elemindex);

	FESolidElement * selem = dynamic_cast<FESolidElement *>(elem);

	// Obtain shape function weights
	double shapeF[FEElement::MAX_NODES];
	assert(selem);
	selem->shape_fnc(shapeF, pt.q.x, pt.q.y, pt.q.z);//check these numbers

	// Determine component of new vessel direction due to nodal collagen fiber orientation
	vec3d coll_angle(0, 0, 0);
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		n = node_map[0][n];//map this value to be in bounds
		coll_angle += GetFiberAtNode(n)*shapeF[i];
	}

	// make unit vector
	lambda = coll_angle.unit();

	return coll_angle;
}
void FiberManager::Update()
{
	// get the domain
	FESolidDomain * sd = dynamic_cast<FESolidDomain*>(material->domainptrs[0]);
	int NN = sd->Nodes();
	int NE = sd->Elements();

	// build the element data array
	for (int n = 0; n < 3; n++)
	{
		fiber_at_int_pts[n].clear();
		fiber_at_int_pts[n].resize(NE);
		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& e = sd->Element(i);
			int nint = e.GaussPoints();
			fiber_at_int_pts[n][i].assign(nint, 0.0);
		}
	}
	

	// this array will store the results
	FESPRProjection map;

	// loop over stress components
	
	// fill the ED array
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd->Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			vec3d fd;
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			assert(angioPt->matPt->Next());
			FEElasticMaterialPoint*  emp = angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();
			for (int n = 0; n < 3; ++n)
			{	
				//fiber_at_int_pts[n][i][j] = emp.m_Q[n][0];
				switch (n)
				{
				case 0:
					fd.x = emp->m_Q[n][0];
					break;
				case 1:
					fd.y = emp->m_Q[n][0];
					break;
				case 2:
					fd.z = emp->m_Q[n][0];
					break;
				}
			}
			//multiply by the defgrad
			fd = emp->m_F * fd;

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
		}
	}
	for (int n = 0; n<3; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], fibers_at_nodes[n], node_map[n]);
	}

	//fill minor axis1
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd->Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			vec3d fd;
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			assert(angioPt->matPt->Next());
			FEElasticMaterialPoint& emp = *angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();
			for (int n = 0; n < 3; ++n)
			{
				//fiber_at_int_pts[n][i][j] = emp.m_Q[n][0];
				switch (n)
				{
				case 0:
					fd.x = emp.m_Q[n][0];
					break;
				case 1:
					fd.y = emp.m_Q[n][0];
					break;
				case 2:
					fd.z = emp.m_Q[n][0];
					break;
				}
			}
			//multiply by the defgrad
			fd = emp.m_F * fd;

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
		}
	}
	for (int n = 0; n<3; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], minoraxis1_at_nodes[n], node_map[n]);
	}

	//fill minor axis1
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd->Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			vec3d fd;
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			assert(angioPt->matPt->Next());
			FEElasticMaterialPoint& emp = *angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();
			for (int n = 0; n < 3; ++n)
			{
				//fiber_at_int_pts[n][i][j] = emp.m_Q[n][0];
				switch (n)
				{
				case 0:
					fd.x = emp.m_Q[n][0];
					break;
				case 1:
					fd.y = emp.m_Q[n][0];
					break;
				case 2:
					fd.z = emp.m_Q[n][0];
					break;
				}
			}
			//multiply by the defgrad
			fd = emp.m_F * fd;

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
		}
	}
	for (int n = 0; n<3; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], minoraxis2_at_nodes[n], node_map[n]);
	}

}

vec3d FiberManager::GetFiberAtNode(int node)
{
	vec3d rv;

	rv.x = fibers_at_nodes[0][node];
	rv.y = fibers_at_nodes[1][node];
	rv.z = fibers_at_nodes[2][node];

	rv.unit();
	return rv;
}

void RandomFiberInitializer::InitializeFibers(FiberManager * fman)
{
	//set the fiber vectors at the nodes of the fiber manager
	//make sure the map has been generated
	for (auto iter = fman->node_map[0].begin(); iter != fman->node_map[0].end(); ++iter)
	{
		int nn = iter->second;
		vec3d rv = fman->material->m_pangio->uniformRandomDirection();

		fman->fibers_at_nodes[0][nn] = rv.x;
		fman->fibers_at_nodes[1][nn] = rv.y;
		fman->fibers_at_nodes[2][nn] = rv.z;
	}

	//now reproject back to integration points
	nodeToInt(fman);
}

void FiberInitializer::nodeToInt(FiberManager * fman)
{
	//project the values from the nodes to integration points
	for (int i = 0; i < fman->material->domainptrs.size(); i++)
	{
		FEDomain * dom = fman->material->domainptrs[i];
		for (int j = 0; j < dom->Elements(); j++)
		{
			FESolidElement *se = dynamic_cast<FESolidElement*>(&dom->ElementRef(j));
			for (int k = 0; k < se->GaussPoints(); k++)
			{
				GridPoint pt;
				GridPointOfIntPoint(se, j, k, pt);
				double lambda;
				vec3d fiber = fman->GetFiberDirection(pt, lambda);
				
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
				FEElasticMaterialPoint *  emp = angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();
				
				emp->m_Q[0][0] = fiber.x;
				emp->m_Q[1][0] = fiber.y;
				emp->m_Q[2][0] = fiber.z;
				//TODO: still need to set minor axes
			}
			
		}
	}
}
void FiberInitializer::GridPointOfIntPoint(FESolidElement * se, int ei, int intp, GridPoint & gp)
{
	gp.ndomain = dynamic_cast<FESolidDomain*>(se->GetDomain());
	gp.q = vec3d(se->gr(intp), se->gs(intp), se->gt(intp));

	//not used in calculations
	//gp.nelem = se->GetID();
	gp.elemindex = ei;

	//the global position could also be set
}
