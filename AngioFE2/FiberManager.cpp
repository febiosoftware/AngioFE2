#include "FiberManager.h"
#include <FEBioMech/FESPRProjection.h>
#include <FEBioMech/FEFiberMaterialPoint.h>
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
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetFiberAtNode(n)*shapeF[i];
	}

	// make unit vector
	lambda = coll_angle.unit();

	return coll_angle;
}

vec3d FiberManager::GetMinorAxisDirection1(GridPoint & pt)
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
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetMinor1AtNode(n)*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}
vec3d FiberManager::GetMinorAxisDirection2(GridPoint & pt)
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
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetMinor2AtNode(n)*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

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
			FEMaterialPoint * mp = el.GetMaterialPoint(j);
			vec3d fd = material->GetFiberVector(*mp);

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
		}
	}
	for (int n = 0; n<3; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], fibers_at_nodes[n]);
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
vec3d FiberManager::GetMinor1AtNode(int node)
{
	vec3d rv;

	rv.x = minoraxis1_at_nodes[0][node];
	rv.y = minoraxis1_at_nodes[1][node];
	rv.z = minoraxis1_at_nodes[2][node];

	rv.unit();
	return rv;
}
vec3d FiberManager::GetMinor2AtNode(int node)
{
	vec3d rv;

	rv.x = minoraxis2_at_nodes[0][node];
	rv.y = minoraxis2_at_nodes[1][node];
	rv.z = minoraxis2_at_nodes[2][node];

	rv.unit();
	return rv;
}

void RandomFiberInitializer::InitializeFibers(FiberManager * fman)
{
	//set the fiber vectors at the nodes of the fiber manager
	//make sure the map has been generated
	for (auto iter = fman->material->node_map.begin(); iter != fman->material->node_map.end(); ++iter)
	{
		size_t nn = iter->second;

		mat3d rm = fman->material->m_pangio->unifromRandomRotationMatrix();

		vec3d q(
			fman->fibers_at_nodes[0][nn] ,
			fman->fibers_at_nodes[1][nn] ,
			fman->fibers_at_nodes[2][nn] 
		);
		vec3d m1(
			fman->minoraxis1_at_nodes[0][nn],
			fman->minoraxis1_at_nodes[1][nn],
			fman->minoraxis1_at_nodes[2][nn]
		);
		vec3d m2(
			fman->minoraxis2_at_nodes[0][nn],
			fman->minoraxis2_at_nodes[1][nn],
			fman->minoraxis2_at_nodes[2][nn]
		);
		vec3d nq = rm*q;
		vec3d nm1 = rm*m1;
		vec3d nm2 = rm*m2;

#ifndef NDEBUG
		double noq = nq.norm();
		double nom1 = nm1.norm();
		double nom2 = nm2.norm();
		double tol = 1.0;
		assert(noq < (1 + tol) && noq >(1 - tol));
		assert(nom1 < (1 + tol) && nom1 >(1 - tol));
		assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif

		fman->fibers_at_nodes[0][nn] = nq.x;
		fman->fibers_at_nodes[1][nn] = nq.y;
		fman->fibers_at_nodes[2][nn] = nq.z;

		fman->minoraxis1_at_nodes[0][nn] = nm1.x;
		fman->minoraxis1_at_nodes[1][nn] = nm1.y;
		fman->minoraxis1_at_nodes[2][nn] = nm1.z;

		fman->minoraxis2_at_nodes[0][nn] = nm2.x;
		fman->minoraxis2_at_nodes[1][nn] = nm2.y;
		fman->minoraxis2_at_nodes[2][nn] = nm2.z;
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
				vec3d m1 = fman->GetMinorAxisDirection1(pt);
				vec3d m2 = fman->GetMinorAxisDirection2(pt);
				
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
				FEElasticMaterialPoint *  emp = angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();
				//fiber
				emp->m_Q[0][0] = fiber.x;
				emp->m_Q[1][0] = fiber.y;
				emp->m_Q[2][0] = fiber.z;
				
				//minor1
				emp->m_Q[0][1] = m1.x;
				emp->m_Q[1][1] = m1.y;
				emp->m_Q[2][1] = m1.z;

				//minor2
				emp->m_Q[0][2] = m2.x;
				emp->m_Q[1][2] = m2.y;
				emp->m_Q[2][2] = m2.z;
			}
			
		}
	}
}

void RandomFiberInitializerNonMangling::InitializeFibers(FiberManager * fman)
{
	//everything stays at the integration points
	for (int i = 0; i < fman->material->domainptrs.size(); i++)
	{
		FEDomain * dom = fman->material->domainptrs[i];
		for (int j = 0; j < dom->Elements(); j++)
		{
			FESolidElement *se = dynamic_cast<FESolidElement*>(&dom->ElementRef(j));
			for (int k = 0; k < se->GaussPoints(); k++)
			{

				

				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEFiberMaterialPoint *  fp = mp->ExtractData<FEFiberMaterialPoint>();
				//FEElasticMaterialPoint *  emp = angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();


				mat3d rm = fman->material->m_pangio->unifromRandomRotationMatrix();
				fp->m_n0 = rm * fp->m_n0;
#ifndef NDEBUG
				double noq = fiber.norm();
				double nom1 = m1.norm();
				double nom2 = m2.norm();
				double tol = 1.0;
				assert(noq < (1 + tol) && noq >(1 - tol));
				assert(nom1 < (1 + tol) && nom1 >(1 - tol));
				assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif
			}

		}
	}
}

void RandomFiberInitializerPE::InitializeFibers(FiberManager * fman)
{
	//everything stays at the integration points
	for (int i = 0; i < fman->material->domainptrs.size(); i++)
	{
		FEDomain * dom = fman->material->domainptrs[i];
		for (int j = 0; j < dom->Elements(); j++)
		{
			mat3d rm = fman->material->m_pangio->unifromRandomRotationMatrix();
			FESolidElement *se = dynamic_cast<FESolidElement*>(&dom->ElementRef(j));
			for (int k = 0; k < se->GaussPoints(); k++)
			{



				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
				FEElasticMaterialPoint *  emp = angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();

				vec3d fiber = vec3d(
					emp->m_Q[0][0],
					emp->m_Q[1][0],
					emp->m_Q[2][0]
				);
				vec3d m1 = vec3d(
					emp->m_Q[0][1],
					emp->m_Q[1][1],
					emp->m_Q[2][1]
				);
				vec3d m2 = vec3d(
					emp->m_Q[0][2],
					emp->m_Q[1][2],
					emp->m_Q[2][2]
				);

				
				fiber = rm*fiber;
				m1 = rm*m1;
				m2 = rm*m2;

				//fiber
				emp->m_Q[0][0] = fiber.x;
				emp->m_Q[1][0] = fiber.y;
				emp->m_Q[2][0] = fiber.z;

				//minor1
				emp->m_Q[0][1] = m1.x;
				emp->m_Q[1][1] = m1.y;
				emp->m_Q[2][1] = m1.z;

				//minor2
				emp->m_Q[0][2] = m2.x;
				emp->m_Q[1][2] = m2.y;
				emp->m_Q[2][2] = m2.z;
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
