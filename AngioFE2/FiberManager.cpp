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
	lambda = 0.0;
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		double lc;//current lambda
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetFiberAtNode(n, lc)*shapeF[i];
		lambda += lc * shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}

vec3d FiberManager::GetMinorAxisDirection1(GridPoint & pt, double  &lambda)
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
	lambda = 0;
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		double lc;
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetMinor1AtNode(n, lc)*shapeF[i];
		lambda += lc*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}
vec3d FiberManager::GetMinorAxisDirection2(GridPoint & pt, double & lambda)
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
	lambda = 0;
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		double lc;
		n = material->node_map[n];//map this value to be in bounds
		coll_angle += GetMinor2AtNode(n, lc)*shapeF[i];
		lambda += lc * shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}
void FiberManager::Update()
{
	// get the domain
	FESolidDomain * sd = dynamic_cast<FESolidDomain*>(material->domainptrs[0]);
	assert(sd);
	int NN = sd->Nodes();
	int NE = sd->Elements();

	// build the element data array
	for (int n = 0; n < 4; n++)
	{
		fiber_at_int_pts[n].clear();
		fiber_at_int_pts[n].resize(NE);
		m1_at_int_pts[n].clear();
		m1_at_int_pts[n].resize(NE);
		m2_at_int_pts[n].clear();
		m2_at_int_pts[n].resize(NE);

		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& e = sd->Element(i);
			int nint = e.GaussPoints();
			fiber_at_int_pts[n][i].assign(nint, 0.0);
			m1_at_int_pts[n][i].assign(nint, 0.0);
			m2_at_int_pts[n][i].assign(nint, 0.0);
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
			FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
			vec3d fd = emp->m_F * emp->m_Q * vec3d(1,0,0);
			double lambda = fd.unit();

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
			fiber_at_int_pts[3][i][j] = lambda;

			//store the data
			fd = emp->m_F * emp->m_Q * vec3d(0, 1, 0);
			lambda = fd.unit();
			m1_at_int_pts[0][i][j] = fd.x;
			m1_at_int_pts[1][i][j] = fd.y;
			m1_at_int_pts[2][i][j] = fd.z;
			m1_at_int_pts[3][i][j] = lambda;

			fd = emp->m_F * emp->m_Q * vec3d(0, 0, 1);
			lambda = fd.unit();
			m2_at_int_pts[0][i][j] = fd.x;
			m2_at_int_pts[1][i][j] = fd.y;
			m2_at_int_pts[2][i][j] = fd.z;
			m2_at_int_pts[3][i][j] = lambda;
		}
	}
	
	for (int n = 0; n<4; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], fibers_at_nodes[n]);
		map.Project(*sd, m1_at_int_pts[n], minoraxis1_at_nodes[n]);
		map.Project(*sd, m2_at_int_pts[n], minoraxis2_at_nodes[n]);
	}
}

vec3d FiberManager::GetFiberAtNode(int node, double & lambda)
{
	vec3d rv;

	rv.x = fibers_at_nodes[0][node];
	rv.y = fibers_at_nodes[1][node];
	rv.z = fibers_at_nodes[2][node];
	lambda = fibers_at_nodes[3][node];
	rv.unit();
	return rv;
}
vec3d FiberManager::GetMinor1AtNode(int node, double  &lambda)
{
	vec3d rv;

	rv.x = minoraxis1_at_nodes[0][node];
	rv.y = minoraxis1_at_nodes[1][node];
	rv.z = minoraxis1_at_nodes[2][node];
	lambda = minoraxis1_at_nodes[3][node];

	rv.unit();
	return rv;
}
vec3d FiberManager::GetMinor2AtNode(int node, double  &lambda)
{
	vec3d rv;

	rv.x = minoraxis2_at_nodes[0][node];
	rv.y = minoraxis2_at_nodes[1][node];
	rv.z = minoraxis2_at_nodes[2][node];
	lambda = minoraxis2_at_nodes[3][node];

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
			assert(se);
			for (int k = 0; k < se->GaussPoints(); k++)
			{
				GridPoint pt;
				GridPointOfIntPoint(se, j, k, pt);
				double lambda;
				vec3d fiber = fman->GetFiberDirection(pt, lambda);
				vec3d m1 = fman->GetMinorAxisDirection1(pt, lambda);
				vec3d m2 = fman->GetMinorAxisDirection2(pt, lambda);
				
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
				FEElasticMaterialPoint *  emp = angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
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
			assert(se);
			for (int k = 0; k < se->GaussPoints(); k++)
			{
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint * angiopt = mp->ExtractData<FEAngioMaterialPoint>();
				FEElasticMaterialPoint *  emp = mp->ExtractData<FEElasticMaterialPoint>();

				FEElasticMaterialPoint *  emp_matrix = angiopt->matPt->ExtractData<FEElasticMaterialPoint>();
				FEElasticMaterialPoint *  emp_vessel = angiopt->vessPt->ExtractData<FEElasticMaterialPoint>();

				mat3d rm = fman->material->m_pangio->unifromRandomRotationMatrix();
				emp->m_Q = emp->m_Q * rm;
				emp_matrix->m_Q = emp->m_Q;
				emp_vessel->m_Q = emp->m_Q;
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
			assert(se);
			for (int k = 0; k < se->GaussPoints(); k++)
			{
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint * angiopt = mp->ExtractData<FEAngioMaterialPoint>();
				FEElasticMaterialPoint *  emp = mp->ExtractData<FEElasticMaterialPoint>();

				FEElasticMaterialPoint *  emp_matrix = angiopt->matPt->ExtractData<FEElasticMaterialPoint>();
				FEElasticMaterialPoint *  emp_vessel = angiopt->vessPt->ExtractData<FEElasticMaterialPoint>();

				emp->m_Q = emp->m_Q * rm;
				emp_matrix->m_Q = emp->m_Q;
				emp_vessel->m_Q = emp->m_Q;
			}
		}
	}
}

void FiberInitializer::GridPointOfIntPoint(FESolidElement * se, int ei, int intp, GridPoint & gp)
{
	gp.ndomain = dynamic_cast<FESolidDomain*>(se->GetDomain());
	gp.q = vec3d(se->gr(intp), se->gs(intp), se->gt(intp));

	//not used in calculations
	//gp.nelem = se->GetID_ang();
	gp.elemindex = ei;

	//the global position could also be set
}
ExplicitDistributionsFiberInitializer::ExplicitDistributionsFiberInitializer(FEModel * model) : FiberInitializer(model)
{
	AddProperty(&alpha, "alpha");
	AddProperty(&alpha, "beta");
	AddProperty(&alpha, "gamma");
}
void ExplicitDistributionsFiberInitializer::Setup()
{
	alpha->StepToTime(0);
	beta->StepToTime(0);
	gamma->StepToTime(0);
}

void ExplicitDistributionsFiberInitializer::InitializeFibers(FiberManager * fman)
{
	//everything stays at the integration points
	for (int i = 0; i < fman->material->domainptrs.size(); i++)
	{
		FEDomain * dom = fman->material->domainptrs[i];
		for (int j = 0; j < dom->Elements(); j++)
		{
			FESolidElement *se = dynamic_cast<FESolidElement*>(&dom->ElementRef(j));
			assert(se);
			for (int k = 0; k < se->GaussPoints(); k++)
			{
				FEMaterialPoint * mp = se->GetMaterialPoint(k);
				FEAngioMaterialPoint * angiopt = mp->ExtractData<FEAngioMaterialPoint>();
				FEElasticMaterialPoint *  emp = mp->ExtractData<FEElasticMaterialPoint>();

				FEElasticMaterialPoint *  emp_matrix = angiopt->matPt->ExtractData<FEElasticMaterialPoint>();
				FEElasticMaterialPoint *  emp_vessel = angiopt->vessPt->ExtractData<FEElasticMaterialPoint>();

				mat3d rm = fman->material->m_pangio->rotationMatrix(alpha->NextValue(fman->material->m_pangio->rengine), beta->NextValue(fman->material->m_pangio->rengine), gamma->NextValue(fman->material->m_pangio->rengine));
				emp->m_Q = emp->m_Q * rm;
				emp_matrix->m_Q = emp->m_Q;
				emp_vessel->m_Q = emp->m_Q;
			}
		}
	}
}
