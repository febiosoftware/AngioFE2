#include "StdAfx.h"
#include "ECMInitializer.h"
#include <FECore/FEMesh.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngioMaterialBase.h"
#include "FEAngio.h"

void ECMInitializer::updateECMdensity(FEAngioMaterialBase* mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	//break even on core in field model
	mat->m_pangio->ForEachElementPar([&](FESolidElement & elem, FESolidDomain & d)
	{
		//these will hold the natural coordinates once the project to nodes is complete 
		double nr[FEElement::MAX_NODES];
		double ns[FEElement::MAX_NODES];
		double nt[FEElement::MAX_NODES];
		//these hold the natural coordinates of the integration points (r,s,t)
		double gr[FEElement::MAX_NODES];
		double gs[FEElement::MAX_NODES];
		double gt[FEElement::MAX_NODES];
		//TODO: if needed get FEBIO to expose the vectors that contain these to avoid this copy
		for (int i = 0; i < elem.Nodes(); i++)
		{
			gr[i] = elem.gr(i);
			gs[i] = elem.gs(i);
			gt[i] = elem.gt(i);
		}

		elem.project_to_nodes(gr, nr);
		elem.project_to_nodes(gs, ns);
		elem.project_to_nodes(gt, nt);

		// For each node in the element...
		for (int j = 0; j<elem.Nodes(); ++j)
		{
			// get the node
			int nnum = elem.m_node[j];
			nnum = mesh->Node(nnum).GetID();
			// get the ecm density and collagen fiber
			double ecm_den = mat->m_pangio->m_fe_node_data[nnum].m_ecm_den0;
			vec3d coll_fib = mat->m_pangio->m_fe_node_data[nnum].m_collfib0;

			/*
			//clamp n* to [1,-1]
			nr[j] = min(max(nr[j], -1), 1);
			ns[j] = min(max(ns[j], -1), 1);
			nt[j] = min(max(nt[j], -1), 1);
			*/

			//round to nearest integer
			nr[j] = round(nr[j]);
			ns[j] = round(ns[j]);
			nt[j] = round(nt[j]);

			// Calculate the deformation gradient tensor and jacobian at the node
			mat3d F;
			double Jacob = d.defgrad(elem, F, nr[j], ns[j], nt[j]);

			//make sure the function is differentiable and preserves orientation
			assert(Jacob > 0.0);

			// Update the collagen fiber orientation vector into the current configuration using F		
			coll_fib = F*coll_fib;
			coll_fib.unit();

			// Update matrix density using the Jacobian
			ecm_den = ecm_den / Jacob;

			// accumulate fiber directions and densities
#pragma omp critical
			{
				mat->m_pangio->m_fe_node_data[nnum].m_collfib += coll_fib;
				mat->m_pangio->m_fe_node_data[nnum].m_ecm_den += ecm_den;

				// increment counter
				mat->m_pangio->m_fe_node_data[nnum].m_ntag++;
			}
		}
	}, matls);
}

void ECMInitializerConstant::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	mat->m_pangio->ForEachNodePar([&](FENode & node)
	{
		mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 = mat->m_cultureParams.m_matrix_density;
		mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den = mat->m_cultureParams.m_matrix_density;
		mat->m_pangio->m_fe_node_data[node.GetID()].m_da = mat->GetAnisotropy();
	}, matls);
}

void ECMInitializerSpecified::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	mat->m_pangio->ForEachElement([this, mat, mesh](FESolidElement & se, FESolidDomain & sd)
	{
		double den[FEElement::MAX_INTPOINTS];
		double anis[FEElement::MAX_INTPOINTS];
		double pden[FEElement::MAX_INTPOINTS];
		double panis[FEElement::MAX_INTPOINTS];
		for (int n = 0; n<se.GaussPoints(); ++n)
		{
			// generate a coordinate transformation at this integration point
			FEMaterialPoint* mpoint = se.GetMaterialPoint(n);
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
			den[n] = angioPt->m_D;
			anis[n] = angioPt->m_DA;
		}
		se.project_to_nodes(den, pden);
		se.project_to_nodes(anis, panis);
		for (int k = 0; k < se.Nodes(); k++)
		{
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ecm_den0 += pden[k];
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ntag++;
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_da += panis[k];
		}
	}, matls);
}

void ECMInitializerNoOverwrite::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());
	mat->m_pangio->ForEachNodePar([&](FENode & node)
	{
		if (mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 == 0.0)
		{
			mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 = mat->m_cultureParams.m_matrix_density;
			mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den = mat->m_cultureParams.m_matrix_density;
			mat->m_pangio->m_fe_node_data[node.GetID()].m_da = mat->GetAnisotropy();
		}
	}, matls);
}