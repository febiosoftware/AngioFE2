#include "GrowDirectionModifier.h"
#include "FEAngio.h"
#include "Culture.h"
#include "angio3d.h"

GrowDirectionModifier::GrowDirectionModifier(FEModel * model) : FEMaterial(model)
{

}

void GrowDirectionModifier::SetCulture(Culture * cp)
{
	culture = cp;
}

GrowDirectionModifiers::GrowDirectionModifiers(FEModel* model) : FEMaterial(model)
{
	AddProperty(&grow_direction_modifiers, "gdm");
}


vec3d GrowDirectionModifiers::ApplyModifiers(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		previous_dir = grow_direction_modifiers[i]->GrowModifyGrowDirection(previous_dir, tip, mat, branch, start_time, grow_time, seg_length);
	}
	return previous_dir;
}

void GrowDirectionModifiers::SetCulture(Culture * c)
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		grow_direction_modifiers[i]->SetCulture(c);
	}
}

GradientGrowDirectionModifier::GradientGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
//begin implementations of grow direction modifiers
vec3d GradientGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	//calculate the density gradinet if above the threshold set the grow direction
	std::vector<double> densities;
	FESolidElement * se = dynamic_cast<FESolidElement*>(&tip.pt.ndomain->ElementRef(tip.pt.elemindex));
	densities = culture->m_pmat->m_pangio->createVectorOfMaterialParameters(se, &FEAngioNodeData::m_ecm_den);
	vec3d gradient = culture->m_pmat->m_pangio->gradient(se, densities, tip.pt.q);
	double gradnorm = gradient.norm();
	Segment seg;
	if (gradnorm > threshold)
	{
		vec3d currentDirection = previous_dir;
		currentDirection.unit();
		vec3d currentDirectionGradientPlane = gradient ^ currentDirection;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ gradient;
		perpendicularToGradient.unit();
		return perpendicularToGradient;
	}
	return previous_dir;
}
BEGIN_PARAMETER_LIST(GradientGrowDirectionModifier, GrowDirectionModifier)
ADD_PARAMETER(threshold, FE_PARAM_DOUBLE, "threshold");
END_PARAMETER_LIST();

AnastamosisGrowDirectionModifier::AnastamosisGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d AnastamosisGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	if (branch)
	{
		return previous_dir;
	}
	Segment * nearest_valid_target = mat->m_cult->tips.nearestCondition(tip.parent, [&tip](Segment * seg){return seg->seed() != tip.parent->seed(); });

	if (nearest_valid_target && (distance(nearest_valid_target->tip(1).pos().x,
			tip.parent->tip(1).pos().x,
		nearest_valid_target->tip(1).pos().y,
		tip.parent->tip(1).pos().y,
		nearest_valid_target->tip(1).pos().z, 
		tip.parent->tip(1).pos().z) > (search_radius + search_multiplier * seg_length )))
	{
		nearest_valid_target = nullptr;
	}

	if (nearest_valid_target)
	{
		//grow towards nearest valid target
		vec3d dir_to_nearest = nearest_valid_target->tip(1).pos() - tip.pos();
		//make this the same length as 
		double new_length = dir_to_nearest.unit();

		
		//reduce the length if too high
		seg_length = std::min(seg_length, new_length);
		if (seg_length == new_length)
		{
			//deactivate the tip
			tip.bactive = false;
			tip.parent->SetFlagOn(Segment::ANAST);
			//increment the anastamosis count of the underlying material
			mat->m_cult->m_num_anastom++;

			//TODO: consider adding connectivity information
			//TODO: consider setting the id of all reachable tips from this network to be the same id
			return dir_to_nearest;
		}

		return dir_to_nearest;
	}

	return previous_dir;
}

BEGIN_PARAMETER_LIST(AnastamosisGrowDirectionModifier, GrowDirectionModifier)
ADD_PARAMETER(search_radius, FE_PARAM_DOUBLE, "search_radius");
ADD_PARAMETER(search_multiplier, FE_PARAM_DOUBLE, "search_multiplier");
END_PARAMETER_LIST();

BranchGrowDirectionModifier::BranchGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d BranchGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	// If new segment is a branch we modify the grow direction a bit
	if (branch)
	{
		// TODO: what's the logic here? Why the 0.5 factor?
		//      If the vessel is aligned with the collagen (and the initial fragments are)
		//      then  the new branch will overlap the old segment.
		vec3d seg_vec = -previous_dir;
		double lambda;
		vec3d coll_fib = culture->m_pmat->CollagenDirection(tip.pt,lambda);

		seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
		seg_vec.unit();
		return seg_vec;
	}
	return previous_dir;
}

DefaultGrowDirectionModifier::DefaultGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{
	
}

vec3d DefaultGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	// Find the component of the new vessel direction determined by collagen fiber orientation    
	double lambda;
	vec3d coll_dir = culture->m_pmat->CollagenDirection(tip.pt,lambda);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = mix(per_dir, coll_dir, culture->m_pmat->m_cultureParams.GetWeightInterpolation(grow_time));
	new_dir.unit();

	return new_dir;
}
