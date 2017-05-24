#include "GrowDirectionModifier.h"
#include "FECore/ElementDataRecord.h"
#include "FEBioPlot/FEBioPlotFile2.h"
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

void GrowDirectionModifiers::Update()
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		grow_direction_modifiers[i]->Update();
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
		vec3d coll_fib = culture->m_pmat->fiber_manager->GetFiberDirection(tip.pt, lambda);

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
	vec3d coll_dir = culture->m_pmat->fiber_manager->GetFiberDirection(tip.pt, lambda);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = mix(per_dir, coll_dir, culture->m_pmat->m_cultureParams.GetWeightInterpolation(grow_time));
	new_dir.unit();

	return new_dir;
}

BaseFiberAwareGrowDirectionModifier::BaseFiberAwareGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}

vec3d BaseFiberAwareGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	// Find the component of the new vessel direction determined by collagen fiber orientation    
	double lambda;
	vec3d coll_dir = culture->m_pmat->fiber_manager->GetFiberDirection(tip.pt, lambda);

	double l_m1;
	double l_m2;
	vec3d m1 = culture->m_pmat->fiber_manager->GetMinorAxisDirection1(tip.pt, l_m1);
	vec3d m2 = culture->m_pmat->fiber_manager->GetMinorAxisDirection2(tip.pt, l_m2);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;
	per_dir.unit();
	coll_dir.unit();
	double theta = per_dir * coll_dir;
	double theta1 = per_dir * m1;
	double theta2 = per_dir * m2;


	double contrib = (theta * lambda + theta1*l_m1 + theta2*l_m2) / 3;
	contrib = abs(contrib);
	contrib = 1 / contrib;
	vec3d new_dir = mix(per_dir, coll_dir, culture->m_pmat->m_cultureParams.GetWeightInterpolation(grow_time)* contrib);
	new_dir.unit();

	return new_dir;
}

UnitLengthGrowDirectionModifier::UnitLengthGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d UnitLengthGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	seg_length = 1.0;
	return previous_dir;
}

DensityScaleGrowDirectionModifier::DensityScaleGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d DensityScaleGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	seg_length *= culture->FindDensityScale(tip.pt);
	return previous_dir;
}
SegmentLengthGrowDirectionModifier::SegmentLengthGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d SegmentLengthGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	seg_length *= culture->SegmentLength(start_time, grow_time);
	return previous_dir;
}

bool DataStoreLengthDoubleGrowDirectionModifier::Init()
{
	if (FEMaterial::Init() == false) return false;
	FEModel * model = GetFEModel();
	FEBioModel * bm = dynamic_cast<FEBioModel*>(model);
	assert(bm);
	DataStore & ds = bm->GetDataStore();
	for(int i =0; i < ds.Size();i++)
	{
		DataRecord * dr = ds.GetDataRecord(i);
		ElementDataRecord *er = dynamic_cast<ElementDataRecord*>(dr);
		if(er)
		{
			if(!strcmp(er->GetName(), field_name))
			{
				record_index = i;
				return true;
			}
		}
	}
	//data record name not found
	return false;
}
DataStoreLengthDoubleGrowDirectionModifier::DataStoreLengthDoubleGrowDirectionModifier(FEModel * model) : GrowDirectionModifier(model)
{

}
vec3d DataStoreLengthDoubleGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	DataStore & ds = culture->m_pmat->m_pangio->m_fem->GetDataStore();
	seg_length *= ds.GetDataRecord(record_index)->Evaluate(tip.pt.nelem, field);
	return previous_dir;
}

BEGIN_PARAMETER_LIST(DataStoreLengthDoubleGrowDirectionModifier, GrowDirectionModifier)
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");
END_PARAMETER_LIST();

void GDMArchive::reset()
{
	fpdata.clear();
}

void GDMArchive::WriteData(int nid, std::vector<float>& data)
{
	if(fpdata.size() <= nid)
	{
		vector<float> temp;
		fpdata.push_back(temp);
	}
		
	fpdata[nid-1].resize(data.size());
	std::copy(data.begin(), data.end(), fpdata[nid-1].begin());
}
mat3dd GDMArchive::GetDataMat3dd(int domain, int element_index)
{
	const int index = 3 * element_index;
	return mat3dd(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2]);
}
mat3ds GDMArchive::GetDataMat3ds(int domain, int element_index)
{
	const int index = 6 * element_index;
	return mat3ds(fpdata[domain][index], fpdata[domain][index+1], fpdata[domain][index+2],
		fpdata[domain][index+3], fpdata[domain][index+4], fpdata[domain][index+5]);

}
mat3d  GDMArchive::GetDataMat3d(int domain, int element_index)
{
	const int index = 9 * element_index;
	return mat3d(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2],
		fpdata[domain][index + 3], fpdata[domain][index + 4], fpdata[domain][index + 5],
		fpdata[domain][index + 6], fpdata[domain][index + 7], fpdata[domain][index + 8]);
}
float  GDMArchive::GetDataFloat(int domain, int element_index)
{
	const int index = element_index;
	return fpdata[domain][index];
}
vec3d  GDMArchive::GetDataVec3d(int domain, int element_index)
{
	const int index = 3 * element_index;
	return vec3d(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2]);
}


bool PlotFile2DoubleGrowDirectionModifier::Init()
{
	if (FEMaterial::Init() == false) return false;
	FEModel * model = GetFEModel();
	FEBioModel * bm = dynamic_cast<FEBioModel*>(model);
	assert(bm);
	FEBioPlotFile2 * pf2 = dynamic_cast<FEBioPlotFile2 *>(bm->GetPlotFile());
	if (!pf2)
		return false;
	const list<FEBioPlotFile2::DICTIONARY_ITEM>& domain_var_list = pf2->GetDictionary().DomainVariableList();


	for (auto i = domain_var_list.begin(); i != domain_var_list.end(); ++i)
	{
		if (!strcmp(i->m_szname, field_name))
		{
			record_index = i;
			return true;
		}
	}
	//data record name not found
	return false;
}

vec3d PlotFile2DoubleGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	switch(datatype)
	{
	case 0:
		seg_length *= archive.GetDataFloat(mat->domains[0], tip.pt.elemindex);
		break;
	case 1:
		seg_length *= archive.GetDataVec3d(mat->domains[0], tip.pt.elemindex).norm();
		break;
	case 2:
		seg_length *= archive.GetDataMat3ds(mat->domains[0], tip.pt.elemindex).det();
		break;
	default:
		assert(false);
	}
	
	return previous_dir;
}

void PlotFile2DoubleGrowDirectionModifier::Update()
{
	record_index->m_psave->Save(*culture->m_pmat->m_pangio->m_fem, archive);
}

BEGIN_PARAMETER_LIST(PlotFile2DoubleGrowDirectionModifier, GrowDirectionModifier)
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");
ADD_PARAMETER(datatype, FE_PARAM_INT, "datatype");
END_PARAMETER_LIST();

bool Plot2GGP::Init()
{
	if (FEMaterial::Init() == false) return false;
	FEModel * model = GetFEModel();
	FEBioModel * bm = dynamic_cast<FEBioModel*>(model);
	assert(bm);
	FEBioPlotFile2 * pf2 = dynamic_cast<FEBioPlotFile2 *>(bm->GetPlotFile());
	if (!pf2)
		return false;
	const list<FEBioPlotFile2::DICTIONARY_ITEM>& domain_var_list = pf2->GetDictionary().DomainVariableList();


	for (auto i = domain_var_list.begin(); i != domain_var_list.end(); ++i)
	{
		if (!strcmp(i->m_szname, field_name))
		{
			record_index = i;
			return true;
		}
	}
	//data record name not found
	return false;
}

mat3d Plot2GGP::Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip)
{
	switch (record_index->m_nfmt)
	{
	case Storage_Fmt::FMT_ITEM:
		{
			switch(record_index->m_ntype)
			{
			case Var_Type::PLT_FLOAT:
				in(0,0) = archive.GetDataFloat(mat->domains[0], tip.pt.elemindex);
				break;
			case Var_Type::PLT_MAT3F:
				{
					//make sure this order is correct
					mat3d temp = archive.GetDataMat3d(mat->domains[0], tip.pt.elemindex);
					in = in * temp;
				}
				break;
			case Var_Type::PLT_MAT3FD:
				{
					mat3dd temp = archive.GetDataMat3dd(mat->domains[0], tip.pt.elemindex);
					in = in *temp;
				}
				break;
			case Var_Type::PLT_MAT3FS:
				{
					mat3ds temp = archive.GetDataMat3ds(mat->domains[0], tip.pt.elemindex);
					in = in *temp;
				}
				break;
			case Var_Type::PLT_VEC3F:
				{
					vec3d col = archive.GetDataVec3d(mat->domains[0], tip.pt.elemindex);
					in(0, 0) = col.x;
					in(1, 0) = col.y;
					in(2, 0) = col.z;
				}
				break;
			default:
				assert(false);
			}
		}
		break;
	default:
		assert(false);
	}

	return GGP::Operation(in, mat, tip);
}

void Plot2GGP::Update()
{
	record_index->m_psave->Save(*culture->m_pmat->m_pangio->m_fem, archive);
	GGP::Update();
}

BEGIN_PARAMETER_LIST(Plot2GGP, GGP)
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");
END_PARAMETER_LIST();

MatrixConverterGGP::MatrixConverterGGP(FEModel * model) : GGP(model)
{
	//initialize the convrsion matrix as an identity matrix
	for(int i=0; i < 9;i++)
	{
		for(int j=0; j < 9; j++)
		{
			m[i][j] = 0.0;
		}
	}

	for(int i=0;i < 9;i++)
	{
		m[i][i] = 1.0;
	}
}

mat3d MatrixConverterGGP::Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip)
{
	double unroll[9];
	double results[9];
	int index = 0;
	for(int i =0;i<3;i++)
	{
		for(int j=0; j < 3;j++)
		{
			unroll[index] = in[i][j];
			index++;
		}
	}

	for(int i =0; i < 9;i++)
	{
		double sum = 0.0;
		for(int j=0; j< 9;j++)
		{
			sum += unroll[j] * m[i][j];
		}
		results[i] = sum;
	}


	index = 0;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			in[i][j] = results[index];
			index++;
		}
	}
	return GGP::Operation(in, mat, tip);
}

BEGIN_PARAMETER_LIST(MatrixConverterGGP, GGP)
ADD_PARAMETER(m[0][0], FE_PARAM_DOUBLE, "m11");
ADD_PARAMETER(m[0][1], FE_PARAM_DOUBLE, "m12");
ADD_PARAMETER(m[0][2], FE_PARAM_DOUBLE, "m13");
ADD_PARAMETER(m[0][3], FE_PARAM_DOUBLE, "m14");
ADD_PARAMETER(m[0][4], FE_PARAM_DOUBLE, "m15");
ADD_PARAMETER(m[0][5], FE_PARAM_DOUBLE, "m16");
ADD_PARAMETER(m[0][6], FE_PARAM_DOUBLE, "m17");
ADD_PARAMETER(m[0][7], FE_PARAM_DOUBLE, "m18");
ADD_PARAMETER(m[0][8], FE_PARAM_DOUBLE, "m19");

ADD_PARAMETER(m[1][0], FE_PARAM_DOUBLE, "m21");
ADD_PARAMETER(m[1][1], FE_PARAM_DOUBLE, "m22");
ADD_PARAMETER(m[1][2], FE_PARAM_DOUBLE, "m23");
ADD_PARAMETER(m[1][3], FE_PARAM_DOUBLE, "m24");
ADD_PARAMETER(m[1][4], FE_PARAM_DOUBLE, "m25");
ADD_PARAMETER(m[1][5], FE_PARAM_DOUBLE, "m26");
ADD_PARAMETER(m[1][6], FE_PARAM_DOUBLE, "m27");
ADD_PARAMETER(m[1][7], FE_PARAM_DOUBLE, "m28");
ADD_PARAMETER(m[1][8], FE_PARAM_DOUBLE, "m29");

ADD_PARAMETER(m[2][0], FE_PARAM_DOUBLE, "m31");
ADD_PARAMETER(m[2][1], FE_PARAM_DOUBLE, "m32");
ADD_PARAMETER(m[2][2], FE_PARAM_DOUBLE, "m33");
ADD_PARAMETER(m[2][3], FE_PARAM_DOUBLE, "m34");
ADD_PARAMETER(m[2][4], FE_PARAM_DOUBLE, "m35");
ADD_PARAMETER(m[2][5], FE_PARAM_DOUBLE, "m36");
ADD_PARAMETER(m[2][6], FE_PARAM_DOUBLE, "m37");
ADD_PARAMETER(m[2][7], FE_PARAM_DOUBLE, "m38");
ADD_PARAMETER(m[2][8], FE_PARAM_DOUBLE, "m39");

ADD_PARAMETER(m[3][0], FE_PARAM_DOUBLE, "m41");
ADD_PARAMETER(m[3][1], FE_PARAM_DOUBLE, "m42");
ADD_PARAMETER(m[3][2], FE_PARAM_DOUBLE, "m43");
ADD_PARAMETER(m[3][3], FE_PARAM_DOUBLE, "m44");
ADD_PARAMETER(m[3][4], FE_PARAM_DOUBLE, "m45");
ADD_PARAMETER(m[3][5], FE_PARAM_DOUBLE, "m46");
ADD_PARAMETER(m[3][6], FE_PARAM_DOUBLE, "m47");
ADD_PARAMETER(m[3][7], FE_PARAM_DOUBLE, "m48");
ADD_PARAMETER(m[3][8], FE_PARAM_DOUBLE, "m49");

ADD_PARAMETER(m[4][0], FE_PARAM_DOUBLE, "m51");
ADD_PARAMETER(m[4][1], FE_PARAM_DOUBLE, "m52");
ADD_PARAMETER(m[4][2], FE_PARAM_DOUBLE, "m53");
ADD_PARAMETER(m[4][3], FE_PARAM_DOUBLE, "m54");
ADD_PARAMETER(m[4][4], FE_PARAM_DOUBLE, "m55");
ADD_PARAMETER(m[4][5], FE_PARAM_DOUBLE, "m56");
ADD_PARAMETER(m[4][6], FE_PARAM_DOUBLE, "m57");
ADD_PARAMETER(m[4][7], FE_PARAM_DOUBLE, "m58");
ADD_PARAMETER(m[4][8], FE_PARAM_DOUBLE, "m59");

ADD_PARAMETER(m[5][0], FE_PARAM_DOUBLE, "m61");
ADD_PARAMETER(m[5][1], FE_PARAM_DOUBLE, "m62");
ADD_PARAMETER(m[5][2], FE_PARAM_DOUBLE, "m63");
ADD_PARAMETER(m[5][3], FE_PARAM_DOUBLE, "m64");
ADD_PARAMETER(m[5][4], FE_PARAM_DOUBLE, "m65");
ADD_PARAMETER(m[5][5], FE_PARAM_DOUBLE, "m66");
ADD_PARAMETER(m[5][6], FE_PARAM_DOUBLE, "m67");
ADD_PARAMETER(m[5][7], FE_PARAM_DOUBLE, "m68");
ADD_PARAMETER(m[5][8], FE_PARAM_DOUBLE, "m69");

ADD_PARAMETER(m[6][0], FE_PARAM_DOUBLE, "m71");
ADD_PARAMETER(m[6][1], FE_PARAM_DOUBLE, "m72");
ADD_PARAMETER(m[6][2], FE_PARAM_DOUBLE, "m73");
ADD_PARAMETER(m[6][3], FE_PARAM_DOUBLE, "m74");
ADD_PARAMETER(m[6][4], FE_PARAM_DOUBLE, "m75");
ADD_PARAMETER(m[6][5], FE_PARAM_DOUBLE, "m76");
ADD_PARAMETER(m[6][6], FE_PARAM_DOUBLE, "m77");
ADD_PARAMETER(m[6][7], FE_PARAM_DOUBLE, "m78");
ADD_PARAMETER(m[6][8], FE_PARAM_DOUBLE, "m79");

ADD_PARAMETER(m[7][0], FE_PARAM_DOUBLE, "m81");
ADD_PARAMETER(m[7][1], FE_PARAM_DOUBLE, "m82");
ADD_PARAMETER(m[7][2], FE_PARAM_DOUBLE, "m83");
ADD_PARAMETER(m[7][3], FE_PARAM_DOUBLE, "m84");
ADD_PARAMETER(m[7][4], FE_PARAM_DOUBLE, "m85");
ADD_PARAMETER(m[7][5], FE_PARAM_DOUBLE, "m86");
ADD_PARAMETER(m[7][6], FE_PARAM_DOUBLE, "m87");
ADD_PARAMETER(m[7][7], FE_PARAM_DOUBLE, "m88");
ADD_PARAMETER(m[7][8], FE_PARAM_DOUBLE, "m89");

ADD_PARAMETER(m[8][0], FE_PARAM_DOUBLE, "m91");
ADD_PARAMETER(m[8][1], FE_PARAM_DOUBLE, "m92");
ADD_PARAMETER(m[8][2], FE_PARAM_DOUBLE, "m93");
ADD_PARAMETER(m[8][3], FE_PARAM_DOUBLE, "m94");
ADD_PARAMETER(m[8][4], FE_PARAM_DOUBLE, "m95");
ADD_PARAMETER(m[8][5], FE_PARAM_DOUBLE, "m96");
ADD_PARAMETER(m[8][6], FE_PARAM_DOUBLE, "m97");
ADD_PARAMETER(m[8][7], FE_PARAM_DOUBLE, "m98");
ADD_PARAMETER(m[8][8], FE_PARAM_DOUBLE, "m99");

END_PARAMETER_LIST();


mat3d ForkedGGP::Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip)
{
	mat3d temp = nest->Operation(in, mat, tip);
	return GGP::Operation(temp, mat, tip);
}

mat3d EigenValuesGGP::Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip)
{
	mat3ds temp(in[0][0], in[1][1], in[2][2], in[0][1], in[1][2], in[0][2]);
	double eig_val[3];
	vec3d eig_vec[3];
	temp.eigen(eig_val, eig_vec);

	mat3d rv(eig_val[0], 0.0,0.0,
		eig_val[1],0.0, 0.0,
		eig_val[2], 0.0, 0.0);
	return rv;
}

mat3d EigenVectorsGGP::Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip)
{
	mat3ds temp(in[0][0], in[1][1], in[2][2], in[0][1], in[1][2], in[0][2]);
	double eig_val[3];
	vec3d eig_vec[3];
	temp.eigen(eig_val, eig_vec);

	mat3d rv(eig_vec[0].x, eig_vec[0].y, eig_vec[0].z,
		eig_vec[1].x, eig_vec[1].y, eig_vec[1].z,
		eig_vec[2].x, eig_vec[2].y, eig_vec[2].z);
	return rv;
}