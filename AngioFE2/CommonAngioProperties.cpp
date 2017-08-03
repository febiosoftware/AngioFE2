#include "CommonAngioProperites.h"

CommonAngioProperties::CommonAngioProperties(FEModel * pfem) : FEMaterial(pfem)
{
	AddProperty(&vessel_material, "vessel");
	
	AddProperty(&fbrancher, "brancher");
	AddProperty(&fiber_initializer, "fiber_initializer");
	AddProperty(&fseeder, "fragment_seeder");
	AddProperty(&bc, "boundary_condition");
	AddProperty(&gdms, "grow_direction_modifiers");
}

void CommonAngioProperties::InitializeFibers(FiberManager * man)
{
	fiber_initializer->InitializeFibers(man);
}

void CommonAngioProperties::UpdateGDMs()
{
	gdms->Update();
}

void CommonAngioProperties::SetCulture(Culture * culture)
{
	fbrancher->SetCulture(culture);
	fseeder->SetCulture(culture);
	gdms->SetCulture(culture);
	bc->SetCulture(culture);
}