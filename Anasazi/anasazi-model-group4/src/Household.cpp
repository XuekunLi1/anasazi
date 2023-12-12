#include "Household.h"
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include <stdio.h>
#include "repast_hpc/Random.h"

Household::Household(repast::AgentId id, int a, int deAge, int mStorage)
{
	householdId = id;
	age = a;
	deathAge = deAge;
	maizeStorage = mStorage;
	assignedField = NULL;
}

Household::~Household()
{

}

int Household::splitMaizeStored(int percentage)
{
	int maizeEndowment;
	maizeEndowment = maizeStorage * percentage;
	maizeStorage = maizeStorage - maizeEndowment;
	return maizeEndowment;
}

int Household::getLoanMaize(int needs)
{
	return (assignedField->getExpectedYield() + maizeStorage - needs);
}

int Household::getMaize(void)
{
	return (assignedField->getExpectedYield() + maizeStorage);
}

int Household::getlackMaize(int needs)
{
	return (assignedField->getExpectedYield() + maizeStorage - needs);
}

bool Household::checkMaize(int needs)
{
	if((assignedField->getExpectedYield() + maizeStorage) > needs)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Household::addMaize(int maize)
{
	maizeStorage = maizeStorage + maize;
}
void Household::removeMaize(int maize)
{
	maizeStorage = maizeStorage - maize;
}

bool Household::death()
{
	if(age>=deathAge)
	{
			return true;
	}
	else
	{
		return false;
	}
}

bool Household::fission(int minFissionAge, int maxFissionAge, double gen, double fProb)
{
	if((age>=minFissionAge && age<=maxFissionAge) && (gen <= fProb))
	{
			return true;
	}
	else
	{
		return false;
	}
}

void Household::nextYear(int needs)
{
	age++;
	maizeStorage = assignedField->getExpectedYield() + maizeStorage - needs;
}

void Household::chooseField(Location* Field)
{
	//set the old location as emtpy
	if (assignedField!=NULL) assignedField->setState(0);

	//set the new location as a field
	Field->setState(2);

	assignedField = Field;
}


/* Setters for Closeness*/
void Household::setCloseness(const repast::AgentId& otherId, double value)
{
    closenessMap[otherId] = value;
}

// Getters for Proximity
double Household::getCloseness(const repast::AgentId& otherId) const
{
    auto it = closenessMap.find(otherId);
    if (it != closenessMap.end()) 
	{
        return it->second;
    }
    return -1.0; // Default value if not found
 }

