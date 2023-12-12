#ifndef HOUSEHOLD
#define HOUSEHOLD

#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/Random.h"
#include "Location.h"
#include <map>
class Household{
private:
	repast::AgentId householdId;
	Location* assignedField;
	int maizeStorage;
	int age;
	int deathAge;
    /* Map to store closeness to other households */
    std::map<repast::AgentId, double > closenessMap;

public:
	Household(repast::AgentId id,int a, int deathAge, int mStorage);
	~Household();

	/* Required Getters */
	virtual repast::AgentId& getId() { return householdId; }
	virtual const repast::AgentId& getId() const { return householdId; }

	void setCloseness(const repast::AgentId& otherId, double value);

    // Getters for Proximity
	double getCloseness(const repast::AgentId& otherId)	const;

	int getLoanMaize(int needs);
	int getMaize(void);
	int getlackMaize(int needs);
	void addMaize(int maize);
	void removeMaize(int maize);

	/* Getters specific to this kind of Agent */
	Location* getAssignedField(){return assignedField; }
	int splitMaizeStored(int percentage);
	
	bool checkMaize(int needs);
	bool death();
	bool fission(int minFissionAge, int maxFissionAge, double gen, double fProb);
	void nextYear(int needs);
	void chooseField(Location* Field);
};

#endif
