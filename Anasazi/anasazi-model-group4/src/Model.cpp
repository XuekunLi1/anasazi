#include <stdio.h>
#include <vector>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/Random.h"
#include "repast_hpc/Schedule.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"
#include <string>
#include <fstream>
#include <stdlib.h>
#include "repast_hpc/Moore2DGridQuery.h"

#include "Model.h"

// substracts b<T> to a<T>
template <typename T>
void substract_vector(std::vector<T>& a, const std::vector<T>& b)
{
	typename std::vector<T>::iterator       it = a.begin();
	typename std::vector<T>::const_iterator it2 = b.begin();

	while (it != a.end())
	{
		while (it2 != b.end() && it != a.end())
		{
			if (*it == *it2)
			{
				it = a.erase(it);
				it2 = b.begin();
			}

			else
				++it2;
		}
		if (it != a.end())
			++it;

		it2 = b.begin();
	}
}

AnasaziModel::AnasaziModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm): context(comm) , locationContext(comm)
{
	props = new repast::Properties(propsFile, argc, argv, comm);
	boardSizeX = repast::strToInt(props->getProperty("board.size.x"));
	boardSizeY = repast::strToInt(props->getProperty("board.size.y"));

	initializeRandom(*props, comm);
	repast::Point<double> origin(0,0);
	repast::Point<double> extent(boardSizeX, boardSizeY);
	repast::GridDimensions gd (origin, extent);

	int procX = repast::strToInt(props->getProperty("proc.per.x"));
	int procY = repast::strToInt(props->getProperty("proc.per.y"));
	int bufferSize = repast::strToInt(props->getProperty("grid.buffer"));

	std::vector<int> processDims;
	processDims.push_back(procX);
	processDims.push_back(procY);
	householdSpace = new repast::SharedDiscreteSpace<Household, repast::StrictBorders, repast::SimpleAdder<Household> >("AgentDiscreteSpace",gd,processDims,bufferSize, comm);
	locationSpace = new repast::SharedDiscreteSpace<Location, repast::StrictBorders, repast::SimpleAdder<Location> >("LocationDiscreteSpace",gd,processDims,bufferSize, comm);

	context.addProjection(householdSpace);
	locationContext.addProjection(locationSpace);

	param.startYear = repast::strToInt(props->getProperty("start.year"));
	param.endYear = repast::strToInt(props->getProperty("end.year"));
	param.maxStorageYear = repast::strToInt(props->getProperty("max.store.year"));
	param.maxStorage = repast::strToInt(props->getProperty("max.storage"));
	param.householdNeed = repast::strToInt(props->getProperty("household.need"));
	param.minFissionAge = repast::strToInt(props->getProperty("min.fission.age"));
	param.maxFissionAge = repast::strToInt(props->getProperty("max.fission.age"));
	param.minDeathAge = repast::strToInt(props->getProperty("min.death.age"));
	param.maxDeathAge = repast::strToInt(props->getProperty("max.death.age"));
	param.maxDistance = repast::strToInt(props->getProperty("max.distance"));
	param.initMinCorn = repast::strToInt(props->getProperty("initial.min.corn"));
	param.initMaxCorn = repast::strToInt(props->getProperty("initial.max.corn"));

	param.annualVariance = repast::strToDouble(props->getProperty("annual.variance"));
	param.spatialVariance = repast::strToDouble(props->getProperty("spatial.variance"));
	param.fertilityProbability = repast::strToDouble(props->getProperty("fertility.prop"));
	param.harvestAdjustment = repast::strToDouble(props->getProperty("harvest.adj"));
	param.maizeStorageRatio = repast::strToDouble(props->getProperty("new.household.ini.maize"));
	param.thresholdSharefood = repast::strToDouble(props->getProperty("threshold.sharefood"));
	param.b1 = repast::strToDouble(props->getProperty("B1"));//new added
	param.b3 = repast::strToDouble(props->getProperty("B3"));//new added
	param.b4 = repast::strToDouble(props->getProperty("B4"));//new added

	year = param.startYear;
	stopAt = param.endYear - param.startYear + 1;
	fissionGen = new repast::DoubleUniformGenerator(repast::Random::instance()->createUniDoubleGenerator(0,1));
	deathAgeGen = new repast::IntUniformGenerator(repast::Random::instance()->createUniIntGenerator(param.minDeathAge,param.maxDeathAge));
	yieldGen = new repast::NormalGenerator(repast::Random::instance()->createNormalGenerator(0,param.annualVariance));
	soilGen = new repast::NormalGenerator(repast::Random::instance()->createNormalGenerator(0,param.spatialVariance));
	initAgeGen = new repast::IntUniformGenerator(repast::Random::instance()->createUniIntGenerator(0,param.minDeathAge));
	initMaizeGen = new repast::IntUniformGenerator(repast::Random::instance()->createUniIntGenerator(param.initMinCorn,param.initMaxCorn));

	Agent_Number = 14; //Number of initial agents
	std::vector<std::vector<int>> adjacencyMatrix(Agent_Number, std::vector<int>(Agent_Number, 0));
	initnetwork();

	string resultFile = props->getProperty("result.file");
	out.open(resultFile);
	out << "Year,Number-of-Households,maxCapacity" << endl;
}

AnasaziModel::~AnasaziModel()
{
	delete props;
	out.close();
}

void AnasaziModel::initAgents()
{
	int rank = repast::RepastProcess::instance()->rank();

	int LocationID = 0;
	for(int i=0; i<boardSizeX; i++ )
	{
		for(int j=0; j<boardSizeY; j++)
		{
			repast::AgentId id(LocationID, rank, 1);
			Location* agent = new Location(id, soilGen->next());
			locationContext.addAgent(agent);
			locationSpace->moveTo(id, repast::Point<int>(i, j));
			LocationID++;
		}
	}

	readCsvMap();
	readCsvWater();
	readCsvPdsi();
	readCsvHydro();
	int noOfAgents  = repast::strToInt(props->getProperty("count.of.agents"));
	repast::IntUniformGenerator xGen = repast::IntUniformGenerator(repast::Random::instance()->createUniIntGenerator(0,boardSizeX-1));
	repast::IntUniformGenerator yGen = repast::IntUniformGenerator(repast::Random::instance()->createUniIntGenerator(0,boardSizeY-1));
	for(int i =0; i< noOfAgents;i++)
	{
		repast::AgentId id(houseID, rank, 2);
		int initAge = initAgeGen->next();
		int mStorage = initMaizeGen->next();
		Household* agent = new Household(id, initAge, deathAgeGen->next(), mStorage);
		context.addAgent(agent);
		std::vector<Location*> locationList;

		newLocation:
		int x = xGen.next();
		int y = yGen.next();
		locationSpace->getObjectsAt(repast::Point<int>(x, y), locationList);

		if(locationList[0]->getState()==2)
		{
			locationList.clear();
			goto newLocation;
		}
		else
		{
			householdSpace->moveTo(id, repast::Point<int>(x, y));
			locationList[0]->setState(1);
		}
		houseID++;
	}
	updateCloseness();
	updateLocationProperties();

	repast::SharedContext<Household>::const_iterator local_agents_iter = context.begin();
	repast::SharedContext<Household>::const_iterator local_agents_end = context.end();

	while(local_agents_iter != local_agents_end)
	{
		Household* household = (&**local_agents_iter);
		if(household->death())
		{
			repast::AgentId id = household->getId();
			local_agents_iter++;

			std::vector<int> loc;
			householdSpace->getLocation(id, loc);

			std::vector<Location*> locationList;
			if(!loc.empty())
			{
				locationSpace->getObjectsAt(repast::Point<int>(loc[0], loc[1]), locationList);
				locationList[0]->setState(0);
			}
			context.removeAgent(id);
		}
		else
		{
			local_agents_iter++;
			fieldSearch(household);
		}
	}

}

void AnasaziModel::doPerTick()
{
	updateLocationProperties();
	writeOutputToFile();
	year++;
	updateHouseholdProperties();
	updateNetwork();//new added
}

void AnasaziModel::initSchedule(repast::ScheduleRunner& runner)
{
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<AnasaziModel> (this, &AnasaziModel::doPerTick)));
	runner.scheduleStop(stopAt);
}

void AnasaziModel::readCsvMap()
{
	int x,y,z , mz;
	string zone, maizeZone, temp;

	std::ifstream file ("data/map.csv");//define file object and open map.csv
	file.ignore(500,'\n');//Ignore first line

	while(1)//read until end of file
	{
		getline(file,temp,',');
		if(!temp.empty())
		{
			x = repast::strToInt(temp); //Read until ',' and convert to int & store in x
			getline(file,temp,',');
			y = repast::strToInt(temp); //Read until ',' and convert to int & store in y
			getline(file,temp,','); //colour
			getline(file,zone,',');// read until ',' and store into zone
			getline(file,maizeZone,'\n');// read until next line and store into maizeZone
			if(zone == "\"Empty\"")
			{
				z = 0;
			}
			else if(zone == "\"Natural\"")
			{
				z = 1;
			}
			else if(zone == "\"Kinbiko\"")
			{
				z = 2;
			}
			else if(zone == "\"Uplands\"")
			{
				z = 3;
			}
			else if(zone == "\"North\"")
			{
				z = 4;
			}
			else if(zone == "\"General\"")
			{
				z = 5;
			}
			else if(zone == "\"North Dunes\"")
			{
				z = 6;
			}
			else if(zone == "\"Mid Dunes\"")
			{
				z = 7;
			}
			else if(zone == "\"Mid\"")
			{
				z = 8;
			}
			else
			{
				z = 99;
			}

			if(maizeZone.find("Empty") != std::string::npos)
			{
				mz = 0;
			}
			else if(maizeZone.find("No_Yield") != std::string::npos)
			{
				mz = 1;
			}
			else if(maizeZone.find("Yield_1") != std::string::npos)
			{
				mz = 2;
			}
			else if(maizeZone.find("Yield_2") != std::string::npos)
			{
				mz = 3;
			}
			else if(maizeZone.find("Yield_3") != std::string::npos)
			{
				mz = 4;
			}
			else if(maizeZone.find("Sand_dune") != std::string::npos)
			{
				mz = 5;
			}
			else
			{
				mz = 99;
			}
			std::vector<Location*> locationList;
			locationSpace->getObjectsAt(repast::Point<int>(x, y), locationList);
			locationList[0]->setZones(z,mz);
		}
		else{
			goto endloop;
		}
	}
	endloop: ;
}

void AnasaziModel::readCsvWater()
{
	//read "type","start date","end date","x","y"
	int type, startYear, endYear, x, y;
	string temp;

	std::ifstream file ("data/water.csv");//define file object and open water.csv
	file.ignore(500,'\n');//Ignore first line
	while(1)//read until end of file
	{
		getline(file,temp,',');
		if(!temp.empty())
		{
			getline(file,temp,',');
			getline(file,temp,',');
			getline(file,temp,',');
			type = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			startYear = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			endYear = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			x = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,'\n');
			y = repast::strToInt(temp); //Read until ',' and convert to int

			std::vector<Location*> locationList;
			locationSpace->getObjectsAt(repast::Point<int>(x, y), locationList);
			locationList[0]->addWaterSource(type,startYear, endYear);
			//locationList[0]->checkWater(existStreams, existAlluvium, x, y, year);
		}
		else
		{
			goto endloop;
		}
	}
	endloop: ;
}

void AnasaziModel::readCsvPdsi()
{
	//read "year","general","north","mid","natural","upland","kinbiko"
	int i=0;
	string temp;

	std::ifstream file ("data/pdsi.csv");//define file object and open pdsi.csv
	file.ignore(500,'\n');//Ignore first line

	while(1)//read until end of file
	{
		getline(file,temp,',');
		if(!temp.empty())
		{
			pdsi[i].year = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			pdsi[i].pdsiGeneral = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			pdsi[i].pdsiNorth = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			pdsi[i].pdsiMid = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			pdsi[i].pdsiNatural = repast::strToDouble(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			pdsi[i].pdsiUpland = repast::strToDouble(temp); //Read until ',' and convert to int
			getline(file,temp,'\n');
			pdsi[i].pdsiKinbiko = repast::strToDouble(temp); //Read until ',' and convert to double
			i++;
		}
		else{
			goto endloop;
		}
	}
	endloop: ;
}

void AnasaziModel::readCsvHydro()
{
	//read "year","general","north","mid","natural","upland","kinbiko"
	string temp;
	int i =0;

	std::ifstream file ("data/hydro.csv");//define file object and open hydro.csv
	file.ignore(500,'\n');//Ignore first line

	while(1)//read until end of file
	{
		getline(file,temp,',');
		if(!temp.empty())
		{
			hydro[i].year = repast::strToInt(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			hydro[i].hydroGeneral = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			hydro[i].hydroNorth = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			hydro[i].hydroMid = repast::strToDouble(temp); //Read until ',' and convert to double
			getline(file,temp,',');
			hydro[i].hydroNatural = repast::strToDouble(temp); //Read until ',' and convert to int
			getline(file,temp,',');
			hydro[i].hydroUpland = repast::strToDouble(temp); //Read until ',' and convert to int
			getline(file,temp,'\n');
			hydro[i].hydroKinbiko = repast::strToDouble(temp); //Read until ',' and convert to double
			i++;
		}
		else
		{
			goto endloop;
		}
	}
	endloop: ;
}

int AnasaziModel::yieldFromPdsi(int zone, int maizeZone)
{
	int pdsiValue, row, col;
	switch(zone)
	{
		case 1:
			pdsiValue = pdsi[year-param.startYear].pdsiNatural;
			break;
		case 2:
			pdsiValue = pdsi[year-param.startYear].pdsiKinbiko;
			break;
		case 3:
			pdsiValue = pdsi[year-param.startYear].pdsiUpland;
			break;
		case 4:
		case 6:
			pdsiValue = pdsi[year-param.startYear].pdsiNorth;
			break;
		case 5:
			pdsiValue = pdsi[year-param.startYear].pdsiGeneral;
			break;
		case 7:
		case 8:
			pdsiValue = pdsi[year-param.startYear].pdsiMid;
			break;
		default:
			return 0;
	}

	/* Rows of pdsi table*/
	if(pdsiValue < -3)
	{
		row = 0;
	}
	else if(pdsiValue >= -3 && pdsiValue < -1)
	{
		row = 1;
	}
	else if(pdsiValue >= -1 && pdsiValue < 1)
	{
		row = 2;
	}
	else if(pdsiValue >= 1 && pdsiValue < 3)
	{
		row = 3;
	}
	else if(pdsiValue >= 3)
	{
		row = 4;
	}
	else
	{
		return 0;
	}

	/* Col of pdsi table*/
	if(maizeZone >= 2)
	{
		col = maizeZone - 2;
	}
	else
	{
		return 0;
	}

	return yieldLevels[row][col];
}

double AnasaziModel::hydroLevel(int zone)
{
	switch(zone)
	{
		case 1:
			return hydro[year-param.startYear].hydroNatural;
		case 2:
			return hydro[year-param.startYear].hydroKinbiko;
		case 3:
			return hydro[year-param.startYear].hydroUpland;
		case 4:
		case 6:
			return hydro[year-param.startYear].hydroNorth;
		case 5:
			return hydro[year-param.startYear].hydroGeneral;
		case 7:
		case 8:
			return hydro[year-param.startYear].hydroMid;
		default:
			return 0;
	}
}

void AnasaziModel::checkWaterConditions()
{
	if ((year >= 280 && year < 360) or (year >= 800 && year < 930) or (year >= 1300 && year < 1450))
	{
		existStreams = true;
	}
	else
	{
		existStreams = false;
	}

	if (((year >= 420) && (year < 560)) or ((year >= 630) && (year < 680)) or	((year >= 980) && (year < 1120)) or ((year >= 1180) && (year < 1230)))
	{
		existAlluvium = true;
	}
	else
	{
		existAlluvium = false;
	}
}

void AnasaziModel::writeOutputToFile()
{
	out << year << "," <<  context.size() << "," << maxCapacity << std::endl;
}

void  AnasaziModel::updateLocationProperties()
{
	checkWaterConditions();
	int x = 0;
	int allHarvest = 0;
	for(int i=0; i<boardSizeX; i++ )
	{
		for(int j=0; j<boardSizeY; j++)
		{
			std::vector<Location*> locationList;
			locationSpace->getObjectsAt(repast::Point<int>(i, j), locationList);
			locationList[0]->checkWater(existStreams,existAlluvium, i, j, year);
			int mz = locationList[0]->getMaizeZone();
			int z = locationList[0]->getZone();
			int y = yieldFromPdsi(z,mz);
			locationList[0]->calculateYield(y, param.harvestAdjustment, yieldGen->next());
			if(locationList[0]->getExpectedYield() >= param.householdNeed)
			{
				allHarvest ++;
			}
		}
	}
	maxCapacity = allHarvest;
}

void AnasaziModel::updateHouseholdProperties()
{
	repast::SharedContext<Household>::const_iterator local_agents_iter = context.begin();
	repast::SharedContext<Household>::const_iterator local_agents_end = context.end();

	while(local_agents_iter != local_agents_end)
	{
		Household* household = (&**local_agents_iter);
		if(household->death())
		{
			local_agents_iter++;
			removeHousehold(household);

		}
		else
		{
			local_agents_iter++;
			if(household->fission(param.minFissionAge,param.maxFissionAge, fissionGen->next(), param.fertilityProbability))
			{
				int rank = repast::RepastProcess::instance()->rank();
				repast::AgentId id(houseID, rank, 2);
				int mStorage = household->splitMaizeStored(param.maizeStorageRatio);
				Household* newAgent = new Household(id, 0, deathAgeGen->next(), mStorage);
				context.addAgent(newAgent);

				std::vector<int> loc;
				householdSpace->getLocation(household->getId(), loc);
				householdSpace->moveTo(id, repast::Point<int>(loc[0], loc[1]));
				fieldSearch(newAgent);
				addNewAgentContacts(houseID);//new added
				houseID++;
			}

			bool fieldFound = true;
			std::vector<int>  locgoal,loctemp;
			if(!(household->checkMaize(param.householdNeed)))
			{
				
				if(!ShareFood(household))
				{
					householdSpace->getLocation(household->getId(), loctemp);
					if(!loctemp.empty())
					{
						locgoal.assign(loctemp.begin(), loctemp.end());
					}
					fieldFound = fieldSearch(household);
				}
				
			}
			if(fieldFound)
			{
				if(Relocateflag)
				{
					MovewithFriends(locgoal, household);
				}
				household->nextYear(param.householdNeed);
			}
			//if(moveoutflag)
			//{
			//	Moveout(locgoal,household);
			//	moveoutflag = false;
			//}
		}
	}
	updateCloseness();
}

bool AnasaziModel::fieldSearch(Household* household)
{
	/******** Choose Field ********/
	std::vector<int> loc;
	householdSpace->getLocation(household->getId(), loc);
	repast::Point<int> center(loc);

	std::vector<Location*> neighbouringLocations;
	std::vector<Location*> checkedLocations;
	repast::Moore2DGridQuery<Location> moore2DQuery(locationSpace);
	int range = 1;
	while(1)
	{
		moore2DQuery.query(loc, range, false, neighbouringLocations);

		for (std::vector<Location*>::iterator it = neighbouringLocations.begin() ; it != neighbouringLocations.end(); ++it)
		{
			Location* tempLoc = (&**it);
			if(tempLoc->getState() == 0)
			{
				if(tempLoc->getExpectedYield() >= param.householdNeed)
				{
					std::vector<int> loc;
					locationSpace->getLocation(tempLoc->getId(), loc);
					household->chooseField(tempLoc);
					goto EndOfLoop;
				}
			}
		}
		range++;
		if(range > boardSizeY)
		{
			removeHousehold(household);
			moveoutflag = true;
			Relocateflag = false;
			return false;
		}
	}
	EndOfLoop:
	if(range >= 10)
	{
		return relocateHousehold(household);
	}
	else
	{
		Relocateflag = false;
		return true;
	}
}

void AnasaziModel::removeHousehold(Household* household)
{
	repast::AgentId id = household->getId();

	std::vector<int> loc;
	householdSpace->getLocation(id, loc);

	std::vector<Location*> locationList;
	std::vector<Household*> householdList;
	if(!loc.empty())
	{
		locationSpace->getObjectsAt(repast::Point<int>(loc[0], loc[1]), locationList);
		householdSpace->getObjectsAt(repast::Point<int>(loc[0], loc[1]), householdList);
		if(householdList.size() == 1)
		{
			locationList[0]->setState(0);
		}
		if(household->getAssignedField()!= NULL)
		{
			std::vector<int> loc;
			locationSpace->getLocation(household->getAssignedField()->getId(), loc);
			locationSpace->getObjectsAt(repast::Point<int>(loc[0], loc[1]), locationList);
			locationList[0]->setState(0);
		}
	}

	context.removeAgent(id);
}

bool AnasaziModel::relocateHousehold(Household* household)
{
	std::vector<Location*> neighbouringLocations;
	std::vector<Location*> suitableLocations;
	std::vector<Location*> waterSources;
	std::vector<Location*> checkedLocations;

	std::vector<int> loc, loc2;
	locationSpace->getLocation(household->getAssignedField()->getId(), loc);
	householdSpace->getLocation(household->getId(),loc2);

	locationSpace->getObjectsAt(repast::Point<int>(loc2[0], loc2[1]), neighbouringLocations);
	Location* householdLocation = neighbouringLocations[0];

	repast::Point<int> center(loc);
	repast::Moore2DGridQuery<Location> moore2DQuery(locationSpace);
	int range = floor(param.maxDistance/100);
	int i = 1;
	bool conditionC = true;

	//get all !Field with 1km
	LocationSearch:
		moore2DQuery.query(loc, range*i, false, neighbouringLocations);
		for (std::vector<Location*>::iterator it = neighbouringLocations.begin() ; it != neighbouringLocations.end(); ++it)
		{
			Location* tempLoc = (&**it);
			if(tempLoc->getState() != 2)
			{
				if(householdLocation->getExpectedYield() < tempLoc->getExpectedYield() && conditionC == true)
				{
					suitableLocations.push_back(tempLoc);
				}
				if(tempLoc->getWater())
				{
					waterSources.push_back(tempLoc);
				}
			}
		}
		if(suitableLocations.size() == 0 || waterSources.size() == 0)
		{
			if(conditionC == true)
			{
				conditionC = false;
			}
			else
			{
				conditionC = true;
				i++;
				if(range*i > boardSizeY)
				{
					removeHousehold(household);
					moveoutflag = true;
					Relocateflag = false;
					return false;
				}
			}
			goto LocationSearch;
		}
		else if(suitableLocations.size() == 1)
		{
			std::vector<int> loc2;
			locationSpace->getLocation(suitableLocations[0]->getId(),loc2);
			householdSpace->moveTo(household->getId(),repast::Point<int>(loc2[0], loc2[1]));
			Relocateflag = true;
			return true;
		}
		else
		{
			std::vector<int> point1, point2;
			std::vector<double> distances;
			for (std::vector<Location*>::iterator it1 = suitableLocations.begin() ; it1 != suitableLocations.end(); ++it1)
			{
				locationSpace->getLocation((&**it1)->getId(),point1);
				for (std::vector<Location*>::iterator it2 = waterSources.begin() ; it2 != waterSources.end(); ++it2)
				{
					locationSpace->getLocation((&**it2)->getId(),point2);
					double distance = sqrt(pow((point1[0]-point2[0]),2) + pow((point1[1]-point2[1]),2));
					distances.push_back(distance);
				}
			}
			int minElementIndex = std::min_element(distances.begin(),distances.end()) - distances.begin();
			minElementIndex = minElementIndex / waterSources.size();
			std::vector<int> loc2;
			locationSpace->getLocation(suitableLocations[minElementIndex]->getId(),loc2);
			householdSpace->moveTo(household->getId(),repast::Point<int>(loc2[0], loc2[1]));
			Relocateflag = true;
			return true;
		}
}

void AnasaziModel::updateCloseness(void)
{
	std::vector<int> loc,loc2;
	double ClosenessTemp;
	repast::SharedContext<Household>::const_iterator local_agents_iter = context.begin();
	repast::SharedContext<Household>::const_iterator local_agents_end = context.end();

    for (auto it1 = local_agents_iter; it1 != local_agents_end; ++it1) 
	{
        Household* household1 = &(**it1);

        // Iterate over all other households
        for (auto it2 = local_agents_iter; it2 != local_agents_end; ++it2) 
		{
            if (it1 == it2) continue; // Skip the same household
            Household* household2 = &(**it2);

            // Check if both households are in a location
            if (household1->getAssignedField() != nullptr && household2->getAssignedField() != nullptr) 
			{
				householdSpace->getLocation(household1->getId(),loc);
				householdSpace->getLocation(household2->getId(),loc2);
				if((loc[0] == loc2[0]) && (loc[1] == loc2[1]))
				{						
					if(household1->getCloseness(household2->getId()) == -1)
					{
						Closeness = new repast::NormalGenerator(repast::Random::instance()->createNormalGenerator(0.5,0.1)); 
						household1->setCloseness(household2->getId(), Closeness->next());
					}
					else
					{
						//closeness is equal to a random number that follows a normal distribution.
						ClosenessTemp = household1->getCloseness(household2->getId());
						ClosenessTemp += 0.01;
						household1->setCloseness(household2->getId(), ClosenessTemp);
					}
				}
				else
				{
					if(household1->getCloseness(household2->getId()) != -1)
					{
						double distance = sqrt(pow((loc[0]-loc2[0]),2) + pow((loc[1]-loc2[1]),2));
						ClosenessTemp = household1->getCloseness(household2->getId());
						ClosenessTemp -= 0.0001 * distance;
						household1->setCloseness(household2->getId(), ClosenessTemp);	
					}
				}
            }
        }
    }
	
}

bool AnasaziModel::ShareFood(Household* household)
{
	std::vector<int> loc;
	std::vector<Household*> householdList;
	std::vector<Household*> tempHouseholdList;
	bdi BDIgiver;
	bdi BDIreciever;
	int loanMaize;
	int addedMaize = 0;
	double ClosenessTemp;
//	int spareMaize = 0;
	bool Shareflag = false;
	bool ShareSucflag = false;
	bool Liveflag = false;

	BDIgiver.w0 = 1;
	BDIreciever.w0 = 1;
	householdSpace->getLocation(household->getId(),loc);
	householdSpace->getObjectsAt(repast::Point<int>(loc[0], loc[1]), householdList);
	
	std::sort(householdList.begin(), householdList.end(), 
        [household](const Household* a, const Household* b) {
            return a->getCloseness(household->getId()) < b->getCloseness(household->getId());
    });

	for (std::vector<Household*>::iterator it = householdList.begin() ; it != householdList.end(); ++it)
	{
		Household* tempHousehold = (&**it);
		repast::AgentId id1 = household->getId();
		repast::AgentId id2 = household->getId();
		int id3 = id1.id();
		int id4 = id2.id();
		if(household->getCloseness(	tempHousehold->getId()) >= param.thresholdSharefood && checkConnection(id3,id4))  //new added
		{
			if(tempHousehold->checkMaize(param.householdNeed))
			{
				loanMaize = tempHousehold->getLoanMaize(param.householdNeed);
				if(household->getlackMaize(param.householdNeed)<=loanMaize)
				{			
					//BDI model
					//TODO: add code of foodshare proposal
					//BDIgiver.ABx = BDIgiver.w0*((-4)*pow((household->getCloseness(tempHousehold->getId()))-0.5, 2)+1);
					//BDIgiver.ADx = BDIgiver.w1*(-exp((-1/400)*(household->getMaize()-param.householdNeed)+1))
					//+ BDIgiver.w2*(-pow(tempHousehold->getCloseness(household->getId())-1, 4)+1);	
					//BDIreciever.ABx = -BDIreciever.w0*(pow((tempHousehold->getCloseness(household->getId()))-1, 4)+1);
					//BDIreciever.ADx = BDIreciever.w1*((-1/pow(800,2))*(tempHousehold->getMaize()-param.householdNeed)+1);
					
					
					household->addMaize(household->getlackMaize(param.householdNeed));
					tempHousehold->removeMaize(household->getlackMaize(param.householdNeed));
					Shareflag = true;
					ClosenessTemp = tempHousehold->getCloseness(household->getId());
					ClosenessTemp += 0.05;
					tempHousehold->setCloseness(household->getId(), ClosenessTemp);
					Liveflag = true;
					return Liveflag;
					break;
				}
				else
				{
					tempHouseholdList.push_back(tempHousehold);
					addedMaize += tempHousehold->getLoanMaize(param.householdNeed);
					if(addedMaize >= household->getlackMaize(param.householdNeed))
					{
						ShareSucflag = true;
						addedMaize = 0;
						Liveflag = true;
						break;
					}
				
				}
			}
		}
	}
	if((ShareSucflag == true )&& (Shareflag == false))
	{
		for(std::vector<Household*>::iterator it = tempHouseholdList.begin() ; it != tempHouseholdList.end(); ++it)
		{
			Household* tempHousehold = (&**it);
			addedMaize += tempHousehold->getLoanMaize(param.householdNeed);

			ClosenessTemp = tempHousehold->getCloseness(household->getId());
			ClosenessTemp += 0.05;
			tempHousehold->setCloseness(household->getId(), ClosenessTemp);

			if(addedMaize < household->getlackMaize(param.householdNeed))
			{
				household->addMaize(tempHousehold->getLoanMaize(param.householdNeed));
				tempHousehold->removeMaize(tempHousehold->getLoanMaize(param.householdNeed));
			}
			else
			{
				household->addMaize(household->getlackMaize(param.householdNeed));
				tempHousehold->removeMaize(household->getlackMaize(param.householdNeed));
				Liveflag = true;
				return Liveflag;
				break;
			}
		}
	}
	return Liveflag;

}

bool AnasaziModel::MovewithFriends(std::vector<int> locgoal, Household* household)
{
	
	std::vector<Household*> householdList;
	std::vector<Location*> neighbouringLocations;
	double pfm = 0;
	repast::DoubleUniformGenerator* pfmTemp;
	repast::Moore2DGridQuery<Location> moore2DQuery(locationSpace);
	Location* tempLocSet;
	
	int range = 1;
	if(!locgoal.empty())
	{
		householdSpace->getObjectsAt(repast::Point<int>(locgoal[0], locgoal[1]), householdList);
	}
	for (std::vector<Household*>::iterator it = householdList.begin() ; it != householdList.end(); ++it)
	{
		Household* tempHousehold = (&**it);
		pfm = ((household->getCloseness(tempHousehold->getId())-(param.thresholdSharefood +0.1))/(1 - (param.thresholdSharefood +0.1)));
		pfmTemp = new repast::DoubleUniformGenerator(repast::Random::instance()->createUniDoubleGenerator(0,1)); 
		if(pfmTemp->next() >= pfm)
		{
			while(1)
			{
				moore2DQuery.query(locgoal, range, false, neighbouringLocations);
				for (std::vector<Location*>::iterator it = neighbouringLocations.begin() ; it != neighbouringLocations.end(); ++it)
				{
					Location* tempLoc = (&**it);
					if(tempLoc->getState() == 0)
					{
						if(tempLoc->getExpectedYield() >= param.householdNeed)
						{
							tempLocSet = tempLoc;
							goto EndOfLoop;
						}
					}
				}
				range++;
				if(range > boardSizeY)
				{
					return false;
				}
			}
			EndOfLoop:			
			if(range >=10)
			{
				return false;
			}
			else
			{
				std::cout << "movewithfriend" << std::endl;	
				tempHousehold->chooseField(tempLocSet);
				householdSpace->moveTo(tempHousehold->getId(), repast::Point<int>(locgoal[0], locgoal[1]));
				tempHousehold->nextYear(param.householdNeed);
				return true;
			}
		}	
	}
}

bool AnasaziModel::Moveout(std::vector<int> locgoal, Household* household)
{
	std::vector<Household*> householdList;
	std::vector<Location*> neighbouringLocations;
	double pfm = 0;
	repast::DoubleUniformGenerator* pfmTemp;
	repast::Moore2DGridQuery<Location> moore2DQuery(locationSpace);
	Location* tempLocSet;
	int range = 1;
	if(!locgoal.empty())
	{
		householdSpace->getObjectsAt(repast::Point<int>(locgoal[0], locgoal[1]), householdList);
	}
    if(householdList.empty())
	{
		return false;
	}
	for (std::vector<Household*>::iterator it = householdList.begin() ; it != householdList.end(); ++it)
	{
		Household* tempHousehold = (&**it);
		pfm = ((household->getCloseness(tempHousehold->getId())-(param.thresholdSharefood +0.1))/(1 - (param.thresholdSharefood +0.1)));
		pfmTemp = new repast::DoubleUniformGenerator(repast::Random::instance()->createUniDoubleGenerator(0,1)); 
		if(pfmTemp->next() >= pfm)
		{	
			if(tempHousehold->getMaize() < 1.2*param.householdNeed)
			{
				std::cout << "moveout" << std::endl;	
				removeHousehold(tempHousehold);
			}
		}	
	}
	return true;
}


void AnasaziModel::addNewAgentContacts(int agentId) {
    int id = agentId;
    int newSize = id;
    adjacencyMatrix.resize(newSize, std::vector<int>(newSize, 0));
    for (auto& row : adjacencyMatrix) {
        row.resize(newSize, 0);
    }
	int i;
	for(i=0;i<id;i++){
		if(rand()<Probability){
			setContact(id, i);
			setContact(i, id);
		}
	}
}

void AnasaziModel::setContact(int agentId1, int agentId2) {
    adjacencyMatrix[agentId1][agentId2] = 1;
}

void AnasaziModel::disContact(int agentId1, int agentId2) {
    adjacencyMatrix[agentId1][agentId2] = 0;
}

int AnasaziModel::getContactStatus(int agentId1, int agentId2){
    return adjacencyMatrix[agentId1][agentId2];
}


void AnasaziModel::initnetwork()
{
	//Initialise the adjacency matrix
	adjacencyMatrix.resize(Agent_Number, std::vector<int>(Agent_Number, 0));
	//Initialise the dynamic network using the probability
	int i;
	int j;
	for (i = 0; i < Agent_Number; i++){
		for (j = i + 1; j < Agent_Number; j++){
			if (rand() < Probability){
				adjacencyMatrix[i][j] = 1;
				adjacencyMatrix[j][i] = 1;
			}
		}
	}
}

void AnasaziModel::updateNetwork()
{
	std::vector<Household*> agents;
    std::vector<Household*>::iterator it3 = agents.begin();
	while(it3 != agents.end()){
		Household* household3;
        	household3 = (*it3);
		repast::AgentId id3 = household3->getId();
		int i = id3.id();
		int kmax = 0;
		std::vector<Household*> agents;
		std::vector<Household*>::iterator it2 = agents.begin();
		while(it2 != agents.end()){
			int sum1 = 0;
			int sum2 = 0;
			int sum3 = 0;
			int sum4 = 0;
			double Fmax = 0;
			double F2 = 0;
			Household* household2;
			household2 = (*it2);
			repast::AgentId id2 = household2->getId();
			int k = id2.id();
			std::vector<Household*> agents;
			std::vector<Household*>::iterator it1 = agents.begin();
			while(it1 != agents.end()){
				Household* household1;
				household1 = (*it1);
				repast::AgentId id1 = household1->getId();
				int j = id1.id();
				if(k!=j){
					sum1 += getContactStatus(j,i);
					sum2 += getContactStatus(j,i)*abs(static_cast<int>(household3->checkMaize(800)-household2->checkMaize(800)));
					std::vector<Household*> agents;
					std::vector<Household*>::iterator it4 = agents.begin();
					while(it4 != agents.end()){
						Household* household4;
						household4 = (*it4);
						repast::AgentId id4 = household4->getId();
						int h = id4.id();
						if(k!=h){
							sum4 += getContactStatus(h,i)*getContactStatus(j,h);
						}
						else{
							if(getContactStatus(h,i)==0){
								sum4 += getContactStatus(j,h);
							}
							else{
								sum4 += 0;
							}
						}
						it4++;
					}
					sum3 += getContactStatus(j,i)*sum4;
				}
				else {
					if(getContactStatus(j,i)==0){
						sum1 += 1;
						sum2 += abs(static_cast<int>(household3->checkMaize(800)-household2->checkMaize(800)));
						std::vector<Household*> agents;
						std::vector<Household*>::iterator it4 = agents.begin();
						while(it4 != agents.end()){
							Household* household4;
							household4 = (*it4);
							repast::AgentId id4 = household4->getId();
							int h = id4.id();
							if(k!=h){
								sum3 += getContactStatus(h,i)*getContactStatus(j,h);
							}
							else{
								if(getContactStatus(h,i)==0){
									sum3 += getContactStatus(j,h);
								}
								else{
									sum3 += 0;
								}
							}
							it4++;
						}
					}
					else{
						sum1 += 0;
						sum2 += 0;
						sum3 += 0;
					}
				}
				it1++;
			}
			F2 = param.b1*sum1+param.b3*sum2+param.b4*sum3;
			if(F2>=Fmax){
				kmax = k;
				Fmax = F2;
			}
			it2++;
		}
		if(kmax != i){
			if(getContactStatus(kmax,i)==0){
				setContact(kmax,i);
			}
			else{
				disContact(kmax,i);
			}
		}
		it3++;
	}
}

bool AnasaziModel::checkConnection(int id1, int id2)
{
	if(getContactStatus(id1,id2)==1&&getContactStatus(id2,id1)==1){
		return true;
	}
	else{
		return false;
	}
}
