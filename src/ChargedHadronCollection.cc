#include "ChargedHadronCollection.hh"

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
// using namespace lcio ;
// using namespace marlin ;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::Vertex;
using IMPL::LCRelationImpl;
using IMPL::ReconstructedParticleImpl;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using std::string;
using std::vector;
using UTIL::LCRelationNavigator;
using HiddenVAnalysis::MathOperator;

ChargedHadronCollection aChargedHadronCollection;

ChargedHadronCollection::ChargedHadronCollection() : Processor("ChargedHadronCollection")
{

	// modify processor description
	_description = "";

	registerProcessorParameter("Charged_or_All",
							"use charged or all particles",
							_charged_bool,
							_charged_bool);
	registerProcessorParameter("momentum_cut",
							"a cut on the momentum",
							_momentum_cut,
							_momentum_cut);
	registerProcessorParameter("costheta_cut",
							"a cut on the costheta",
							_costheta_cut,
							_costheta_cut);
	// input collections
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
							"PFOCollection",
							"PFO collection name",
							_colName,
							string("PandoraPFOs"));
	registerInputCollection(LCIO::MCPARTICLE,
							"MCCollectionName",
							"Name of the MC collection",
							_MCColName,
							std::string("MCParticlesSkimmed"));
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
							"MCOutputName",
							"MC stable particles",
							_newMCColName,
							std::string("StableParticles"));
}

ChargedHadronCollection::~ChargedHadronCollection() {}

void ChargedHadronCollection::init()
{

	printParameters();
}

ReconstructedParticle *ChargedHadronCollection::getRecoParticle(MCParticle *particle)
{

	float momentum[3];
	momentum[0] = particle->getMomentum()[0];
	momentum[1] = particle->getMomentum()[1];
	momentum[2] = particle->getMomentum()[2];
	float mass = particle->getMass();
	float energy = particle->getEnergy();
	float charge = particle->getCharge();

	ReconstructedParticleImpl *result = new ReconstructedParticleImpl();
	result->setCharge(charge);
	result->setEnergy(energy);
	result->setMomentum(momentum);
	return result;
}

std::vector<MCParticle *> ChargedHadronCollection::GetQQbarStables(EVENT::LCCollection *myCollection)
{
	// std::cout << "##################################" << std::endl;
	// std::cout << "##         Stable level         ##" << std::endl;
	// std::cout << "##################################" << std::endl;

	int number = myCollection->getNumberOfElements();
	vector<MCParticle *> stables;
	if (number < 3)
		return stables; // return empty

	vector<MCParticle *> stable_stables;
	vector<MCParticle *> isr_stables;
	vector<MCParticle *> overlay_stables;

	int istable; // the border of stableization

	int isrphotons[2];
	isrphotons[0] = -1;
	isrphotons[1] = -1;

	vector<MCParticle *> pairISR;
	for (int i = 2; i < number; i++)
	{
		MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
		if (particle->getPDG() == 22 && pairISR.size() < 2 && particle->getGeneratorStatus() == 1)
		{
			pairISR.push_back(particle);
			isrphotons[pairISR.size() - 1] = i;
		}
	}
	for (int i = isrphotons[pairISR.size() - 1] + 1; i < number; i++)
	{

		MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
		vector<MCParticle *> daughters = particle->getDaughters();

		streamlog_out(DEBUG) << "\n MCCollection, particle:" << i;
		streamlog_out(DEBUG) << " pdg=" << particle->getPDG();
		streamlog_out(DEBUG) << " satus=" << particle->getGeneratorStatus();
		streamlog_out(DEBUG) << " Ndaughters=" << daughters.size();
		streamlog_out(DEBUG) << " E=" << particle->getEnergy();
		streamlog_out(DEBUG) << " px=" << particle->getMomentum()[0];
		streamlog_out(DEBUG) << " py=" << particle->getMomentum()[1];
		streamlog_out(DEBUG) << " pz=" << particle->getMomentum()[2];
		streamlog_out(DEBUG) << " m=" << particle->getMass();
		streamlog_out(DEBUG) << " charge=" << particle->getCharge();
		if (daughters.size() == 0 && particle->isOverlay() == true)
		{
			overlay_stables.push_back(particle);
			stables.push_back(particle);
			streamlog_out(DEBUG) << " ----> IS OVERLAY";
		}

		if (daughters.size() == 0 && particle->isOverlay() == false)
		{

			int ISRstable = 0;

			vector<MCParticle *> parents = particle->getParents();
			for (int j = 0; j < parents.size(); j++)
			{
				if (parents.at(j) == pairISR.at(0) || parents.at(j) == pairISR.at(1))
				{
					ISRstable = 1;
				}
				else
				{
					vector<MCParticle *> parents2 = parents.at(j)->getParents();
					for (int j2 = 0; j2 < parents2.size(); j2++)
					{
						if (parents2.at(j2) == pairISR.at(0) || parents2.at(j2) == pairISR.at(1))
						{
							ISRstable = 1;
						}
						else
						{
							vector<MCParticle *> parents3 = parents2.at(j2)->getParents();
							for (int j3 = 0; j3 < parents3.size(); j3++)
							{
								if (parents3.at(j3) == pairISR.at(0) || parents3.at(j3) == pairISR.at(1))
								{
									ISRstable = 1;
								}
							}
						}
					}
				}
			}

			if (ISRstable == 0)
			{
				stable_stables.push_back(particle);
				stables.push_back(particle);
				streamlog_out(DEBUG) << " ----> IS STABLE";
			}

			if (ISRstable == 1)
			{
				isr_stables.push_back(particle);
				stables.push_back(particle);
				streamlog_out(DEBUG) << " ----> IS ISR";
			}
		}
	}

	/*stables.push_back(stable_stables);
	stables.push_back(isr_stables);
	stables.push_back(overlay_stables);*/
	return stables;
} // GetQQbarStables()

void ChargedHadronCollection::processRunHeader(LCRunHeader *run)
{
}

void ChargedHadronCollection::processEvent(LCEvent *evt)
{

	try
	{
		LCCollection *mccol = evt->getCollection(_MCColName);
		LCCollection *pfocol = evt->getCollection(_colName);

		vector<MCParticle *> stable_particles = GetQQbarStables(mccol);

		IMPL::LCCollectionVec *newMCCol= new IMPL::LCCollectionVec(EVENT::LCIO::RECONSTRUCTEDPARTICLE);

		for (unsigned int i = 0; i < stable_particles.size(); i++)
		{
			ReconstructedParticle *recoparticle = getRecoParticle(stable_particles.at(i));
			vector<float> direction = MathOperator::getDirection(recoparticle->getMomentum());
			float angle = fabs(MathOperator::getAngles(direction)[1]);
			float costheta = std::cos(angle);
			float mom=MathOperator::getModule(recoparticle->getMomentum());
			if(fabs(costheta)>_costheta_cut) continue;
			if(mom>_momentum_cut) continue;
			if (stable_particles.at(i)->getCharge() != 0)
			{
				newMCCol->addElement(recoparticle);
			} else {
				if(_charged_bool==0) newMCCol->addElement(recoparticle);
			}
		}
		
		evt->addCollection(newMCCol, _newMCColName);
	
	}catch (DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "Whoops!....\n";
			streamlog_out(DEBUG) << e.what();
		}
	}

	void ChargedHadronCollection::check(LCEvent * evt)
	{
		// nothing to check here - could be used to fill checkplots in reconstruction processor
	}

	void ChargedHadronCollection::end()
	{
	}
