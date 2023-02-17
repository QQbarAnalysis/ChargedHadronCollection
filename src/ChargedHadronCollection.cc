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

	for (int i = 0 i < number; i++)
	{

		MCParticle *particle = dynamic_cast<MCParticle *>(myCollection->getElementAt(i));
		vector<MCParticle *> daughters = particle->getDaughters();
		int status=particle->getGeneratorStatus();

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
		if (status==1 && particle->isOverlay() == true)
		{
			overlay_stables.push_back(particle);
			//stables.push_back(particle);
			streamlog_out(DEBUG) << " ----> IS OVERLAY";
		}

		if ( status==1 && particle->isOverlay() == false)
		{
			stables.push_back(particle);	
		}
	}


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
