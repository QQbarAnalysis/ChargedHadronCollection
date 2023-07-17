#ifndef ChargedHadronCollection_h
#define ChargedHadronCollection_h 1
#include <iomanip>
#include <EVENT/LCRelation.h>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <cmath>
#include <string>
#include <vector>

#include <TF1.h>
#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include "MathOperator.hh"
// #include <TLorentzVector.h>
using namespace lcio;
using namespace marlin;

class ChargedHadronCollection : public Processor
{

public:
  virtual Processor *newProcessor() { return new ChargedHadronCollection; }

  ChargedHadronCollection();
  virtual ~ChargedHadronCollection();

  /** Called at the begin of the job before anything is read.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  virtual void end();

private:
  EVENT::ReconstructedParticle *getRecoParticle(MCParticle *particle);
  std::vector<EVENT::MCParticle *> GetQQbarStables(EVENT::LCCollection *myCollection);

  int _onlycharged_bool=0;
  float _momentum_cut=1;
  float _pt_cut=0.5;
  float _costheta_cut=0.99;
  std::string _newMCColName;
  std::string _MCColName;
  std::string _colName;
};

#endif
