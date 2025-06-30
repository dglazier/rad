#pragma once

#include "ParticleCreator.h"
#include "ConfigReaction.h"

namespace rad{
  namespace generator{
    //using rad::names::data_type::Rec;
    //using rad::names::data_type::MC;
    using ROOT::RVecD;
    
    template<typename Tp, typename Tm>
	int ParticleCreateTwoBody(const int &id, const RVec<RVecD> components, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec <Tm> &m, const RVecI& iafter){
	
	auto idx = px.size();
	
	px.push_back(components[id][0]);
	py.push_back(components[id][1]);
	pz.push_back(components[id][2]);
	m.push_back(components[id][3]);
      
	return idx;
      }
    
    template<typename Tp, typename Tm>
      RVec<RVecD> CalculateTwoBody(const int &pidx, const RVec<Tm> &masses, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec <Tm> &m){
	
	//now to actually get the parent info
	//auto parent4 = FourVector(pidx);
	
	//for testing functions all play nice together
	//and actually return two particles in the end
	RVec<double> part1 = {0,0,0,0};
	RVec<double> part2 = {0,0,0,0};
	RVec<RVecD> result;
	result.push_back(part1);
	result.push_back(part2);
	return result;
      }
      
	  
      
    class ParticleGenerator : public rad::config::ParticleCreator{

    public:
      
      ParticleGenerator() = default;
    ParticleGenerator(rad::config::ConfigReaction &cr):rad::config::ParticleCreator{cr}{};
      
      template<typename Tp, typename Tm>
	void TwoBody(const std::vector<std::string> &names, const string &meson_parent, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec <Tm> &m){
	  
	int pidx=4;
	double emass=0.000511;
	RVec<double> masses={emass,emass};
	RVec<RVecD> two_body_moms = CalculateTwoBody(pidx,masses,px,py,pz,m);
	//cout << names[0] << " " << names[1] << endl;
	Reaction()->Define("two_body_first", [two_body_moms]() { return two_body_moms[0]; }, {});
	Reaction()->Define("two_body_second", [two_body_moms]() { return two_body_moms[1]; }, {});
	TString expr0 = Form("rad::generator::ParticleCreateTwoBody(%s,%s", "0", "two_body_first");
	TString expr1 = Form("rad::generator::ParticleCreateTwoBody(%s,%s", "1", "two_body_second");
	DefineParticle(names[0],std::vector<string>(),expr0.Data()); 
	DefineParticle(names[1],std::vector<string>(),expr1.Data()); 
	
      }
	  
    };// end Class ParticleGenerator

  }//end namespace generator
}//end namespace rad
