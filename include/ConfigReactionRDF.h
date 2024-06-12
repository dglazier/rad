#pragma once

//!  Base configuration class for electron scattering reactions

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
*/
#include "DefineNames.h"

#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

namespace rdf{

  //! Code simplifications
  
  using ROOT::RVecI;
  
  /*! RDFstep is used to update the dataframe handle after each filter or define step */
  using RDFstep = ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> ;

  /*! RVecIndexMap is gives access to specific particle indices in action functions.
    Use reaction part index getters e.g. ScatEleIdx() to get the index from the map.
    e.g. rdf::RVecIndexMap& react; ... ; auto index = react[ScatEleIdx()]; 
  */
  using RVecIndexMap = ROOT::RVec<ROOT::RVecI>;


  class ConfigReaction {


  public:

  ConfigReaction(const std::string_view treeName, const std::string_view fileNameGlob) : _orig_df{treeName,fileNameGlob} ,_curr_df{_orig_df}{

    }

    /** 
     * Get the current dataframe to add further actions
     */
    RDFstep CurrFrame(){return _curr_df;}
    /** 
     * Set the current dataframe after adding further actions
     */
    void setCurrFrame(RDFstep  df){ _curr_df = df ;}
    
    /** 
     * Add an alias for a branch and update the current frame to the aliased one
     */
     void setBranchAlias(const string_view name,const string_view alias){

      setCurrFrame(CurrFrame().Alias(alias,name));
    }
    
     /** 
     * Set constant index in collection for particle
     * This assumes constant position in collection (e.g in some HepMC3 files)
     * and update the current frame to the aliased one
     */
     void setParticleIndex(const string_view particle, const int idx){
       setCurrFrame(CurrFrame().Define(particle,[idx](){return idx;}));
     }
     
     /** 
     * Set function, func,  which defines variable index in collection for particle
     * The given function func is responsible for giving the correct position of particle in the collection
     * Names of required branches or identifiers given as vector of column names : columns
     * and update the current frame to the aliased one
     */
     template<typename Lambda>
       void setParticleIndex(const string_view particle, Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
       setCurrFrame(CurrFrame().Define(particle,func,columns));
     }

     /**
      * Make map that links particle names to indices in user functions
      * in C++ functions you can use the RVecIndexMap object indexed by 
      * name of the reaction component you need
      */
     void makeParticleMap(){
       //note, ordering in arguments, map and names must be maintained
       setCurrFrame(CurrFrame().Define(nameReactionMap(),
		     [](const int& beamion, const int& beamel,const int& scatel,
			const RVecI& baryons,const RVecI& mesons){
					 return RVecIndexMap{{beamion},{beamel},{scatel},baryons,mesons};},
				       {nameBeamIon().data(),nameBeamEle().data(),
					nameScatEle().data(),nameBaryons().data(),nameMesons().data()}));
       
     }


    /** 
     * Set constant index for beam electron, scattered electron and beam ion
     * This assumes constant position in collection (e.g in some HepMC3 files)
     * and update the current frame to the aliased one
     */
     void setBeamElectronIndex(const int idx){
       setParticleIndex(nameBeamEle(),idx);
     }
     void setScatElectronIndex(const int idx){
       setParticleIndex(nameScatEle(),idx);
     }
     void setBeamIonIndex(const int idx){
       setParticleIndex(nameBeamIon(),idx);
     }
     /**
      * Collect constant indices for final state mesons and baryons
      */
     void setMesonIndices(const RVecI& indices){
       setCurrFrame(CurrFrame().Define(nameMesons(),[indices](){return indices;}));
     } 
     void setBaryonIndices(RDFstep& df,const RVecI& indices){
       setCurrFrame(CurrFrame().Define(nameBaryons(),[indices](){return indices;}));
     }
     /**
      * Collect variable indices for final state mesons and baryons
      */
     void setMesonIndices(const  ROOT::RDF::ColumnNames_t& particles){
       setGroupParticles(nameMesons(),particles);
     }
     void setBaryonIndices(const  ROOT::RDF::ColumnNames_t& particles){
       setGroupParticles(nameMesons(),particles);
     }
     
     template <typename... Args>
       void setGroupParticles(const string_view name,const  ROOT::RDF::ColumnNames_t& particles){
       auto func = [](auto&&... args) {  RVecI vec; (vec.emplace_back(std::forward<Args>(args)), ...); return vec;};
        
     }
     /*
     void setGroupParticles(const string_view name,const ColumnNames_t &particles){
       //should be able to do this with variadic args, but don't know how to pass as argument to lambda
       switch( particles.size() ){
	   case 1 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0){return RVecI{p0};},particles));
	     break;
 	   case 2 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0,const int p1){return RVecI{p0,p1};},particles));
	     break;
 	   case 3 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0,const int p1,const int p2){return RVecI{p0,p1,p2};},particles));
	     break;
 	   case 4 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0,const int p1,const int p2,const int p3){return RVecI{p0,p1,p2,p3};},particles));
	     break;
 	   case 5 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4){return RVecI{p0,01,p2,p3,p4};},particles));
	     break;
 	   case 6 : 
	     setCurrFrame(CurrFrame().Define(name,[](const int p0,const int p1,const int p2,const int p3,const int p4,const int p5){return RVecI{p0,01,p2,p3,p4,p5};},particles));
	     break;
       default:
	 std::cerr<<"setGroupParticles only defined up to 6 particles! "<<std::endl;
	 exit(0);
	 break;
       }
	     
     } 
     */
     

  private :
    
     /**
     * _bare_df
     * Base dataframe constructed prior to adding actions
     * It should not be necessary to use this object for anything
     */
     ROOT::RDataFrame _orig_df;
   /** 
     * _curr_df 
     * Handle for the current dataframe state. 
     * This includes all added defines and filters so far.
     * Additional actions will be applied to this.
     */
    RDFstep _curr_df;
    
  };//ConfigReaction

}//rdf
