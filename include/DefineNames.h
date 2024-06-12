#pragma once


//!  Helper functions for names of reaction parts

namespace rad{
  namespace names{

  /**
   * Names used to identify reaction components
   */
  constexpr const std::string_view  ReactionMap()  {return "reaction_map"; }
  constexpr const std::string_view  Mesons()  {return "meson"; }
  constexpr const std::string_view  Baryons() {return "baryon";}
  constexpr const std::string_view  ScatEle() {return "scat_ele";}
  constexpr const std::string_view  BeamEle() {return "beam_ele";}
  constexpr const std::string_view  BeamIon() {return "beam_ion";}
  constexpr const std::string_view  Beams() {return "{beam_ele,beam_ion}";}

  /**
   * Links to reaction component links for users c++ functions
   */
  enum class PartGroup{ BeamIon,BeamEle,ScatEle,Baryons,Mesons}; //ordering below must match this
  constexpr uint  BeamIonIdx() {return static_cast<uint>(PartGroup::BeamIon);}
  constexpr uint  BeamEleIdx() {return static_cast<uint>(PartGroup::BeamEle);}
  constexpr uint  ScatEleIdx() {return static_cast<uint>(PartGroup::ScatEle);}
  constexpr uint  BaryonsIdx() {return static_cast<uint>(PartGroup::Baryons);}
  constexpr uint  MesonsIdx() {return static_cast<uint>(PartGroup::Mesons);}

  }//names
}//rad
