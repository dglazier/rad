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
  constexpr const std::string_view  TargetIon() {return "tar_ion";}
  constexpr const std::string_view  BeamGamma() {return "beam_gamma";}
    //constexpr const std::string_view  Beams() {return "{beam_ele,beam_ion}";}

  /**
   * Links to reaction component links for users c++ functions
   */
    enum class InitGroup{Bot,Top}; //ordering below must match this
    enum class ElectroGroup{ BeamIon,BeamEle,ScatEle=4}; //ordering below must match this
    enum class PhotoGroup{ TarIon,BeamGam}; //ordering below must match this
    enum class FinalGroup{ Baryons=2,Mesons}; //ordering below must match this

    
    constexpr uint  InitialTopIdx() {return static_cast<uint>(InitGroup::Top);}
    constexpr uint  InitialBotIdx() {return static_cast<uint>(InitGroup::Bot);}

    constexpr uint  ElectroIonIdx() {return InitialBotIdx();}
    constexpr uint  ElectroEleIdx() {return InitialTopIdx();}
    
    constexpr uint  ScatEleIdx() {return static_cast<uint>(ElectroGroup::ScatEle);}
    
    constexpr uint  BaryonsIdx() {return static_cast<uint>(FinalGroup::Baryons);}
    constexpr uint  MesonsIdx() {return static_cast<uint>(FinalGroup::Mesons);}
    
    constexpr uint  PhotoGammaIdx() {return InitialTopIdx();}
    constexpr uint  PhotoIonIdx() {return InitialBotIdx();}

  }//names
}//rad
