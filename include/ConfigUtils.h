#pragma once


namespace rad{
  namespace config{

   //using config::RVecIndexMap;
 
  void PrintDefinedColumnNames(RDFstep  df){
      std::cout<<"Print Column Names : ";
      auto cols =  df.GetDefinedColumnNames();
      for(auto& col:cols){
	std::cout<<col<<", ";
      }
      cout<<"\n";
    }
   void PrintAllColumnNames(RDFstep  df){
      std::cout<<"Print Column Names : ";
      auto cols =  df.GetColumnNames();
      for(auto& col:cols){
	std::cout<<col<<", ";
      }
      cout<<"\n";
    }

    bool ColumnExists(const string& col,RDFstep  df){
      auto cols =  df.GetDefinedColumnNames();
      if(std::find(cols.begin(),cols.end(),col)==cols.end()) return false;
      return true;
    }
    
    
  }
}
