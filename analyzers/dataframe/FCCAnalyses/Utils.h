#ifndef  UTILS_ANALYZERS_H
#define  UTILS_ANALYZERS_H

#include <cmath>
#include <algorithm>
#include "ROOT/RVec.hxx"

namespace FCCAnalyses {
  namespace Utils {

    template<typename T> inline auto getsize( T& vec){ return vec.size();};
    template<typename T> inline ROOT::VecOps::RVec<ROOT::VecOps::RVec<T> >  as_vector(const ROOT::VecOps::RVec<T>& in){return ROOT::VecOps::RVec<ROOT::VecOps::RVec<T> >(1, in);};     
    template <typename T> inline ROOT::VecOps::RVec<int> index_range(const ROOT::VecOps::RVec<T> & in){
      ROOT::VecOps::RVec<int> indices(in.size()); 
      std::iota(indices.begin(),indices.end(), 0);
      return indices; 
    };  
    /// @brief count the number of valid (>=0, < size of collection) indices in an index list 
    /// @param in: index list
    /// @param ref: particle vector to which the indices refer. 
    /// @return integer count of valid indices 
    template <typename T> inline int count_valid_indices(const ROOT::VecOps::RVec<int> & in,
                                                         const ROOT::VecOps::RVec<T> & ref){
      int maxSize = ref.size(); 
      return std::count_if(in.begin(),in.end(),[maxSize](const int & i){return (i >=0 && i < maxSize);}); 
    };

    // @brief for a given list of indices, returns a list of particle (copies) 
    /// @param idx : Indices of desired particles within the full set ("in")
    /// @param in : Full set of particles 
    /// @return A vector of particles with the desired indices. 
    template <typename T> inline ROOT::VecOps::RVec<T> sel_byIndex( const ROOT::VecOps::RVec<int> & idx, const ROOT::VecOps::RVec<T> & in){
      ROOT::VecOps::RVec<T> found; 
      for (int index : idx){
        if (index < 0 || index >= in.size()) continue; 
        found.push_back(in.at(index)); 
      }
      return found; 
    }
    // @brief merge (concatenate) two collections of arbitrary content 
    /// @param x : first collection - entries will be copied in-order
    /// @param y : second collection - entries will be copied in-order after the last element of the first
    /// @return A combined collection of size (x.size()+y.size()), containing the content of x followed by that of y 
    template <typename T> inline ROOT::VecOps::RVec<T> merge( const ROOT::VecOps::RVec<T> & x, const ROOT::VecOps::RVec<T> & y){
      ROOT::VecOps::RVec<T> merged;
      merged.reserve(x.size()+y.size()); 
      merged.insert(merged.end(), x.begin(), x.end());  
      merged.insert(merged.end(), y.begin(), y.end());  
      return merged; 
    }

    template <typename outType, typename inType> inline ROOT::VecOps::RVec<outType> convertVec ( const ROOT::RVec<inType> & in, 
                 std::function<outType(const inType&)> convertor){
      ROOT::VecOps::RVec<outType> out;
      out.reserve(in.size());
      for (const inType & input : in){
        out.push_back(std::move(convertor(input))); 
      }
      return out; 
    }

    /// @brief Helper struct to select entries matching a certain predicate. 
    /// Supports two signatures - either a list of candidates is passed and a list of accepted candidates returned,
    /// Or a list of indices in a vector of candidates is passed and a list of accepted indices returned. 
    /// The latter is more compatible with index-based selection logic.  
    template <class content> struct selByPredicate{
      selByPredicate(std::function<bool(const content &)> thePredicate):m_predicate(thePredicate){}
      std::function<bool(const content &)> m_predicate; 
      ROOT::VecOps::RVec<content>  operator() (const ROOT::VecOps::RVec<content> & in){
          ROOT::VecOps::RVec<content> result;
          result.reserve(in.size());
          for (auto & p : in) {
            if (m_predicate(p)) result.emplace_back(p);
          } 
          return result; 
      }
      ROOT::VecOps::RVec<int> operator() (const ROOT::VecOps::RVec<int> & indices, const ROOT::VecOps::RVec<content> & in){
        ROOT::VecOps::RVec<int> result;
        result.reserve(in.size());
        for (int index : indices) {
          if (index < 0 || index >= in.size()) continue; 
          if (m_predicate(in[index])) result.emplace_back(index);
        } 
        return result; 
      }
      ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> operator() (
        const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> & setsOfIndices, 
        const ROOT::VecOps::RVec<content> & in){
        ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result(setsOfIndices.size());
        for (int elem = 0; elem <  setsOfIndices.size(); ++elem){
          result[elem] = this->operator()(setsOfIndices[elem], in); 
        }
        return result; 
      }
    };
    }
}

#endif
