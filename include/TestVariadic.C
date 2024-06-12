#include <array>
#include <iostream>
#include <tuple>

using namespace std;

template<typename... Args> 
void mc_vector_(vector<int>& vec, Args&&... args) {
  (vec.emplace_back(std::forward<Args>(args)), ...); // fold-expression
}
template<typename... Args> 
vector<int> mc_vector_2(Args&&... args) {
  vector<int> vec;
  (vec.emplace_back(std::forward<Args>(args)), ...); // fold-expression
  return vec;
}
// template< typename... Args> 
// vector<int>  mc_vector_3(Args&&... args) {
  
//   auto func  = [](auto&&... param){
//     vector<int> vec;
//     (vec.emplace_back(std::forward<Args>(param)), ...);
//     return vec;};
  
//   return func(args);
// }
// template< typename... Args> 
// vector<int>  mc_vector_3(Args&&... args) {
  
//   auto func  = [](auto&&... param){
//     vector<int> vec;
//     (vec.emplace_back(std::forward<Args>(param)), ...);
//     return vec;};
  
//   return func(args);
// }
template<int N>
struct A {
    constexpr A() : arr() {
        for (auto i = 0; i != N; ++i)
            arr[i] = i; 
    }
  size_t arr[N];
};
   constexpr std::array<uint,4> arr = {0,1,2,3}; 

constexpr size_t number() {return 0;}
constexpr size_t numberi() {return arr[0];}

template <typename... Args>
vector<int> myMap(Args&&... args)
{
    constexpr std::size_t n = sizeof...(Args);
    constexpr auto a = A<n>();
    static_assert(n >= 1, "must pass at least one argument");
    std::cout << n << std::endl;
    auto tuple = std::tie(args...);
    vector<int> vec;
    uint i=0;
    //  std::apply([&vec](auto&&... args) {((std::cout << args << '\n'), ...);}, tuple);
    auto func = [](auto&&... args) {  vector<int> vec; (vec.emplace_back(std::forward<Args>(args)), ...); return vec;};
    vec = std::apply(func , tuple);
    // for (auto i = 0; i < n ; ++i){
    //   vec[i] = std::get<numberi()>(tuple);
    // }
    return vec;
}  
//template <typename ... Args>
void TestVariadic(){
  vector<int> myvec;
  mc_vector_(myvec,1,2,3,4,5,6);
  for(auto i:myvec) cout<<" "<<i;
  cout<<endl;
  auto myvec2=mc_vector_2(1,2,3,4,5,6);
  for(auto i:myvec2) cout<<" "<<i;
  cout<<endl;


  auto myvec3 = myMap(1,2,3);
  for(auto i:myvec3) cout<<" "<<i;
  cout<<endl;

}
