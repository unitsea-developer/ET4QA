/* 
 * File:   main.cpp
 * Author: matthewsupernaw
 *
 * Created on May 20, 2014, 8:02 AM
 */
//#define GOOGLE

#ifdef GOOGLE
#include <iostream>
#include <google/dense_hash_set>
using namespace std;
//using google;      // namespace where class lives by default

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
      if (s1 && s2)
          return strcmp(s1, s2) == 0;
      else
          return s1 == NULL && s2 == NULL;
  }
};


void lookup(const google::dense_hash_set<const char*, std::tr1::hash<const char*>, eqstr>& Set,
            const char* word)
{
  google::dense_hash_set<const char*,  std::tr1::hash<const char*>, eqstr>::const_iterator it
    = Set.find(word);
  cout << word << ": "
       << (it != Set.end() ? "present" : "not present")
       << endl;
}

int main()
{
  google::dense_hash_set<uint32_t,  std::tr1::hash<uint32_t> > Set;
  Set.set_empty_key(NULL);
//  Set.insert(1);
  for(uint32_t i =1; i< 1000000; i++){
      Set.insert(i);
  }
  
//  for(uint32_t i =1; i< 1000000; i+=10000){
//  if(Set.find(i)!=Set.end()){
//      std::cout<<i<<" was found\n";
//  }
//  }
  
  google::dense_hash_set<uint32_t,  std::tr1::hash<uint32_t> > Set2;
  
  Set2.set_empty_key(NULL);
  Set2.insert(Set.begin(),Set.end());
  
  
//  Set.insert("plum");
//  Set.insert("apple");
//  Set.insert("mango");
//  Set.insert("apricot");
//  Set.insert("banana");

//  lookup(Set, "mango");
//  lookup(Set, "apple");
//  lookup(Set, "durian");
}

#else

#include <cstdlib>

#include "ET4AD.hpp"
#include "CatchAtAge.hpp"
#include <string>
#include "BigFloat.hpp"
#include "support/cache-table-0.2/mm/cache_map.hpp"
using namespace std;

template<char name >
class TestBase {
};

class V : public TestBase < 'y' > {
};

class IV : public TestBase < 'y' > {
};

/*
 * 
 */
int main(int argc, char** argv) {




//                ad::Variable<double>::SetSupportArbitraryOrder(true);
//                while (1) {
                    CatchAtAge<double > ca;
                    ca.SetPrintInterval(5);
                    ca.Run();
                    ////            }
    //               std::cout<<ad::Variable<double>::misses_g<<" map misses!";
                    exit(0);
    typedef ad::Variable<double> var;
    //    //
    //    //
    var x(3.1459, true);
    var y = 10.0;
    var z(0, true);
    var f = x + y;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << " exp size = " << f.Size() << "\n";
    //    
    f = y + x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = x + 10.0;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = 10.0 + x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";


    f = x * y;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = y * x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = x * 10.0;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = 10.0 * x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = x / y;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = y / x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = x / 10.0;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f = 10.0 / x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f += y;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";

    f += x;
    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";


//    for (int i = 0; i < 100000000; i++) {
//        f += f;
//        f += x;
//        f += x;
//        f += x;
//    }


    //    f = std::log(f);
    //    std::cout << "f = " << f << ", df/dx = "  <<" == "<<f.WRT(x)<< "\n";


    //    f = f;//std::exp(f);
    //    var ff =f;
    //    std::cout << "f = " << f << ", df/dx = " << ff.Diff(x) <<" == "<<ff.WRT(x)<< "\n";
    //    //
    //    //    for(int i =0; i < 10000; i++){
    //    //        
    //    //    }
    //

    //    


    return 0;
}

#endif