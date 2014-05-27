/* 
 * File:   main.cpp
 * Author: matthewsupernaw
 *
 * Created on May 20, 2014, 8:02 AM
 */

#include <cstdlib>

#include "ET4AD.hpp"
#include "CatchAtAge.hpp"
#include <string>
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
                    CatchAtAge<double> ca;
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

