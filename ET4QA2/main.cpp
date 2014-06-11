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

struct eqstr {

    bool operator()(const char* s1, const char* s2) const {
        if (s1 && s2)
            return strcmp(s1, s2) == 0;
        else
            return s1 == NULL && s2 == NULL;
    }
};

void lookup(const google::dense_hash_set<const char*, std::tr1::hash<const char*>, eqstr>& Set,
        const char* word) {
    google::dense_hash_set<const char*, std::tr1::hash<const char*>, eqstr>::const_iterator it
            = Set.find(word);
    cout << word << ": "
            << (it != Set.end() ? "present" : "not present")
            << endl;
}

int main() {
    google::dense_hash_set<uint32_t, std::tr1::hash<uint32_t> > Set;
    Set.set_empty_key(NULL);
    //  Set.insert(1);
    for (uint32_t i = 1; i < 1000000; i++) {
        Set.insert(i);
    }

    //  for(uint32_t i =1; i< 1000000; i+=10000){
    //  if(Set.find(i)!=Set.end()){
    //      std::cout<<i<<" was found\n";
    //  }
    //  }

    google::dense_hash_set<uint32_t, std::tr1::hash<uint32_t> > Set2;

    Set2.set_empty_key(NULL);
    Set2.insert(Set.begin(), Set.end());


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
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h> 
#include <sstream>


#include <curl/curl.h>

std::stringstream datastream;

size_t writeCallback(char* buf, size_t size, size_t nmemb, void* up) { //callback must have this declaration
    //buf is a pointer to the data that curl has for us
    //size*nmemb is the size of the buffer
    //        std::cout<<
    for (int c = 0; c < size * nmemb; c++) {
        //            sdata.push_back(buf[c]);
        datastream << buf[c];
    }
    return size*nmemb; //tell curl how many bytes we handled
}

int GetData(const std::string &name) {
    datastream.str("");

    CURL *curl;
    CURLcode res;

    curl_global_init(CURL_GLOBAL_DEFAULT);

    curl = curl_easy_init();

    std::stringstream ss;
    //
    ss << "http://ichart.finance.yahoo.com/table.csv?s=" << name << "&ignore=.csv";
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, ss.str().c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &writeCallback);

#ifdef SKIP_PEER_VERIFICATION
        /*
         * If you want to connect to a site who isn't using a certificate that is
         * signed by one of the certs in the CA bundle you have, you can skip the
         * verification of the server's certificate. This makes the connection
         * A LOT LESS SECURE.
         *
         * If you have a CA cert for the server stored someplace else than in the
         * default bundle, then the CURLOPT_CAPATH option might come handy for
         * you.
         */
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
#endif

#ifdef SKIP_HOSTNAME_VERIFICATION
        /*
         * If the site you're connecting to uses a different host name that what
         * they have mentioned in their server certificate's commonName (or
         * subjectAltName) fields, libcurl will refuse to connect. You can skip
         * this check, but this will make the connection less secure.
         */
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);
#endif

        /* Perform the request, res will get the return code */
        res = curl_easy_perform(curl);
        /* Check for errors */
        if (res != CURLE_OK)
            fprintf(stderr, "curl_easy_perform() failed: %s\n",
                curl_easy_strerror(res));

        /* always cleanup */
        curl_easy_cleanup(curl);
    }

    curl_global_cleanup();



}

int testcurl() {
    CURL *curl;
    CURLcode res;

    curl_global_init(CURL_GLOBAL_DEFAULT);

    curl = curl_easy_init();
    if (curl) {
        curl_easy_setopt(curl, CURLOPT_URL, "http://ichart.finance.yahoo.com/table.csv?s=aapl&ignore=.csv");

#ifdef SKIP_PEER_VERIFICATION
        /*
         * If you want to connect to a site who isn't using a certificate that is
         * signed by one of the certs in the CA bundle you have, you can skip the
         * verification of the server's certificate. This makes the connection
         * A LOT LESS SECURE.
         *
         * If you have a CA cert for the server stored someplace else than in the
         * default bundle, then the CURLOPT_CAPATH option might come handy for
         * you.
         */
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
#endif

#ifdef SKIP_HOSTNAME_VERIFICATION
        /*
         * If the site you're connecting to uses a different host name that what
         * they have mentioned in their server certificate's commonName (or
         * subjectAltName) fields, libcurl will refuse to connect. You can skip
         * this check, but this will make the connection less secure.
         */
        curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);
#endif

        /* Perform the request, res will get the return code */
        res = curl_easy_perform(curl);
        /* Check for errors */
        if (res != CURLE_OK)
            fprintf(stderr, "curl_easy_perform() failed: %s\n",
                curl_easy_strerror(res));

        /* always cleanup */
        curl_easy_cleanup(curl);
    }

    curl_global_cleanup();



}

void testurl() {
    std::stringstream ss;
    //    ss << "POST " << "" << " HTTP/1.0\r\n";
    ss << "Host: " << "http://ichart.finance.yahoo.com/" << "\r\n";
    ss << "GET: " << "/table.csv?s=MSFT&a=0&b=1&c=2010&d=0&e=1&f=2013&g=w&q=q&y=0&x=.csv\n\n";
    //https://ichart.finance.yahoo.com/table.csv?s=%5EGSPC&d=5&e=2&f=2014&g=d&a=0&b=3&c=1950&ignore=.csv
    int status;
    struct addrinfo host_info; // The struct that getaddrinfo() fills up with data.
    struct addrinfo *host_info_list; // Pointer to the to the linked list of host_info's.

    // The MAN page of getaddrinfo() states "All  the other fields in the structure pointed
    // to by hints must contain either 0 or a null pointer, as appropriate." When a struct 
    // is created in C++, it will be given a block of memory. This memory is not necessary
    // empty. Therefor we use the memset function to make sure all fields are NULL.     
    memset(&host_info, 0, sizeof host_info);

    std::cout << "Setting up the structs..." << std::endl;

    host_info.ai_family = AF_UNSPEC; // IP version not specified. Can be both.
    host_info.ai_socktype = SOCK_STREAM; // Use SOCK_STREAM for TCP or SOCK_DGRAM for UDP.

    // Now fill up the linked list of host_info structs with google's address information.
    status = getaddrinfo("http://ichart.finance.yahoo.com/", "80", &host_info, &host_info_list);
    // getaddrinfo returns 0 on succes, or some other value when an error occured.
    // (translated into human readable text by the gai_gai_strerror function).
    if (status != 0) std::cout << "getaddrinfo error" << gai_strerror(status);

    std::cout << "Creating a socket..." << std::endl;
    int socketfd; // The socket descripter
    socketfd = socket(host_info_list->ai_family, host_info_list->ai_socktype,
            host_info_list->ai_protocol);
    if (socketfd == -1) std::cout << "socket error ";

    std::cout << "Connect()ing..." << std::endl;
    status = connect(socketfd, host_info_list->ai_addr, host_info_list->ai_addrlen);
    if (status == -1) std::cout << "connect error";

    std::cout << "send()ing message..." << std::endl;
    char *msg = "GET / HTTP/1.1\nhost: /table.csv?s=aapl&ignore=.csv\n\n";
    int len;
    ssize_t bytes_sent;
    len = strlen(msg);
    bytes_sent = send(socketfd, ss.str().c_str(), ss.str().size(), 0);


    std::cout << "Waiting to recieve data..." << std::endl;
    ssize_t bytes_recieved;
    char incoming_data_buffer[1000];
    bytes_recieved = recv(socketfd, incoming_data_buffer, 1000, 0);
    // If no data arrives, the program will just wait here until some data arrives.
    if (bytes_recieved == 0) std::cout << "host shut down." << std::endl;
    if (bytes_recieved == -1)std::cout << "recieve error!" << std::endl;
    std::cout << bytes_recieved << " bytes recieved :" << std::endl;
    std::cout << incoming_data_buffer << std::endl;

    std::cout << "Receiving complete. Closing socket..." << std::endl;
    freeaddrinfo(host_info_list);
    close(socketfd);
}

template<class T>
class Simple : public ad::FunctionMinimizer<T> {
public:
    unsigned int numberOfObservations;

    ad::Variable<T> m;
    ad::Variable<T> b;

    std::vector<T> x;
    std::vector<T> y;

    std::vector<ad::Variable<T> > predictedY;

    void Initialize() {
        this->numberOfObservations = 1000003;
        this->x = std::vector<double>(this->numberOfObservations);
        this->y = std::vector<double>(this->numberOfObservations);
        //        this->predictedY = std::vector<ad::Variable<T> >(this->numberOfObservations);
        T init = T(.0001);
        for (int i = 0; i < this->numberOfObservations; i++) {
            T r = ((double) rand() / (RAND_MAX));
            T rv = -15.0 + (30.0 * r);
            init += .0001;
            x[i] = init;
            y[i] = T(((1.9) * 4.8 * x[i] + (rv)));
        }
        //        slope = 1.9;
        //        intercept = 
        m.SetName("m");
        b.SetName("b");
        this->Register(m);
        this->Register(b);

    }

    void ObjectiveFunction(ad::Variable<T>& f) {
        f = 0.0;
        ad::Variable<T> temp;
        ad::Variable<T> sum = 0;
        for (int i = 0; i< this->y.size(); i++) {
            //            this->predictedY[i] = this->slope * x[i] + this->intercept;
            temp = this->m * x[i] + this->b;
            sum += (temp - y[i])*(temp - y[i]);

        }

        //        for (int i = 0; i< this->predictedY.size(); i++) {
        //            sum += (this->predictedY[i] - y[i])*(this->predictedY[i] - y[i]);
        //        }

        f = sum; //this->numberOfObservations / 2.0 * std::log(sum);


    }

    void Finalize() {
        //        std::cout << "x,observed y,predicted y" << std::endl;
        //        for (int i = 0; i < this->numberOfObservations; i++) {
        //                        std::cout << x[i] << "," << y[i] << "," << (slope * x[i] + intercept) << std::endl;
        //        }
    }

};

#include "Portfolio.hpp"
#include "Statistics.hpp"
#include "Regression.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    //    Simple<double> s;
    //    s.Run();
    //    exit(0);
    int SIZE = 10000;
    //    typedef ad::BigFloat<double> real;
    typedef double real;
    std::vector<real > x(SIZE);
    std::vector<real > y(SIZE);
    ////    double init = double(.0001);
    //    for (int i = 0; i < SIZE; i++) {
    //      
    //        x[i] = i+1;
    //        
    //        y[i] = x[i]*x[i];
    //        
    ////        std::cout<<y[i]<<"\n";
    //    }


    real init = real(.0001);
    for (int i = 0; i < SIZE; i++) {
        real r = ((double) rand() / (RAND_MAX));
        real rv = -0.025+ (.05 * r);
        init += .0001;
        x[i] = init;
        y[i] = (22.1234 * std::pow(x[i], 2.0) + 2.00*x[i] + 1.234567)+rv;
         //3.041 real(((1.9) * 4.8 * x[i] + (rv)));
    }
    std::cout << ad::Mean(y) << "\n";
    ad::LogrithmicRegression<real> lm(x, y);
    lm.SetVerbose(true);
    lm.Run();
    std::cout << lm.ToString() << "\n";

            for (int i = 0; i < SIZE; i++) {
                std::cout<<y[i]<<","<<lm.Evaluate(x[i])<<"\n";
            }
        sleep(1);
    //    std::cout<<lm.Evaluate(1000003);
    exit(0);

    std::setprecision(50);
    double data[] = {6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49};
    double data2[] = {-3, -7, -15, -4, -39, -40 / 2.0, -8, -42, -22, -100, -49};
    std::vector<ad::BigFloat<double> > xv;
    std::vector<ad::BigFloat<double> > yv;
    for (int i = 0; i < 11; i++) {
        xv.push_back(data[i]);
        yv.push_back(data2[i]);
    }
    //    xv.push_back(1.0);
    // 
    //       xv.push_back(4);
    //       xv.push_back(2);
    //       
    //    std::copy(&data, &data + 11, xv.begin());
    std::cout << ad::CorrelationCoefficient<ad::BigFloat<double> >(xv, yv).ToString() << " == " << (12.0 / 7.0);

    exit(0);

    int i = 0;
    std::ifstream in;
    in.open("sp500");
    std::vector<Asset<float> > assets;
    while (in.good()) {
        std::string line;
        std::getline(in, line);
        std::cout << "fetching data for " << line << "...";
        GetData(line);
        if (datastream.str().at(0) == 'D') {//if it starts with Date...otherwise theres an error!
            assets.push_back(Asset<float>(line, datastream.str()));
        } else {
            std::cout << "error fetching\n";
        }
    }



    //    std::cout<<datastream.str()<<"\n";
    //testcurl();

    exit(0);

    typedef ad::Variable<double> var;
    typedef ad::Variable<double, 1 > var1;
    //    var::SetRecording(false);
    var x11 = 10.0;
    var1 x22 = 10.0;
    var z11 = x11 + x22;

    std::cout << var::IsRecording() << " " << var1::IsRecording() << "  " << z11;
    //      exit(0);

    //                    ad::Variable<double>::SetSupportArbitraryOrder(true);
    //    //                while (1) {
    CatchAtAge<double> ca;
    //    ca.SetTolerance(1e-4); 
    //    //                    ca.SetPrintInterval(1);
    ca.Run();
    std::stringstream serialized;
    ca.GetFunctionValue().Serialize(serialized);
    var g = ca.GetFunctionValue().Deserialize(serialized);

    std::cout << "g = " << g << "\n";
    std::cout.precision(50);
    std::cout << std::scientific << "df/deffort_devs[575] = " << ca.GetFunctionValue().WRT(ca.effort_devs[575]) << "\n";
    std::cout << std::scientific << "dg/deffort_devs[575] = " << g.WRT(ca.effort_devs[575]) << "\n";

    //    ////            }
    //    //               std::cout<<ad::Variable<double>::misses_g<<" map misses!";
    exit(0);
    //
    //    //    //
    //    //    //
    //    var x(3.1459, true);
    //    var y(3.1459, true);
    //    var z(0, true);
    //    var f = x + y;
    //    f.SetName("f");
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << " exp size = " << f.Size() << "\n";
    //    //    
    //    f = y + x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = x + 10.0;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = 10.0 + x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //    std::stringstream ss;
    //
    //    f.Serialize(ss);
    //
    //    var sv = f.Deserialize(ss);
    //
    //
    //
    //
    //
    //    f = x * y;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = y * x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = x * 10.0;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = 10.0 * x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = x / y;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = y / x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = x / 10.0;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f = 10.0 / x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f += y;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //    f += x;
    //    std::cout << "f = " << f << ", df/dx = " << " == " << f.WRT(x) << "\n";
    //
    //
    //    //    for (int i = 0; i < 100000000; i++) {
    //    //        f += f;
    //    //        f += x;
    //    //        f += x;
    //    //        f += x;
    //    //    }
    //
    //
    //    //    f = std::log(f);
    //    //    std::cout << "f = " << f << ", df/dx = "  <<" == "<<f.WRT(x)<< "\n";
    //
    //
    //    //    f = f;//std::exp(f);
    //    //    var ff =f;
    //    //    std::cout << "f = " << f << ", df/dx = " << ff.Diff(x) <<" == "<<ff.WRT(x)<< "\n";
    //    //    //
    //    //    //    for(int i =0; i < 10000; i++){
    //    //    //        
    //    //    //    }
    //    //
    //
    //    //    


    return 0;
}

#endif