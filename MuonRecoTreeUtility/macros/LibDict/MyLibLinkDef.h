#include <vector>
#include <map>
#include <utility>
#include <string>
//#include <TH1>
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<double>+;
//#pragma link C++ class std::pair<int,int >+;
#pragma link C++ class std::pair<int,std::vector<int> >+;
#pragma link C++ class std::pair<int,std::vector<double> >+;
#pragma link C++ class std::map<std::string,double>;
//#pragma link C++ class std::map<std::string,TH1*>+;
#pragma link C++ class std::map<int,std::vector<int> >+;
#pragma link C++ class std::map<int,int >+;
#pragma link C++ class std::map<int,std::vector<double> >+;
#pragma link C++ class std::map<int,std::vector<int> >::iterator+;
#pragma link C++ class std::map<int,std::vector<double> >::iterator+;
