#include "chamber.hpp"
#include <fstream>
#include <string>
using namespace std;

chamber::chamber(std::ifstream& chamberFile)
{
string temp;

chamberFile>>temp>>plateDistance>>temp;
chamberFile>>temp>>diameter>>temp;
chamberFile>>temp>>Nw>>temp;
chamberFile>>temp>>kQ>>temp;
}

