%module dynmeans

%include "std_vector.i"

%{
#include "dynmeans.hpp"
%}

namespace std{
    %template(dV) vector<double>;
    %template(iV) vector<int>;
}


%include "dynmeans.hpp"
