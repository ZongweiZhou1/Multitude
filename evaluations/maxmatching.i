%module maxmatching

%include "std_vector.i"
%include "std_map.i"

%{

#include "maxmatching.hpp"
%}


namespace std{
    %template (vect_int) vector<int>;
    %template (vect_double) vector<double>;
    %template (map_int2int) map<int, int>;
}


%include "maxmatching.hpp"