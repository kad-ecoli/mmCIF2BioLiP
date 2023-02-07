#ifndef cssr_struct_H
#define cssr_struct_H 1
#include <cstring>
#include <set>
#include <map>
#include <algorithm>
#include "PDBParser.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 
using namespace std;

const float PI=3.141592653589793;
const float Extra=1.0e-4;
const float UpMax=1.0e+10;

inline float rad2deg(float rad);
inline void subtract(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);
inline void crossproduct(const vector<float> &c1, const vector<float> &c2, vector<float> &cc);
inline float innerproduct(const vector<float> &c1, const vector<float> &c2);
inline bool norm(const vector<float> &c, vector<float> &cc);
inline float Points2Distance(const vector<float> &c1, const vector<float> &c2);
inline float Points2Distance2(const vector<float> &c1, const vector<float> &c2);
inline float Points4Angle(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3, const vector<float> &c4);
inline float Points2Dihedral(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3, const vector<float> &c4);


bool check_balance_bracket(const vector<char>&dot_bracket, const size_t r1, const size_t r2);
size_t check_type_bracket(const vector<char>&dot_bracket, const size_t r1, const size_t r2);
void cssr_bp2dot(const vector<string>&res_str_vec, 
    const vector<size_t> filtered_bp_vec,
    const vector<pair<float,vector<string> > >&bp_vec, vector<char>&dot_bracket);
void cssr_bpseq(const vector<string>&res_str_vec,
    const vector<size_t>&filtered_bp_vec,
    const vector<pair<float,vector<string> > >&bp_vec,
    vector<vector<int> >&bpseq_vec);
void filter_bp(const vector<pair<float,vector<string> > >&bp_vec,
    vector<size_t> &filtered_bp_vec);
void removeUnspecifiedAtom(ModelUnit &pdb_entry, const string atom);
inline bool bp_tor_score(
    const vector<float> &c1, const vector<float> &c2,
    const vector<float> &c3, const vector<float> &c4,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator);
inline bool bp_ang_score(
    const vector<float> &c1, const vector<float> &c2,
    const vector<float> &c3, const vector<float> &c4,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator);
inline bool bp_len_score(const vector<float>&c1, const vector<float>&c2,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator);
inline bool bp_nn_score(const bool previnextj, const bool nextiprevj,
    const bool has_prev_ci,       const bool has_next_ci,
    const bool has_prev_cj,       const bool has_next_cj,
    const vector<float> &prev_ci, const vector<float> &next_ci,
    const vector<float> &prev_cj, const vector<float> &next_cj,
    const float mu, const float sd, const float tol, const float weight,
    float &nominator, float &denominator);
void cssr(const ModelUnit &pdb_entry, vector<string>&res_str_vec,
    vector<pair<float,vector<string> > >&bp_vec, const int fastOpt, 
    vector<pair<float,pair<int,int> > >&RR_list, const bool interchain);
#endif
