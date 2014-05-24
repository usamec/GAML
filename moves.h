#ifndef MOVES_H__
#define MOVES_H__

#include "graph.h"
#include "prob_calculator.h"

bool BreakPath(vector<vector<int> >&new_paths, Graph& gr, int threshold);
bool LocalChange2(vector<vector<int> >& new_paths, Graph& gr, int threshold,
                  int path_id, int ps, int pt, ProbCalculator& prob_calc);
bool FixMultiLocal(vector<vector<int> >& new_paths, Graph& gr, int threshold);
bool FixRep(vector<vector<int> >& new_paths, Graph& gr, int threshold);
bool FixSelfLoops(vector<vector<int> >& new_paths, Graph& gr, int threshold);
bool LocalChange(vector<vector<int> >& new_paths, Graph& gr, int threshold, int &path_id,
    int &xx, int &yy, ProbCalculator& prob_calc);
void ReversePath(vector<int>& path);
bool ExtendPathsAlt(vector<vector<int> >& paths, Graph& gr, int threshold);
bool ExtendPaths(vector<vector<int> >& new_paths, Graph& gr, int threshold,
    ProbCalculator& prob_calc);
int SamplePathByLength(vector<vector<int> >& paths, Graph& gr);
bool FixGapLength(vector<vector<int> >& paths, int path_id, int gap_pos,
                  ProbCalculator& prob_calc, int prev_len);
bool ExtendPathsAdv(vector<vector<int> >& paths, Graph&gr, int threshold,
                    PacbioReadSet& rs, int kmer, ProbCalculator& prob_calc);
bool ExtendPathsAdv(vector<vector<int> >& paths, Graph&gr, int threshold,
                    ReadSet& rs1, ReadSet& rs2, int kmer, ProbCalculator& prob_calc);
bool FixGapLength(vector<vector<int> >& paths, ProbCalculator& prob_calc);
bool SplitOnNode(int node, vector<vector<int>>& paths);
void FixRepForNode2(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                   int node, ProbCalculator& prob_calc);
bool FixBigReps(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                ProbCalculator& prob_calc);
bool FixSomeBigReps(vector<vector<int>>& paths, Graph&gr, int threshold, bool disjoin_similar,
                ProbCalculator& prob_calc);
bool FixRepForNode(int node, vector<vector<int>>& paths, int threshold, Graph& gr,
                   ProbCalculator& prob_calc);
#endif
