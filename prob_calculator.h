#ifndef PROB_CALCULATOR_H__
#define PROB_CALCULATOR_H__

#include "graph.h"
#include "utility.h"

struct SingleReadConfig {
  SingleReadConfig() {}
  SingleReadConfig(double pc, double s, double mp, double mps, double w, bool a) :
      penalty_constant(pc), step(s),
      min_prob_per_base(mp), min_prob_start(mps), weight(w), advice(a) {}
  double penalty_constant;
  double step;
  double min_prob_per_base;
  double min_prob_start;
  double weight;
  bool advice;
};

struct PairedReadConfig {
  PairedReadConfig() {}
  PairedReadConfig(double pc, double s, double im, double is, double mp,
                   double mps, double w, bool a) :
      penalty_constant(pc), step(s), insert_mean(im), insert_std(is),
      min_prob_per_base(mp), min_prob_start(mps), weight(w), advice(a) {}
  
  double penalty_constant;
  double step;
  double insert_mean;
  double insert_std;
  double min_prob_per_base;
  double min_prob_start;
  double weight;
  bool advice;
};

class ProbCalculator {
 public:
  ProbCalculator(
      const vector<pair<SingleReadConfig, ReadSet*>>& single_reads,
      const vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>>& paired_reads,
      const vector<pair<SingleReadConfig, PacbioReadSet*>>& pacbio_reads,
      Graph& gr) :
        single_reads(single_reads), paired_reads(paired_reads),
        pacbio_reads(pacbio_reads), gr(gr) {
    paired_scoring_states.resize(paired_reads.size());
  }

  vector<vector<int>> NormalizePaths(vector<vector<int>>& paths) {
    vector<vector<int>> ret;
    for (auto &p: paths) {
      vector<int> r = InvertPath(p);
      if (r < p) {
        ret.push_back(r);
      } else {
        ret.push_back(p);
      }
    }
    return ret;
  }


  double CalcProb(vector<vector<int>>& pathso,
                  vector<pair<int, int>>& zeros,
                  int& total_len) {
//    vector<vector<int>> paths = NormalizePaths(pathso);
    vector<vector<int>> paths=pathso;
    zeros.clear();
    double prob = 0;
    for (auto &e: single_reads) {
      int zero = 0;
      prob += CalcScoreForPaths(
          gr, paths, *e.second, zero, total_len,
          true, e.first.penalty_constant, e.first.step,
          e.first.min_prob_per_base, e.first.min_prob_start) * e.first.weight;
      zeros.push_back(make_pair(zero, e.second->GetNumberOfReads()));
    }
    int ind = 0;
    for (auto &e: paired_reads) {
/*      double score_slow = CalcScoreForPaths(
          gr, paths, *e.second.first, *e.second.second, 
          e.first.insert_mean, e.first.insert_std, zero,
          total_len, true, e.first.penalty_constant,
          e.first.step, true,
          e.first.min_prob_per_base, e.first.min_prob_start) * e.first.weight;
      int zero2, t2;*/
      int zero = 0;
      double score_fast = CalcScoreForPathsNew(
          gr, paths, *e.second.first, *e.second.second,
          e.first.insert_mean, e.first.insert_std,
          zero, total_len, paired_scoring_states[ind],
          true, e.first.penalty_constant,
          e.first.step, true, e.first.min_prob_per_base,
          e.first.min_prob_start) * e.first.weight;
//      printf("cmp %lf %lf\n", score_slow, score_fast);
      zeros.push_back(make_pair(zero, e.second.first->GetNumberOfReads()));
      prob += score_fast;
      ind++;
    }
    for (auto &e: pacbio_reads) {
      int zero = 0;
      prob += CalcScoreForPacbio(
          gr, paths, *e.second, zero, total_len, true,
          e.first.penalty_constant, e.first.step,
          e.first.min_prob_per_base, e.first.min_prob_start) * e.first.weight;
      zeros.push_back(make_pair(zero, e.second->GetNumberOfReads()));
    }
    return prob;
  }
  double CalcProb(vector<vector<int> >& paths,
                  int& total_len) {
    vector<pair<int, int>> zeros;
    return CalcProb(paths, zeros, total_len);
  }
  double CalcProb(vector<vector<int> >& paths) {
    int tl;
    return CalcProb(paths, tl);
  }
  vector<pair<SingleReadConfig, ReadSet*>> single_reads;
  vector<pair<PairedReadConfig, pair<ReadSet*, ReadSet*>>> paired_reads;
  vector<pair<SingleReadConfig, PacbioReadSet*>> pacbio_reads;
  vector<ScoringState> paired_scoring_states;
  Graph& gr;
};


#endif
