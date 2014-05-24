#ifndef PROB_CALCULATOR_H__
#define PROB_CALCULATOR_H__

#include "graph.h"

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
  }


  double CalcProb(vector<vector<int> >& paths,
                  vector<pair<int, int>>& zeros,
                  int& total_len) {
    double prob = 0;
    int zero;
    for (auto &e: single_reads) {
      prob += CalcScoreForPaths(
          gr, paths, *e.second, zero, total_len,
          true, e.first.penalty_constant, e.first.step,
          e.first.min_prob_per_base, e.first.min_prob_start) * e.first.weight;
      zeros.push_back(make_pair(zero, e.second->GetNumberOfReads()));
    }
    for (auto &e: paired_reads) {
      prob += CalcScoreForPaths(
          gr, paths, *e.second.first, *e.second.second, 
          e.first.insert_mean, e.first.insert_std, zero,
          total_len, true, e.first.penalty_constant,
          e.first.step, true,
          e.first.min_prob_per_base, e.first.min_prob_start) * e.first.weight;
      zeros.push_back(make_pair(zero, e.second.first->GetNumberOfReads()));
    }
    for (auto &e: pacbio_reads) {
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
  Graph& gr;

/*  double CalcProb(vector<vector<int> >& paths,
                  int& zero1, int& zero2, int& total_len) {
    return CalcScoreForPacbio(
        gr_, paths, kmer_, rs11_,
        zero1, total_len, use_caching_, no_cov_1_,
        exp_cov_move_1_)*0.01 +
        CalcScoreForPaths(
        gr_, paths, kmer_, rs21_, rs22_,
        insert_mean2_, insert_std2_, zero2, total_len, use_caching_, no_cov_2_,
        exp_cov_move_2_, true)+
        CalcScoreForPaths(
        gr_, paths, kmer_, rs31_, rs32_,
        insert_mean3_, insert_std3_, zero3, total_len, use_caching_, 0,
        500, true);
  }
  double CalcProb(vector<vector<int> >& paths,
                  int& zero1, int& zero2, int& zero3, int& total_len) {
    return CalcScoreForPacbio(
        gr_, paths, kmer_, rs11_,
        zero1, total_len, use_caching_, no_cov_1_,
        exp_cov_move_1_)*0.01 +
        CalcScoreForPaths(
        gr_, paths, kmer_, rs21_, rs22_,
        insert_mean2_, insert_std2_, zero2, total_len, use_caching_, no_cov_2_,
        exp_cov_move_2_, true)+
        CalcScoreForPaths(
        gr_, paths, kmer_, rs31_, rs32_,
        insert_mean3_, insert_std3_, zero3, total_len, use_caching_, 0,
        500, true);
  }*/
};


#endif
