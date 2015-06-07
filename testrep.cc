#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cassert>
#include "graph.h"
#include "utility.h"

using namespace std;

struct Scaffold {
  vector<string> contigs;
  vector<vector<int>> contig_paths;
  vector<int> gaps;

  Scaffold() {}

  Scaffold(const string& scf) {
    assert(scf[0] != 'N' && scf[0] != 'n');
    int gap_len = 0;
    string ctg_buf = "";
    for (int i = 0; i < scf.size(); i++) {
      if (scf[i] == 'n' || scf[i] == 'N') {
        if (gap_len == 0) {
          assert(ctg_buf.size() > 0);
          contigs.push_back(ctg_buf);
/*          if (ctg_buf.size() > 100000)
            printf("c %d\n", ctg_buf.size());*/
          ctg_buf = "";
        }
        gap_len++;
      } else {
        if (gap_len > 0) {
          assert(contigs.size() > 0);
          gaps.push_back(gap_len);
//          printf("g %d\n", gap_len);
          gap_len = 0;
        }
        ctg_buf += scf[i];
      }
    }
    assert(ctg_buf.size() > 0);
    contigs.push_back(ctg_buf);
    assert(gaps.size() + 1 == contigs.size());
    contig_paths.resize(contigs.size());
//    printf("e\n");
  }
};

struct Coord {
  int scf, ctg, pos;
  Coord() {}
  Coord(int s, int c, int p) : scf(s), ctg(c), pos(p) {}
};

struct KmerDB {
  unordered_map<string, int> db;
  unordered_map<int, Coord> coords;
  vector<vector<int>> cons;
  vector<vector<int>> big_cons;
  int Get(const string& x, Coord c) {
    if (db.count(x) == 0) {
      int id = db.size();
      db[x] = id;
      coords[id] = c;
      db[ReverseSeq(x)] = id+1;
    }
    return db[x];
  }

  int Get(const string& x) {
    assert(db.count(x));
    return db[x];
  }

  void AddCon(int from, int to) {
    if (cons.size() <= from + 1) {
      cons.resize(from + 2);
    }
    if (cons.size() <= to + 1) {
      cons.resize(to + 2);
    }
    cons[from].push_back(to);
    cons[to^1].push_back(from^1);
  }

  void AddBigCon(int from, int to) {
    if (big_cons.size() <= from + 1) {
      big_cons.resize(from + 2);
    }
    big_cons[from].push_back(to);
  }

};

int main(int argc, char** argv) {
  int k = 101;
  ifstream input(argv[1]);
  vector<string> scfs;
  string buf, l;
  while (getline(input, l)) {
    if (l[0] == '>') {
      if (buf.length() > 0) {
        scfs.push_back(buf);
      }
      buf = "";
    } else {
      buf += l;
    }
  }
  if (buf.length() > 0) {
    scfs.push_back(buf);
  }

  vector<Scaffold> scaffolds;
  for (auto &e: scfs) {
    scaffolds.push_back(Scaffold(e));
  }

  unordered_map<size_t, vector<Coord>> poses;
  int k2 = 500;
  for (int i = 0; i < scaffolds.size(); i++) {
    for (int j = 0; j < scaffolds[i].contigs.size(); j++) {
      for (int p = 0; p + k2 <= scaffolds[i].contigs[j].size(); p++) {
        size_t h = hash<string>()(scaffolds[i].contigs[j].substr(p, k2));
        poses[h].push_back(Coord(i, j, p));
      }
    }
  }
  printf("table done\n");
  for (auto &e: poses) {
    if (e.second.size() > 1) {
      printf("hash hit\n");
      if (scaffolds[e.second[0].scf].contigs[e.second[0].ctg].substr(e.second[0].pos, k2) ==
          scaffolds[e.second[1].scf].contigs[e.second[1].ctg].substr(e.second[1].pos, k2)) {
        printf("normal hit (%d %d %d) (%d %d %d) %s %s\n", e.second[0].scf, e.second[0].ctg,
            e.second[0].pos, e.second[1].scf, e.second[1].ctg, e.second[1].pos,
           scaffolds[e.second[0].scf].contigs[e.second[0].ctg].substr(e.second[0].pos, 10).c_str(),
           scaffolds[e.second[1].scf].contigs[e.second[1].ctg].substr(e.second[1].pos, 10).c_str() );
      }
    }
  }

  return 0;

  KmerDB kmerdb;
  
  for (int si = 0; si < scaffolds.size(); si++) {
    for (int ci = 0; ci < scaffolds[si].contigs.size(); ci++) {
      auto &c = scaffolds[si].contigs[ci];
      int prev = -1;
      for (int i = 0; i + k <= c.size(); i++) {
        int id = kmerdb.Get(c.substr(i, k), Coord(si, ci, i));
        if (prev != -1) {
          kmerdb.AddCon(prev, id);
        }
        prev = id;
      }
    }
  }
  printf("%d\n", kmerdb.db.size());

  unordered_set<int> ignored;
  for (int i = 0; i < kmerdb.db.size(); i++) {
    if (kmerdb.cons[i].size() == 1) {
      int next = kmerdb.cons[i][0];
      if (kmerdb.cons[next^1].size() == 1) {
        ignored.insert(next);
      }
    }
  }
  printf("k %d\n", ignored.size());

  unordered_map<int, vector<int>> intervals;
  for (auto &s: scaffolds) {
    int contig_id = 0;
    for (auto &c: s.contigs) {
      vector<int> cur_int;
      for (int i = 0; i + k <= c.size(); i++) {
        int id = kmerdb.Get(c.substr(i, k));
        if (ignored.count(id)) {
          if (cur_int.size() > 0) {
            cur_int.push_back(id);
          }
        } else {
          if (cur_int.size() > 0) {
//            printf("add %d %d\n", cur_int[0], cur_int.back());
            if (intervals.count(cur_int[0]) == 0 || cur_int.size() > intervals[cur_int[0]].size()) {
              intervals[cur_int[0]] = cur_int;
            }
            kmerdb.AddBigCon(cur_int[0], id);
            s.contig_paths[contig_id].push_back(cur_int[0]);
          }
          cur_int.clear();
          cur_int.push_back(id);
        }
      }
      if (cur_int.size() > 0) {
//        printf("add %d %d\n", cur_int[0], cur_int.back());
        if (intervals.count(cur_int[0]) == 0 || cur_int.size() > intervals[cur_int[0]].size()) {
          intervals[cur_int[0]] = cur_int;
        }
      }
      cur_int.clear();
      auto cr = ReverseSeq(c);
      for (int i = 0; i + k <= cr.size(); i++) {
        int id = kmerdb.Get(cr.substr(i, k));
        if (ignored.count(id)) {
          if (cur_int.size() > 0) {
            cur_int.push_back(id);
          }
        } else {
          if (cur_int.size() > 0) {
//            printf("addr %d %d\n", cur_int[0], cur_int.back());
            if (intervals.count(cur_int[0]) == 0 || cur_int.size() > intervals[cur_int[0]].size()) {
              intervals[cur_int[0]] = cur_int;
            }
            kmerdb.AddBigCon(cur_int[0], id);
          }
          cur_int.clear();
          cur_int.push_back(id);
        }
      }
      if (cur_int.size() > 0) {
//        printf("addr %d %d\n", cur_int[0], cur_int.back());
        if (intervals.count(cur_int[0]) == 0 || cur_int.size() > intervals[cur_int[0]].size()) {
          intervals[cur_int[0]] = cur_int;
        }
      }
      contig_id++;
    }
  }
  printf("i %d\n", intervals.size());
  for (auto &inter: intervals) {
    auto inter_rev = InvertPath(inter.second);
    assert(intervals.count(inter_rev[0]));
    assert(intervals[inter_rev[0]] == inter_rev);
  }

  unordered_map<int, int> renumber;

  for (auto &e: intervals) {
    if (renumber.count(e.second[0]) == 0) {
      int id = renumber.size();
      renumber[e.second[0]] = id;
      renumber[e.second.back()^1] = id+1;
    }
  }
  printf("renumber done\n");

  Graph gr;
  gr.nodes.resize(renumber.size());
  
  for (auto &e: intervals) {
    Node *node = new Node;
    node->id = renumber[e.second[0]];
//    node->inv = node->id^1;
    node->s = "";
    for (auto &x: e.second) {
      if (x % 2 == 0) {
        auto &c = kmerdb.coords[x];
        node->s += scaffolds[c.scf].contigs[c.ctg][c.pos+k-1];
      } else {
        auto &c = kmerdb.coords[x^1];
        node->s += ReverseBase(scaffolds[c.scf].contigs[c.ctg][c.pos]); 
      }
    }
    gr.nodes[renumber[e.second[0]]] = node;
  }
  printf("gr done %d\n", gr.nodes.size());
  vector<vector<int>> paths;
  for (auto &s: scaffolds) {
    vector<int> path;
    for (int i = 0; i < s.contigs.size(); i++) {
      for (int j = 0; j < s.contig_paths[i].size(); j++) {
        path.push_back(renumber[s.contig_paths[i][j]]);
      }
      if (i + 1 < s.contigs.size()) {
        path.push_back(-s.gaps[i]);
      }
    }
    paths.push_back(path);
  }
  printf("paths done %d\n", paths.size());
}
