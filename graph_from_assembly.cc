#include "graph_from_assembly.h"
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
  int sc_size;

  Scaffold() {}

  Scaffold(const string& scf) {
    sc_size = scf.size();
    assert(scf[0] != 'N' && scf[0] != 'n');
    int gap_len = 0;
    string ctg_buf = "";
    for (int i = 0; i < scf.size(); i++) {
      if (scf[i] == 'n' || scf[i] == 'N') {
        if (gap_len == 0) {
          assert(ctg_buf.size() > 0);
          contigs.push_back(ctg_buf);
          ctg_buf = "";
        }
        gap_len++;
      } else {
        if (gap_len > 0) {
          assert(contigs.size() > 0);
          gaps.push_back(gap_len);
          gap_len = 0;
        }
        switch (scf[i]) {
          case 'A':
          case 'C':
          case 'G':
          case 'T': 
            ctg_buf += scf[i];
            break;
          case 'R':
          case 'M':
            ctg_buf += 'A';
            break;
          case 'Y':
          case 'S':
            ctg_buf += 'C';
            break;
          case 'K':
            ctg_buf += 'G';
            break;
          case 'W':
            ctg_buf += 'T';
            break;
        }
      }
    }
    assert(ctg_buf.size() > 0);
    contigs.push_back(ctg_buf);
    assert(gaps.size() + 1 == contigs.size());
//    printf("ct %d %d\n", contigs.size(), scf.size());
    contig_paths.resize(contigs.size());
    int sc_len = 0;
    int c_len = 0;
    for (auto &c: contigs) sc_len += c.size();
    c_len = sc_len;
    for (auto &e: gaps) sc_len += e;
    printf("%d/%d ", sc_len, c_len);
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

  void AddConAndCheck(int from, int to) {
    if (cons.size() <= from + 1) {
      cons.resize(from + 2);
    }
    for (auto &t: cons[from]) {
      if (t == to)
        return;
    }
    cons[from].push_back(to);
  }

  void AddCon(int from, int to) {
    AddConAndCheck(from, to);
    AddConAndCheck(to^1, from^1);
  }

  void AddBigCon(int from, int to) {
    if (big_cons.size() <= from + 1) {
      big_cons.resize(from + 2);
    }
    big_cons[from].push_back(to);
  }

};

void GetGraphFromAssembly(string filename, Graph& gr, vector<vector<int>>& paths) {
  int k = 101;
  ifstream input(filename.c_str());
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

  printf("sc len ");
  vector<Scaffold> scaffolds;
  for (auto &e: scfs) {
    scaffolds.push_back(Scaffold(e));
  }
  printf("\n");

/*  for (auto &s: scaffolds) {
    vector<int> sc_path;
    for (int i = 0; i < s.contigs.size(); i++) {
      int node_id = gr.nodes.size();
      Node *n1 = new Node;
      Node *n2 = new Node;
      n1->id = node_id;
      n2->id = node_id+1;
      n1->s = s.contigs[i];
      n2->s = ReverseSeq(s.contigs[i]);
      gr.nodes.push_back(n1);
      gr.nodes.push_back(n2);
      sc_path.push_back(node_id);
      if (i + 1 < s.contigs.size()) {
        sc_path.push_back(-s.gaps[i]);
      }
    }

    paths.push_back(sc_path);
  }
  printf("done\n");
  return;*/

  KmerDB kmerdb;
  unordered_set<int> end_markers;
  
  for (int si = 0; si < scaffolds.size(); si++) {
    for (int ci = 0; ci < scaffolds[si].contigs.size(); ci++) {
      auto &c = scaffolds[si].contigs[ci];
      int prev = -1;
      for (int i = 0; i + k <= c.size(); i++) {
        int id = kmerdb.Get(c.substr(i, k), Coord(si, ci, i));
        if (prev != -1) {
          kmerdb.AddCon(prev, id);
        }
        if (i == 0) {
          end_markers.insert(id);
          end_markers.insert(id^1);
        }
        if (i + k == c.size()) {
          end_markers.insert(id);
          end_markers.insert(id^1);
        }
        prev = id;
      }
    }
  }
  printf("kmerdb size %d\n", kmerdb.db.size());

  unordered_set<int> ignored;
  for (int i = 0; i < kmerdb.db.size(); i++) {
    if (kmerdb.cons[i].size() == 1 && end_markers.count(i) == 0) {
      int next = kmerdb.cons[i][0];
      if (next == (i^1)) {
        continue;
      }
      if (kmerdb.cons[next^1].size() == 1 && end_markers.count(next) == 0) {
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
  int total_int_size = 0;
  for (auto &inter: intervals) {
    auto inter_rev = InvertPath(inter.second);
    assert(intervals.count(inter_rev[0]));
    assert(intervals[inter_rev[0]] == inter_rev);
    total_int_size += inter.second.size();
  }
  printf("ti %d\n", total_int_size);

  unordered_map<int, int> renumber;

  for (auto &e: intervals) {
    if (renumber.count(e.second[0]) == 0) {
      assert(renumber.count(e.second.back()^1) == 0);
      assert(e.second[0] != (e.second.back()^1));
      int id = renumber.size();
      renumber[e.second[0]] = id;
      renumber[e.second.back()^1] = id+1;
    }
  }
  printf("renumber done\n");

  gr.nodes.resize(renumber.size());

  assert(gr.nodes.size() == intervals.size());
  for (auto &e: intervals) {
    Node *node = new Node;
    assert(e.first == e.second[0]);
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
    assert(gr.nodes[renumber[e.second[0]]] == NULL);
    gr.nodes[renumber[e.second[0]]] = node;
  }
  printf("gr done %d\n", gr.nodes.size());
  for (auto &s: scaffolds) {
    vector<int> path;
    for (int i = 0; i < s.contigs.size(); i++) {
      for (int j = 0; j < s.contig_paths[i].size(); j++) {
        assert(renumber.count(s.contig_paths[i][j]));
        assert(gr.nodes[renumber[s.contig_paths[i][j]]]->s.size() > 0);
        path.push_back(renumber[s.contig_paths[i][j]]);
      }
      if (i + 1 < s.contigs.size()) {
        path.push_back(-(s.gaps[i] + k - 1));
      }
    }
    paths.push_back(path);
  }
  printf("paths done %d\n", paths.size());
  int ind = 0;
  printf("pl ");
  for (auto &p: paths) {
    int len = 0;
    for (auto &e: p) {
      if (e < 0) len += -e;
      else len += gr.nodes[e]->s.size();
    }
    printf("%d(%d/%d) ", len, p.size(), scaffolds[ind].sc_size);
    ind++;
  }
  printf("\n");
}
