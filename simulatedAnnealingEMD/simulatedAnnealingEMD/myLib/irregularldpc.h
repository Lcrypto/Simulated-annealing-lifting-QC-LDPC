/*
Copyright(c) 2012, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#pragma once
#include"..\myLib\regularLDPC.h"
#include"..\myLib\CycleEnum.h"


bool girthAtLeast6Manual(const vector<vector<vector<int> > >& a, ll circulant);
void print(const vector<vector<vector<int> > >& a);
void eprint(const vector<vector<vector<int> > >& a);
bool findBalancedCycle(int g, const vector<vector<int> >& protograph, int circulant);
bool isPossible(ll targetGirth, const vector<vector<int> >& protograph);
bool noCycles(int g, const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant);
int getGirth(const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant);
ll getRand(ll mod);
void readME(vector<vector<vector<int> > >& mtr, istream& in);


void readME(vector<vector<vector<int> > >& mtr, istream& in = cin) {
    int m = mtr.size();
    int n = mtr[0].size();
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            string s;
            in >> s;
            int cur = 0;
            for (int k = 0; k < s.size(); ++k) {
                if (s[k] == '&') {
                    mtr[i][j].push_back(cur);
                    cur = 0;
                    continue;
                }
                if ((s[k] >= '0') && (s[k] <= '9')) {
                    cur = cur * 10 + s[k] - '0';
                    continue;
                }
                break;
            }
            mtr[i][j].push_back(cur);
        }
    }
}

ll getRand(ll mod) {
    while (true) {
        ll q = RAND_MAX + 1;
        ll x = rand();
        while (q < mod) {
            q *= (RAND_MAX + 1);
            x = x * (RAND_MAX + 1) + rand();
        }
        if (x <= q - 1 - (q % mod))
            return x % mod;
    }

}


bool findBalancedCycle(int g, const vector<vector<int> >& protograph, ll circulant) {
    CycleEnum cycle(g, protograph);
    if (!cycle.init()) {
        return 0;
    }
    while (true) {
        bool balanced = 1;
        map<pair<pii, int>, int> m;
        for (int i = 0; i < g; i += 2) {
            ++m[make_pair(pii(cycle.cycle[i].r, cycle.cycle[i].c), cycle.cycle[i].id)];
        }
        for (int i = 1; i < g; i += 2) {
            --m[make_pair(pii(cycle.cycle[i].r, cycle.cycle[i].c), cycle.cycle[i].id)];
            if (m[make_pair(pii(cycle.cycle[i].r, cycle.cycle[i].c), cycle.cycle[i].id)] < 0) {
                balanced = 0;
                break;
            }
        }
        if (balanced) {
            cerr << "Can't achieve girth more than " << g << " because of the balanced cycle in protograph:\n";
            for (int i = 0; i < g; ++i) {
                cerr << "(" << cycle.cycle[i].r << ", " << cycle.cycle[i].c << ", " << cycle.cycle[i].id << ")->";
            }
            cerr << endl;
            return 1;
        }
        if (!cycle.next())
            return 0;
    }
    return 0;
}


//cycle counted multiple times
int getCycles(int g, const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant) {
    CycleEnum cycle(g, protograph);
    if (!cycle.init())
        return 0;

    int count = 0;
    while (true) {
        ll res = 0;
        for (int i = 0; i * 2 < g; ++i) {
            res = res + circulant + a[cycle.cycle[2 * i].r][cycle.cycle[2 * i].c][cycle.cycle[2 * i].id] - a[cycle.cycle[2 * i + 1].r][cycle.cycle[2 * i + 1].c][cycle.cycle[2 * i + 1].id];
        }
        if (res % circulant == 0) {
            /*for (int i = 0; i < 4; ++i) {
            cout << cycle.cycle[i].r << " " << cycle.cycle[i].c << " " << cycle.cycle[i].id << endl;
            }*/
            ++count;

        }
        if (!cycle.next())
            return count;
    }

}

bool noCycles(int g, const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant) {
    return (getCycles(g, a, protograph, circulant) == 0);
}

int getGirth(const vector<vector<vector<int> > >& a, const vector<vector<int> >& protograph, ll circulant) {
    for (int g = 4; g < 1000; g += 2) {
        if (!noCycles(g, a, protograph, circulant))
            return g;
    }
    cerr << "girth >= 1000\n";
    return (1 << 31) - 1;
}



bool isPossible(ll targetGirth, const vector<vector<int> >& protograph, ll circulant) {
    for (int g = 6; g < targetGirth; g += 2) {
        if (findBalancedCycle(g, protograph, circulant)) {
            return 0;
        }
    }
    return 1;
}

bool girthAtLeast6Manual(const vector<vector<vector<int> > >& a, ll circulant) {
    //cycles inside one entry
    int checkNodes = a.size();
    int variableNodes = a[0].size();
    for (int i = 0; i < checkNodes; ++i) {
        for (int j = 0; j < variableNodes; ++j) {
            for (size_t i1 = 0; i1 < a[i][j].size(); ++i1) {
                for (size_t i2 = i1 + 1; i2 < a[i][j].size(); ++i2) {
                    for (size_t i3 = i1; i3 < a[i][j].size(); ++i3) {
                        if (i3 == i2)
                            continue;
                        for (size_t i4 = i1 + 1; i4 < a[i][j].size(); ++i4) {
                            if (i3 == i4)
                                continue;
                            if ((a[i][j][i1] + a[i][j][i3] - a[i][j][i2] - a[i][j][i4] + 2 * circulant) % circulant == 0)
                                return 0;
                        }
                    }
                }
            }
        }
    }
    //cycles between two entries in one row
    for (int i = 0; i < checkNodes; ++i) {
        for (int i1 = 0; i1 < variableNodes; ++i1) {
            if (a[i][i1].size() < 2)
                continue;
            set<ll> dif;
            for (size_t i11 = 0; i11 < a[i][i1].size(); ++i11) {
                for (size_t i12 = i11 + 1; i12 < a[i][i1].size(); ++i12) {
                    dif.insert(a[i][i1][i12] - a[i][i1][i11]);
                    dif.insert(circulant + a[i][i1][i11] - a[i][i1][i12]);
                }
            }
            for (int i2 = i1 + 1; i2 < variableNodes; ++i2) {
                if (a[i][i2].size() < 2)
                    continue;
                for (size_t i21 = 0; i21 < a[i][i2].size(); ++i21) {
                    for (size_t i22 = i21 + 1; i22 < a[i][i2].size(); ++i22) {
                        if (dif.find(a[i][i2][i22] - a[i][i2][i21]) != dif.end())
                            return 0;
                    }
                }
            }
        }
    }
    //cycles between two entries in one column
    for (int i = 0; i < variableNodes; ++i) {
        for (int i1 = 0; i1 < checkNodes; ++i1) {
            if (a[i1][i].size() < 2)
                continue;
            set<ll> dif;
            for (size_t i11 = 0; i11 < a[i1][i].size(); ++i11) {
                for (size_t i12 = i11 + 1; i12 < a[i1][i].size(); ++i12) {
                    dif.insert(a[i1][i][i12] - a[i1][i][i11]);
                    dif.insert(circulant + a[i1][i][i11] - a[i1][i][i12]);
                }
            }
            for (int i2 = i1 + 1; i2 < checkNodes; ++i2) {
                if (a[i2][i].size() < 2)
                    continue;
                for (size_t i21 = 0; i21 < a[i2][i].size(); ++i21) {
                    for (size_t i22 = i21 + 1; i22 < a[i2][i].size(); ++i22) {
                        if (dif.find(a[i2][i][i22] - a[i2][i][i21]) != dif.end())
                            return 0;
                    }
                }
            }
        }
    }
    //cycles between four entries
    for (int r1 = 0; r1 < checkNodes; ++r1) {
        for (int c1 = 0; c1 < variableNodes; ++c1) {
            if (a[r1][c1].empty())
                continue;
            for (int r2 = r1 + 1; r2 < checkNodes; ++r2) {
                if (a[r2][c1].empty())
                    continue;
                for (int c2 = c1 + 1; c2 < variableNodes; ++c2) {
                    if ((a[r1][c2].empty()) || (a[r2][c2].empty()))
                        continue;
                    for (size_t i11 = 0; i11 < a[r1][c1].size(); ++i11) {
                        for (size_t i12 = 0; i12 < a[r1][c2].size(); ++i12) {
                            for (size_t i21 = 0; i21 < a[r2][c1].size(); ++i21) {
                                for (size_t i22 = 0; i22 < a[r2][c2].size(); ++i22) {
                                    if ((a[r1][c1][i11] + a[r2][c2][i22] - a[r1][c2][i12] - a[r2][c1][i21] + 2 * circulant) % circulant == 0)
                                        return 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 1;
}

void print(const vector<vector<vector<int> > >& a) {
    if (a.empty())
        return;
    int CHECK_NODES = a.size();
    int VARIABLE_NODES = a[0].size();
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            for (size_t k = 0; k + 1 < a[i][j].size(); ++k)
                cout << a[i][j][k] << "&";
            if (a[i][j].empty())
                cout << -1 << "\t";
            else
                cout << a[i][j][a[i][j].size() - 1] << "\t";
        }
        cout << endl;
    }
}

void eprint(const vector<vector<vector<int> > >& a) {
    if (a.empty())
        return;
    int CHECK_NODES = a.size();
    int VARIABLE_NODES = a[0].size();
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            for (size_t k = 0; k + 1 < a[i][j].size(); ++k)
                cerr << a[i][j][k] << "&";
            if (a[i][j].empty())
                cerr << -1 << "\t";
            else
                cerr << a[i][j][a[i][j].size() - 1] << "\t";
        }
        cerr << endl;
    }
}

