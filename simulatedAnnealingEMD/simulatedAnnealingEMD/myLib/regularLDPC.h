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
#include<iostream>
#include<vector>
#include<algorithm>
#include<set>
#include<sstream>
#include<ctime>
#include<cassert>
#include<map>
#include<queue>

using namespace std;

#define ll long long
#define pii pair<int, int>

bool nextCombination(vector<int> & a, int n);
ll getBinomial(ll n, ll k);
double getBigBinomial(ll n, ll k);
void print(const vector<int>& a);
void print(const vector<vector<int> >& a);
void eprint(const vector<int>& a);
void eprint(const vector<vector<int> >& a);
bool noCyclesOfLength4(const vector<vector<int> >& a, int m);
bool noCyclesOfLength6(const vector<vector<int> >& a, int m);
bool girthAtLeast6(const vector<vector<int> >& a, int m);
bool girthAtLeast8(const vector<vector<int> >& a, int m);
bool noCyclesofLength4ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m);
bool noCyclesofLength6ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m);
bool girthAtLeast6ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m);
bool girthAtLeast8ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m);
void printGapMatrix(const vector<vector<int> >& a);
bool isZeroes(const vector<int>& a);
string toStr(ll x);
bool toUnsignedInt(string s, ll& x);

string toStr(ll x) {
    stringstream s;
    s << x;
    return s.str();
}

bool toUnsignedInt(string s, ll& x) {
    for (size_t i = 0; i < s.size(); ++i)
        if ((s[i] < '0') || (s[i] > '9'))
            return 0;
    x = 0;
    for (size_t i = 0; i < s.size(); ++i)
        x = x * 10 + s[i] - '0';
    return 1;
}

bool isZeroes(const vector<int>& a) {
	for (size_t i = 0; i < a.size(); ++i)
		if (a[i] != 0)
			return 0;
	return 1;
}

void printGapMatrix(const vector<vector<int> >& a) {
	cout << "M:=[";
	for (size_t i = 0; i + 1 < a.size(); ++i) {
		cout << "[";
		for (size_t j = 0; j + 1 < a[i].size(); ++j)
			cout << a[i][j] << ", ";
		cout << a[i][a[i].size() - 1] << "],";
	}
	cout << "[";
	for (size_t j = 0; j + 1 < a[a.size() - 1].size(); ++j)
		cout << a[a.size() - 1][j] << ", ";
	cout << a[a.size() - 1][a[a.size() - 1].size() - 1] << "]";
	cout << "];";
}


bool noCyclesofLength4ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m) {
	int j = a.size();
	int l = a[0].size();
	for (int r1 = 1; r1 < j; ++r1) {
		for (int r2 = r1 + 1; r2 < j; ++r2) {
			for (int c1 = 0; c1 < l; ++c1) {
				if ((a[r1][c1] == -1) || (a[r2][c1] == -1))
					continue;
				for (int c2 = c1 + 1; c2 < l; ++c2) {
					if ((a[r1][c2] == -1) || (a[r2][c2] == -1))
						continue;
					if ((a[r1][c1] + a[r2][c2] - a[r1][c2] - a[r2][c1] + 2 * m) % m == 0)
						return 0;
				}
			}
		}
	}
	for (int r = 1; r < j; ++r) {
		for (int c1 = 0; c1 < l; ++c1) {
			if (a[r][c1] == -1)
				continue;
			for (int c2 = c1 + 1; c2 < l; ++c2) {
				if (a[r][c2] == -1)
					continue;
				if (a[r][c1] == a[r][c2])
					return 0;
			}
		}
	}
	return 1;
}


bool noCyclesofLength6ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m) {
	int j = a.size();
	int l = a[0].size();
	for (int r1 = 1; r1 < j; ++r1) {
		for (int r2 = r1 + 1; r2 < j; ++r2) {
			for (int r3 = r2 + 1; r3 < j; ++r3) {
				for (int c1 = 0; c1 < l; ++c1) {
					if (a[r1][c1] == -1)
						continue;
					for (int c2 = 0; c2 < l; ++c2) {
						if (c1 == c2)
							continue;
						if (a[r2][c2] == -1)
							continue;
						for (int c3 = 0; c3 < l; ++c3) {
							if (a[r3][c3] == -1)
								continue;
							if ((c3 == c1) || (c3 == c2))
								continue;
							if ((a[r1][c2] != -1) && (a[r2][c3] != -1) && (a[r3][c1] != -1) &&
								(((a[r1][c1] + a[r2][c2] + a[r3][c3] - a[r1][c2] - a[r2][c3] - a[r3][c1] + 3 * m) % m) == 0))
								return 0;

							if ((a[r1][c3] != -1) && (a[r2][c1] != -1) && (a[r3][c2] != -1) &&
								(((a[r1][c1] + a[r2][c2] + a[r3][c3] - a[r1][c3] - a[r2][c1] - a[r3][c2] + 3 * m) % m) == 0))
								return 0;
						}
					}
				}
			}
		}
	}
	for (int r1 = 1; r1 < j; ++r1) {
		for (int r2 = r1 + 1; r2 < j; ++r2) {
			for (int c1 = 0; c1 < l; ++c1) {
				if ((a[r1][c1] == -1) || (a[r2][c1] == -1))
					continue;
				for (int c2 = 0; c2 < l; ++c2) {
					if (a[r2][c2] == -1)
						continue;
					if (c1 == c2)
						continue;
					for (int c3 = 0; c3 < l; ++c3) {
						if (a[r1][c3] == -1)
							continue;
						if (c3 == c1)
							continue;
						if ((a[r1][c1] + a[r2][c2] - a[r2][c1] - a[r1][c3] + 2 * m) % m == 0)
							return 0;
					}
				}
			}
		}
	}
	return 1;
}


bool girthAtLeast6ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m) {
	if (!isZeroes(a[0])) {
		cerr << "NONZERO FIRST ROW\n";
		assert(0);
	}
	return noCyclesofLength4ForMatrixWithFirstRowOfZeroes(a, m);
}

bool girthAtLeast8ForMatrixWithFirstRowOfZeroes(const vector<vector<int> >& a, int m) {
	if (!isZeroes(a[0])) {
		cerr << "NONZERO FIRST ROW\n";
		assert(0);
	}
	return noCyclesofLength4ForMatrixWithFirstRowOfZeroes(a, m) && noCyclesofLength6ForMatrixWithFirstRowOfZeroes(a, m);
}


bool girthAtLeast6(const vector<vector<int> >& a, int m) {
	return noCyclesOfLength4(a, m);
}

bool girthAtLeast8(const vector<vector<int> >& a, int m) {
	return noCyclesOfLength4(a, m) && noCyclesOfLength6(a, m);
}

bool noCyclesOfLength6(const vector<vector<int> >& a, int m) {
	int j = a.size();
	int l = a[0].size();
	for (int r1 = 0; r1 < j; ++r1) {
		for (int r2 = r1 + 1; r2 < j; ++r2) {
			for (int r3 = r2 + 1; r3 < j; ++r3) {
				for (int c1 = 0; c1 < l; ++c1) {
					if (a[r1][c1] == -1)
						continue;
					for (int c2 = 0; c2 < l; ++c2) {
						if (c1 == c2)
							continue;
						if (a[r2][c2] == -1)
							continue;
						for (int c3 = 0; c3 < l; ++c3) {
							if (a[r3][c3] == -1)
								continue;
							if ((c3 == c1) || (c3 == c2))
								continue;
							if ((a[r1][c2] != -1) && (a[r2][c3] != -1) && (a[r3][c1] != -1) &&
								(((a[r1][c1] + a[r2][c2] + a[r3][c3] - a[r1][c2] - a[r2][c3] - a[r3][c1] + 3 * m) % m) == 0))
								return 0;

							if ((a[r1][c3] != -1) && (a[r2][c1] != -1) && (a[r3][c2] != -1) &&
								(((a[r1][c1] + a[r2][c2] + a[r3][c3] - a[r1][c3] - a[r2][c1] - a[r3][c2] + 3 * m) % m) == 0))
								return 0;
						}
					}
				}
			}
		}
	}
	return 1;
}

bool noCyclesOfLength4(const vector<vector<int> >& a, int m) {
	int j = a.size();
	int l = a[0].size();
	for (int r1 = 0; r1 < j; ++r1) {
		for (int r2 = r1 + 1; r2 < j; ++r2) {
			for (int c1 = 0; c1 < l; ++c1) {
				if ((a[r1][c1] == -1) || (a[r2][c1] == -1))
					continue;
				for (int c2 = c1 + 1; c2 < l; ++c2) {
					if ((a[r1][c2] == -1) || (a[r2][c2] == -1))
						continue;
					if ((a[r1][c1] + a[r2][c2] - a[r1][c2] - a[r2][c1] + 2 * m) % m == 0)
						return 0;
				}
			}
		}
	}
	return 1;
}

bool nextCombination(vector<int>& a, int n) {
	int k = (int)a.size();
	for (int i = k - 1; i >= 0; --i)
		if (a[i] < n - k + i + 1) {
			++a[i];
			for (int j = i + 1; j < k; ++j)
				a[j] = a[j - 1] + 1;
			return true;
		}
	return false;
}

ll getBinomial(ll n, ll k) {
	ll res = 1;
	for (int i = 1; i <= k; ++i) {
		res = res * (n - i + 1);
		res /= i;
	}
	return res;
}

double getBigBinomial(ll n, ll k) {
	double res = 1.0;
	for (int i = 1; i <= k; ++i) {
		res = res * (n - i + 1);
		res /= i;
	}
	return res;
}

void print(const vector<int>& a) {
	for (size_t i = 0; i < a.size(); ++i) {
		cout << a[i] << "\t";
	}
	cout << endl;
}

void print(const vector<vector<int> >& a) {
	for (size_t i = 0; i < a.size(); ++i) {
        print(a[i]);
	}
}

void eprint(const vector<vector<int> >& a) {
	for (size_t i = 0; i < a.size(); ++i) {
        eprint(a[i]);
	}
}

void eprint(const vector<int>& a) {
	for (size_t i = 0; i < a.size(); ++i) {
		cerr << a[i] << "\t";
	}
	cerr << endl;
}
