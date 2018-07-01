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
#include"..\myLib\CycleEnum.h"
#include<iostream>
#include<ctime>
#include<algorithm>
#include<map>

#define pii pair<int, int>
const double eps = 1e-9;

static struct Tiii {
    int first, second, third;
    Tiii() {}
    Tiii(int _first, int _second, int _third) : first(_first), second(_second), third(_third) {}
    Tiii(const entry& e) : first(e.r), second(e.c), third(e.id) {}
    bool operator==(const Tiii& rhs) const {
        return (first == rhs.first) && (second == rhs.second) && (third == rhs.third);
    }
    bool operator!=(const Tiii& rhs) const {
        return (first != rhs.first) || (second != rhs.second) || (third != rhs.third);
    }
    bool operator<(const Tiii& rhs) const {
        return (first < rhs.first) ||
            ((first == rhs.first) && (second < rhs.second)) ||
            ((first == rhs.first) && (second == rhs.second) && (third == rhs.third));
    }
};

class emdOpt {
    struct Cycle {
        vector<Tiii> cycle;
        int uniqueNodes;
        bool oneEntry;
        Cycle(const vector<entry>& _cycle) {
			uniqueNodes = 1;
			oneEntry = true;
            cycle.resize(_cycle.size());
            for (size_t i = 0; i < cycle.size(); ++i)
                cycle[i] = _cycle[i];
            for (int i = 1; i < cycle.size(); ++i) {
                if (cycle[i] != cycle[i - 1]) {
                    oneEntry = false;
                    break;
                }
            }
            for (int i = 1; i < cycle.size(); ++i) {
                ++uniqueNodes;
                for (int j = 0; j < i; ++j)
                    if (cycle[i] == cycle[j]) {
                        --uniqueNodes;
                        break;
                    }
            }
        }
    };
    int checkNodes, variableNodes;
    vector<vector<int> > protograph;
    long long circulant, targetEmd, upGirth, targetGirth;
    vector<vector<vector<int> > > mtr;
    vector<vector<int> > sparseChecks, sparseVars;
    vector<int> posByLiftVals;
    bool noCycles;
public:
    emdOpt(long long _circulant, long long _upGirth, const vector<vector<vector<int> > >& _mtr) {
		noCycles = false;
        circulant = _circulant, upGirth = _upGirth, mtr = _mtr;
        checkNodes = mtr.size(), variableNodes = mtr[0].size();
        protograph.assign(checkNodes, vector<int>(variableNodes));
        sparseChecks.resize(checkNodes);
        sparseVars.resize(variableNodes);
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                protograph[r][c] = mtr[r][c].size();
                for (int id = 0; id < mtr[r][c].size(); ++id) {
                    sparseChecks[r].push_back(c);
                    sparseVars[c].push_back(r);
                }
            }
        }

    }
    emdOpt(long long _circulant, long long _upGirth, long long _targetGirth, long long _targetEmd, const vector<vector<vector<int> > >& _mtr) {
		noCycles = false;
        circulant = _circulant, upGirth = _upGirth, targetGirth = _targetGirth, targetEmd = _targetEmd, mtr = _mtr;
        checkNodes = mtr.size(), variableNodes = mtr[0].size();
        protograph.assign(checkNodes, vector<int>(variableNodes));
        sparseChecks.resize(checkNodes);
        sparseVars.resize(variableNodes);
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                protograph[r][c] = mtr[r][c].size();
                for (int id = 0; id < mtr[r][c].size(); ++id) {
                    sparseChecks[r].push_back(c);
                    sparseVars[c].push_back(r);
                }
            }
        }
    }
    emdOpt(long long _circulant, long long _upGirth, long long _targetGirth, long long _targetEmd, const vector<vector<int> >& _mtr) {
		noCycles = false;
        circulant = _circulant, upGirth = _upGirth, targetEmd = _targetEmd, targetGirth = _targetGirth;
        checkNodes = _mtr.size(), variableNodes = _mtr[0].size();
        mtr.assign(checkNodes, vector<vector<int> >(variableNodes));
        protograph.assign(checkNodes, vector<int>(variableNodes));
        sparseChecks.resize(checkNodes);
        sparseVars.resize(variableNodes);
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                if (_mtr[r][c] != -1)
                    mtr[r][c].push_back(_mtr[r][c]);
                protograph[r][c] = mtr[r][c].size();
                for (int id = 0; id < mtr[r][c].size(); ++id) {
                    sparseChecks[r].push_back(c);
                    sparseVars[c].push_back(r);

                }

            }
        }
    }
    bool optimizeEmd() {
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));

        long long numberOfValuesToAssign = 0;
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));

                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));

            }
        }
        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth))
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];

            } while (enumerator.next());
        }

        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        long long moves = 0;
        while (true) {
            for (int r = 0; r < checkNodes; ++r) {
                for (int c = 0; c < variableNodes; ++c) {
                    for (int id = 0; id < protograph[r][c]; ++id) {
                        ++movesWithoutChange;
                        ++moves;
                        if (movesWithoutChange > numberOfValuesToAssign) {
                            noCycles = !cycleWithLowEmdOrGirthExists;
                            return noCycles;
                        }
                        numberOfCyclesWithTheseValues[r][c][id].assign(circulant, 0);
                        numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0);

                        for (int girth = 4; girth <= upGirth; girth += 2) {
                            int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];//numberOfCyclesOfFixedLenContainsThisEntry[r][c][id][girth];
                            if (numberOfCyclesToProcess == 0)
                                continue;
                            CycleEnum enumerator(girth, protograph);
                            enumerator.init(r, c, id);
                            for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                                while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                                    enumerator.next();
                                Cycle cycle(enumerator.cycle);
                                processCycle(cycle, numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant);
                                enumerator.next();
                            }
                        }
                        for (int i = 0; i < circulant; ++i) {
                            if (i == mtr[r][c][id])
                                continue;
                            //if (numberOfCyclesWithTheseValues[r][c][id][i] < numberOfCyclesWithTheseValues[r][c][id][mtr[r][c][id]]) {
                            if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                                movesWithoutChange = 0;
                                mtr[r][c][id] = i;
                                cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                            }
                        }
                        cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }
                }
            }
        }
    }
    vector<vector<vector<int> > > getMatrix() {
        return mtr;
    }
    vector<vector<int> > getRegMatrix() {
        vector<vector<int> > res(mtr.size(), vector<int>(mtr[0].size()));
        for (int i = 0; i < res.size(); ++i) {
            for (int j = 0; j < res[i].size(); ++j) {
                if (mtr[i][j].empty())
                    res[i][j] = -1;
                else
                    res[i][j] = mtr[i][j][0];
            }
        }
        return res;
    }
    bool annealEmd(vector<Tiii> order = vector<Tiii>()) {
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<int> > > isInQueue(checkNodes, vector<vector<int> >(variableNodes));
        long long numberOfValuesToAssign = 0;
        long long totalNumberOfCyclesWithLowEmdOrGirth = 0;
        long long moves = numberOfValuesToAssign - order.size();
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                isInQueue[r][c].assign(protograph[r][c], 0);
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));

                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));

            }
        }
        for (int i = 0; i < order.size(); ++i) {
            isInQueue[order[i].first][order[i].second][order[i].third] = 1;
        }
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (isInQueue[r][c][i])
                        continue;
                    order.push_back(Tiii(r, c, i));
                }
            }
        }

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                /*++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                ++totalNumberOfCycles;*/
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth)) {
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                    ++totalNumberOfCyclesWithLowEmdOrGirth;
                }

                /*if (curEmd < targetEmd) {
                ++numberOfCyclesOfFixedLenWithLowEmdContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                }*/
            } while (enumerator.next());
        }
        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        while (true) {
            for (int it = 0; it < order.size(); ++it) {
                ++moves;
                int r = order[it].first, c = order[it].second, id = order[it].third;
                double temperature = 1.0 * totalNumberOfCyclesWithLowEmdOrGirth / moves / moves;
                ++movesWithoutChange;
                if (movesWithoutChange > 20 * numberOfValuesToAssign) {
                    return !cycleWithLowEmdOrGirthExists;
                }
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0);
                for (int girth = 4; girth <= upGirth; girth += 2) {
                    //int numberOfCyclesToProcess = numberOfCyclesOfFixedLenContainsThisEntry[r][c][id][girth];
                    int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];
                    if (numberOfCyclesToProcess == 0)
                        continue;
                    CycleEnum enumerator(girth, protograph);
                    enumerator.init(r, c, id);
                    for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                        while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                            enumerator.next();
                        Cycle cycle(enumerator.cycle);
                        processCycle(cycle, numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant);
                        enumerator.next();
                    }
                }
                vector<int> used(circulant, 0);
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (i == id)
                        continue;
                    used[mtr[r][c][i]] = 1;
                }
                for (int i = 0; i < circulant; ++i) {
                    if (i == mtr[r][c][id])
                        continue;
                    if (used[i])
                        continue;
                    if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                        movesWithoutChange = 0;
                        //print(mtr);
                        /*freopen("err_log.txt", "a", stderr);
                        eprint(mtr);
                        cerr << endl;
                        fclose(stderr);*/
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }
                }
                cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] == 0)
                    continue;
                vector<double> prob(circulant);
                double sumProb = 0;
                for (int i = 0; i < circulant; ++i) {
                    if (used[i])
                        continue;
                    prob[i] = exp(-(numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] - numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) / temperature);
                    /*if (i != mtr[r][c][id])
                    prob[i] /= circulant;*/
                    /*if ((numberOfCyclesWithTheseValues[r][c][id][i] - numberOfCyclesWithTheseValues[r][c][id][mtr[r][c][id]] == 0) && (i != mtr[r][c][id]))
                    prob[i] = exp(-0.5 / temperature);*/
                    sumProb += prob[i];
                }
                double randMove = sumProb * rand() / RAND_MAX;
                double sum = 0;
                for (int i = 0; i < circulant; ++i) {
                    sum += prob[i];
                    if (sum > randMove) {
                        if (mtr[r][c][id] == i)
                            break;
                        /*freopen("err_log.txt", "a", stderr);
                        eprint(mtr);
                        cerr << endl;
                        fclose(stderr);*/
                        if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] != numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]])
                            movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                        break;
                    }
                }
            }
        }
    }


    bool annealEmdWithFixed(const vector<vector<vector<int> > >& fixed, double time = -1, vector<Tiii> order = vector<Tiii>()) {
        clock_t startTime = clock();
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<int> > > isInQueue(checkNodes, vector<vector<int> >(variableNodes));
        long long numberOfValuesToAssign = 0;
        long long totalNumberOfCyclesWithLowEmdOrGirth = 0;
        long long moves = numberOfValuesToAssign - order.size();
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                isInQueue[r][c].assign(protograph[r][c], 0);
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));

                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));

            }
        }
        for (int i = 0; i < order.size(); ++i) {
            isInQueue[order[i].first][order[i].second][order[i].third] = 1;
        }
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (isInQueue[r][c][i])
                        continue;
                    order.push_back(Tiii(r, c, i));
                }
            }
        }

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth)) {
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                    ++totalNumberOfCyclesWithLowEmdOrGirth;
                }

            } while (enumerator.next());
        }
        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        while (true) {
            for (int it = 0; it < order.size(); ++it) {
                int r = order[it].first, c = order[it].second, id = order[it].third;
                ++moves;

                //                 if (fixed[r][c][id])
                //                     continue;

                clock_t curTime = clock();
                if ((time > -0.5) && (curTime - startTime > CLOCKS_PER_SEC * time))
                    return 0;

                double temperature = 1.0 * totalNumberOfCyclesWithLowEmdOrGirth / moves / moves;
                ++movesWithoutChange;
                if (movesWithoutChange > 20 * numberOfValuesToAssign) {
                    return !cycleWithLowEmdOrGirthExists;
                }
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0);
                for (int girth = 4; girth <= upGirth; girth += 2) {
                    int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];
                    if (numberOfCyclesToProcess == 0)
                        continue;
                    CycleEnum enumerator(girth, protograph);
                    enumerator.init(r, c, id);
                    for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                        while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                            enumerator.next();
                        Cycle cycle(enumerator.cycle);
                        processCycle(cycle, numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant);
                        enumerator.next();
                    }
                }
                vector<int> used(circulant, 0);
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (i == id)
                        continue;
                    used[mtr[r][c][i]] = 1;
                }
                for (int i = 0; i < circulant; ++i) {
                    if (i == mtr[r][c][id])
                        continue;
                    if ((used[i]) || (fixed[r][c][id]))
                        continue;
                    if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                        movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }
                }
                cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] == 0)
                    continue;
                if (fixed[r][c][id])
                    continue;
                vector<double> prob(circulant);
                double sumProb = 0;
                for (int i = 0; i < circulant; ++i) {
                    if (used[i])
                        continue;
                    prob[i] = exp(-(numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] - numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) / temperature);
                    sumProb += prob[i];
                }
                double randMove = sumProb * rand() / RAND_MAX;
                double sum = 0;
                for (int i = 0; i < circulant; ++i) {
                    sum += prob[i];
                    if (sum > randMove) {
                        if (mtr[r][c][id] == i)
                            break;
                        if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] != numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]])
                            movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                        break;
                    }
                }
            }
        }
    }


    bool annealEmdWithFixedAndCycleCost(const vector<vector<vector<int> > >& fixed, double time = -1, vector<Tiii> order = vector<Tiii>()) {
        clock_t startTime = clock();
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));


        vector<vector<vector<vector<double> > > > costOfCyclesWithTheseValues(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<double> > >(variableNodes));

        vector<vector<vector<int> > > isInQueue(checkNodes, vector<vector<int> >(variableNodes));
        long long numberOfValuesToAssign = 0;
        long long totalNumberOfCyclesWithLowEmdOrGirth = 0;
        long long moves = numberOfValuesToAssign - order.size();
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                isInQueue[r][c].assign(protograph[r][c], 0);
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));


                costOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<double>(circulant));
                costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<double>(circulant));
                costOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<double>(upGirth + 1, 0));
                costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<double>(upGirth + 1, 0));


            }
        }
        for (int i = 0; i < order.size(); ++i) {
            isInQueue[order[i].first][order[i].second][order[i].third] = 1;
        }
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (isInQueue[r][c][i])
                        continue;
                    order.push_back(Tiii(r, c, i));
                }
            }
        }

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                double cost = getCost(girth, curEmd);
                costOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth] += cost;

                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth)) {
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                    costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth] += cost;

                    ++totalNumberOfCyclesWithLowEmdOrGirth;
                }

            } while (enumerator.next());
        }
        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        while (true) {
            for (int it = 0; it < order.size(); ++it) {
                int r = order[it].first, c = order[it].second, id = order[it].third;
                ++moves;

                //                 if (fixed[r][c][id])
                //                     continue;

                clock_t curTime = clock();
                if ((time > -0.5) && (curTime - startTime > CLOCKS_PER_SEC * time))
                    return 0;

                double temperature = 1.0 * totalNumberOfCyclesWithLowEmdOrGirth / moves / moves;
                ++movesWithoutChange;
                if (movesWithoutChange > 20 * numberOfValuesToAssign) {
                    return !cycleWithLowEmdOrGirthExists;
                }
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0);
                costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0.0);

                for (int girth = 4; girth <= upGirth; girth += 2) {
                    int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];
                    if (numberOfCyclesToProcess == 0)
                        continue;
                    CycleEnum enumerator(girth, protograph);
                    enumerator.init(r, c, id);
                    for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                        while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                            enumerator.next();
                        Cycle cycle(enumerator.cycle);
                        int curEmd = getEmd(cycle);
                        processCycle(cycle, getCost(girth, curEmd), costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant);
                        enumerator.next();
                    }
                }
                vector<int> used(circulant, 0);
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (i == id)
                        continue;
                    used[mtr[r][c][i]] = 1;
                }
                for (int i = 0; i < circulant; ++i) {
                    if (i == mtr[r][c][id])
                        continue;
                    if ((used[i]) || (fixed[r][c][id]))
                        continue;
                    /*if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                        movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }*/
                    if (costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] + eps < costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                        movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }
                }
                cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] == 0)
                    continue;
                if (fixed[r][c][id])
                    continue;
                vector<double> prob(circulant);
                double sumProb = 0;
                for (int i = 0; i < circulant; ++i) {
                    if (used[i])
                        continue;
                    prob[i] = exp(-(costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] - costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) / temperature);
                    sumProb += prob[i];
                }
                double randMove = sumProb * rand() / RAND_MAX;
                double sum = 0;
                for (int i = 0; i < circulant; ++i) {
                    sum += prob[i];
                    if (sum > randMove) {
                        if (mtr[r][c][id] == i)
                            break;
                        if (costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] + eps < costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]])
                            movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                        break;
                    }
                }
            }
        }
    }


    bool annealEmdWithFixedAndCycleCostAndInnerCode(const vector<vector<vector<int> > >& fixed, double time = -1, vector<Tiii> order = vector<Tiii>()) {
        clock_t startTime = clock();
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));


        vector<vector<vector<vector<double> > > > costOfCyclesWithTheseValues(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<double> > >(variableNodes));
        vector<vector<vector<vector<double> > > > costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<double> > >(variableNodes));

        vector<vector<vector<int> > > isInQueue(checkNodes, vector<vector<int> >(variableNodes));
        long long numberOfValuesToAssign = 0;
        long long totalNumberOfCyclesWithLowEmdOrGirth = 0;
        long long moves = numberOfValuesToAssign - order.size();
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                isInQueue[r][c].assign(protograph[r][c], 0);
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(circulant));
                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));


                costOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<double>(circulant));
                costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<double>(circulant));
                costOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<double>(upGirth + 1, 0));
                costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<double>(upGirth + 1, 0));


            }
        }
        for (int i = 0; i < order.size(); ++i) {
            isInQueue[order[i].first][order[i].second][order[i].third] = 1;
        }
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (isInQueue[r][c][i])
                        continue;
                    order.push_back(Tiii(r, c, i));
                }
            }
        }

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                double cost = getCostInnerCode(enumerator.cycle, girth, curEmd);
                costOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth] += cost;

                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth)) {
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                    costOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth] += cost;

                    ++totalNumberOfCyclesWithLowEmdOrGirth;
                }

            } while (enumerator.next());
        }
        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        while (true) {
            for (int it = 0; it < order.size(); ++it) {
                int r = order[it].first, c = order[it].second, id = order[it].third;
                ++moves;

                //                 if (fixed[r][c][id])
                //                     continue;

                clock_t curTime = clock();
                if ((time > -0.5) && (curTime - startTime > CLOCKS_PER_SEC * time))
                    return 1;

                double temperature = 1.0 * totalNumberOfCyclesWithLowEmdOrGirth / moves / moves;
                ++movesWithoutChange;
                if (movesWithoutChange > 20 * numberOfValuesToAssign) {
                    return 1;// !cycleWithLowEmdOrGirthExists;
                }
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0);
                costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(circulant, 0.0);

                for (int girth = 4; girth <= upGirth; girth += 2) {
                    int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];
                    if (numberOfCyclesToProcess == 0)
                        continue;
                    CycleEnum enumerator(girth, protograph);
                    enumerator.init(r, c, id);
                    for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                        while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                            enumerator.next();
                        Cycle cycle(enumerator.cycle);
                        int curEmd = getEmd(cycle);
                        processCycle(cycle, getCost(girth, curEmd), costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant);
                        enumerator.next();
                    }
                }
                vector<int> used(circulant, 0);
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (i == id)
                        continue;
                    used[mtr[r][c][i]] = 1;
                }
                for (int i = 0; i < circulant; ++i) {
                    if (i == mtr[r][c][id])
                        continue;
                    if ((used[i]) || (fixed[r][c][id]))
                        continue;
                    /*if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                    movesWithoutChange = 0;
                    mtr[r][c][id] = i;
                    cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }*/
                    if (costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] + eps < costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) {
                        cerr << costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] - costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] << endl;
                        movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                    }
                }
                cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] == 0)
                    continue;
                if (fixed[r][c][id])
                    continue;
                vector<double> prob(circulant);
                double sumProb = 0;
                for (int i = 0; i < circulant; ++i) {
                    if (used[i])
                        continue;
                    prob[i] = exp(-(costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] - costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]]) / temperature);
                    sumProb += prob[i];
                }
                double randMove = sumProb * rand() / RAND_MAX;
                double sum = 0;
                for (int i = 0; i < circulant; ++i) {
                    sum += prob[i];
                    if (sum > randMove) {
                        if (mtr[r][c][id] == i)
                            break;
                        if (costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] + eps < costOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]])
                            movesWithoutChange = 0;
                        mtr[r][c][id] = i;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][mtr[r][c][id]] > 0);
                        break;
                    }
                }
            }
        }
    }


    bool annealEmdWithFixed(const vector<vector<vector<int> > >& fixed, const vector<int>& liftVals, double time = -1, vector<Tiii> order = vector<Tiii>()) {
        //matrix values are supposed to be from liftVals
        int numOfLifts = liftVals.size();
        posByLiftVals.assign(circulant, -1);
        for (int i = 0; i < liftVals.size(); ++i)
            posByLiftVals[liftVals[i]] = i;
        clock_t startTime = clock();
        vector<vector<vector<vector<int> > > > numberOfCyclesWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesWithLowEmdOrGirthWithTheseValues(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));
        vector<vector<vector<vector<int> > > > numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry(checkNodes, vector<vector<vector<int> > >(variableNodes));

        vector<vector<vector<int> > > isInQueue(checkNodes, vector<vector<int> >(variableNodes));
        long long numberOfValuesToAssign = 0;
        long long totalNumberOfCyclesWithLowEmdOrGirth = 0;
        long long moves = numberOfValuesToAssign - order.size();
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                isInQueue[r][c].assign(protograph[r][c], 0);
                numberOfValuesToAssign += protograph[r][c];
                numberOfCyclesWithTheseValues[r][c].assign(protograph[r][c], vector<int>(numOfLifts));
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c].assign(protograph[r][c], vector<int>(numOfLifts));

                numberOfCyclesOfFixedLenContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));
                numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c].assign(protograph[r][c], vector<int>(upGirth + 1, 0));

            }
        }
        for (int i = 0; i < order.size(); ++i) {
            isInQueue[order[i].first][order[i].second][order[i].third] = 1;
        }
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (isInQueue[r][c][i])
                        continue;
                    order.push_back(Tiii(r, c, i));
                }
            }
        }

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int curEmd = getEmd(enumerator.cycle);
                ++numberOfCyclesOfFixedLenContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                if ((curEmd < targetEmd) || (girth + 2 <= targetGirth)) {
                    ++numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[enumerator.cycle[0].r][enumerator.cycle[0].c][enumerator.cycle[0].id][girth];
                    ++totalNumberOfCyclesWithLowEmdOrGirth;
                }

            } while (enumerator.next());
        }
        int movesWithoutChange = 0;
        bool cycleWithLowEmdOrGirthExists = 0;
        while (true) {
            for (int it = 0; it < order.size(); ++it) {
                int r = order[it].first, c = order[it].second, id = order[it].third;
                ++moves;

                //                 if (fixed[r][c][id])
                //                     continue;

                clock_t curTime = clock();
                if ((time > -0.5) && (curTime - startTime > CLOCKS_PER_SEC * time))
                    return 0;

                double temperature = 1.0 * totalNumberOfCyclesWithLowEmdOrGirth / moves / moves;
                ++movesWithoutChange;
                if (movesWithoutChange > 20 * numberOfValuesToAssign) {
                    return !cycleWithLowEmdOrGirthExists;
                }
                numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id].assign(numOfLifts, 0);
                for (int girth = 4; girth <= upGirth; girth += 2) {
                    int numberOfCyclesToProcess = numberOfCyclesOfFixedLenWithLowEmdOrGirthContainsThisEntry[r][c][id][girth];
                    if (numberOfCyclesToProcess == 0)
                        continue;
                    CycleEnum enumerator(girth, protograph);
                    enumerator.init(r, c, id);
                    for (int cyclesCounter = 0; cyclesCounter < numberOfCyclesToProcess; ++cyclesCounter) {
                        while ((girth + 2 > targetGirth) && (getEmd(enumerator.cycle) >= targetEmd))
                            enumerator.next();
                        Cycle cycle(enumerator.cycle);
                        processCycle(cycle, numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id], mtr, circulant, 0);
                        enumerator.next();
                    }
                }
                vector<int> used(numOfLifts, 0);
                for (int i = 0; i < protograph[r][c]; ++i) {
                    if (i == id)
                        continue;
                    used[posByLiftVals[mtr[r][c][i]]] = 1;
                }
                for (int i = 0; i < numOfLifts; ++i) {
                    int curLiftVal = liftVals[i];
                    if (curLiftVal == mtr[r][c][id])
                        continue;
                    if ((used[i]) || (fixed[r][c][id]))
                        continue;
                    if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] < numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][posByLiftVals[mtr[r][c][id]]]) {
                        movesWithoutChange = 0;
                        mtr[r][c][id] = curLiftVal;
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] > 0);
                    }
                }
                cycleWithLowEmdOrGirthExists = cycleWithLowEmdOrGirthExists || (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][posByLiftVals[mtr[r][c][id]]] > 0);
                if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][posByLiftVals[mtr[r][c][id]]] == 0)
                    continue;
                if (fixed[r][c][id])
                    continue;
                vector<double> prob(numOfLifts);
                double sumProb = 0;
                for (int i = 0; i < numOfLifts; ++i) {
                    if (used[i])
                        continue;
                    prob[i] = exp(-(numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] - 
                        numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][posByLiftVals[mtr[r][c][id]]]) / temperature);
                    sumProb += prob[i];
                }
                double randMove = sumProb * rand() / RAND_MAX;
                double sum = 0;
                for (int i = 0; i < numOfLifts; ++i) {
                    sum += prob[i];
                    if (sum > randMove) {
                        if (mtr[r][c][id] == liftVals[i])
                            break;
                        if (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] != numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][posByLiftVals[mtr[r][c][id]]])
                            movesWithoutChange = 0;
                        mtr[r][c][id] = liftVals[i];
                        cycleWithLowEmdOrGirthExists = (numberOfCyclesWithLowEmdOrGirthWithTheseValues[r][c][id][i] > 0);
                        break;
                    }
                }
            }
        }
    }



    pair<int, int> getGirthAndEmd() {
        int resGirth = upGirth + 2;
        int resEmd = 1000000;

        for (int girth = 4; girth <= upGirth; girth += 2) {
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int sum = 0;
                for (int i = 0; i < enumerator.cycle.size(); ++i) {
                    if (i & 1)
                        sum += mtr[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id];
                    else
                        sum -= mtr[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id];
                }
                sum = ((sum % circulant) + circulant) % circulant;
                if (sum != 0)
                    continue;
                int curEmd = getEmd(enumerator.cycle);
                resGirth = min(resGirth, girth);
                resEmd = min(resEmd, curEmd);
            } while (enumerator.next());
        }
        return make_pair(resGirth, resEmd);
    }

    double getCost(int girth, int emd) {
        return 1000 / pow(emd + 1, 2) / pow(10, girth);
    }

    double getCostInnerCode(const Cycle& cycle, int girth, int curEmd) {
        double res = getCost(girth, curEmd);
        bool oneGroup = 1;
        for (int i = 1; i < cycle.cycle.size(); ++i) {
            if (cycle.cycle[i].second % 3 != cycle.cycle[i - 1].second % 3) {
                oneGroup = 0;
                break;
            }
        }
        if (oneGroup)
            res *= 10;
        return res;
    }

    vector<vector<pair<int, int> > > getEMDDistr() {
        vector<vector<pair<int, int> > > res(upGirth + 1);
        for (int girth = 4; girth <= upGirth; girth += 2) {
            map<int, int> emdCounter;
            CycleEnum enumerator(girth, protograph);
            if (!enumerator.init()) {
                continue;
            }
            do {
                int sum = 0;
                for (int i = 0; i < enumerator.cycle.size(); ++i) {
                    if (i & 1)
                        sum += mtr[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id];
                    else
                        sum -= mtr[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id];
                }
                sum = ((sum % circulant) + circulant) % circulant;
                if (sum != 0)
                    continue;
                int curEmd = getLiftedEmd(enumerator.cycle);
                emdCounter[curEmd] += 1;
            } while (enumerator.next());
			for (map<int, int>::iterator it = emdCounter.begin(); it != emdCounter.end(); ++it) {
                res[girth].push_back(*it);
            }
        }
        return res;

    }


private:

    int getEmd(const Cycle& cycle) {
        int res = 0;
        vector<int> cnt(checkNodes, 0);
        vector<int> var(variableNodes, 0);
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            var[cycle.cycle[i].second] = 1;
        }
        for (int i = 0; i < var.size(); ++i) {
            if (!var[i])
                continue;
            for (int j = 0; j < sparseVars[i].size(); ++j) {
                ++cnt[sparseVars[i][j]];
            }
        }
        return count(cnt.begin(), cnt.end(), 1);
    }
    int getLiftedEmd(const Cycle& cycle) {
        vector<pair<pii, int> > varNodes(1, make_pair(pii(cycle.cycle[0].second, cycle.cycle[0].third), 0));
        for (int i = 1; i < cycle.cycle.size(); i += 2) {
            varNodes.push_back(make_pair(pii(cycle.cycle[i].second, cycle.cycle[i].third), 
                (varNodes[i / 2].second + mtr[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third] -
                mtr[cycle.cycle[i - 1].first][cycle.cycle[i - 1].second][cycle.cycle[i - 1].third] + circulant) % circulant));
        }
        sort(varNodes.begin(), varNodes.end());
        varNodes.resize(unique(varNodes.begin(), varNodes.end()) - varNodes.begin());
        int res = 0;
        for (int r = 0; r < checkNodes; ++r) {
            vector<int> ls;
            for (int i = 0; i < varNodes.size(); ++i) {
                if (mtr[r][varNodes[i].first.first].empty())
                    continue;
                ls.push_back((mtr[r][varNodes[i].first.first][varNodes[i].first.second] - varNodes[i].second + circulant) % circulant);
            }
            sort(ls.begin(), ls.end());
            for (int i = 0; i < ls.size(); ++i) {
                if ((i > 0) && (ls[i] == ls[i - 1]))
                    continue;
                if ((i + 1 < ls.size()) && (ls[i] == ls[i + 1]))
                    continue;
                ++res;
            }
        }
        return res;
    }
    

    long long gcd(long long a, long long b) {
        return b ? gcd(b, a % b) : a;
    }

    void gcd(long long a, long long b, long long& x, long long& y) {
        if (b == 0) {
            x = 1, y = 0;
            return;
        }
        gcd(b, a % b, x, y);
        long long q = a / b;
        long long xx = y;
        y = x - y * q;
        x = xx;
    }

    long long inverse(long long a, long long m) {
        long long x, y;
        gcd(a, m, x, y);
        return ((x % m) + m) % m;
    }

    bool lexMin(const Cycle& cycle) {
        int len = cycle.cycle.size();
        for (int i = 1; i < len; ++i) {
            if (cycle.cycle[i] != cycle.cycle[0])
                continue;
            if ((i & 1) == 0) {
                for (int j = i + 1, k = 1; k < len; ++j, ++k) {
                    if (j == len)
                        j = 0;
                    if (cycle.cycle[k] < cycle.cycle[j])
                        break;
                    if (cycle.cycle[j] < cycle.cycle[k])
                        return false;

                }
            }
            else {
                for (int j = i - 1, k = 1; k < len; --j, ++k) {
                    if (j == -1)
                        j = len - 1;
                    if (cycle.cycle[j] < cycle.cycle[k])
                        return false;
                    if (cycle.cycle[k] < cycle.cycle[j])
                        break;
                }
            }
        }
        return true;
    }

    void processCycle(const Cycle& cycle, vector<int>& numberOfCycles, const vector<vector<vector<int> > >& a, long long circulant, bool allLiftVals = 1, long long mult = 1) {
        mult = cycle.uniqueNodes;
        Tiii cur = cycle.cycle[0];
        long long c = 0, b = 0;//cx=b mod circulant
        bool one = true;
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            if (cur == cycle.cycle[i]) {
                if (i & 1)
                    ++c;
                else
                    --c;
                if (i)
                    one = false;
            }
            else {
                if (i & 1)
                    b -= a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
                else
                    b += a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
            }
        }
        if (!one) {
            if (!lexMin(cycle))
                return;
        }
        c = ((c % circulant) + circulant) % circulant;
        b = ((b % circulant) + circulant) % circulant;
        if (c == 0) {
            if (b != 0)
                return;
            for (int i = 0; i < numberOfCycles.size(); ++i)
                numberOfCycles[i] += mult;
            return;
        }
        long long d1 = gcd(c, circulant);
        if ((b % d1) != 0)
            return;
        c /= d1, b /= d1, circulant /= d1;
        long long x = (inverse(c, circulant) * b) % circulant;
        if (allLiftVals) {
            for (int i = 0; i < d1; ++i) {
                numberOfCycles[x + circulant * i] += mult;
            }
        }
        else {
            for (int i = 0; i < d1; ++i) {
                int pos = posByLiftVals[x + circulant * i];
                if (pos != -1)
                    numberOfCycles[pos] += mult;
            }
        }
    }

    void processCycle(const Cycle& cycle, double cycleCost, vector<double>& costOfCycles, const vector<vector<vector<int> > >& a, long long circulant, bool allLiftVals = 1, long long mult = 1) {
        mult = cycle.uniqueNodes;
        Tiii cur = cycle.cycle[0];
        long long c = 0, b = 0;//cx=b mod circulant
        bool one = true;
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            if (cur == cycle.cycle[i]) {
                if (i & 1)
                    ++c;
                else
                    --c;
                if (i)
                    one = false;
            }
            else {
                if (i & 1)
                    b -= a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
                else
                    b += a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third];
            }
        }
        if (!one) {
            if (!lexMin(cycle))
                return;
        }
        c = ((c % circulant) + circulant) % circulant;
        b = ((b % circulant) + circulant) % circulant;
        if (c == 0) {
            if (b != 0)
                return;
            for (int i = 0; i < costOfCycles.size(); ++i)
                costOfCycles[i] += mult * cycleCost;
            return;
        }
        long long d1 = gcd(c, circulant);
        if ((b % d1) != 0)
            return;
        c /= d1, b /= d1, circulant /= d1;
        long long x = (inverse(c, circulant) * b) % circulant;
        if (allLiftVals) {
            for (int i = 0; i < d1; ++i) {
                costOfCycles[x + circulant * i] += mult * cycleCost;
            }
        }
        else {
            for (int i = 0; i < d1; ++i) {
                int pos = posByLiftVals[x + circulant * i];
                if (pos != -1)
                    costOfCycles[pos] += mult * cycleCost;
            }
        }
    }
};






