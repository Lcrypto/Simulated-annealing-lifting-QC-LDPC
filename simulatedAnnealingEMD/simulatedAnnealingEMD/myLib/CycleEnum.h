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
#include<vector>
using namespace std;
struct entry {
    int r, c, id;
    int index, rowIndex, colIndex;
    entry() {}
    entry(int _r, int _c, int _id) : r(_r), c(_c), id(_id) {}
    entry(int _r, int _c, int _id, int _index, int _rowIndex, int _colIndex) : r(_r), c(_c), id(_id), index(_index), rowIndex(_rowIndex), colIndex(_colIndex) {}

};
class CycleEnum {
    

    vector<entry> entries;
    vector<vector<entry> > entriesByRow, entriesByCol;
    vector<vector<vector<entry> > > entriesByRowAndCol;

    void enumerateEntries() {
        entriesByRow.resize(checkNodes);
        entriesByCol.resize(variableNodes);
        entriesByRowAndCol.resize(checkNodes);
        for (int i = 0; i < checkNodes; ++i)
            entriesByRowAndCol[i].resize(variableNodes);
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                for (int i = 0; i < protograph[r][c]; ++i) {
                    entry x(r, c, i, entries.size(), entriesByRow[r].size(), entriesByCol[c].size());
                    entries.push_back(x);
                    entriesByRow[r].push_back(x);
                    entriesByCol[c].push_back(x);
                    entriesByRowAndCol[r][c].push_back(x);
                }
            }
        }
    }
public:
    vector<entry> cycle;
    int sizeOfCycle, checkNodes, variableNodes;
    vector<vector<int> > protograph;

    CycleEnum(int _sizeOfCycle, const vector<vector<int> >& _protograph) {
        sizeOfCycle = _sizeOfCycle;
        protograph = _protograph;
        checkNodes = protograph.size();
        variableNodes = protograph[0].size();
        cycle.resize(sizeOfCycle);
        enumerateEntries();
    }
    bool next(int step = -1) {
        if (step == -1)
            step += sizeOfCycle;
        if (step == sizeOfCycle - 1) {
            int r = cycle[step - 1].r;
            int c = cycle[0].c;
            for (int i = cycle[step].id + 1; i < protograph[r][c]; ++i) {
                if ((c == cycle[step - 1].c) && (i == cycle[step - 1].id))
                    continue;
                if ((r == cycle[0].r) && (i == cycle[0].id))
                    continue;
                cycle[step] = entriesByRowAndCol[r][c][i];
                return 1;
            }
            return next(step - 1);
        }
        if (step == 0) {
            for (size_t i = cycle[0].index + 1; i < entries.size(); ++i) {
                cycle[0] = entries[i];
                if (init(1))
                    return 1;
            }
            return 0;
        }
        if (step & 1) {
            int r = cycle[step - 1].r;
            for (size_t i = cycle[step].rowIndex + 1; i < entriesByRow[r].size(); ++i) {
                if (cycle[step - 1].rowIndex == i)
                    continue;
                cycle[step] = entriesByRow[r][i];
                if (init(step + 1))
                    return 1;
            }
            return next(step - 1);
        }
        int c = cycle[step - 1].c;
        for (size_t i = cycle[step].colIndex + 1; i < entriesByCol[c].size(); ++i) {
            if (cycle[step - 1].colIndex == i)
                continue;
            cycle[step] = entriesByCol[c][i];
            if (init(step + 1))
                return 1;
        }
        return next(step - 1);
    }
    bool init(int step) {
        if (step == sizeOfCycle)
            return 1;
        if (step + 1 == sizeOfCycle) {
            int r = cycle[step - 1].r;
            int c = cycle[0].c;
            for (int i = 0; i < protograph[r][c]; ++i) {
                if ((c == cycle[step - 1].c) && (i == cycle[step - 1].id))
                    continue;
                if ((r == cycle[0].r) && (i == cycle[0].id))
                    continue;
                cycle[step] = entriesByRowAndCol[r][c][i];
                return 1;
            }
            return 0;
        }
        if (step & 1) {
            int r = cycle[step - 1].r;
            for (size_t i = 0; i < entriesByRow[r].size(); ++i) {
                if (i == cycle[step - 1].rowIndex)
                    continue;
                cycle[step] = entriesByRow[r][i];
                if (init(step + 1))
                    return 1;
            }
            return 0;
        }
        int c = cycle[step - 1].c;
        for (size_t i = 0; i < entriesByCol[c].size(); ++i) {
            if (i == cycle[step - 1].colIndex)
                continue;
            cycle[step] = entriesByCol[c][i];
            if (init(step + 1))
                return 1;
        }
        return 0;
    }
    bool init() {//gen lex-min cycle
        size_t id = 0;
        while (id < entries.size()) {
            cycle[0] = entries[id];
            if (init(1))
                break;
            ++id;
        }
        return id < entries.size();
    }
    bool init(int rowId, int colId, int id) {//gen lex-min cycle with this entry
        size_t index = 0;
        while (index < entriesByRowAndCol[rowId][colId].size()) {
            if (entriesByRowAndCol[rowId][colId][index].id == id) {
                cycle[0] = entriesByRowAndCol[rowId][colId][index];
                return init(1);
            }
            ++index;
        }
        return 0;
    }
};
