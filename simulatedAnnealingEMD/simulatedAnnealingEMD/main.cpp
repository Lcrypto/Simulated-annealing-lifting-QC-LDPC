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


#include".\myLib\emdOptimization.h"
#include".\myLib\irregularLDPC.h"




int preprocess(vector<vector<int> >& protograph) {
    int m = protograph.size();
    int n = protograph[0].size();
    if (n < m)
        return -1;
    vector<int> degs(m - 1, 0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m - 1; ++j)
            degs[j] += protograph[i][j];
    int ones = count(degs.begin(), degs.end(), 1), twos = count(degs.begin(), degs.end(), 2);
    if (ones + twos != m - 1)
        return -1;
    int res = twos;
    vector<pii> a(twos);
    for (int colId1 = ones; colId1 < m - 1; ++colId1) {
        if (degs[colId1] != 2) {
            return -1;
        }
    }
    int first = -1;
    int firstCol = -1;
    for (int r = 0; r < twos; ++r) {
        int d = 0;
        for (int j = ones; j < m - 1; ++j) {
            d += protograph[r][j];
            if (protograph[r][j])
                firstCol = j;
        }
        if (d == 1) {
            first = r;
            break;
        }
    }
    if (first == -1)
        return -1;
    swap(protograph[first], protograph[0]);
    for (int i = 0; i < m; ++i) {
        swap(protograph[i][ones], protograph[i][firstCol]);
    }
    for (int c = ones; c < m - 2; ++c) {
        int newr = -1;
        for (int r = c - ones + 1; r < m; ++r) {
            if (protograph[r][c]) {
                newr = r;
                break;
            }
        }
        /*for (int rr = 0; rr < res; ++rr) {
            for (int cc = ones; cc < m - 1; ++cc) {
                cerr << protograph[rr][cc] << " ";
            }
            cerr << endl;
        }
        cerr << endl;*/

        if (newr == -1)
            return -1;
        swap(protograph[newr], protograph[c - ones + 1]);
        int colIdSw = -1;
        for (int j = c + 1; j < m - 1; ++j) {
            if (protograph[c - ones + 1][j]) {
                colIdSw = j;
                break;
            }
        }
        if (colIdSw == -1)
            return -1;
        for (int i = 0; i < m; ++i) {
            swap(protograph[i][colIdSw], protograph[i][c + 1]);
        }

    }
    return res;
}

int gen(int checkNodes, int variableNodes, vector<vector<int> >& protograph, ll circulant, ll targetGirth, ll targetEmd, ll upGirth,
    vector<int> liftVals, vector<vector<vector<int> > >& mtr, const vector<vector<int> >& fixed2d, 
    const vector<vector<int> >& fixedVals, bool cycleCost = 0) {

    vector<vector<vector<int> > > fixed(protograph.size(), vector<vector<int> > (protograph[0].size()));
    for (int r = 0; r < protograph.size(); ++r) {
        for (int c = 0; c < protograph[r].size(); ++c) {
            fixed[r][c].assign(protograph[r][c], fixed2d[r][c]);
        }
    }

    
    mtr.assign(checkNodes, vector<vector<int> >(variableNodes));
    for (int c = 0; c < variableNodes; ++c) {
        for (int r = 0; r < checkNodes; ++r) {
            mtr[r][c].resize(protograph[r][c]);
            while (true) {
                for (int id = 0; id < protograph[r][c]; ++id) {
                    if (fixed2d[r][c])
                        mtr[r][c][id] = fixedVals[r][c];
                    else
                        mtr[r][c][id] = liftVals[getRand(liftVals.size())];
                }
                bool ok = 1;
                for (int i = 0; i < protograph[r][c]; ++i) {
                    for (int j = 0; j < i; ++j) {
                        if (mtr[r][c][i] == mtr[r][c][j])
                            ok = 0;
                    }
                }
                if (ok)
                    break;
            }
        }
    }
    
    

    emdOpt opt(circulant, upGirth, targetGirth, targetEmd, mtr);

    if (cycleCost) {
        opt.annealEmdWithFixedAndCycleCostAndInnerCode(fixed);
        mtr = opt.getMatrix();
        return 1;
        
    }
    if (liftVals.size() == circulant) {
        if (opt.annealEmdWithFixed(fixed)) {
            mtr = opt.getMatrix();
            return 1;
        }
    }
    else {
        if (opt.annealEmdWithFixed(fixed, liftVals)) {
            mtr = opt.getMatrix();
            return 1;
        }
    }
    return 0;
}

bool checkInput(ll GIRTH, ll UP_GIRTH, ll EMD, ll SEED, ll CIRCULANT_SIZE, ll DESIRED_NUMBER_OF_MATRICES, string INPUT_FILENAME, bool regular) {
    bool validInput = 1;
    if (UP_GIRTH < 4) {
        validInput = 0;
        cerr << "wrong upGirth\n";
        cerr << endl;
    }

    if ((GIRTH < 4) || ((GIRTH > 0) && (GIRTH & 1))) {
        validInput = 0;
        cerr << "girth must be >= 4 and even\n";
        cerr << endl;
    }
    if (SEED < 0) {
        validInput = 0;
        cerr << "seed must be non-negative\n";
        cerr << "example: -seed 123\n";
        cerr << endl;
    }
    if (CIRCULANT_SIZE < 0) {
        validInput = 0;
        cerr << "circulant must be non-negative\n";
        cerr << "example: -circulant 666\n";
        cerr << endl;
    }
    if (DESIRED_NUMBER_OF_MATRICES < 0) {
        validInput = 0;
        cerr << "numberOfMatrices must be non-negative\n";
        cerr << "example: -numberOfMatrices 15\n";
        cerr << endl;
    }
    if ((INPUT_FILENAME == "") ^ (regular)) {
        validInput = 0;
        cerr << "wrong protograph\n";
        cerr << "example1: -regular 12 4\n";
        cerr << "example2: -file input.txt\n";
        cerr << endl;
    }
    return validInput;
}

int main(int argc, char* argv[]) {
    //Initialization
    bool validInput = 1;
    bool regular = 0;
    ll SEED = -1;
    vector<vector<int> > PROTOGRAPH, fixed;
    ll CIRCULANT_SIZE = -1;
    ll VARIABLE_NODES;
    ll CHECK_NODES;
    ll DESIRED_NUMBER_OF_MATRICES = -1;
    string INPUT_FILENAME = "";
    string LIFTING_CONSTR_FILENAME = "";
    ll UP_GIRTH = 0;
    ll EMD = 0;
    ll GIRTH = 0;
    bool LIFTING_MASK = 0;
    bool LIFTING_VAL = 0;
    bool CYCLE_COST = 0;
    bool FIXED = 0;
    bool LIFTED_INPUT = 0;

    //reading parameters from cmd
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "-seed") {
            validInput = validInput && toUnsignedInt(argv[i + 1], SEED);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-upGirth") {
            validInput = validInput && toUnsignedInt(argv[i + 1], UP_GIRTH);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-girth") {
            validInput = validInput && toUnsignedInt(argv[i + 1], GIRTH);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-circulant") {
            validInput = validInput && toUnsignedInt(argv[i + 1], CIRCULANT_SIZE);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-numberOfMatrices") {
            validInput = validInput && toUnsignedInt(argv[i + 1], DESIRED_NUMBER_OF_MATRICES);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-file") {
            INPUT_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-emd") {
            validInput = validInput && toUnsignedInt(argv[i + 1], EMD);
            ++i;
            continue;
        }
        if (string(argv[i]) == "-regular") {
            validInput = validInput && toUnsignedInt(argv[i + 1], VARIABLE_NODES) && (i + 2 < argc) && toUnsignedInt(argv[i + 2], CHECK_NODES);
            i += 2;
            regular = 1;
            continue;
        }
        if (string(argv[i]) == "-CYCLE_COST") {
            CYCLE_COST = 1;
            continue;
        }
        if (string(argv[i]) == "-LIFTING_MASK") {
            LIFTING_MASK = 1;
            LIFTING_CONSTR_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-LIFTING_VAL") {
            LIFTING_VAL = 1;
            LIFTING_CONSTR_FILENAME = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-FIXED") {
            FIXED = 1;
            continue;
        }
        if (string(argv[i]) == "-LIFTED_INPUT") {
            LIFTED_INPUT = 1;
            continue;
        }
    }
    if ((!checkInput(GIRTH, UP_GIRTH, EMD, SEED, CIRCULANT_SIZE, DESIRED_NUMBER_OF_MATRICES, INPUT_FILENAME, regular)) || (!validInput))
        return 0;

    if (UP_GIRTH == 0) {
        UP_GIRTH = GIRTH;
        if (EMD > 0)
            UP_GIRTH += 2;
    }

    srand(SEED);


    //protograph initialization
    if (!regular) {
        freopen(INPUT_FILENAME.c_str(), "r", stdin);
        cin >> VARIABLE_NODES >> CHECK_NODES;
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES));
        for (int i = 0; i < CHECK_NODES; ++i) {
            for (int j = 0; j < VARIABLE_NODES; ++j) {
                cin >> PROTOGRAPH[i][j];
            }
        }
        fixed.assign(CHECK_NODES, vector<int>(VARIABLE_NODES, 0));
        if (FIXED) {
            for (int i = 0; i < CHECK_NODES; ++i) {
                for (int j = 0; j < VARIABLE_NODES; ++j) {
                    cin >> fixed[i][j];
                }
            }
        }

        fclose(stdin);
    }
    else {
        PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES, 1));
    }

    vector<vector<int> > fixedVals(CHECK_NODES, vector<int>(VARIABLE_NODES, 0));
    if (LIFTED_INPUT) {
        for (int i = 0; i < CHECK_NODES; ++i) {
            for (int j = 0; j < VARIABLE_NODES; ++j) {
                fixedVals[i][j] = PROTOGRAPH[i][j];
                PROTOGRAPH[i][j] = (PROTOGRAPH[i][j] >= 0);
            }
        }
    }
    vector<int> liftMask(CIRCULANT_SIZE, 1);
    vector<int> liftVals;
    if (LIFTING_MASK) {
        freopen(LIFTING_CONSTR_FILENAME.c_str(), "r", stdin);
        for (int i = 0; i < CIRCULANT_SIZE; ++i) {
            cin >> liftMask[i];
        }
        fclose(stdin);
        
    }
    if (LIFTING_VAL) {
        freopen(LIFTING_CONSTR_FILENAME.c_str(), "r", stdin);
        int numVals;
        cin >> numVals;
        liftVals.resize(numVals);
        for (int i = 0; i < numVals; ++i) {
            cin >> liftVals[i];
        }
        sort(liftVals.begin(), liftVals.end());
        fclose(stdin);
    }
    else {
        for (int i = 0; i < liftMask.size(); ++i) {
            if (liftMask[i])
                liftVals.push_back(i);
        }
    }
    

    //folder and filenames generation
    string folderName = toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH) + "upGirth" + toStr(UP_GIRTH) + "emd" + toStr(EMD);
    string outputFilename = folderName + "/" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE) + "girth" + toStr(GIRTH) + "upGirth" + toStr(UP_GIRTH)
        + "emd" + toStr(EMD) + "seed" + toStr(SEED);
    if (regular)
        outputFilename += "regular_protograph_matrix";
    else
        outputFilename += "protograph_from_" + INPUT_FILENAME + "_matrix";





    time_t start = time(NULL);
    ll iterationCount = 0;
    ll successCount = 0;
    ll power = 1;


    //matrix generation for fixed emd -- START
    while (successCount < DESIRED_NUMBER_OF_MATRICES) {
        ++iterationCount;
        vector<vector<vector<int> > > a;
        int resGen = gen(CHECK_NODES, VARIABLE_NODES, PROTOGRAPH, CIRCULANT_SIZE, GIRTH, EMD, UP_GIRTH, liftVals, a, fixed, fixedVals, CYCLE_COST);
        if (resGen == -1) {
            return 0;
        }
        if (resGen) {
            if (successCount == 0)
                system(("mkdir " + folderName).c_str());
            freopen((outputFilename + toStr(iterationCount) + ".txt").c_str(), "w", stdout);
            ++successCount;
            cout << VARIABLE_NODES << "\t" << CHECK_NODES << "\t" << CIRCULANT_SIZE << endl;
            print(a);
            cout << endl;
            fclose(stdout);
            cerr << "girth = " << getGirth(a, PROTOGRAPH, CIRCULANT_SIZE) << endl;
            eprint(a);
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
        if (iterationCount == power) {
            power *= 2;
            cerr << iterationCount << " iterations\n";
            cerr << time(NULL) - start << " seconds\n";
            cerr << endl;
        }
    }

    //matrix generation for fixed emd -- END


    return 0;
}