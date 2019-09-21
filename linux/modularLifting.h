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
#include<string>
#include<sstream>
#include<iostream>
#include"regularLDPC.h"

using std::string;
using std::vector;
using std::cout;
using std::cin;
using std::endl;

//std::string toStr(long long x) {
//    std::stringstream s;
//    s << x;
//    return s.str();
//}


vector<vector<int> > modularLifting(vector<vector<int> > protograph, long long start, long long end, long long numberOfIterationsForFirstMatrix, long long numberOfIterationsForOneStep, 
    long long seed, bool RA, bool leftRA, bool EIRA, bool leftEIRA) {
    //freopen("in.txt", "w", stdout);
    string inputFile = "inputForModLifting.txt";
    string outputFile = "out.txt";
    FILE * in = fopen(inputFile.c_str(), "w");
    fprintf(in, "%d\t%d\n", protograph[0].size(), protograph.size());
    for (int i = 0; i < protograph.size(); ++i) {
        for (int j = 0; j < protograph[i].size(); ++j)
            fprintf(in, "%d\t", protograph[i][j]);
        fprintf(in, "\n");
    }
    fclose(in);
    string command = "modularLifting.exe -seed " + toStr(seed) + " -limits " + toStr(start) + " " + toStr(end) 
        + " -numberOfAttemptsForFirstMatrix " + toStr(numberOfIterationsForFirstMatrix) + " -numberOfAttemptsForOneStep " + toStr(numberOfIterationsForOneStep) + 
        " -file inputForModLifting.txt";

    if (EIRA)
        command += " -EIRA";
    if (RA)
        command += " -RA";
    if (leftEIRA)
        command += " -leftEIRA";
    if (leftRA)
        command += " -leftRA";
    //command += "\"";
    std::cerr << command << endl;
    system(command.c_str());
    int var, check, circ;
    FILE * out = fopen(outputFile.c_str(), "r");
    fscanf(out, "%d %d %d", &var, &check, &circ);
    vector<vector<int> > mtr(check, vector<int>(var));
    for (int i = 0; i < check; ++i) {
        for (int j = 0; j < var; ++j) {
            fscanf(out, "%d", &mtr[i][j]);
        }
    }
    fclose(out);
    system("del inputForModLifting.txt");
    //std::cerr << "del inputForModLifting.txt" << endl;
    system("del out.txt");

    return mtr;

    

}
