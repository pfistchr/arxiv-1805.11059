pxy = {{1 / 10000, 9997 / 60000, 9997 / 60000},
       {9997 / 60000, 1 / 10000, 9997 / 60000},
       {9997 / 60000, 9997 / 60000, 1 / 10000}};

rateA = 3898 / 2^17;
rateB = 3984 / 2^17;

upperboundA = 58593464420737815 / 2^56;
upperboundB = 58382556630811219 / 2^56;

rxyA = {{16^^00d18e2d53dba4 / 2^56, 16^^6c6ebcb6c6ea40 / 2^56, 16^^6c6ebcb6c6ea40 / 2^56},
        {16^^006ff71d804e2a / 2^56, 16^^03d405476786bd / 2^56, 16^^0ee47fcda75307 / 2^56},
        {16^^006ff71d804e2a / 2^56, 16^^0ee47fcda75307 / 2^56, 16^^03d405476786bd / 2^56}};

rxyB = {{16^^0184ae0a6be14a / 2^56, 16^^35ba25f4e1fd7f / 2^56, 16^^870eb8aa072ec5 / 2^56},
        {16^^02598735ff8940 / 2^56, 16^^057e6f74c876f3 / 2^56, 16^^35ba25f4e1fd7f / 2^56},
        {16^^00422176958a36 / 2^56, 16^^02598735ff8940 / 2^56, 16^^0184ae0a6be14a / 2^56}};

Print[Min[pxy] > 0];
Print[Total[pxy, 2] == 1];

Print[Min[rxyA] > 0];
Print[Total[rxyA, 2] == 1];

Print[Min[rxyB] > 0];
Print[Total[rxyB, 2] == 1];

rxryA = rxyA . {{1}, {1}, {1}} . {{1, 1, 1}} . rxyA;
rxryB = rxyB . {{1}, {1}, {1}} . {{1, 1, 1}} . rxyB;

Print[Total[rxyA * Log[rxyA / rxryA], 2] <= rateA];
Print[Total[rxyB * Log[rxyB / rxryB], 2] <= rateB];

Print[Total[rxyA * Log[rxyA / pxy], 2] <= upperboundA];
Print[Total[rxyB * Log[rxyB / pxy], 2] <= upperboundB];

Print["finish"];
