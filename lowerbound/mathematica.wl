pxy = {{1 / 10000}, {9997 / 60000}, {9997 / 60000},
       {9997 / 60000}, {1 / 10000}, {9997 / 60000},
       {9997 / 60000}, {9997 / 60000}, {1 / 10000}};

rate = 3941 / 2^17;
lowerbound = 58488010525784883 / 2^56;

verify[alpha_, beta_, qxy_] := Print[-(Log[Norm[Flatten[(pxy ^ alpha) + beta], 1 / alpha] - Min[(qxy ^ (1 - alpha)) . beta]] + (1 - alpha) * rate) / alpha >= lowerbound];

Print[Min[pxy] > 0];
Print[Total[pxy, 2] == 1];

If[10*29^(7151557967478901/9007199254740992) < 145, <<"input.wl", Print["please use Mathematica >= 12.0.0"]];
