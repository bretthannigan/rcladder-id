% [Calzavara2021] Example 1

S_true = [-4 1 1 2; 1 -1 0 0; 1 0 -2 1; 2 0 1 -3];
C = [1 0 0 0; 0 1 0 0; 0 0 1 0];
Ahat = [-10 -4 -23 4; 1 -1 3 -1; 3 1 7 -2; 1 -1 3 -6];
Chat = [1 0 3 -1; 0 1 0 0; 0 -1 1 0];
sys = ss(Ahat, [], Chat, []);

[G] = Calzavara2021(sys, C)