syms l1 l2 l3 d1 d2 d3;
syms y1 y2  r1 r2 r3;
syms alpha1 n1 n2 alpha2; 
syms x1 x2 x22 x3 x4;
syms dl1 dl2 dl3 dl4; %Dicken der Linsen
syms Rl1 Rl2 Rl3 Rl4; %Radien der Linsen

x1 = 80 - 2 - 1.5;
x2 = 20 * (1 + 1/sqrt(2)) - 3;
x22 = 20 * (1 + sqrt(2)) - 3;
x31x32 = 20 * (1 + sqrt(2)) - 20 * (1 + 1/sqrt(2)); %Distanz zwischen x'3 und x''3
dl1 = 4;
dl2 = 3;
dl3 = 3;
dl4 = 3;
Rl1 = 78.672364;
Rl2 = 39.603940;
Rl3 = 19.013503;
Rl4 = 7.804226;
n1 = 1;
n2 = 1.5;
% Transfermatrix Linse1
BL1e = [1, 0;(1/Rl1)*((n1/n2)-1), n1/n2];
BL1a = [1, 0;(-1/Rl1)*((n2/n1)-1), n2/n1];
TlL1 = [1 ,dl1;0,1];
Tl1 = BL1a * TlL1 * BL1e 

% Transfermatrix Linse 2
BL2e = [1, 0;(1/Rl2)*((n1/n2)-1), n1/n2];
BL2a = [1, 0;(-1/Rl2)*((n2/n1)-1), n2/n1];
TlL2 = [1 ,dl2;0,1];
Tl2 = BL2a * TlL2 * BL2e 

% Transfermatrix Linse 3
BL3e = [1, 0;(1/Rl3)*((n1/n2)-1), n1/n2];
BL3a = [1, 0;(-1/Rl3)*((n2/n1)-1), n2/n1];
TlL3 = [1 ,dl3;0,1];
TlL3B = [1,dl3/2;0,1];
Tl3 = BL3a * TlL3 * BL3e
Tl3B = TlL3B * BL3e; % Transfermatrix fuer die Halbe Linse 3 zur berechnung der Vergroesserung

% Transfermatrix Linse 4
BL4e = [1, 0;(1/Rl4)*((n1/n2)-1), n1/n2];
BL4a = [1, 0;(-1/Rl4)*((n2/n1)-1), n2/n1];
TlL4 = [1 ,dl4;0,1];
Tl4 = BL4a * TlL4 * BL4e

% Transfermatrizen fuer Distanzen zwischen den Linsen
Tx1 = [1,x1;0,1]
Tx2 = [1,x2;0,1]
Tx22 = [1,x22;0,1]
Tx3 = [1,x3;0,1]

% Berechnung des Abstandes x3 (Zwischen Linse 3 und 4), damit system scharf abbildet.
T = Tl4 * Tx3 * Tl3 * Tx2 * Tl2 * Tx1 * Tl1;
T2 = Tl4 * Tx3 * Tl3 * Tx22 * Tl2 * Tx1 * Tl1;
vpa(T);
vpa(T2);
eqd1 = T(2,1);
eqd2 = T2(2,1);
Dx2 = solve(eqd1, x3);
Dx3 = solve(eqd2, x3);

% Ausrechnen der Absoluten Position der Linse 4
vpa(Dx2);
vpa(Dx3);
Dtotal1 = Dx2 + x1 + x2 + 9.5;
Dtotal2 = Dx3 + x1 + x22 + 9.5;
pos1 = vpa(Dtotal1)
Pos2 = vpa(Dtotal2)

% Transfermatrix des gesmaten Systems
x3 = Dx2;
Tx3 = [1,x3;0,1];
T = Tl4 * Tx3 * Tl3 * Tx2 * Tl2 * Tx1 * Tl1;

%Berechung der Vergroesserung 1
syms Ytest;
Tx31x32 = [1,x31x32-1.5;0,1];
vpa(T);
Tv1 = Tx31x32 * Tl3 * Tx2 * Tl2 * Tx1 * Tl1;
vpa(Tv1);
light = [Ytest, alpha1]';
lightOut = Tv1 * light
eq = lightOut(1);
winkl = solve(eq,alpha1) / Ytest;
vpa(winkl);
light2 = [alpha2/winkl,alpha2]';
vpa(light2);
V1 = T * light2;
vpa(V1);
VergrFrac1 = V1(2)/alpha2';
Vergr1 = vpa(VergrFrac1)

% Berechnung der Vergroesserung 2
x3 = Dx3;
Tx3 = [1,x3;0,1];

T2 = Tl4 * Tx3 * Tl3 * Tx22 * Tl2 * Tx1 * Tl1;
vpa(T2);
Tv2 = Tl3B * Tx22 * Tl2 * Tx1 * Tl1;
vpa(Tv2);
light = [Ytest, alpha1]';
lightOut = Tv2 * light;
eq = lightOut(1);
winkl = solve(eq,alpha1) / Ytest;
vpa(winkl);
light2 = [alpha2/winkl,alpha2]';
vpa(light2)
V2 = T2 * light2;
vpa(V2)
VergrFrac2 = V2(2)/alpha2';
Vergr2 = vpa(VergrFrac2)

% Lichtsrathl zum zeigen, dass System das Bild nicht "Verkehrt" abbildet.
licht = [1,0]'
outlicht = T * licht
