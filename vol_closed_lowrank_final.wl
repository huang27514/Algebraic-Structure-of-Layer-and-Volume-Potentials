(* ::Package:: *)

(* ======================================================================
   vol_closed_lowrank.wl
   Volume potential below the curve by the EXPLICIT formulas of Sec. 3
   (vol_section_vertical_cut.tex; branch cut of the logarithm on the negative
   imaginary axis), specialized to
     (A1) b <= s(a) (targets below Gamma)  -> N_top = 0, c(w) = b,
     (A2) real f,  (A3) only the real part is needed:
   No particular solution and no Green's identity are used.

     D = {-1 <= x <= 1,  yb <= y <= s(x)},  yb = -4/3  (bottom side y = yb),
     vol(w) = Sum_{n,m} c_nm H_nm(w),  H_nm = Int_D e^{i(nx+my)} L(x+iy-w) dA,
     V(w)   = Re vol(w) = Re[ vol_SLP,top(w) + vol_branch(a,b) + vol_low-rank(w) ].

   vol_SLP,top (eq:volslp, closed form eq:volslpclosed):
     Int_{-1}^{1} Qt(x) L(P(x)-w) dx + i c00 w Int_{-1}^{1} L(P(x)-w) dx,
     Qt(x) = (1+i s'(x)) Sum' c_nm/(n+im) e^{i(n x + m s(x))} - i c00 (x + i s(x)),
     Int Qt L(P-w) = Qcal(1) Log(P(1)-w) - Q_Qcal(w) - Sum_j Qcal(w_j) Lambda(w_j)
                   = Qcal(1) Log(P(1)-w) - Sum_j Int_{-1}^{1} Qcal(x)/(x - w_j) dx
       (partial fractions P'/(P-w) = Sum_j 1/(x-w_j)). DOUBLE PRECISION version:
       Qcal is kept in the Chebyshev basis (Clenshaw); for the root w_j near [-1,1]
       the term is Int (Qcal(x)-Qcal(w_j))/(x-w_j) dx + Qcal(w_j) Lambda(w_j), for the
       far roots it is integrated directly by Gauss-Legendre. (The monomial form
       Q_Qcal(w) + Sum Qcal(w_j) Lambda(w_j) is unusable in double: at the far roots,
       |w_j| ~ 1.4-2.4, it loses 9 (degree 20) to 20+ (degree 60) digits.)
     Qcal(x) = Int_{-1}^{x} Qt;  same formula with Qcal = x + 1 for the c00 term.
   vol_branch (eq:volbranch, eq:Cnm, eq:C00), explicit in (a, b):
     c00 [pi (a+1)^2 + 2 pi i (a+1)(b-yb)]
     + Sum' c_nm/(n+im) [ 2 pi e^{ina} Int_{yb}^{b} e^{imy} dy
                          - 2 pi i e^{i m yb} Int_{-1}^{a} e^{inx} dx ].
   vol_low-rank (eq:vollowrank):
     Sum' c_nm/(n+im) [ -e^{i m yb} Phi_bot(n,w) - i e^{in} Phi_R(m,w) + i e^{-in} Phi_L(m,w) ]
     + i c00 Int_{-1}^{1} (x+i yb-w) Log(x+i yb-w) dx - c00 Int_{-1}^{1} (s(x)-yb) dx,
     Phi_bot = Int_{-1}^{1} e^{inx} Log(x+i yb-w) dx,
     Phi_R   = Int_{yb}^{s(1)} e^{imy} Log(1+iy-w) dy,
     Phi_L   = Int_{yb}^{s(-1)} e^{imy} [Log(1-iy+w) + i pi] dy.
   Taylor expansion about ONE center w0 = -i/3 with the moments
     A^(bot)_{n,0} = Int e^{inx} Log(x+i yb-w0), A^(bot)_{n,k} = Int e^{inx} (x+i yb-w0)^-k,
     A^(R)_{m,0}  = Int e^{imy} Log(1+iy-w0),   A^(R)_{m,k}  = Int e^{imy} (1+iy-w0)^-k,
     A^(L)_{m,0}  = Int e^{imy}[Log(1-iy+w0)+i pi], A^(L)_{m,k} = Int e^{imy}(-1+iy-w0)^-k,
     Phi(w) = A_0 - Sum_{k>=1} A_k (w-w0)^k / k,
   and, for the c00 bottom term with d = x+i yb-w0, L0 = Log(d):
     (d - dw) Log(d - dw) = d L0 - (L0 + 1) dw + Sum_{k>=2} d^(1-k) dw^k /(k(k-1)).

   Test problem:  f = 1 + cos(2 pi x) cos(pi y),  s = x^2/3 + x^3/8 - x^4/6,
   c00 = 1, c_{+-2pi,+-pi} = 1/4 (non-integer frequencies are allowed).
   Targets: the 101 x 100 grid of the layer-potential tests (Figs. 5/7);
   errors are evaluated and reported ONLY on the region
     R = [-1/3, 1/3] x [-2/3, min(0, s(x))]   (9963 of the 10100 nodes),
   nodes outside R are skipped (exported with inRegion = 0).
   On R, max|w - w0| = sqrt(2)/3 = 0.471; dist(w0, bottom side) = 1. Phi_*
   continue analytically across the straight sides and are singular only at
   the corners -1+i yb, 1+i yb, 1+i s(1), -1+i s(-1); the nearest is
   -1 + i s(-1) = -1 + i/24 (distance 1.068 from w0); asymptotic rate
   rho = 0.471/1.068 = 0.441.

   Reference ("direct"): at every node,
     Re vol_SLP,top by adaptive quadrature with the exact density (no closed form)
     + Re vol_branch + Re vol_low-rank by direct summation of Phi_* over
     Gauss-Legendre nodes on the flat sides (no Taylor expansion).
   vol_branch cancels in V_K - V_ref; it is verified by a full 2D volume
   quadrature at a few nodes (Step 7).

   ALL computations are in machine double precision (closed form, Chebyshev
   interpolation, moments, Taylor sums, branch term, and the reference quadratures).
   Expected (Python double-precision prototype of exactly these formulas):
     top closed form vs adaptive quadrature: <= 3.4e-15 (degree 40)
     max |error| on R, Taylor about -i/3 (yb = -4/3):  K = 6: 3.1e-6,
       K = 12: 5.0e-9,  K = 20: 1.3e-12,  K = 25: 1.2e-14,
       K >= 30: ~4.4e-15 (double-precision floor; median ~4e-16)
     total vs 2D volume quadrature (double): <= 9.1e-15; V at node (i,j):
       (51,1) -0.6885360837  (11,51) -1.1558870554  (91,100) -1.3919234676
       (1,40) -0.9605829786  (76,31) -1.0845949197
   Not yet run in Mathematica.

   Outputs (vol_closed_out/), MAT files:
     zdataxy_vol.mat : X, Y (nx x ny), Nx, Ny, w0 ([Re Im]), rhoEst, inRegion (nx x ny, 0/1),
                       xReg, yRegMin
     zerror_vol.mat  : logErrAll (nx x ny x (Kmax+1)), Klist, maxErrK, medErrK,
                       logErrK6, logErrK12, logErrK20, logErrTop, maxErrTop,
                       Vref, Vfast20 (nx x ny), inRegion, check2D (rows: i, j, V_fast(K=20), V_2D)
     (maxErrK, medErrK, maxErrTop are over R only; map entries outside R are 0)
   ====================================================================== *)

ClearAll["Global`*"];
$HistoryLength = 0;
Needs["NumericalDifferentialEquationAnalysis`"];   (* GaussianQuadratureWeights *)

(* ---------- 0. Parameters (ALL computations in machine double precision) ---------- *)
nCheb     = 40;      (* Chebyshev degree of Qt(x); tail at round-off in double *)
nGT       = 200;     (* Gauss-Legendre nodes on [-1,1] for the root terms     *)
rCut      = 3/2;     (* Bernstein-ellipse parameter separating near/far roots *)
nGL       = 300;     (* Gauss-Legendre nodes per flat side                     *)
Kmax      = 40;
w0c       = -I/3;    (* single expansion center                                *)
yBot      = -4/3;    (* bottom side of D: y in [yBot, s(x)]                    *)
pgRef = 14; agRef = 15; maxRec = 30;      (* reference quadrature, machine precision *)
nintMethod = {"GlobalAdaptive", Method -> {"GaussKronrodRule", "Points" -> 15},
   "MaxErrorIncreases" -> 4000, "SymbolicProcessing" -> 0};
stride    = 1;       (* 5 for a quick test *)
do2DCheck = True;
idx2D     = {{51, 1}, {11, 51}, {91, 100}, {1, 40}, {76, 31}};   (* all in R *)
xReg = {-1/3, 1/3}; yRegMin = -2/3;      (* error region R: y in [yRegMin, min(0, s(x))] *)
errFloor  = 10^-17;

(* ---------- 1. Geometry ---------- *)
sc[x_]  := x^2/3 + x^3/8 - x^4/6;
scp[x_] := 2 x/3 + 3 x^2/8 - 2 x^3/3;
P[x_]   := x + I sc[x];
dP[x_]  := 1 + I scp[x];
pCoef = CoefficientList[P[xx], xx];          (* {0, 1, 2 I, 0, 5 I} *)
Jdeg  = Length[pCoef] - 1;

(* ---------- 2. Density and Fourier modes ---------- *)
fxy[x_, y_] := 1 + Cos[2 Pi x] Cos[Pi y];
c00   = 1;
modes = {{2 Pi, Pi, 1/4}, {2 Pi, -Pi, 1/4}, {-2 Pi, Pi, 1/4}, {-2 Pi, -Pi, 1/4}};  (* {n, m, c_nm} *)
Print["[Check] Fourier series - f = ",
  Simplify[ExpToTrig[c00 + Total[Table[md[[3]] Exp[I (md[[1]] xx + md[[2]] yy)], {md, modes}]]] - fxy[xx, yy]]];
Print["[Check] real f: c_{-n,-m} = conj(c_{n,m}) : ",
  And @@ Table[MemberQ[modes, {-md[[1]], -md[[2]], Conjugate[md[[3]]]}], {md, modes}]];

(* ---------- 3. Branch-cut part: explicit function of (a, b) ---------- *)
intExp[k_, lo_, hi_] := If[k === 0, hi - lo, (Exp[I k hi] - Exp[I k lo])/(I k)];
C00[a_, b_] := Pi (a + 1)^2 + 2 Pi I (a + 1) (b - yBot);
Cnm[n_, m_, a_, b_] := 2 Pi Exp[I n a] intExp[m, yBot, b] - 2 Pi I Exp[I m yBot] intExp[n, -1, a];
volBranch[a_, b_] = c00 C00[a, b] + Total[Table[md[[3]]/(md[[1]] + I md[[2]]) Cnm[md[[1]], md[[2]], a, b], {md, modes}]];
Print["[Branch] Re vol_branch(a,b) = ",
  FullSimplify[ComplexExpand[Re[volBranch[a, b]]], Element[{a, b}, Reals]]];

(* ---------- 4. Top single layer potential: closed form in double precision ---------- *)
QtExact[x_] := (1 + I scp[x]) Total[Table[md[[3]]/(md[[1]] + I md[[2]]) Exp[I (md[[1]] x + md[[2]] sc[x])], {md, modes}]] -
   I c00 P[x];

(* 4a. Chebyshev interpolant of Qt (machine precision), Qcal = Int_{-1}^x Qt in the
       Chebyshev basis, Clenshaw evaluation. The monomial basis is NOT used: in double
       precision its evaluation at the far roots of P(z) = w (|z| ~ 1.4-2.4) loses
       up to 20+ digits. *)
chebEval[c_, x_] := Module[{b1 = 0, b2 = 0, t},                  (* Clenshaw; x scalar or list *)
   Do[t = 2 x b1 - b2 + c[[k]]; b2 = b1; b1 = t, {k, Length[c], 2, -1}];
   x b1 - b2 + c[[1]]];
chebIntegrate[c_] := Module[{n = Length[c], b},                   (* Int_{-1}^x, Chebyshev coeffs *)
   b = ConstantArray[0. + 0. I, n + 1];
   b[[2]] += c[[1]];                                              (* Int T0 = T1      *)
   If[n >= 2, b[[3]] += c[[2]]/4];                                (* Int T1 = T2/4 + c *)
   Do[b[[k + 2]] += c[[k + 1]]/(2 (k + 1)); b[[k]] -= c[[k + 1]]/(2 (k - 1)), {k, 2, n - 1}];
   b[[1]] -= chebEval[b, -1.];                                    (* Qcal(-1) = 0 *)
   b];
chebX  = N[Table[Cos[Pi k/nCheb], {k, 0, nCheb}]];
fvals  = QtExact /@ chebX;
cosTab = N[Table[Cos[Pi j k/nCheb], {j, 0, nCheb}, {k, 0, nCheb}]];
endW   = Table[If[k == 0 || k == nCheb, 1/2, 1], {k, 0, nCheb}];
chebC  = N[(2/nCheb) endW (cosTab . (endW fvals))];
cQ     = chebIntegrate[chebC];                                    (* Qcal *)
cOne   = {1., 1.};                                                (* x + 1 = T0 + T1 *)
Print["[Cheb] degree ", nCheb, ", |c_{N-2..N}| = ", Abs[chebC[[-3 ;;]]]];
Print["[Cheb] max |Qt - interpolant| at 200 random x: ",
  Max[Table[With[{xr = RandomReal[{-1, 1}]}, Abs[QtExact[xr] - chebEval[chebC, xr]]], {200}]]];

(* 4b. Root terms. With P'/(P-w) = Sum_j 1/(x - w_j),
         Int Qcal P'/(P-w) dx = Sum_j I_j,  I_j = Int_{-1}^{1} Qcal(x)/(x - w_j) dx,
       which is Q_Qcal(w) + Sum_j Qcal(w_j) Lambda(w_j) redistributed over the roots.
       near root (Bernstein parameter <= rCut): singularity subtraction
         I_j = Int (Qcal(x) - Qcal(w_j))/(x - w_j) dx + Qcal(w_j) Lambda(w_j)
         (the difference quotient is a polynomial of degree nCheb: Gauss-Legendre
          with nGT >= nCheb/2 + 1 nodes integrates it exactly);
       far root: I_j by Gauss-Legendre directly (smooth integrand). *)
{xGT, wGT} = Transpose[GaussianQuadratureWeights[nGT, -1, 1]];
QxGT  = chebEval[cQ, xGT];
OneGT = xGT + 1.;
QcalAt1 = chebEval[cQ, 1.];
Lambda[a_] := -2 ArcTanh[1/a];                 (* = Log((a-1)/(a+1)), cut a in [-1,1] *)
bernstein[a_] := With[{r = a + Sqrt[a - 1] Sqrt[a + 1]}, Max[Abs[r], 1/Abs[r]]];
rootTerm[c_, fxG_, a_] := If[bernstein[a] > rCut,
   Total[wGT fxG/(xGT - a)],
   With[{fa = chebEval[c, a]}, Total[wGT (fxG - fa)/(xGT - a)] + fa Lambda[a]]];

(* Int Qt L(P-w) dx = Qcal(1) Log(P(1)-w) - Sum_j I_j;  L(P(1)-w) = Log(P(1)-w)
   since Re(P(1)-w) = 1-a > 0; no crossing term by (A1). *)
topClosed[wv_] := Module[{ww = N[wv], den, rts, zz, S1, S2},
   den = N[pCoef]; den[[1]] = den[[1]] - ww;
   rts = zz /. NSolve[Sum[den[[k]] zz^(k - 1), {k, Jdeg + 1}] == 0, zz];
   If[Length[rts] =!= Jdeg, Print["[topClosed] wrong root count at w = ", N[wv]]; Abort[]];
   S1 = Total[rootTerm[cQ, QxGT, #] & /@ rts];
   S2 = Total[rootTerm[cOne, OneGT, #] & /@ rts];
   (QcalAt1 Log[P[1] - ww] - S1) + I c00 ww (2 Log[P[1] - ww] - S2)];       (* complex *)

(* 4c. Reference: adaptive quadrature (machine precision) with the exact density *)
Ldown[z_?NumericQ] := Module[{t = Arg[z]}, Log[Abs[z]] + I If[t <= -Pi/2, t + 2 Pi, t]];
topIntegrand[t_?NumericQ, wv_] := Re[(QtExact[t] + I c00 wv) Ldown[P[t] - wv]];
nearestParam[a_, b_] := Module[{s0, sv},
   s0 = sv /. FindRoot[(sv - a) + (sc[sv] - b) scp[sv] == 0, {sv, N[a]},
       AccuracyGoal -> 10, PrecisionGoal -> 10];
   Clip[s0, {-1 + 10^-6, 1 - 10^-6}]];
topRef[wv_] := With[{ss = nearestParam[Re[wv], Im[wv]]},
   NIntegrate[topIntegrand[t, N[wv]], {t, -1, ss, 1},
     PrecisionGoal -> pgRef, AccuracyGoal -> agRef,
     MaxRecursion -> maxRec, Method -> nintMethod]];             (* real part *)

(* ---------- 5. Low-rank part ---------- *)
glB = GaussianQuadratureWeights[nGL, -1, 1];                     (* machine precision *)
glR = GaussianQuadratureWeights[nGL, N[yBot], N[sc[1]]];
glL = GaussianQuadratureWeights[nGL, N[yBot], N[sc[-1]]];
{xB, wB} = Transpose[glB]; {yR, wR} = Transpose[glR]; {yL, wL} = Transpose[glL];
w0n = N[w0c];
cS  = N[c00 Integrate[sc[x] - yBot, {x, -1, 1}]];                   (* c00 Int (s-yb) *)

(* 5a. moments A_{.,k}(w0), k = 0..Kmax, for each mode *)
momentList[nodes_, wts_, zeta_, freq_, logPart_] := Module[{e, d},
   e = wts Exp[I freq nodes]; d = zeta - w0n;
   Prepend[Table[Total[e d^(-k)], {k, 1, Kmax}], Total[e logPart]]];
alphaW0 = ConstantArray[0. + 0. I, Kmax + 1];
Do[Module[{n = md[[1]], m = md[[2]], kc = md[[3]]/(md[[1]] + I md[[2]]), Ab, AR, AL, pre, comb},
    Ab = momentList[xB, wB, xB + I yBot, n, Log[xB + I yBot - w0n]];
    AR = momentList[yR, wR, 1 + I yR, m, Log[1 + I yR - w0n]];
    AL = momentList[yL, wL, -1 + I yL, m, Log[1 - I yL + w0n] + I Pi];
    pre  = {-Exp[I m yBot], -I Exp[I n], I Exp[-I n]};
    comb = pre[[1]] Ab + pre[[2]] AR + pre[[3]] AL;
    alphaW0 += kc Prepend[-comb[[2 ;;]]/Range[Kmax], comb[[1]]]], {md, modes}];
(* c00 bottom term i c00 Int (x+i yb-w) Log(x+i yb-w) dx *)
Module[{d = xB + I yBot - w0n, L0},
  L0 = Log[d];
  alphaW0[[1]] += I c00 Total[wB d L0] - cS;
  alphaW0[[2]] += I c00 Total[wB (-L0 - 1)];
  Do[alphaW0[[k + 1]] += I c00 Total[wB d^(1 - k)]/(k (k - 1)), {k, 2, Kmax}]];

cornerPts = {-1 + I yBot, 1 + I yBot, 1 + I sc[1], -1 + I sc[-1]};
Print["[Center] w0 = ", N[w0c], ",  dist to bottom side = ", N[Im[w0c] - yBot],
  ",  dist to nearest corner = ", N[Min[Abs[cornerPts - w0c]]]];
Print["[Center] |alpha_k|^(1/k), k = 10, 20, ", Kmax, ": ",
  N[Abs[alphaW0[[{11, 21, Kmax + 1}]]]^(1/{10, 20, Kmax}), 4]];

LRtaylorAll[wv_] := Module[{dw = N[wv - w0c]},
   Accumulate[alphaW0 Prepend[dw^Range[1, Kmax], 1]]];       (* complex, K = 0..Kmax *)

(* 5b. Direct summation of vol_low-rank (Phi_* over the side nodes) *)
rhoB = -Total[Table[md[[3]] Exp[I md[[2]] yBot]/(md[[1]] + I md[[2]]) Exp[I md[[1]] xB], {md, modes}]];
rhoR = -I Total[Table[md[[3]] Exp[I md[[1]]]/(md[[1]] + I md[[2]]) Exp[I md[[2]] yR], {md, modes}]];
rhoL =  I Total[Table[md[[3]] Exp[-I md[[1]]]/(md[[1]] + I md[[2]]) Exp[I md[[2]] yL], {md, modes}]];
LRdirect[wv_] := Module[{ww = N[wv]},
   Total[wB rhoB Log[xB + I yBot - ww]] + Total[wR rhoR Log[1 + I yR - ww]] +
   Total[wL rhoL (Log[1 - I yL + ww] + I Pi)] +
   I c00 Total[wB (xB + I yBot - ww) Log[xB + I yBot - ww]] - cS];  (* complex *)
Print["[Check] |Taylor(K=Kmax) - direct| at w0 + 0.1 = ",
  Abs[LRtaylorAll[w0c + 1/10][[-1]] - LRdirect[w0c + 1/10]]];

(* ---------- 6. Target grid and evaluation ---------- *)
nx = 101; ny = 100;
grid = Table[With[{x0 = i/150}, x0 + I (sc[x0] - (1 + 3/2 sc[x0]) j/150)],
   {i, -50, 50}, {j, 1, ny}];
iList = Range[1, nx, stride]; jList = Range[1, ny, stride];
Print["[Grid] ", Length[iList], " x ", Length[jList], " nodes"];
inRegion[wv_] := With[{a = Re[wv], b = Im[wv]},
   xReg[[1]] <= a <= xReg[[2]] && yRegMin <= b <= Min[0, sc[a]]];
maskR  = Map[Boole[inRegion[#]] &, grid, {2}];
gridR  = Pick[Flatten[grid], Flatten[maskR], 1];
rhoEst = N[Max[Abs[gridR - w0c]]/Min[Abs[cornerPts - w0c]]];
Print["[Grid] region R: ", Length[gridR], " of ", nx ny, " nodes;  max|w - w0| on R = ",
  N[Max[Abs[gridR - w0c]]], ",  rho = ", rhoEst];

DistributeDefinitions[sc, scp, P, dP, pCoef, Jdeg, c00, modes, yBot, volBranch, C00, Cnm, intExp,
  QtExact, chebEval, cQ, cOne, xGT, wGT, QxGT, OneGT, QcalAt1, Lambda, bernstein, rCut, rootTerm,
  topClosed, Ldown, topIntegrand, nearestParam, topRef, pgRef, agRef, maxRec, nintMethod,
  xB, wB, yR, wR, yL, wL, rhoB, rhoR, rhoL, cS, LRdirect, w0c, alphaW0, LRtaylorAll,
  Kmax, grid, inRegion, xReg, yRegMin];

Print["[Eval] closed form + reference on the nodes of R ..."];
res = ParallelTable[
   If[! inRegion[grid[[i, j]]], ConstantArray[0., Kmax + 5],          (* skipped *)
   Module[{wv = grid[[i, j]], tc, tr, br, lrd, lrt},
     tc  = topClosed[wv];                         (* complex *)
     tr  = topRef[wv];                            (* real part, by quadrature *)
     br  = volBranch[N[Re[wv]], N[Im[wv]]];       (* complex, explicit, double *)
     lrd = LRdirect[wv];                          (* complex, direct summation *)
     lrt = LRtaylorAll[wv];                       (* complex, K = 0..Kmax *)
     Join[{N[Abs[Re[tc] - tr]], N[tr + Re[br + lrd]]},
       N[Abs[Re[tc + lrt] - (tr + Re[lrd])]],
       {N[Re[tc + br + lrt[[21]]]], 1.}]]],
   {i, iList}, {j, jList}, Method -> "CoarsestGrained"];

errTop  = res[[All, All, 1]];
Vref    = res[[All, All, 2]];
errK    = res[[All, All, 3 ;; Kmax + 3]];
Vfast20 = res[[All, All, -2]];
inR     = Round[res[[All, All, -1]]];          (* 1 on R, 0 outside *)
onR[e_] := Pick[Flatten[e], Flatten[inR], 1];
maxErrK = Table[Max[onR[errK[[All, All, k]]]], {k, Kmax + 1}];
medErrK = Table[Median[onR[errK[[All, All, k]]]], {k, Kmax + 1}];
Print["[Eval] errors over R only (", Total[inR, 2], " nodes)"];
Print["[Eval] max |Re top closed - Re top quadrature| = ", ScientificForm[Max[onR[errTop]], 3]];
Do[Print["[Eval] K = ", K, " : max |V_K - V_ref| = ", ScientificForm[maxErrK[[K + 1]], 3],
   ",  median = ", ScientificForm[medErrK[[K + 1]], 3]], {K, {0, 2, 4, 6, 8, 10, 12, 16, 20, 25, 30, 35, 40}}];
(* expected on R (yb = -4/3): K = 6: 3.1e-6, K = 12: 5.0e-9, K = 20: 1.3e-12, K = 30: 1.8e-16, K = 40: ~2e-20 *)

(* ---------- 7. Independent 2D check of the full decomposition ---------- *)
vol2D[wv_] := Module[{a = N[Re[wv]], b = N[Im[wv]], inner},
   inner[x_?NumericQ] := NIntegrate[fxy[x, y] Log[Abs[x + I y - wv]],
      Evaluate[If[yBot < b < sc[x], {y, yBot, b, sc[x]}, {y, yBot, sc[x]}]],
      PrecisionGoal -> 14, AccuracyGoal -> 15, MaxRecursion -> 30];      (* machine precision *)
   NIntegrate[inner[x], {x, -1, a, 1},
     PrecisionGoal -> 13, AccuracyGoal -> 14, MaxRecursion -> 30]];
check2D = {};
If[do2DCheck,
  check2D = Table[Module[{wv = grid[[p[[1]], p[[2]]]], vf, v2},
      vf = Re[topClosed[wv] + volBranch[N[Re[wv]], N[Im[wv]]] + LRtaylorAll[wv][[21]]];
      v2 = vol2D[wv];
      Print["[2D] node ", p, "  w = ", N[wv, 6], "  V_fast(K=20) = ", NumberForm[N[vf, 16], 16],
        "  V_2D = ", NumberForm[N[v2, 16], 16], "  |diff| = ", ScientificForm[N[Abs[vf - v2]], 2]];
      N[{p[[1]], p[[2]], vf, v2}]], {p, idx2D}]];
(* expected V (yb = -4/3): (51,1) -0.6885360837, (11,51) -1.1558870554, (91,100) -1.3919234676,
   (1,40) -0.9605829786, (76,31) -1.0845949197; |diff| <= ~1e-12 at K = 20 (Taylor) *)

(* ---------- 8. Export ---------- *)
logE[e_] := N[Log10[Map[Max[#, errFloor] &, e, {Depth[e] - 1}]]];
outDir = FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName],
     Quiet@Check[NotebookDirectory[], Directory[]]], "vol_closed_out"}];
Quiet@CreateDirectory[outDir, CreateIntermediateDirectories -> True];
gsub = grid[[iList, jList]];
Export[FileNameJoin[{outDir, "zdataxy_vol.mat"}],
  <|"X" -> N[Re[gsub]], "Y" -> N[Im[gsub]], "Nx" -> Length[iList], "Ny" -> Length[jList],
    "w0" -> N[{Re[w0c], Im[w0c]}], "rhoEst" -> rhoEst,
    "inRegion" -> inR, "xReg" -> N[xReg], "yRegMin" -> N[yRegMin]|>, "MAT"];
Export[FileNameJoin[{outDir, "zerror_vol.mat"}],
  <|"logErrAll" -> logE[errK], "Klist" -> Range[0, Kmax],
    "maxErrK" -> maxErrK, "medErrK" -> medErrK,
    "logErrK6" -> logE[errK[[All, All, 7]]], "logErrK12" -> logE[errK[[All, All, 13]]],
    "logErrK20" -> logE[errK[[All, All, 21]]],
    "logErrTop" -> logE[errTop], "maxErrTop" -> Max[onR[errTop]], "inRegion" -> inR,
    "Vref" -> Vref, "Vfast20" -> Vfast20,
    "check2D" -> If[check2D === {}, {{0, 0, 0, 0}}, check2D],
    "nCheb" -> nCheb, "precision" -> "MachinePrecision"|>, "MAT"];
Print["[Done] output in ", outDir];



