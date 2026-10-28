(* ::Package:: *)
(* ======================================================================
   poly_closed_form_LaurentNewton.wl
   Closed-form (exact) evaluation of the standardized layer potential
       phi(w) = Int_{-1}^{1} f(x) / (P(x) - w) dx,   P(x) = x + i s(x),
   and its derivative phi'(w) = Int f / (P - w)^2 dx, using Sec. 2.3
   (label sec:polydiv) of the paper:

     f(x) = q(x,w) (P(x) - w) + r(x,w),   deg_x r < J = deg s,
     r(w_j, w) = f(w_j)  at the J roots w_j of P(z) = w, hence
     phi(w)  = Q(w) + Sum_j f(w_j)/P'(w_j) Lambda(w_j),
     Lambda(a) = Int_{-1}^{1} dx/(x - a) = Log((a-1)/(a+1)) = -2 ArcTanh(1/a),
     Q(w)    = Int_{-1}^{1} q(x,w) dx    (polynomial in w, degree floor(N/J)-1),
     phi'(w) = Q'(w) + Sum_j [ (f' Lambda + 2 f/(z^2-1)) P' - f Lambda P'' ] / P'^3
               evaluated at z = w_j.

   This replaces poly_final_version.nb, which used residues in the upper
   half plane plus far-field Taylor expansions of the tail integrals
   Int_{|x|>1} (and an arc term for l = 3) truncated at order n. The closed
   form has no truncation parameter: its error is set only by the
   arithmetic precision, so the "n" axis of the old Fig. 5 / Fig. 7
   disappears.

   Test problem (as in the paper's Fig. 5 / Fig. 7):
     s(x) = 5 x^4 + 2 x^2  (J = 4),  f(x) = T_5(x)  (N = 5),
     targets: 101 x 100 curve-following grid below Gamma (exact rationals).

   Outputs (poly_closed_out/), MAT files with named variables:
     zdataxy_closed.mat : X, Y (nx x ny), dist (distance to Gamma), Nx, Ny
     zerror_closed.mat  : logErrMach, logErrHP, logErrDMach, logErrDHP (nx x ny),
                          maxErrMach, maxErrHP, maxErrDMach, maxErrDHP,
                          evalPrecHP, targetDigits
     ("Mach" = closed form in machine precision, which is what a C/Fortran/
      MATLAB implementation achieves; "HP" = closed form at evalPrecHP digits,
      which verifies exactness against the reference.)

   Evaluation (machine precision and evalPrecHP digits):
     - all roots polished by 2 Newton steps on P(z) - w;
     - the closest root pair is summed by the Laurent/Taylor residue series
       (Sec. 3b) when (max|delta|/rho)^2 < ratioMax = 1/4, for phi and phi';
       clustered roots are not Newton-polished.

   Independent checks (Python/mpmath, same grid, 30-digit closed form as
   reference; max / median absolute error in double precision):
                            phi                 phi'
     plain root sum       1.07e-13 / 7.2e-15   4.03e-12 / 3.0e-14
     plain + Newton       3.36e-14 / 1.7e-15   1.52e-12 / 5.6e-15
     pair-Taylor          3.73e-14 / 5.3e-15   3.49e-13 / 1.5e-14
     pair-Taylor + Newton 7.83e-15 / 1.5e-15   3.28e-14 / 3.5e-15   (quadratic factor)
     Laurent + Newton     7.16e-15 / 1.6e-15   3.28e-14 / 3.5e-15   <- this script
   On this grid |phi| <= 2.61 and |phi'| <= 23.4, so the last row is at
   round-off (relative 3e-15 and 1.4e-15 of the maxima). The cluster series is
   active at 3982 of 10100 points. Note: direct quadrature in double
   precision (scipy) is itself off by up to 2.1e-12 for phi' near Gamma.
   ====================================================================== *)

ClearAll["Global`*"];
$HistoryLength = 0;

(* ---------- 0. Precision budget (target: 16 correct digits) ----------
   Reference integral: the near-singular integrand peaks at ~1/gap ~ 150-500
   for values O(1), so ~3 digits are lost. Goal 16 + 2 guard digits;
   WorkingPrecision = PrecisionGoal + 3 + 8 (the 8 guard digits let the
   error estimator certify the goal).
   Closed form, HP mode: 30 digits leaves >10 digits of margin over the
   ~3 digits lost near the critical value. *)
targetDigits = 16;
pgRef  = targetDigits + 2;                  (* 18 *)
wpRef  = pgRef + 3 + 8;                     (* 29 *)
evalPrecHP = 30;
maxRec = 20;
nintMethod = {"GlobalAdaptive", Method -> {"GaussKronrodRule", "Points" -> 15},
   "MaxErrorIncreases" -> 2000, "SymbolicProcessing" -> 0};
errFloor = 10^-(targetDigits + 2);          (* floor for log10 output *)
doDerivative = True;                         (* also evaluate phi' (Fig. 7) *)

(* ---------- 1. Geometry and density (exact) ---------- *)
sc[x_]  := 5 x^4 + 2 x^2;
scp[x_] := 20 x^3 + 4 x;
P[x_]   := x + I sc[x];
dP[x_]  := 1 + I scp[x];
d2P[x_] := I (60 x^2 + 4);
Jdeg = Exponent[sc[xx], xx];                 (* 4 *)
cJ   = Coefficient[P[xx], xx, Jdeg];         (* leading coeff of P = 5 I (NOT of s) *)

nCheb = 5;
f[x_]  := ChebyshevT[nCheb, x];
fp[x_] := Evaluate[D[ChebyshevT[nCheb, t], t] /. t -> x];

Print["[Init] J = ", Jdeg, ", cJ = ", cJ, ", f = T_", nCheb, ", targetDigits = ", targetDigits];

(* ---------- 2. Division f = q (P - w) + r, and Q(w) ---------- *)
(* I is temporarily replaced by a symbol so that PolynomialQuotientRemainder
   treats the coefficients as polynomials (as in the original notebook). *)
divideByPminusW[fx_, x_, w_] := Module[{ii, Pii, q, r},
   Pii = Expand[(P[x] /. I -> ii) - w];
   {q, r} = PolynomialQuotientRemainder[Expand[fx /. I -> ii], Pii, x];
   Expand[{q, r} /. ii -> I]];

{qxw, rxw} = divideByPminusW[f[x], x, w];
Print["[Div] identity f - (q(P-w)+r) = ", Simplify[Expand[f[x] - (qxw (P[x] - w) + rxw)]],
  ";  deg_x r = ", Exponent[rxw, x], " (< J = ", Jdeg, ")"];
Qw  = Expand[Integrate[qxw, {x, -1, 1}]];    (* exact polynomial in w *)
dQw = D[Qw, w];
Print["[Div] Q(w) = ", Qw];

(* ---------- 3. Closed form ---------- *)
Lambda[a_] := -2 ArcTanh[1/a];   (* = Log((a-1)/(a+1)); cut a in [-1,1]; no
                                     cancellation for large |a| *)
setP[v_, prec_] := If[prec === MachinePrecision, N[v], SetPrecision[v, prec]];

rootsP[wv_, prec_] := Module[{ww = setP[wv, prec], z, rts},
   rts = NSolveValues[P[z] == ww, z, WorkingPrecision -> prec];
   If[Length[rts] =!= Jdeg, Print["[rootsP] wrong root count at w = ", N[wv]]; Abort[]];
   rts];

(* ---------- 3a. Newton polishing of the roots ----------
   Machine-precision roots (companion/eigenvalue solvers) carry errors of
   ~1e-15; with |d(term)/dz| ~ 30 this produced errors of ~3e-14 in phi.
   Two Newton steps on p(z) = P(z) - w reduce the root errors to ~1e-16. *)
newtonSteps = 2;
polishRoots[rts_, ww_] := Nest[# - (P[#] - ww)/dP[#] &, rts, newtonSteps];

(* ---------- 3b. Laurent/Taylor series for clustered roots ----------
   Let C be the cluster (here: the closest pair, p = 2), w0 = Mean of the
   clustered roots (z-plane, ~ z_c), delta_j = w_j - w0.
     P(z) - w = cJ Prod_k (z - w_k),   G_q(z) = F(z) / Prod_{k not in C}(z - w_k)^q.
   Expanding Prod_{j in C}(z - w_j)^(-q) in powers of (z - w0)^(-1)
   (geometric series, valid for |z - w0| > max|delta_j|) and G_q(z) =
   Sum_n gamma^(q)_n (z - w0)^n (Taylor, valid for |z - w0| < rho), the
   residue on the annulus max|delta| < |z - w0| < rho gives
     phi  : Sum_{j in C} F(w_j)/P'(w_j)
              = cJ^-1 Sum_n gamma^(1)_n h_(n-p+1)(delta_1..delta_p),
     phi' : Sum_{j in C} Res F/(P-w)^2
              = cJ^-2 Sum_n gamma^(2)_n h_(n-2p+1)(delta_1,delta_1,..,delta_p,delta_p),
   h_k = complete homogeneous symmetric polynomial of degree k. For a pair
   with delta = +-h only odd n (phi) contribute: cJ^-1 Sum_k gamma_(2k+1) h^(2k).
   This equals the quadratic-factor pair-Taylor series.
   rho = dist(w0, [-1,1] U {other roots}); terms shrink by (max|delta|/rho)^2.
   The clustered roots are NOT Newton-polished: P' ~ 0 there and individual
   Newton steps degrade them (Python: errors up to 4.5e-11). The other
   (simple) roots are polished.
   Python check (f = T5, |w - w_c| = 1e-1..1e-12, 50-digit reference):
     w_c = -0.6981 i          : phi 4.5e-16..2.8e-15
     w_c = 0.0291 + 0.1490 i  : phi 3.4e-15..1.8e-14,  phi' <= 1.6e-13
                                (|phi'| ~ 3.9; quadratic factor gives the same) *)
ratioMax = 1/4;
tolTaylor[prec_] := If[prec === MachinePrecision, 10^-17, 10^-(prec + 1)];
distToSegment[z_] := If[-1 <= Re[z] <= 1, Abs[Im[z]], Min[Abs[z - 1], Abs[z + 1]]];
closestPair[rts_] := First@MinimalBy[Subsets[Range[Length[rts]], {2}],
    Abs[rts[[#[[1]]]] - rts[[#[[2]]]]] &];

shiftCoeffs[poly_, zv_, m_, nmax_] := Module[{tt},
   PadRight[CoefficientList[Expand[poly /. zv -> m + tt], tt], nmax + 1]];
cauchy[u_, v_, nmax_] := Table[Sum[u[[k + 1]] v[[n - k + 1]], {k, 0, n}], {n, 0, nmax}];
seriesRecip[a_, nmax_] := Module[{inv = ConstantArray[0, nmax + 1]},
   inv[[1]] = 1/a[[1]];
   Do[inv[[n + 1]] = -Sum[a[[k + 1]] inv[[n - k + 1]], {k, 1, n}]/a[[1]], {n, 1, nmax}];
   inv];

(* Complete homogeneous symmetric polynomials h_0..h_kmax of the offsets:
   coefficients of Prod_j 1/(1 - delta_j t). *)
completeH[deltas_, kmax_] := Module[{h = PadRight[{1}, kmax + 1]},
   Do[Do[h[[k + 1]] += d h[[k]], {k, 1, kmax}], {d, deltas}];
   h];

(* Taylor coefficients 0..nmax of G_q = f Lambda / Prod_{k in oth}(z - w_k)^q about w0. *)
taylorGq[w0_, oth_, zv_, nmax_, q_] := Module[{fs, lam, Rpoly, Rs},
   fs  = shiftCoeffs[f[zv], zv, w0, nmax];
   lam = Table[If[k == 0, Lambda[w0], (-1)^(k - 1) ((w0 - 1)^-k - (w0 + 1)^-k)/k], {k, 0, nmax}];
   Rpoly = Expand[Times @@ ((zv - #) & /@ oth)];          (* = 1 if oth = {} *)
   Rs  = shiftCoeffs[Rpoly, zv, w0, nmax];
   If[q == 2, Rs = cauchy[Rs, Rs, nmax]];
   cauchy[cauchy[fs, lam, nmax], seriesRecip[Rs, nmax], nmax]];

Fphi[z_] := f[z] Lambda[z];
dTerm[z_] := ((fp[z] Lambda[z] + 2 f[z]/(z^2 - 1)) dP[z] - f[z] Lambda[z] d2P[z])/dP[z]^3;

(* Returns {phi, phi'} and whether the cluster branch was used. *)
evalClosed[wv_, prec_] := Module[
   {ww = setP[wv, prec], rts, pr, cl, oth, w0, dl, p, rho, ratio, nt, nmax, zv,
    g1, g2, h1, h2, phiC, dphiC},
   rts = rootsP[wv, prec];
   pr = closestPair[rts]; cl = rts[[pr]]; oth = Delete[rts, List /@ pr];
   w0 = Mean[cl]; dl = cl - w0; p = Length[cl];
   rho = Min[If[oth === {}, Infinity, Min[Abs[oth - w0]]], distToSegment[w0]];
   ratio = (Max[Abs[dl]]/rho)^2;
   If[ratio >= ratioMax,
     (* well-separated roots: Newton-polished direct sums *)
     rts = polishRoots[rts, ww];
     Return[<|"phi" -> (Qw /. w -> ww) + Total[(Fphi /@ rts)/(dP /@ rts)],
              "dphi" -> (dQw /. w -> ww) + Total[dTerm /@ rts], "paired" -> False|>]];
   (* cluster: polish only the simple roots; series for the cluster *)
   oth = polishRoots[oth, ww];
   rho = Min[If[oth === {}, Infinity, Min[Abs[oth - w0]]], distToSegment[w0]];
   ratio = (Max[Abs[dl]]/rho)^2;
   nt = If[ratio == 0, 1, Max[1, Ceiling[Log[tolTaylor[prec]]/Log[ratio]] + 2]];
   nmax = 2 nt + 2 p + 1;
   g1 = taylorGq[w0, oth, zv, nmax, 1];
   g2 = taylorGq[w0, oth, zv, nmax, 2];
   h1 = completeH[dl, nmax];
   h2 = completeH[Join[dl, dl], nmax];
   phiC  = Sum[g1[[n + 1]] h1[[n - p + 2]], {n, p - 1, nmax}]/cJ;          (* h_(n-p+1) *)
   dphiC = Sum[g2[[n + 1]] h2[[n - 2 p + 2]], {n, 2 p - 1, nmax}]/cJ^2;    (* h_(n-2p+1) *)
   <|"phi"  -> (Qw /. w -> ww) + Total[(Fphi /@ oth)/(dP /@ oth)] + phiC,
     "dphi" -> (dQw /. w -> ww) + Total[dTerm /@ oth] + dphiC, "paired" -> True|>];

phiClosed[wv_, prec_]  := evalClosed[wv, prec]["phi"];
dphiClosed[wv_, prec_] := evalClosed[wv, prec]["dphi"];

(* ---------- 4. Reference (direct quadrature, once per point) ---------- *)
(* Split at s*, the parameter of the boundary point nearest to w. *)
nearestParam[a_, b_] := Module[{s0},
   s0 = s /. FindRoot[(s - a) + (sc[s] - b) scp[s] == 0, {s, a},
       WorkingPrecision -> wpRef, AccuracyGoal -> 10, PrecisionGoal -> 10];
   Clip[s0, {-1 + 10^-6, 1 - 10^-6}]];

refInt[wv_, ss_, pw_] := NIntegrate[f[t]/(P[t] - wv)^pw, {t, -1, ss, 1},
   WorkingPrecision -> wpRef, PrecisionGoal -> pgRef, AccuracyGoal -> pgRef,
   MaxRecursion -> maxRec, Method -> nintMethod];

(* ---------- 5. Target grid (exact rationals, same as Figs. 3-7) ---------- *)
nx = 101; ny = 100;
grid = Table[With[{x0 = i/150}, x0 + I (sc[x0] - (1 + 3/2 sc[x0]) j/150)],
   {i, -50, 50}, {j, 1, ny}];
Print["[Grid] ", Dimensions[grid], " points"];

(* ---------- 6. Evaluate ---------- *)
DistributeDefinitions[sc, scp, P, dP, d2P, Jdeg, cJ, nCheb, f, fp, Qw, dQw, w,
  Lambda, setP, rootsP, newtonSteps, polishRoots, ratioMax, tolTaylor, distToSegment,
  closestPair, shiftCoeffs, cauchy, seriesRecip, completeH, taylorGq,
  Fphi, dTerm, evalClosed, phiClosed, dphiClosed, nearestParam, refInt, wpRef,
  pgRef, maxRec, nintMethod, evalPrecHP, doDerivative];

Print["[Eval] reference + closed form on ", nx ny, " points ..."];
res = ParallelTable[
   Module[{wv = grid[[i, j]], ss, ref, refD, cM, cH},
     ss  = nearestParam[Re[wv], Im[wv]];
     ref = refInt[wv, ss, 1];
     cM  = evalClosed[wv, MachinePrecision];
     cH  = evalClosed[wv, evalPrecHP];
     refD = If[doDerivative, refInt[wv, ss, 2], 0];
     N[{Abs[P[ss] - wv],
        Abs[cM["phi"] - ref], Abs[cH["phi"] - ref],
        If[doDerivative, Abs[cM["dphi"] - refD], 0],
        If[doDerivative, Abs[cH["dphi"] - refD], 0],
        Boole[cM["paired"]]}]],
   {i, nx}, {j, ny}, Method -> "CoarsestGrained"];

distG = res[[All, All, 1]];
{eM, eH, eDM, eDH} = Table[res[[All, All, k]], {k, 2, 5}];
pairedG = res[[All, All, 6]];
Print["[Eval] cluster (Laurent/Taylor) series used at ", Total[pairedG, 2], " of ", nx ny, " points"];
logE[e_] := N[Log10[Map[Max[#, errFloor] &, e, {2}]]];
Print["[Eval] max |err phi|  : machine ", ScientificForm[Max[eM], 3], ",  HP ", ScientificForm[Max[eH], 3]];
If[doDerivative, Print["[Eval] max |err phi'| : machine ", ScientificForm[Max[eDM], 3], ",  HP ", ScientificForm[Max[eDH], 3]]];
(* Expected (Python check, double precision): max |err phi| ~ 7.2e-15,
   max |err phi'| ~ 3.3e-14 (Python: 7.16e-15 / 3.28e-14); cluster series at
   ~3982 points. HP should be
   limited by the reference goal (~1e-18). Not yet run in Mathematica. *)

(* ---------- 7. Export ---------- *)
outDir = FileNameJoin[{If[$InputFileName =!= "", DirectoryName[$InputFileName],
     Quiet@Check[NotebookDirectory[], Directory[]]], "poly_closed_out"}];
Quiet@CreateDirectory[outDir, CreateIntermediateDirectories -> True];
Export[FileNameJoin[{outDir, "zdataxy_closed.mat"}],
  <|"X" -> N[Re[grid]], "Y" -> N[Im[grid]], "dist" -> distG, "Nx" -> nx, "Ny" -> ny|>, "MAT"];
Export[FileNameJoin[{outDir, "zerror_closed.mat"}],
  <|"logErrMach" -> logE[eM], "logErrHP" -> logE[eH],
    "logErrDMach" -> logE[eDM], "logErrDHP" -> logE[eDH],
    "maxErrMach" -> Max[eM], "maxErrHP" -> Max[eH],
    "maxErrDMach" -> Max[eDM], "maxErrDHP" -> Max[eDH],
    "paired" -> pairedG, "ratioMax" -> N[ratioMax], "newtonSteps" -> newtonSteps,
    "evalPrecHP" -> evalPrecHP, "targetDigits" -> targetDigits|>, "MAT"];
Print["[Done] output in ", outDir];
