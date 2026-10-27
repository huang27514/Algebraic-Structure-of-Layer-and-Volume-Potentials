(* ::Package:: *)

(* ======================================================================
   poly_closed_form_PairTaylorNewton.wl
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
     - the closest root pair is summed by the pair-Taylor expansion when
       (|h|/rho)^2 < ratioMax = 1/4, for both phi and phi'.

   Independent checks (Python/mpmath, same grid, 30-digit closed form as
   reference; max / median absolute error in double precision):
                            phi                 phi'
     plain root sum       1.07e-13 / 7.2e-15   4.03e-12 / 3.0e-14
     plain + Newton       3.36e-14 / 1.7e-15   1.52e-12 / 5.6e-15
     pair-Taylor          3.73e-14 / 5.3e-15   3.49e-13 / 1.5e-14
     pair-Taylor + Newton 7.83e-15 / 1.5e-15   3.28e-14 / 3.5e-15   <- this script
   On this grid |phi| <= 2.61 and |phi'| <= 23.4, so the last row is at
   round-off (relative 3e-15 and 1.4e-15 of the maxima). Pair-Taylor is
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

(* ---------- 3b. Pair-Taylor expansion for the closest pair ----------
   For a pair a, b of roots of P(z) = w:
     P(z) - w = D(z) R(z),   D(z) = z^2 - s z + t,
       (s, t) refined by Bairstow-Newton from s = a + b, t = a b;
       the quadratic factor is well conditioned even when a = b.
     m = s/2,  h^2 = s^2/4 - t.
   phi  : Res_{a} + Res_{b} of F/(P-w)   = (F/R)[a,b]
                                         = Sum_{k>=0} g_(2k+1) h^(2k),
          g_n = Taylor coefficients of F/R about m;
   phi' : Res_{a} + Res_{b} of F/(P-w)^2 = (F/R^2)[a,a,b,b]
                                         = Sum_{k>=0} (k+1) G_(2k+3) h^(2k),
          G_n = Taylor coefficients of F/R^2 about m.
   (Sum of residues of g/Prod(z-x_i) is the divided difference g[x_1..x_n];
    the node offsets (h,h,-h,-h) give the generating function
    1/(1 - h^2 t^2)^2 = Sum (k+1) h^(2k) t^(2k).)
   Coefficients are exact: f(m+t) by polynomial shift,
     Lambda^(k)(m)/k! = (-1)^(k-1) ((m-1)^-k - (m+1)^-k)/k,
     1/R(m+t)^p by series reciprocal.
   Only even powers of h occur: no small denominator. Terms shrink by
   (|h|/rho)^2, rho = dist(m, [-1,1] U {other roots}). The pair branch is
   used when (|h|/rho)^2 < ratioMax; otherwise all roots are summed directly. *)
ratioMax = 1/4;
tolTaylor[prec_] := If[prec === MachinePrecision, 10^-17, 10^-(prec + 1)];
distToSegment[z_] := If[-1 <= Re[z] <= 1, Abs[Im[z]], Min[Abs[z - 1], Abs[z + 1]]];
closestPair[rts_] := First@MinimalBy[Subsets[Range[Length[rts]], {2}],
    Abs[rts[[#[[1]]]] - rts[[#[[2]]]]] &];

(* Bairstow-Newton refinement of the quadratic factor of p (variable zv):
     p = q Dq + r,  q = q2 Dq + r2,  dr/dt = -r2,  dr/ds = (zv r2) mod Dq. *)
coeffs2[e_, zv_] := PadRight[CoefficientList[e, zv], 2];
quadFactor[p_, zv_, a_, b_, prec_, maxIt_: 10] :=
  Module[{s = a + b, t = a b, eps, Dq, q, r, r2, rs, ds, dt},
   eps = If[prec === MachinePrecision, 4 $MachineEpsilon, 10^-prec];
   Do[
     Dq = zv^2 - s zv + t;
     {q, r} = PolynomialQuotientRemainder[p, Dq, zv];
     r2 = PolynomialRemainder[q, Dq, zv];
     rs = PolynomialRemainder[Expand[zv r2], Dq, zv];
     {ds, dt} = LinearSolve[Transpose[{coeffs2[rs, zv], -coeffs2[r2, zv]}], -coeffs2[r, zv]];
     s += ds; t += dt;
     If[Abs[ds] + Abs[dt] < eps (1 + Abs[s] + Abs[t]), Break[]],
     {maxIt}];
   {s, t, PolynomialQuotient[p, zv^2 - s zv + t, zv]}];

(* Taylor coefficients 0..nmax of F/R^pw about m, F = f Lambda. *)
shiftCoeffs[poly_, zv_, m_, nmax_] := Module[{tt},
   PadRight[CoefficientList[Expand[poly /. zv -> m + tt], tt], nmax + 1]];
cauchy[u_, v_, nmax_] := Table[Sum[u[[k + 1]] v[[n - k + 1]], {k, 0, n}], {n, 0, nmax}];
seriesRecip[a_, nmax_] := Module[{inv = ConstantArray[0, nmax + 1]},
   inv[[1]] = 1/a[[1]];
   Do[inv[[n + 1]] = -Sum[a[[k + 1]] inv[[n - k + 1]], {k, 1, n}]/a[[1]], {n, 1, nmax}];
   inv];
taylorFR[m_, Rpoly_, zv_, nmax_, pw_] := Module[{fs, lam, Rs},
   fs  = shiftCoeffs[f[zv], zv, m, nmax];
   lam = Table[If[k == 0, Lambda[m], (-1)^(k - 1) ((m - 1)^-k - (m + 1)^-k)/k], {k, 0, nmax}];
   Rs  = shiftCoeffs[Rpoly, zv, m, nmax];
   If[pw == 2, Rs = cauchy[Rs, Rs, nmax]];
   cauchy[cauchy[fs, lam, nmax], seriesRecip[Rs, nmax], nmax]];

Fphi[z_] := f[z] Lambda[z];
dTerm[z_] := ((fp[z] Lambda[z] + 2 f[z]/(z^2 - 1)) dP[z] - f[z] Lambda[z] d2P[z])/dP[z]^3;

(* Returns {phi, phi'} and whether the pair branch was used. *)
evalClosed[wv_, prec_] := Module[
   {ww = setP[wv, prec], rts, pr, a, b, m0, h20, oth0, rho, ratio, zv, p, s, t, R,
    m, h2, oth, nt, g1, g2, phiPair, dphiPair},
   rts = rootsP[wv, prec];
   pr = closestPair[rts]; {a, b} = rts[[pr]]; oth0 = Delete[rts, List /@ pr];
   m0 = (a + b)/2; h20 = ((a - b)/2)^2;
   rho = Min[If[oth0 === {}, Infinity, Min[Abs[oth0 - m0]]], distToSegment[m0]];
   ratio = Abs[h20]/rho^2;
   If[ratio >= ratioMax,
     (* well-separated roots: Newton-polished direct sums *)
     rts = polishRoots[rts, ww];
     Return[<|"phi" -> (Qw /. w -> ww) + Total[(Fphi /@ rts)/(dP /@ rts)],
              "dphi" -> (dQw /. w -> ww) + Total[dTerm /@ rts], "paired" -> False|>]];
   (* clustered pair: quadratic factor + Taylor sums; others polished *)
   p = Expand[P[zv] - ww];
   {s, t, R} = quadFactor[p, zv, a, b, prec];
   m = s/2; h2 = s^2/4 - t;
   oth = polishRoots[NSolveValues[R == 0, zv, WorkingPrecision -> prec], ww];
   rho = Min[If[oth === {}, Infinity, Min[Abs[oth - m]]], distToSegment[m]];
   ratio = Abs[h2]/rho^2;
   nt = If[ratio == 0, 1, Max[1, Ceiling[Log[tolTaylor[prec]]/Log[ratio]] + 2]];
   g1 = taylorFR[m, R, zv, 2 nt + 3, 1];
   g2 = taylorFR[m, R, zv, 2 nt + 3, 2];
   phiPair  = Sum[g1[[2 k + 2]] h2^k, {k, 0, nt}];             (* g_(2k+1) h^(2k) *)
   dphiPair = Sum[(k + 1) g2[[2 k + 4]] h2^k, {k, 0, nt}];     (* (k+1) G_(2k+3) h^(2k) *)
   <|"phi"  -> (Qw /. w -> ww) + Total[(Fphi /@ oth)/(dP /@ oth)] + phiPair,
     "dphi" -> (dQw /. w -> ww) + Total[dTerm /@ oth] + dphiPair, "paired" -> True|>];

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
  closestPair, coeffs2, quadFactor, shiftCoeffs, cauchy, seriesRecip, taylorFR,
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
Print["[Eval] pair-Taylor used at ", Total[pairedG, 2], " of ", nx ny, " points"];
logE[e_] := N[Log10[Map[Max[#, errFloor] &, e, {2}]]];
Print["[Eval] max |err phi|  : machine ", ScientificForm[Max[eM], 3], ",  HP ", ScientificForm[Max[eH], 3]];
If[doDerivative, Print["[Eval] max |err phi'| : machine ", ScientificForm[Max[eDM], 3], ",  HP ", ScientificForm[Max[eDH], 3]]];
(* Expected (Python check, double precision): max |err phi| ~ 7.8e-15,
   max |err phi'| ~ 3.3e-14; pair-Taylor at ~3982 points. HP should be
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



