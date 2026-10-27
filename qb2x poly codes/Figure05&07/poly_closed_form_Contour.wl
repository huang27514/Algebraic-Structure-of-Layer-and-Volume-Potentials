(* ::Package:: *)
(* ======================================================================
   poly_closed_form.wl
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

   Independent checks (Python/mpmath, same grid):
     - closed form vs reference at 40 digits: 5.5e-40 at the worst point;
     - machine precision: max |err phi| = 1.07e-13, median 7.2e-15;
       max |err phi'| = 4.1e-12, median 3.1e-14;
     - both maxima occur next to the critical value w_c = -0.6981 i
       (P'(z_c) = 0), 0.031 below the deepest grid point. There two roots
       approach each other (min separation 0.137 on this grid) and
       round-off in the roots is amplified by 1/P'.
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

(* Clustered roots (|w_i - w_j| < clusterTol) are summed by the contour
   divided difference of Sec. 2.3, eq. (cluster):
     Sum_{j in C} F(w_j)/P'(w_j) = (1/cJ) (1/2 pi i) Oint G(z)/Prod_{j in C}(z - w_j) dz,
     G(z) = F(z)/Prod_{k not in C}(z - w_k),
   on a circle that encloses the cluster and excludes [-1,1] and the other
   roots. The trapezoidal rule converges geometrically. On the Fig. 5 grid
   the minimum separation is 0.137 > clusterTol, so this branch is not
   triggered; it protects targets closer to a critical value. *)
clusterTol = 10^-2;  nTrap = 64;
distToSegment[z_] := If[-1 <= Re[z] <= 1, Abs[Im[z]], Min[Abs[z - 1], Abs[z + 1]]];

clusterSum[Fz_, idx_, rts_, prec_] := Module[{cl, oth, zc, rad, th, zz, G},
   cl = rts[[idx]]; oth = Delete[rts, List /@ idx];
   If[Length[cl] == 1, Return[Fz[cl[[1]]]/dP[cl[[1]]]]];
   zc = Mean[cl]; rad = 2 Max[Abs[cl - zc]] + clusterTol;
   If[distToSegment[zc] <= 2 rad || (oth =!= {} && Min[Abs[oth - zc]] <= 2 rad),
     (* circle not admissible: fall back to the direct sum *)
     Return[Total[(Fz /@ cl)/(dP /@ cl)]]];
   th = N[2 Pi Range[0, nTrap - 1]/nTrap, prec];
   zz = zc + rad Exp[I th];
   G  = Fz[zz]/(Times @@ ((zz - #) & /@ oth))/(Times @@ ((zz - #) & /@ cl));
   Mean[G rad Exp[I th]]/cJ];

Fphi[z_] := f[z] Lambda[z];

phiClosed[wv_, prec_] := Module[{rts, groups},
   rts = rootsP[wv, prec];
   groups = Gather[Range[Jdeg], Abs[rts[[#1]] - rts[[#2]]] < clusterTol &];
   (Qw /. w -> setP[wv, prec]) + Total[clusterSum[Fphi, #, rts, prec] & /@ groups]];

(* Derivative: direct root sum of eq. (dclosed). Clustered roots would need
   the divided-difference form of the derivative; not needed on this grid. *)
dphiClosed[wv_, prec_] := Module[{rts = rootsP[wv, prec]},
   (dQw /. w -> setP[wv, prec]) + Total[Function[z,
       ((fp[z] Lambda[z] + 2 f[z]/(z^2 - 1)) dP[z] - f[z] Lambda[z] d2P[z])/dP[z]^3] /@ rts]];

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
  Lambda, setP, rootsP, clusterTol, nTrap, distToSegment, clusterSum, Fphi,
  phiClosed, dphiClosed, nearestParam, refInt, wpRef, pgRef, maxRec, nintMethod,
  evalPrecHP, doDerivative];

Print["[Eval] reference + closed form on ", nx ny, " points ..."];
res = ParallelTable[
   Module[{wv = grid[[i, j]], ss, ref, refD, out},
     ss  = nearestParam[Re[wv], Im[wv]];
     ref = refInt[wv, ss, 1];
     out = {Abs[P[ss] - wv],
            Abs[phiClosed[wv, MachinePrecision] - ref],
            Abs[phiClosed[wv, evalPrecHP] - ref]};
     If[doDerivative,
       refD = refInt[wv, ss, 2];
       out = Join[out, {Abs[dphiClosed[wv, MachinePrecision] - refD],
                        Abs[dphiClosed[wv, evalPrecHP] - refD]}],
       out = Join[out, {0, 0}]];
     N[out]],
   {i, nx}, {j, ny}, Method -> "CoarsestGrained"];

distG = res[[All, All, 1]];
{eM, eH, eDM, eDH} = Table[res[[All, All, k]], {k, 2, 5}];
logE[e_] := N[Log10[Map[Max[#, errFloor] &, e, {2}]]];
Print["[Eval] max |err phi|  : machine ", ScientificForm[Max[eM], 3], ",  HP ", ScientificForm[Max[eH], 3]];
If[doDerivative, Print["[Eval] max |err phi'| : machine ", ScientificForm[Max[eDM], 3], ",  HP ", ScientificForm[Max[eDH], 3]]];
(* Expected (Python check): machine 1.07e-13 (phi), 4.1e-12 (phi').
   HP should be limited by the reference goal (~1e-18); not yet run. *)

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
    "evalPrecHP" -> evalPrecHP, "targetDigits" -> targetDigits|>, "MAT"];
Print["[Done] output in ", outDir];
