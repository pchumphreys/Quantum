(* ::Package:: *)

BeginPackage[ "CVtoolbox`"]

\[Kappa]1::usage = "\[Kappa]1[] is a constant to define the chosen quantum optics convention";
\[Kappa]2::usage = "\[Kappa]1[] is a constant to define the chosen quantum optics convention";
\[Kappa]3::usage = "\[Kappa]1[] is a constant to define the chosen quantum optics convention";

P::usage = "P[n] is as defined by Monras";
J::usage = "J[n] is as defined by Monras";
\[Omega]::usage = "\[Omega][n] is as defined by Monras";
\[Lambda]To\[CapitalLambda]::usage = "\[Lambda]To\[CapitalLambda][\[Lambda]] is as defined in Monras, gets between different coordinate systems";

MatDirSum::usage = "MatDirSum[{sqMatrix1, sqMatrix2,...}] calculates the direct sum of n square matrices";
singleModeToMultiMode::usage = "singleModeToMultiMode[\[CapitalSigma]singleModeMat,i,n] takes a single mode symplectic matrix, and embeds it in an identity matrix of size 2n at posn 2i-1:2i";
twoModeToMultiMode::usage = "twoModeToMultiMode[\[CapitalSigma]twoModeMat,i,j,n] takes a two mode symplectic matrix, and embeds it in an identity matrix of size 2n at posns 2i-1:2i, 2j-1:2j";

UnitaryToSimplect::usage = "UnitaryToSimplect[UMat] converts a unitary matrix to symplectic matrix";
ActOnCovMatWithSympMat::usage = "ActOnCovMatWithSympMat[SympMat,CovMat] gives the transformed covariance matrix resulting from acting SympMat on CovMat.";
ActOnRMatWithSympMat::usage = "ActOnRMatWithSympMat[SympMat,RMat] gives the transformed covariance matrix resulting from acting SympMat on CovMat.";

rotMat::usage = "rotMat[\[Theta]] gives a simple single mode symplectic rotation matrix";
Displacement::usage = "Displacement[R,\[Lambda]] displaces the first order moments R by \[Lambda]";
SymBSUnitary::usage = "SymBSUnitary[\[Eta]] gives the symplectic form for the BS unitary defined in terms of \[Eta]";
SymBSUnitaryPhase::usage = "SymBSUnitaryPhase[\[Phi]] gives the symplectic form for the BS unitary defined in terms of \[Phi]";
RotateSingleModeQuadrature::usage = "Rotate the quadratures of one mode in n modes";
WignerFunction::usage = "Calculate the Wigner function (for a Gaussian state) given quadratures, covariance matrix and dimension";
WignerFunctionMax::usage = "The max value of the Wigner function";

\[Mu]::usage = "\[Mu][r] gives a component of the squeezing matrices";
\[Nu]::usage = "\[Nu][r,\[Phi]] gives a component of the squeezing matrices";
Rsq::usage = "Rsq[r,\[Phi]] gives a component of the squeezing matrices";

\[CapitalSigma]forSingleModeSqueezing::usage = "\[CapitalSigma]forSingleModeSqueezing[r,\[Phi]] gives the single mode squeezing symplectic matrix";
\[CapitalSigma]forTwoModeSqueezing::usage = "\[CapitalSigma]forTwoModeSqueezing[r,\[Phi]] gives the two mode squeezing symplectic matrices";

\[Sigma]ThermalState::usage = "\[Sigma]ThermalState[nth] gives the covariance matrix of a thermal state";
\[Sigma]TwoModeThermalState::usage = "\[Sigma]TwoModeThermalState[nth] gives the covariance matrix of a two mode thermal state.";
quadThermalState::usage = "quadThermalState[] gives the trivial quadrature vector for a thermal state";
quadTwoModeThermalState::usage = "quadTwoModeThermalState[] gives the trivial quadrature vector for a two mode thermal state";

\[CapitalSigma]forSingleModeSqueezing::usage = "Symplectic matrix for single mode squeezing";
\[CapitalSigma]forTwoModeSqueezing::usage = "Symplectic matrix for two mode squeezing";

Nsq::usage = "Mean photon number in a squeezed state";

DBtoR::usage = "Helper function to convert dB of squeezing to r";
RtoDB::usage = "Helper function to convert dB of squeezing to r";

quadSingleModeSqueezedState::usage = "quadSingleModeSqueezedState[] gives the trivial quadrature vector for a vacuum state";
\[Sigma]SingleModeSqueezedState::usage = "\[Sigma]SingleModeSqueezedState[r,\[Phi]] gives the covariance matrix of a single mode squeezed vacuum state";

quadTwoModeSqueezedState::usage = "quadTwoModeSqueezedState gives the trivial quadrature vector for a two mode squeezed vacuum state";
\[Sigma]TwoModeSqueezedState::usage = "\[Sigma]TwoModeSqueezedState[r,\[Phi]] gives the covariance matrix of a two mode squeezed vacuum state";


\[Sigma]Loss::usage = "\[Sigma]Loss[\[Eta],\[Sigma]State,m,n] calculates covariance matrix after loss in one (specified) mode";
quadLoss::usage = "quadLoss[\[Eta],quadState,m,n] calculates quadrature vector after loss in one (specified) mode";

\[CapitalSigma]forAmplifying::usage = "The covariance matrix for amplification";
\[Sigma]Amp::usage = "\[Sigma]Amp[G,\[Sigma]State,m,n] calculates covariance matrix after amplification in one (specified) mode";
quadAmp::usage = "\[Sigma]Amp[G,\[Sigma]State,m,n] calculates covariance matrix after amplification in one (specified) mode";

\[Sigma]SecondModeGain::usage = "Calculated covariance matrix after phase insensitive amplification in second mode";
quadSecondModeGain::usage = "Calculated quadrature vector after phase insensitive amplification in second mode";


\[Delta]\[Sigma]::usage = "Helper function for QFI";
\[Delta]Quadrature::usage = "Helper function for QFI";
QFI::usage = "Calculate quantum fisher information from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";
QFIDisplacementOnly::usage = "Calculate displacement contribution to quantum fisher information from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";
QFICovarianceOnly::usage = "Calculate covariance contribution to quantum fisher information from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";
QFICovarianceOnlyEvaluated::usage = "Calculate covariance contribution to quantum fisher information from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";

QFIDisplacementOnlyMulti::usage = "QFIDisplacementOnlyMulti[quadratureQFI, \[Sigma]QFI, n, var1, var2] is my quick gueess at the form of the displacement QFI for multi parameter estimation"

QFIIsothermal::usage = "Calculate covariance contribution to quantum fisher information FOR ISOTHERMAL STATES (INC PURE STATES) from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";
QFIIsothermalCovarianceOnly::usage = "Calculate covariance contribution to quantum fisher information FOR ISOTHERMAL STATES (INC PURE STATES) from the quadrature, covariance, the dimension. Var is the variable to differentiate wrt.";
IsothermalCheck::usage = "Check if covariance matrix is singular (need to check, think this means pure state)";
IsothermalComponent::usage = "Gives \[Nu] for the QFI Isothermal model";

PhotonVariance::usage = "Photon number variance for a single mode gaussian state";
MeanPhotonNum::usage = "Mean photon number for a single mode gaussian state";


Begin[ "Private`"]

(* Conventions used here (see Monras (Phase space formalism for quantum estimation of Gaussian states)
 / Ferraro (Gaussian states in quantum information) *)

(*\[Kappa]3= 1/Sqrt[2];
\[Kappa]1= 1/Sqrt[2];
\[Kappa]2= \[Kappa]1;*)

(* NOTE I HAVE CHANGED THE CONVENTION TO THE QUANTUM OPTICS CONVENTION *)
\[Kappa]1[] := 1;
\[Kappa]2[] := \[Kappa]1[];
\[Kappa]3[] := 1/(2 \[Kappa]1[]);

(*qk = (1/(2\[Kappa]1[]))(ak +ak\[Dagger])*)

P[n_]:= Table[If[k<=n,KroneckerDelta[2k-1,l],0]+If[k>n,KroneckerDelta[2(k-n),l],0],{k,1,2n},{l,1,2n}];
J[n_]:=Table[-KroneckerDelta[k+n,l] +KroneckerDelta[l+n,k],{k,1,2n},{l,1,2n}];
\[Omega][n_]:=-Transpose[P[n]].J[n].P[n];
\[Lambda]To\[CapitalLambda][\[Lambda]_]:=\[Lambda]/(\[Kappa]1[]); (*This is slightly different from Ferraro, but is the only way that makes sense to me *)

(* Helper functions*)

MatDirSum[sqMatrices_]:=Module[{t,dims,rowPaddings},dims=Prepend[Length/@sqMatrices,0];t=Total[dims];
rowPaddings=Rest[FoldList[#1+{#2[[1]],-#2[[2]]}&,{0,t},Partition[dims,2,1]]];
Join@@MapThread[ArrayPad[#1,{{0},#2}]&,{sqMatrices,rowPaddings}]];

ActOnCovMatWithSympMat[SympMat_,CovMat_]:=SympMat.CovMat.Transpose[SympMat];
ActOnRMatWithSympMat[SympMat_,RMat_]:=SympMat.RMat;

UnitaryToSimplect[Unitary_]:=ComplexExpand[Transpose[P[Dimensions[Unitary][[1]]]]. ArrayFlatten[{{Re[Unitary],-Im[Unitary]},{Im[Unitary],Re[Unitary]}}].P[Dimensions[Unitary][[1]]]]

singleModeToMultiMode[\[CapitalSigma]singleModeMat_,i_,n_]:=Module[{m},m=IdentityMatrix[2n];m[[(2i-1);;(2i),(2i-1);;(2i)]]=\[CapitalSigma]singleModeMat[[1;;2,1;;2]];m];
twoModeToMultiMode[\[CapitalSigma]twoModeMat_,i_,j_,n_]:=Module[{m},m=IdentityMatrix[2n];m[[(2i-1);;(2i),(2i-1);;(2i)]]=\[CapitalSigma]twoModeMat[[1;;2,1;;2]];m[[(2i-1);;(2i),(2j-1);;(2j)]]=\[CapitalSigma]twoModeMat[[1;;2,3;;4]];m[[(2j-1);;(2j),(2i-1);;(2i)]]=\[CapitalSigma]twoModeMat[[3;;4,1;;2]];m[[(2j-1);;(2j),(2j-1);;(2j)]]=\[CapitalSigma]twoModeMat[[3;;4,3;;4]];m];

(* Unitary operations *)

rotMat[\[Theta]_] :={{Cos[\[Theta]], Sin[\[Theta]]},{-Sin[\[Theta]],Cos[\[Theta]]}};
Displacement[R_,\[Lambda]_]:=R+ \[Lambda]To\[CapitalLambda][\[Lambda]];

SymBSUnitary[\[Eta]_]:={{Sqrt[\[Eta]],0,Sqrt[1-\[Eta]],0},{0,Sqrt[\[Eta]],0,Sqrt[1-\[Eta]]},{-Sqrt[1-\[Eta]],0,Sqrt[\[Eta]],0},{0,-Sqrt[1-\[Eta]],0,Sqrt[\[Eta]]}};
SymBSUnitaryPhase[\[Phi]_]:={{Cos[\[Phi]],0,Sin[\[Phi]],0},{0,Cos[\[Phi]],0,Sin[\[Phi]]},{-Sin[\[Phi]],0,Cos[\[Phi]],0},{0,-Sin[\[Phi]],0,Cos[\[Phi]]}};
RotateSingleModeQuadrature[\[Theta]_,i_,n_]:=singleModeToMultiMode[rotMat[\[Theta]],i,n];

(* Wigner functions *)

WignerFunction[phaseSpacePosn_,quad_,\[Sigma]_,n_] := 1/((2\[Pi])^n \[Kappa]2[]^(2n) Sqrt[Det[\[Sigma]]]) Exp[-(1/2)Transpose[(phaseSpacePosn-quad)].Inverse[\[Sigma]].(phaseSpacePosn-quad)]
WignerFunctionMax[\[Sigma]_,n_] := 1/((2\[Pi])^n \[Kappa]2[]^(2n) Sqrt[Det[\[Sigma]]]);

(* Thermal states and coherent states *)

\[Sigma]ThermalState[nth_] :=1/(4\[Kappa]1[]^2) (2nth +1) IdentityMatrix[2];
quadThermalState[] := {{0},{0}};
\[Sigma]TwoModeThermalState[nth_] := MatDirSum[{\[Sigma]ThermalState[nth],\[Sigma]ThermalState[nth]}];
quadTwoModeThermalState[] := {{0},{0},{0},{0}};

(* Squeezed states *)

\[Mu][r_]:= Cosh[r];
\[Nu][r_,\[Phi]_]:= E^(I \[Phi]) Sinh[r];
Rsq[r_,\[Phi]_] := ComplexExpand[{{Re[\[Nu][r,\[Phi]]],Im[\[Nu][r,\[Phi]]]},{Im[\[Nu][r,\[Phi]]],-Re[\[Nu][r,\[Phi]]]}}];

\[CapitalSigma]forSingleModeSqueezing[r_,\[Phi]_] := \[Mu][r] IdentityMatrix[2] + Rsq[r,\[Phi]];
\[CapitalSigma]forTwoModeSqueezing[r_,\[Phi]_] := ArrayFlatten[{{\[Mu][r] IdentityMatrix[2],Rsq[r,\[Phi]]},{Rsq[r,\[Phi]],\[Mu][r] IdentityMatrix[2]}}];

Nsq[\[Alpha]_,nth_,r_]:=Abs[\[Alpha]]^2+nth*Cosh[2r]+Sinh[r]^2;(*Mean photon number in a squeezed state*)

DBtoR[db_]:=-0.5*Log[10^(-0.1*db)];
RtoDB[r_]:=-10*Log[10,E^(-2*r)];

(* Single mode squeezed state *)

quadSingleModeSqueezedState[] := quadThermalState[];
\[Sigma]SingleModeSqueezedState[r_,\[Phi]_] :=FullSimplify[ActOnCovMatWithSympMat[\[CapitalSigma]forSingleModeSqueezing[r,\[Phi]],\[Sigma]ThermalState[0]]];

(* Two mode squeezed state *)

quadTwoModeSqueezedState[r_,\[Phi]_] := quadTwoModeThermalState[];
\[Sigma]TwoModeSqueezedState[r_,\[Phi]_] := Simplify[ActOnCovMatWithSympMat[\[CapitalSigma]forTwoModeSqueezing[r,\[Phi]],\[Sigma]TwoModeThermalState[0]]];

(* Loss *)

\[Sigma]Loss[\[Eta]_,\[Sigma]State_,m_,n_]:= (ActOnCovMatWithSympMat[twoModeToMultiMode[SymBSUnitary[\[Eta]],m,n+1,n+1],MatDirSum[{\[Sigma]State,\[Sigma]ThermalState[0]}]])[[1;;2n,1;;2n]];
quadLoss[\[Eta]_,quadState_,m_,n_]:=(ActOnRMatWithSympMat[twoModeToMultiMode[SymBSUnitary[\[Eta]],m,n+1,n+1],Join[quadState,quadThermalState[]]])[[1;;2n]];

(* Amplification *)
(*\[CapitalSigma] = Simplify[\[CapitalSigma]forTwoModeSqueezing[ ArcCosh[Sqrt[G]],0],Assumptions\[Rule]{G\[GreaterEqual]1}]*)
\[CapitalSigma]forAmplifying[G_] :={{Sqrt[G],0,Sqrt[-1+G],0},{0,Sqrt[G],0,-Sqrt[-1+G]},{Sqrt[-1+G],0,Sqrt[G],0},{0,-Sqrt[-1+G],0,Sqrt[G]}};

\[Sigma]Amp[G_,\[Sigma]State_,m_,n_]:=(ActOnCovMatWithSympMat[twoModeToMultiMode[\[CapitalSigma]forAmplifying[G],m,n+1,n+1],MatDirSum[{\[Sigma]State,\[Sigma]ThermalState[0]}]])[[1;;2n,1;;2n]];
quadAmp[G_,quadState_,m_,n_]:=(ActOnRMatWithSympMat[twoModeToMultiMode[\[CapitalSigma]forAmplifying[G],m,n+1,n+1],Join[quadState,quadThermalState[]]])[[1;;2n]];

\[Sigma]SecondModeGain[G_,nth_,\[Sigma]State_]:=(ActOnCovMatWithSympMat[\[CapitalSigma]forAmplifying[G],MatDirSum[{\[Sigma]State,\[Sigma]ThermalState[nth]}]])[[1;;4,1;;4]];
quadSecondModeGain[G_,quadState_] := {1,1,Sqrt[G],Sqrt[G]}*quadState;



(* QFI stuff NOTE THAT HAVENT PROPERLY INCORPORATED THE KAPPAS*)

\[Delta]\[Sigma][\[Sigma]Mat_,var_] := D[Transpose[{Flatten[\[Sigma]Mat]}], var];
\[Delta]Quadrature[quadrature_,var_] := D[quadrature, var];

(* MONRAS DEFINITION - different to what we need since the covariance matrix is defined differently by a factor of 1/2. \[Sigma]Alex = 2 \[Sigma]Matteo.
 Monras: QFI[quadratureQFI_,\[Sigma]QFI_]:=1/2Tr[Transpose[\[Delta]\[Sigma][\[Sigma]QFI]].(Inverse[KroneckerProduct[ \[Sigma]QFI,\[Sigma]QFI] - KroneckerProduct[ \[Omega][1],\[Omega][1]] ]).\[Delta]\[Sigma][\[Sigma]QFI]]+ 2Transpose[\[Delta]Quadrature[quadratureQFI]].(Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI])//FullSimplify; *)
(* THESE SEEM WRONG (PH JAN15) QFI[quadratureQFI_, \[Sigma]QFI_, n_,var_] :=(4\[Kappa]1[]^2)^2/2 Tr[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]].(Inverse[4 KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).\[Delta]\[Sigma][\[Sigma]QFI,var]] + 2/(4\[Kappa]1[]^2) Transpose[\[Delta]Quadrature[quadratureQFI,var]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var];

QFIDisplacementOnly[quadratureQFI_, \[Sigma]QFI_, n_,var_] := 2/(4\[Kappa]1[]^2) Transpose[\[Delta]Quadrature[quadratureQFI,var]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var];

QFICovarianceOnly[\[Sigma]QFI_, n_,var_] := (4\[Kappa]1[]^2)^2/2 Tr[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]].(Inverse[4 KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).\[Delta]\[Sigma][\[Sigma]QFI,var]];

QFICovarianceOnlyEvaluated[\[Sigma]QFI_, n_,var_,varVal_] := ((4\[Kappa]1[]^2)^2/2 Tr[N[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]]/.var->varVal].(Inverse[4 N[KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI]/.var->varVal] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).N[\[Delta]\[Sigma][\[Sigma]QFI,var]/.var->varVal]]);
*)
QFI[quadratureQFI_, \[Sigma]QFI_, n_,var_] :=8 (\[Kappa]1[]^4) Tr[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]].(Inverse[16 (\[Kappa]1[]^4) KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).\[Delta]\[Sigma][\[Sigma]QFI,var]] + (Transpose[\[Delta]Quadrature[quadratureQFI,var]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var])[[1]];

QFIDisplacementOnly[quadratureQFI_, \[Sigma]QFI_, n_,var_] := (Transpose[\[Delta]Quadrature[quadratureQFI,var]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var])[[1]];

QFICovarianceOnly[\[Sigma]QFI_, n_,var_] := 8 (\[Kappa]1[]^4) Tr[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]].(Inverse[16(\[Kappa]1[]^4) KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).\[Delta]\[Sigma][\[Sigma]QFI,var]];

QFICovarianceOnlyEvaluated[\[Sigma]QFI_, n_,var_,varVal_] :=  8 (\[Kappa]1[]^4)  Tr[N[Transpose[\[Delta]\[Sigma][\[Sigma]QFI,var]]/.var->varVal].(Inverse[16 (\[Kappa]1[]^4) N[KroneckerProduct[ \[Sigma]QFI, \[Sigma]QFI]/.var->varVal] - KroneckerProduct[ \[Omega][n], \[Omega][n]] ]).N[\[Delta]\[Sigma][\[Sigma]QFI,var]/.var->varVal]];

QFIDisplacementOnlyMulti[quadratureQFI_, \[Sigma]QFI_, n_,var1_,var2_] := (Transpose[\[Delta]Quadrature[quadratureQFI,var1]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var2])[[1]];

IsothermalCheck[\[Sigma]QFI_, n_] := (4\[Kappa]1[]^2)^2 Simplify[(\[Sigma]QFI.\[Omega][n]).(\[Sigma]QFI.\[Omega][n])];
IsothermalComponent[\[Sigma]QFI_, n_] := Sqrt[-(4\[Kappa]1[]^2)^2 Simplify[(\[Sigma]QFI.\[Omega][n]).(\[Sigma]QFI.\[Omega][n])][[1,1]]];
(* Isothermal models (see Monras), includes pure states *)
QFIIsothermal[quadratureQFI_,\[Sigma]QFI_,var_,\[Nu]_]:= 1/2 \[Nu]^2/(1+\[Nu]^2) Tr[Simplify[D[\[Sigma]QFI,var].Inverse[\[Sigma]QFI].D[\[Sigma]QFI,var].Inverse[\[Sigma]QFI]]] + (Transpose[\[Delta]Quadrature[quadratureQFI,var]].Inverse[\[Sigma]QFI].\[Delta]Quadrature[quadratureQFI,var])[[1]];
QFIIsothermalCovarianceOnly[\[Sigma]QFI_,var_,\[Nu]_]:= 1/2 \[Nu]^2/(1+\[Nu]^2) Tr[Simplify[D[\[Sigma]QFI,var].Inverse[\[Sigma]QFI].D[\[Sigma]QFI,var].Inverse[\[Sigma]QFI]]];

(* Photon number and variance for Gaussian states*)

(*A useful paper is http://pra.aps.org/pdf/PRA/v49/i4/p2993_1. 
It has expressions for the photon number and photon number variance for general squeezed states
 and gaussian states.*)

(*For a SINGLE MODE general Gaussian State:*)
T[Cov_]:=Tr[Cov];
d[Cov_]:=Det[Cov];
PhotonVariance[Cov_,Quad_]:= ( (4 \[Kappa]1[]^4) T[Cov]^2/2 - (4 \[Kappa]1[]^4) d[Cov]-1/4 + (4 \[Kappa]1[]^4)Transpose[Quad].Cov.Quad)[[1,1]]
MeanPhotonNum[Cov_,Quad_]:=(((2 \[Kappa]1[]^2) T[Cov]-1)/2 + 1/2 (2 \[Kappa]1[]^2) Transpose[Quad].Quad)[[1,1]]

End[]

EndPackage[]









