#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Basic Functionality Package*)


BeginPackage["BasicFunctions`"];


(* ::Chapter:: *)
(*Function Usage*)


USolve::usage="Solve for Time evolution operator given Hamiltonian";
DensityMatrix::usage="Outputs density matrix given a wave function(Pure State for now)";
DensityMatrixQ::usage="Returns true if input is a density matrix";
WaveFunction::usage="Creates a wave function from phase and azimuthal angle.";
BlochVector::usage="Creates the Bloch vector for a Density matrix";
EvolveDensityMatrix::usage="Evaluates the inputed density matrix by the inpute time evolution operator.";
BlochPlot::usage="Creates a plot of the bloch vector for the given density matrix under given time evolution operator.";
RandomDensityMatrix::usage="Generates a random density matrix.";
PureStateQ::usage="Outputs true if a given density matrix represents a pure state.";
BlochListPlot::usage="Plots a list of points on the bloch sphere given a list of density matricies.";
SquarePulseHamiltonian::usage="Generates the hamiltonian of a square pulse.";
SquarePulseSoln::usage="Numerically solves the hamiltonian for a given square pulse.";
TransitionProbibility::usage="Calculates the probability of a transition between two states."
HGate::usage="";
SGate::usage="";
TGate::usage="";


Begin["`Private`"];


\[HBar]=1;
(* Common Gates *)
HGate[t_]:=1/Sqrt[2] ({
 {1, 1},
 {1, -1}
});
SGate[t_]:=({
 {1, 0},
 {0, I}
});
TGate[t_]:=({
 {1, 0},
 {0, Exp[I Pi/4]}
});

(* How to approximately check PositiveSemidefinite *)
DensityMatrixQ[A_?SquareMatrixQ]:=Block[{\[Epsilon]=10^(-4)},(1-\[Epsilon]<Tr[A]<1+\[Epsilon])&&(PositiveSemidefiniteMatrixQ[A])&&(HermitianMatrixQ[A])];

WaveFunction::imaginaryMessage="Non-trivial imaginary part";
WaveFunction[\[Theta]_,\[Phi]_]:=Block[{},
If[ Im[\[Theta]]>10^(-6)||  Im[\[Phi]]>10^(-6),Message[WaveFunction::imaginaryMessage]];
Cos[\[Theta]/2]{{1},{0}}+Sin[\[Theta]/2]Exp[I \[Phi]]{{0},{1}}];

USolve[Ham_,UInit_,tMin_?NumericQ,tMax_?NumericQ,\[HBar]_?NumericQ,opts___]:=
	Block[{USol,tSol,intdF,t},
        intdF=USol/.NDSolve[{I \[HBar] USol'[tSol]==Ham[tSol].USol[tSol],USol[tMin]==UInit},USol,{tSol,tMin,tMax},opts]//First;
        Function[{t},intdF[t]//Evaluate]
    ];
    
DensityMatrix[\[Psi]_]:=\[Psi].ConjugateTranspose[\[Psi]];

BlochVector[\[Rho]_?DensityMatrixQ]:={
 \[Rho][[2,1]]+\[Rho][[1,2]],(\[Rho][[2,1]]-\[Rho][[1,2]])/I,2\[Rho][[1,1]]-1
}//Chop;

EvolveDensityMatrix[\[Rho]_?DensityMatrixQ,U_,t_]:=U[t].\[Rho].ConjugateTranspose[U[t]]//Chop;

BlochPlot[\[Rho]_?DensityMatrixQ,U_,tMin_,tMax_]:=Block[{a},
a[t_]:=BlochVector[EvolveDensityMatrix[\[Rho],U,t]];
Show[SphericalPlot3D[1,{\[Theta],0,Pi},{\[Phi],0,2Pi},PlotStyle->{Opacity[.2],Gray},MeshStyle->{Opacity[.1]}],ParametricPlot3D[{a[t][[1]],a[t][[2]],a[t][[3]]},{t,tMin,tMax},PlotStyle->{Thick,Red}]]];

RandomDensityMatrix[dim_]:=Block[{m=RandomComplex[{-1-I,1+I},{dim,dim}],h=(m+ConjugateTranspose[m])/2,\[Rho]},
 \[Rho]=h.ConjugateTranspose[h]/Tr[h.ConjugateTranspose[h]]
 ];

PureStateQ[\[Rho]_?DensityMatrixQ]:=\[Rho]==\[Rho].\[Rho];

BlochListPlot[list__,args___]:=Block[{points},
points=BlochVector/@list;
ListPointPlot3D[points,BoxRatios->{1, 2, 1}]];

SquarePulseHamiltonian[\[Tau]_,\[CapitalOmega]_,\[Omega]_,E0_,E1_,t_]:=Block[{H0=E0 {{1},{0}}.{{1,0}}+E1 {{0},{1}}.{{0,1}},Hp=\[CapitalOmega] Cos[\[Omega] t]({{1},{0}}.{{0,1}}+{{0},{1}}.{{1,0}})},H0+Hp( UnitStep[t]- UnitStep[t-\[Tau]])];

SquarePulseSoln[\[Tau]_?NumericQ,\[CapitalOmega]_?NumericQ,\[Omega]_?NumericQ,E0_?NumericQ,E1_?NumericQ]:=Block[{H},
H[t_]:=SquarePulseHamiltonian[\[Tau],\[CapitalOmega],\[Omega],E0,E1,t];
USolve[H,IdentityMatrix[2],0,20,\[HBar]]];

TransitionProbibility[t_,U_,i_,j_]:=Abs[U[t][[i,j]]]^2;



(* ::Section:: *)
(*End*)


End[]
EndPackage[]
