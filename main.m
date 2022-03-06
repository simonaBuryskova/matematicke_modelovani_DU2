(* ::Package:: *)

(* ::Section:: *)
(*Kvadratick\[YAcute] \[CHacek]len rozvoje*)


(* ::Input:: *)
(*parametry = {g , l , \[Theta]init, \[Omega]init };*)
(*k = Sin[\[Theta]init/2];*)


(* ::Input:: *)
(*(* Explicitni rozvoj periody *)*)
(*(* period = 2*\[Pi]*Sqrt[l/g]*(1+a*\[Theta]init^2+b*\[Theta]init^4+O[\[Theta]init]^6) *)*)


(* ::Input:: *)
(*(* Tayloruv rozvoj integrandu *)*)
(*Integrand = 1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2];*)
(*seriesIntegrand = Series[Integrand, {\[Theta]init, 0,2}]*)


(* ::Input:: *)
(*(* Perioda (integr\[AAcute]l) p\[RHacek]es Taylor\[URing]v rozvoj -- vychazi stejne jako ve skriptech! *)*)


(* ::Input:: *)
(*TseriesAprox = 4*Sqrt[l/g]*Integrate[seriesIntegrand, {\[CurlyPhi], 0, \[Pi]/2}]*)


(* Koeficient u \[Theta]^2 *)
(* SeriesCoefficient[TseriesAprox,{\[Theta]init,0,2}] *)


(* Explicitne spocteny elipticky integral *)
TIntegralExplicit = 4*Sqrt[l/g]Integrate[1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2], {\[CurlyPhi], 0 , \[Pi]/2}]


(* ::Section:: *)
(*Numericky*)


g = 9.81;
l = 1;
\[Theta]init = 10/(2\[Pi]);
\[Omega]init = 0;
Sol = Reap[ (* "Reap" posb\[IAcute]r\[AAcute] hodnoty z\[IAcute]skan\[EAcute] v "Sow" *)
	NDSolve[ (* \[CapitalRHacek]e\[SHacek]en\[IAcute] ODR *)
		{
			\[Theta]''[t] + g/l * Sin[\[Theta][t]] == 0,
			\[Theta][0] == \[Theta]init,
			\[Theta]'[0] == \[Omega]init
		},
		\[Theta][t],
		{t, 0, 4\[Pi]},
		Method->{
			"EventLocator", 
			"Event"-> \[Theta][t],
			"EventAction" :> Sow[t]
		}
	]
];
NumSol = Sol[[1]]
ZeroPoints = Sol[[2]]


Plot[Evaluate[\[Theta][t] /. NumSol],{t,0,4\[Pi]}]


(* Vzdalenosti mezi nulovymi body jsou vsude stejne! *)
DifsZeroPoints = Differences[ZeroPoints[[1]]]
(* Vzdalenosti odpovidaji ctvrtine periody *)
TNumSol = 4*Mean[DifsZeroPoints]



