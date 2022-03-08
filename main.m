(* ::Package:: *)

(* ::Section::Closed:: *)
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
(*seriesIntegrand = Series[Integrand, {\[Theta]init, 0,4}]   (*Dle toho \[UAcute]kolu m\[AAcute]me rov\[IAcute]jet do jednoho \[RHacek]\[AAcute]du nav\[IAcute]c.*)*)


(* ::Input:: *)
(*(* Perioda (integr\[AAcute]l) p\[RHacek]es Taylor\[URing]v rozvoj -- vychazi stejne jako ve skriptech! *)*)


(* ::Input:: *)
(*TseriesAprox = 4*Sqrt[l/g]*Integrate[seriesIntegrand, {\[CurlyPhi], 0, \[Pi]/2}]*)


TseriesAprox /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1


(* Koeficient u \[Theta]^2 *)
(* SeriesCoefficient[TseriesAprox,{\[Theta]init,0,2}] *)


(* Explicitne spocteny elipticky integral *)
TIntegralExplicit = 4*Sqrt[l/g]Integrate[1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2], {\[CurlyPhi], 0 , \[Pi]/2}]


TIntegralExplicit /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1


(* ::Section:: *)
(*Numericky*)


konstg = 9.81;
konstl = 1;
konst\[Theta]init = 10/(2\[Pi]);    (*Tady se dosazovalo do glob\[AAcute]ln\[IAcute] prom\[EHacek]nn\[EAcute], dle kter\[EAcute] se neho\[RHacek]e rozv\[IAcute]j\[IAcute] a pak to d\[EHacek]lalo blbosti, je lep\[SHacek]\[IAcute] to pojmenovat jinak*)
\[Omega]init = 0;
sol = Reap[ (* "Reap" posb\[IAcute]r\[AAcute] hodnoty z\[IAcute]skan\[EAcute] v "Sow" *)(*Poj\[DHacek]me na\[SHacek]e nov\[EAcute] funkce zna\[CHacek]it s mal\[YAcute]mi p\[IAcute]smeny*)
			NDSolve[ (* \[CapitalRHacek]e\[SHacek]en\[IAcute] ODR *)
					{
							\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0,
							\[Theta][0] == konst\[Theta]init,
							\[Theta]'[0] == \[Omega]init
					},
							\[Theta][t],
							{t, 0, 4\[Pi]},
						Method-> {
								"EventLocator", 
								"Event"-> \[Theta][t],
								"EventAction" :> Sow[t]
								}
					]
		];
NumSol = sol[[1]]
ZeroPoints = sol[[2]]


Plot[Evaluate[\[Theta][t] /. NumSol],{t,0,4\[Pi]}]


(* Vzdalenosti mezi nulovymi body jsou vsude stejne! *)
DifsZeroPoints = Differences[ZeroPoints[[1]]]
(* Vzdalenosti odpovidaji ctvrtine periody *)(*Polovin\[EHacek] ne? To pak vyjde i stejn\[AAcute] hodnota*)
TNumSol = 2*Mean[DifsZeroPoints]











(* ::Section:: *)
(*Numericky - jin\[EAcute] metody (\[UAcute]kol 3)*)


(*pouzijeme metodu LinearlyImplicitEuler*)
solLinearlyEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == \[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1.1895988355758171/100, 
							Method->{"TimeIntegration"->"LinearlyImplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 
printResults = False; 


solLinearlyEuler
Plot[{Evaluate[\[Theta][t] /. solLinearlyEuler],  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"solLinearlyEuler, LinearlyImplicitEuler", "NumSol, Default settings"}]


(*pouzijeme metodu ImplicitRungeKutta -> vyjde stejne jako metoda whenEvent*)
solImpRungeKutta:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == \[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1/100,
							Method->"ImplicitRungeKutta",  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False; 
solImpRungeKutta
Plot[{Evaluate[\[Theta][t] /. solImpRungeKutta] ,  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"ImpRungeKutta", "NumSol, Default settings"}]


(*pouzijeme metodu ExplicitEuler*)
solExpEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == \[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->2*\[Pi]/1000, 
							Method->{"TimeIntegration"->"ExplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False
solExpEuler 
Plot[{Evaluate[\[Theta][t] /. solExpEuler],  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"solExpEuler, ExplicitEuler", "NumSol, Default settings"}]



