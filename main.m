(* ::Package:: *)

(* ::Section:: *)
(*Kvadratick\[YAcute] \[CHacek]len rozvoje*)


(* ::Input:: *)
(*parametry = {g, l, \[Theta]init, \[Omega]init};*)
(*k = Sin[\[Theta]init / 2];*)


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
konst\[Omega]init = 0;
sol = Reap[ (* "Reap" posb\[IAcute]r\[AAcute] hodnoty z\[IAcute]skan\[EAcute] v "Sow" *)(*Poj\[DHacek]me na\[SHacek]e nov\[EAcute] funkce zna\[CHacek]it s mal\[YAcute]mi p\[IAcute]smeny*)
			NDSolve[ (* \[CapitalRHacek]e\[SHacek]en\[IAcute] ODR *)
					{
							\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0,
							\[Theta][0] == konst\[Theta]init,
							\[Theta]'[0] == konst\[Omega]init
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
solLinearlyEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1.1895988355758171/100, 
							Method->{"TimeIntegration"->"LinearlyImplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 
printResults = False; 


solLinearlyEuler
Plot[{Evaluate[\[Theta][t] /. solLinearlyEuler],  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"solLinearlyEuler, LinearlyImplicitEuler", "NumSol, Default settings"}]


(*pouzijeme metodu ImplicitRungeKutta -> vyjde stejne jako metoda whenEvent*)
solImpRungeKutta:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1/100,
							Method->"ImplicitRungeKutta",  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False; 
solImpRungeKutta
Plot[{Evaluate[\[Theta][t] /. solImpRungeKutta] ,  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"ImpRungeKutta", "NumSol, Default settings"}]


(*pouzijeme metodu ExplicitEuler*)
solExpEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->2*\[Pi]/1000, 
							Method->{"TimeIntegration"->"ExplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False
solExpEuler 
Plot[{Evaluate[\[Theta][t] /. solExpEuler],  Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi), -10/(2*Pi)}, {t, 0, 4\[Pi]},
	  PlotLegends->{"solExpEuler, ExplicitEuler", "NumSol, Default settings"}]





(* ::Section:: *)
(*Porovn\[AAcute]n\[IAcute] analytick\[EAcute] aproximace a numerick\[EAcute]ho \[RHacek]e\[SHacek]en\[IAcute] (\[UAcute]kol 5)*)


\[Omega]nula=Sqrt[konstg/konstl];
\[Omega]PL = Sqrt[\[Omega]nula^2/(1+((konst\[Theta]init^2)/8))];
solPoincareLindstedt=konst\[Theta]init*Cos[\[Omega]PL*t]+konst\[Theta]init^3/192*(Cos[3*\[Omega]PL*t]-Cos[\[Omega]PL*t]);


Plot[{solPoincareLindstedt,Evaluate[\[Theta][t] /. NumSol]},{t,0,4\[Pi]}]


TseriesAprox /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1   (*Analytick\[AAcute] aproximace*)
TIntegralExplicit /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1 (*P\[RHacek]es eliptick\[YAcute] integr\[AAcute]l*)
TNumSol (*Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute]*)
TPoincareLindstedt=2\[Pi]/\[Omega]PL  (*Aproximace Poincar\[EAcute]-Lindstedt*)


(* ::Section:: *)
(*Dal\[SHacek]\[IAcute] \[CHacek]len v rozvoji (\[UAcute]kol 6)*)


\[Theta]series[t_] := \[Epsilon]*\[Theta]A[t] + \[Epsilon]^3*\[Theta]B[t] + \[Epsilon]^5*\[Theta]C[t] (* Potrebujeme nejvyssi clen eps^5 *)
\[Omega]02 = \[Omega]2 - \[Epsilon]^2 * \[Alpha] - \[Epsilon]^4*\[Beta] (* eps^6 uz nepotrebujeme a bereme to rovnou s minusem*)
sinSeries = Normal[Series[Sin[\[Theta]],{\[Theta],0,5}]] /. \[Theta] -> \[Theta]series[t] (* Normal tu je z nejakych obskurnich duvodu ze mocniny nemuzou byt variable *)


\[Theta]series''[t] + \[Omega]02 * sinSeries


koef = CoefficientList[ExpandAll[\[Theta]series''[t] + \[Omega]02 * sinSeries],\[Epsilon]]; (* Indexovani je o 1 posunute, protoze prvni odpovida nulte mocnine *)
koef[[2]] == 0 // TraditionalForm
koef[[4]] == 0 // TraditionalForm
koef[[6]] == 0 // TraditionalForm


(* Zbab\[EHacek]le dosad\[IAcute]me \[RHacek]e\[SHacek]en\[IAcute] z dokumentu *)
rceProThetaC = koef[[6]] \
	/. \[Theta]A[t] -> \[Theta]tI * Cos[Sqrt[\[Omega]2]*t] \
	/. \[Theta]B[t] -> \[Theta]tI^3 / 192 * { Cos[3*Sqrt[\[Omega]2]*t] - Cos[Sqrt[\[Omega]2]*t] } \
	/. \[Alpha] -> - \[Omega]2 * \[Theta]tI^2 / 8 ;
rceProThetaC // Simplify // TraditionalForm


thetaCsol = DSolve[
{
	rceProThetaC[[1]] == 0,
	\[Theta]C[0] == 0, (* Pocatecni podminka *)
	\[Theta]C'[0] == 0
},\[Theta]C[t],t] // Simplify // TraditionalForm


(* Vid\[IAcute]me, \[ZHacek]e se tu vyskytuje op\[EHacek]t diverguj\[IAcute]c\[IAcute] \[CHacek]len  *)
90 t (512 \[Beta]+5 \[Theta]tI^4 \[Omega]2) // TraditionalForm
(* Situaci op\[EHacek]t zachr\[AAcute]n\[IAcute]me polo\[ZHacek]en\[IAcute]m z\[AAcute]vorky rovn\[EAcute] nule a \[RHacek]e\[SHacek]en\[IAcute] je *)
Style[\[Beta] == -(5/512)HoldForm[\!\(TraditionalForm\`\*
SubsuperscriptBox[
OverscriptBox["\[Theta]", "~"], 
StyleBox[
RowBox[{"i", "n", "i", "t"}], "TI"], "4"] 
\*SuperscriptBox[\(\[Omega]\), \(2\)]\)] ,FontSize->30]// TraditionalForm


(* Finalni reseni *)
thetaCsol /. \[Beta] -> -5/512 *\[Theta]tI^4 *\[Omega]2 /. Sqrt[\[Omega]2] -> \[Omega] /. 1/Sqrt[\[Omega]2] -> 1/\[Omega] /. \[Theta]tI -> HoldForm[\!\(TraditionalForm\`SubscriptBox[OverscriptBox["\<\[Theta]\>", "\<~\>"], StyleBox[RowBox[{"\<i\>", "\<n\>", "\<i\>", "\<t\>"}], "\<TI\>"]]\)]


(* ::Section:: *)
(*Marn\[EAcute] pokusy o \[RHacek]e\[SHacek]en\[IAcute] 6. \[UAcute]kolu p\[RHacek]es mathematicu*)


sol61 = DSolve[{koef[[2]]==0,
	\[Theta]A[0] == \[Theta]tI, (* Pocatecni podminka *)
	\[Theta]A'[0] == 0
},\[Theta]A,t]


\[Theta]A /. sol61


f[1]



DSolve[{koef[[4]]==0,
	\[Theta]A[t] == \[Theta]A /. sol61,
	\[Theta]B[0] == 0, (* Pocatecni podminka *)
	\[Theta]B'[0] == 0
},
\[Theta]B[t],t]


sol = DSolve[
{
	koef[[2]] == 0,
	koef[[4]] == 0,
	\[Theta]A[0] == \[Theta]tI, (* Pocatecni podminka *)
	\[Theta]A'[0] == 0,
	\[Theta]B[0] == 0, (* Pocatecni podminka *)
	\[Theta]B'[0] == 0
},
{\[Theta]A[t],\[Theta]B[t]},t]


\[Theta]B[t] /. sol /. \[Alpha] -> - \[Omega]2*\[Theta]tI^2 /8
