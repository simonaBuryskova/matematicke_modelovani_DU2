(* ::Package:: *)

(* ::Section::Closed:: *)
(*Kyvadlo*)
(*team Pendulum: Benda, Bury\[SHacek]kov\[AAcute], Hrube\[SHacek], Van\[EHacek]k*)


ClearAll[]
(*nacteme si balik MaTeX*)
ResourceFunction["MaTeXInstall"][] 
<<MaTeX`


(*font jako v LaTeXu*)
styleLaTeX = Directive[FontFamily->"Latin Modern Roman 10",FontSize->10, FontColor->Black];


(* ::Section::Closed:: *)
(*1 Dal\[SHacek]\[IAcute] \[CHacek]len rozvoje*)


(* ::Input:: *)
(*parametry = {g, l, \[Theta]init, \[Omega]init};*)
(*k = Sin[\[Theta]init / 2];*)


(* ::Input:: *)
(*(* Explicitni rozvoj periody *)*)
(*(* period = 2*\[Pi]*Sqrt[l/g]*(1+a*\[Theta]init^2+b*\[Theta]init^4+O[\[Theta]init]^6) *)*)


(* ::Input:: *)
(*(* Tayloruv rozvoj integrandu *)*)
(*Integrand = 1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2];*)
(*seriesIntegrand = Series[Integrand, {\[Theta]init, 0,4}]  *)


(* ::Input:: *)
(*(* Perioda (integr\[AAcute]l) p\[RHacek]es Taylor\[URing]v rozvoj -- vychazi stejne jako ve skriptech! *)*)


(* ::Input:: *)
(*TseriesAprox = 4*Sqrt[l/g]*Integrate[seriesIntegrand, {\[CurlyPhi], 0, \[Pi]/2}];*)


TseriesAprox /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1


(* Koeficient u \[Theta]^2 *)
(* SeriesCoefficient[TseriesAprox,{\[Theta]init,0,2}] *)


(* Explicitne spocteny elipticky integral *)
TIntegralExplicit = 4*Sqrt[l/g]Integrate[1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2], {\[CurlyPhi], 0 , \[Pi]/2}]


TIntegralExplicit /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1


(* ::Section:: *)
(*2 Perioda numericky*)


konstg = 9.81;
konstl = 1;
konst\[Theta]init = 10/(2\[Pi]); 
konst\[Omega]init = 0;
sol = Reap[ (* "Reap" posb\[IAcute]r\[AAcute] hodnoty z\[IAcute]skan\[EAcute] v "Sow" *)
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



Plot[{Evaluate[\[Theta][t] /. NumSol],10/(2*Pi),-10/(2*Pi)},{t,0,4\[Pi]},PlotLabel->"Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute], default. metoda", AxesLabel-> MaTeX/@ {t, \[Theta][t]}, 
PlotLegends-> {MaTeX["\\text{Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute]}"],MaTeX["\\theta_{max}"]},BaseStyle->styleLaTeX,Frame->False,FrameStyle->BlackFrame, PlotStyle->{Blue, {Blue,Dashed},{Blue,Dashed}}]


(* Vzdalenosti mezi nulovymi body jsou vsude stejne! *)
DifsZeroPoints = Differences[ZeroPoints[[1]]]
(* Vzdalenosti odpovidaji polovine periody *)
TNumSol = 2*Mean[DifsZeroPoints]


(* ::Section::Closed:: *)
(*3 Perioda numericky - jin\[EAcute] metody*)


(*pouzijeme metodu LinearlyImplicitEuler*)
solLinearlyEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1.1895988355758171/100, 
							Method->{"TimeIntegration"->"LinearlyImplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 
printResults = False; 


solLinearlyEuler
plotLE = Plot[Evaluate[\[Theta][t] /. solLinearlyEuler], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->Red];

plotAmp = Plot[{10/(2*Pi),-10/(2*Pi)}, {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->{{Black, Dotted}, {Black, Dotted}}];

plotNumsol = Plot[Evaluate[\[Theta][t] /. NumSol], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, PlotStyle->Blue];
Show[{plotLE, plotNumsol, plotAmp},PlotLabel->"Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute], default. metoda vs. LinearlyImplicitEuler",BaseStyle->styleLaTeX, AxesLabel->{t,\[Theta][t]}]


(*pouzijeme metodu ImplicitRungeKutta -> vyjde stejne jako metoda whenEvent*)
solImpRungeKutta:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->1/100,
							Method->"ImplicitRungeKutta",  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False; 
solImpRungeKutta
legendsSRK = {Style["ImpRungeKutta, ImplicitRungeKutta", styleLaTeX], Style["NumSol, Default settings", styleLaTeX]};

plotIRK = Plot[Evaluate[\[Theta][t] /. solImpRungeKutta], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->{Blue,Dashed}];

plotAmp = Plot[{10/(2*Pi),-10/(2*Pi)}, {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->{{Black, Dotted}, {Black, Dotted}}];

plotNumsol = Plot[Evaluate[\[Theta][t] /. NumSol], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, PlotStyle->Green];
Show[{plotNumsol, plotIRK,plotAmp},PlotLabel->"Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute], default. metoda vs. ImplicitRungeKutta",BaseStyle->styleLaTeX, AxesLabel->{t,\[Theta][t]}]





(*pouzijeme metodu ExplicitEuler*)
solExpEuler:= NDSolve[{\[Theta]''[t] + konstg/konstl * Sin[\[Theta][t]] == 0, \[Theta][0] == konst\[Theta]init, \[Theta]'[0] == konst\[Omega]init}, \[Theta][t], {t, 0, 4\[Pi]}, 
							StartingStepSize->2*\[Pi]/1000, 
							Method->{"TimeIntegration"->"ExplicitEuler"},  
							StepMonitor:>If[printResults, Print["Step number [", nSteps++,"]: \[Theta][t] = ",\[Theta][t],", t =",t], {}]]
nSteps = 1; 


printResults = False;
solExpEuler 
legendsIE = {Style["solExpEuler, ExplicitEuler", styleLaTeX], Style["NumSol, Default settings", styleLaTeX]};

plotIE = Plot[Evaluate[\[Theta][t] /. solExpEuler], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->Red];

plotAmp = Plot[{10/(2*Pi),-10/(2*Pi)}, {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->{{Black, Dotted}, {Black, Dotted}}];

plotNumsol = Plot[Evaluate[\[Theta][t] /. NumSol], {t, 0, 4\[Pi]},BaseStyle->styleLaTeX, PlotStyle->Blue];
Show[{plotIE, plotNumsol, plotAmp},PlotLabel->"Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute], default. metoda vs. ExplicitEuler",BaseStyle->styleLaTeX, AxesLabel->{t,\[Theta][t]}]



(* ::Section::Closed:: *)
(*4 Z\[AAcute]kon zachov\[AAcute]n\[IAcute] energie a numerick\[EAcute] metody*)


(*ClearAll["Global'*"]*)
h=0.01;
tMax = 100;
mtd ="ExplicitEuler";
\[Omega]0 = 1.7;

(*Metoda navr\[ZHacek]en\[AAcute] Mathematicou*)
goodSolution = NDSolve[{\[Omega]'[t]== -Sin[\[Theta][t]],\[Theta]'[t]==\[Omega][t],\[Omega][0]==\[Omega]0,\[Theta][0]==0 },{\[Theta][t], \[Omega][t]},{t,0,tMax}];
goodPlot = ParametricPlot[Evaluate[{\[Theta][t], \[Omega][t]} /. goodSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Blue,Dashed}];
(*Vl\[AAcute]stn\[IAcute] volba metody*)
badSolution =NDSolve[{\[Omega]'[t]== -Sin[\[Theta][t]],\[Theta]'[t]==\[Omega][t],\[Omega][0]==\[Omega]0,\[Theta][0]==0 },{\[Theta][t], \[Omega][t]},{t,0,tMax},StartingStepSize->h,Method->{"TimeIntegration"->mtd}];
badPlot =ParametricPlot[Evaluate[{\[Theta][t], \[Omega][t]} /. badSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->Red];
(*Porovn\[AAcute]n\[IAcute] v\[YAcute]voje energie*)
goodEnergyPlot = Plot[Evaluate[0.5*\[Omega][t]^2-Cos[\[Theta][t]]]/. goodSolution, {t, 0 ,tMax}, BaseStyle->styleLaTeX,PlotStyle->{Blue,Dashed}, PlotRange-> {{0,tMax},{0,1}}];
badEnergyPlot =Plot[Evaluate[0.5*\[Omega][t]^2-Cos[\[Theta][t]]]/. badSolution, {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->Red, PlotRange-> {{0,tMax},{0.5*\[Omega]0^2-1-0.5,0.5*\[Omega]0^2-1+0.5}}];

GraphicsRow[{Show[{badEnergyPlot,goodEnergyPlot},PlotLabel->"Celkov\[AAcute] energie syst\[EAcute]mu",AxesLabel->{t,"E"}],Show[{badPlot,goodPlot},PlotLabel->"F\[AAcute]zov\[YAcute] prostor",AxesLabel->{\[Theta],\[Omega]}]}]



(* ::InheritFromParent:: *)
(**)


h=0.5;
tMax = 100;
mtd ="ExplicitEuler";
\[Omega]0 = 1.7;
(*Pokus napravit chybu pomoc\[IAcute] projekce*)
goodSolution = NDSolve[{\[Omega]'[t]== -Sin[\[Theta][t]],\[Theta]'[t]==\[Omega][t],\[Omega][0]==\[Omega]0,\[Theta][0]==0 },{\[Theta][t], \[Omega][t]},{t,0,tMax}];
projSolution =NDSolve[{\[Omega]'[t]== -Sin[\[Theta][t]],\[Theta]'[t]==\[Omega][t],\[Omega][0]==\[Omega]0,\[Theta][0]==0 },{\[Theta][t], \[Omega][t]},{t,0,tMax},StartingStepSize->h,
						Method->{"Projection",Method->mtd, "Invariants"->{0.5*\[Omega][t]^2-Cos[\[Theta][t]]}}];
(*Porovn\[AAcute]n\[IAcute] v\[YAcute]voj\[URing] v \[CHacek]ase*)					
projPlot = ParametricPlot[Evaluate[{\[Theta][t], \[Omega][t]} /. projSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Darker[Green],Dashed},PlotLabel->"F\[AAcute]zov\[YAcute] prostor s projek\[CHacek]n\[IAcute] metodou",AxesLabel->{\[Theta],\[Omega]}];
goodPlot = ParametricPlot[Evaluate[{\[Theta][t], \[Omega][t]} /. goodSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Blue,Dashed}];
goodTimePlot= Plot[Evaluate[{\[Theta][t]} /. goodSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Blue},PlotLabel->"V\[YAcute]chylka kyvadla v \[CHacek]ase",AxesLabel->{t,TraditionalForm[HoldForm[\[Theta][t]]]}];
projTimePlot= Plot[Evaluate[{\[Theta][t]} /. projSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Darker[Green]}];

phasePlot = Show[{projPlot,goodPlot}];
thetaPlot =Show[{goodTimePlot,projTimePlot}];
GraphicsRow[{projPlot,thetaPlot}]


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


h=0.01;

SPRKSolution = NDSolve[{\[Omega]'[t]== -Sin[\[Theta][t]],\[Theta]'[t]==\[Omega][t],\[Omega][0]==\[Omega]0,\[Theta][0]==0 },{\[Theta][t], \[Omega][t]},{t,0,tMax}, 
					Method -> {"SymplecticPartitionedRungeKutta", "DifferenceOrder" -> 8,"PositionVariables" -> {\[Theta][t]}},  StartingStepSize -> h];
SPRKPlot = ParametricPlot[Evaluate[{\[Theta][t], \[Omega][t]} /. SPRKSolution], {t, 0 ,tMax}, BaseStyle->styleLaTeX, PlotStyle->{Brown,Dashed}]

SPRKTimePlot= Plot[Evaluate[{\[Theta][t]} /. projSolution], {t, 0 ,tMax}, PlotStyle->{Brown}, BaseStyle->styleLaTeX];
Show[{goodTimePlot,SPRKTimePlot}]
     


(* ::Section::Closed:: *)
(*5 Porovn\[AAcute]n\[IAcute] analytick\[EAcute] aproximace a numerick\[EAcute]ho \[RHacek]e\[SHacek]en\[IAcute]*)


\[Omega]nula=Sqrt[konstg/konstl];
\[Omega]PL = Sqrt[\[Omega]nula^2/(1+((konst\[Theta]init^2)/8))];
solPoincareLindstedt=konst\[Theta]init*Cos[\[Omega]PL*t]+konst\[Theta]init^3/192*(Cos[3*\[Omega]PL*t]-Cos[\[Omega]PL*t]);


Plot[{solPoincareLindstedt,Evaluate[\[Theta][t] /. NumSol], 10/(2*Pi),-10/(2*Pi)},{t,0,4\[Pi]},BaseStyle->styleLaTeX, 
PlotStyle->{Red,Blue, {Black, Dotted},{Black,Dotted}}, PlotLabel->"Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute] vs. Poincar\[EAcute]-Lindstedt", 
PlotLegends->{MaTeX /@ solPoincareLindstedt,  MaTeX["\\text{Numerick\[EAcute] \[RHacek].,default}"],MaTeX["\\theta_{max}"]   }]


(* ::Section:: *)
(*6.1 Dal\[SHacek]\[IAcute] \[CHacek]len v rozvoji*)


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
thetaCsol /. \[Beta] -> -5/512 *\[Theta]tI^4 *\[Omega]2 /. Sqrt[\[Omega]2] -> \[Omega] /. 1/Sqrt[\[Omega]2] -> 1/\[Omega] 


TPLctvrtyRad = 2*\[Pi]*Sqrt[l/g] * ( 1 + \[Theta]^2/16 + 3/1024*\[Theta]^4 ) // TraditionalForm
TPLctvrtyRad /. l -> 1 /. g->9.81 /. \[Theta]->10/(2\[Pi])


(* ::Section::Closed:: *)
(*6.2 Dodatek*)


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


\[Omega]PL4 = Sqrt[ g/l /(1+\[Theta]^2/16 + 5/512*\[Theta]^4 )]
Sqrt[g/(l (1+\[Theta]^2/16+(5 \[Theta]^4)/512))]
TPL4radDleMathematicy = Simplify[Series[2\[Pi]/\[Omega]PL4,{\[Theta],0,5}]]
Normal[TPL4radDleMathematicy] /. l -> 1 /. g->9.81 /. \[Theta]->10/(2\[Pi])
SeriesData[\[Theta], 0, {2 (g/l)^Rational[-1, 2] Pi, 0, Rational[1, 16] (g/l)^Rational[-1, 2] Pi, 0, Rational[9, 1024] (g/l)^Rational[-1, 2] Pi}, 0, 6, 1]
2.221425034572327`


(* ::Section::Closed:: *)
(*7 Shrnut\[IAcute] v\[YAcute]sledk\[URing]*)


MaTeX["\\text{1) Analytick\[AAcute] aproximace:}"]
MaTeX["T_{analytic}"]== TseriesAprox /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1 //TraditionalForm

MaTeX["\\text{2) P\[RHacek]es eliptick\[YAcute] integr\[AAcute]l:}"] 
MaTeX["T_{elliptic}"]==TIntegralExplicit /. \[Theta]init->10/(2\[Pi]) /. g->9.81 /. l->1 //TraditionalForm 

MaTeX["\\text{3) Numerick\[EAcute] \[RHacek]e\[SHacek]en\[IAcute]:}"]
MaTeX["T_{numeric}"]==TNumSol //TraditionalForm 

MaTeX["\\text{4) Aproximace Poincar\[EAcute]-Lindstedt:}"]
MaTeX["T_{PL}"]== 2\[Pi]/\[Omega]PL  //TraditionalForm 



