(* ::Package:: *)

(* ::Section:: *)
(*Kvadratick\[YAcute] \[CHacek]len rozvoje*)


(* ::Input:: *)
(*parametry = {g , l , \[Theta]init };*)


(* ::Input:: *)
(*(* Explicitni rozvoj periody *)*)
(*period = 2*\[Pi]*Sqrt[l/g]*(1+a*\[Theta]init^2+b*\[Theta]init^4+O[\[Theta]init]^6)*)


(* ::Input:: *)
(*(* Explicitne spocteny elipticky integral *)*)
(*TIntegralExplicit = 4*Sqrt[l/g]Integrate[1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2], {\[CurlyPhi], 0 , \[Pi]/2}]*)


(* ::Input:: *)
(*(* Tayloruv rozvoj integrandu *)*)
(*argument = 1/Sqrt[1-k^2*Sin[\[CurlyPhi]]^2];*)
(*k = Sin[\[Theta]init/2];*)
(*seriesIntegrand = Series[argument, {\[Theta]init, 0,2}]*)


(* ::Input:: *)
(*(* Perioda (integr\[AAcute]l) p\[RHacek]es Taylor\[URing]v rozvoj -- vychazi stejne jako ve skriptech! *)*)


(* ::Input:: *)
(*TseriesAprox = 4*Sqrt[l/g]*Integrate[seriesIntegrand, {\[CurlyPhi], 0, \[Pi]/2}]*)


(* ::Input:: *)
(**)
