(* ::Package:: *)

(* ::Title:: *)
(*Cyclic Voltammetry Modeling, version 2*)


(* ::Text:: *)
(*by Mike Gerhardt*)
(**)
(*Changes since v1:*)
(*	All model functions now output a list of points in the form {voltage, current}*)
(*	Many variables, like "electrons" and "scanrate" are now defined as options with default values. For example, modelEr[0] assumes a two-electron reversible CV with a reduction potential of 0 V. modelEr[0, electrons -> 1] does a one-electron reversible CV.*)


(* ::Chapter:: *)
(*Under the Hood - function definitions etc*)


(* ::Subchapter:: *)
(*Default data import*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*aqds=Import["aqds.csv"];*)
(*dhaq1mm=Import["dhaq1mM.csv"];*)
(*aqs = Import["aqs.csv"];*)
(*aqdsph14 = Import["aqds_pH_14_50mVs.csv"];*)
(*(*The data in the .csv files in this notebook's directory are in the form (voltage in volts vs Ag/AgCl, current in A) so there are some helper functions to massage it as needed*)*)
(*dhaqds=Import["dhaqds_V_vs_SHE_25mVs_from_Nature_paper.csv"];*)
(*(*Note that the DHAQDS data is in volts vs SHE, current in mA/cm2. *)*)


(* ::Subchapter:: *)
(*voltage[], weights[], findpeakvoltages[], voltadjust[], currentdensityconvert[]*)


(* ::Input:: *)
(*voltage[t_,erev_,scanrate_,start_,direction_]:=erev-direction Abs[scanrate*t-direction*erev+direction*start];(*Equation 2:2 in Oldham. direction is 1 for a positive sweep and -1 for a negative sweep*)*)
(**)
(*weights[steps_]:=RecurrenceTable[{w[n]==(2n-1)w[n-1]/(2n),w[0]==1.},w,{n,1,steps}];(*weightlist=weights[1200]*);*)
(**)
(*(*These are some helper functions for later analysis. They are copied verbatim from V1 and may require some fine-tuning.*)*)
(**)
(*findpeakvoltages[data_]:=*)
(*With[{currentdata = Thread[data][[2]]},*)
(*With[{maxpos = Ordering[currentdata,-1], minpos = Ordering[currentdata, 1]},*)
(*{data[[maxpos]], data[[minpos]]}]];*)
(**)
(*voltadjust[data_,adjustment_]:= With[{voltagedata = Thread[data][[1]], currentdata = Thread[data][[2]]},Thread[{voltagedata + adjustment, currentdata}]];*)
(**)
(*currentdensityconvert[data_,mfactor_]:= With[{voltagedata = Thread[data][[1]], currentdata = Thread[data][[2]]},Thread[{voltagedata, currentdata*mfactor}]];*)


(* ::Subchapter:: *)
(*xi[] - massaging voltage into the proper form*)


(* ::Input:: *)
(*Options[xi] = {electrons -> 2,  estart -> 0.2,erev ->-0.2,f->96485,r->8.314,temp->298,scanrate->0.025};*)
(*xi[t_,ezero_,OptionsPattern[]]:=With[{direction =(OptionValue[erev] - OptionValue[estart])/Abs[(OptionValue[erev] - OptionValue[estart])], electrons = OptionValue[electrons], f = OptionValue[f], r=OptionValue[r], temp=OptionValue[temp], erev=OptionValue[erev], scanrate=OptionValue[scanrate],estart = OptionValue[estart]},Exp[direction*electrons*f/(r*temp)*(voltage[t,erev,scanrate,estart,direction]-ezero)]]*)


(* ::Input:: *)
(*Options[makevoltagetable]={estart ->0.2, erev ->-0.2, temp->298,f->96485,r->8.314,scanrate->0.025,steps->1200};*)
(*makevoltagetable[OptionsPattern[]]:=With[{erev=OptionValue[erev],estart=OptionValue[estart],temp=OptionValue[temp],f=OptionValue[f],r=OptionValue[r],scanrate=OptionValue[scanrate],direction=(OptionValue[erev] - OptionValue[estart])/Abs[(OptionValue[erev] - OptionValue[estart])],*)
(*delta=2*Abs[OptionValue[erev]-OptionValue[estart]]/(OptionValue[scanrate]*OptionValue[steps]),*)
(*steps= OptionValue[steps]},*)
(*Table[voltage[t*delta,erev,scanrate,estart,direction],{t,1,steps}]]*)


(* ::Subchapter:: *)
(*modelEr[] - a reversible two-electron wave*)


(* ::Input:: *)
(*Options[modelEr]={electrons->2,estart->0.2, erev->-0.2,f->96485,r->8.314,temp->298,scanrate->0.025,steps->1200,area->0.07065*10^-4,cbulk -> 1,dreactant->3.8*10^-10,dproduct->3.8*10^-10};*)
(**)
(*modelEr[ezero_,opts:OptionsPattern[]]:=*)
(*Module[{currentvalue,currents,weightlist=weights[OptionValue[steps]],voltagetable=makevoltagetable[FilterRules[{opts},Options[makevoltagetable]]],direction=(OptionValue[erev] - OptionValue[estart])/Abs[(OptionValue[erev] - OptionValue[estart])],f=OptionValue[f],area=OptionValue[area],cb=OptionValue[cbulk],dr=OptionValue[dreactant],dp=OptionValue[dproduct], delta=2*Abs[OptionValue[erev]-OptionValue[estart]]/(OptionValue[scanrate]*OptionValue[steps]),steps=OptionValue[steps],electrons=OptionValue[electrons]},*)
(**)
(*For[i=2;currents={direction*electrons*f*area*cb*Sqrt[dr/delta]/(1+Sqrt[dr/dp]/xi[delta, ezero, FilterRules[{opts},Options[xi]]])};,i<=steps, i++,*)
(*currentvalue = direction*electrons*f*area*cb*Sqrt[dr/delta]/(1+Sqrt[dr/dp]/xi[i*delta, ezero, FilterRules[{opts},Options[xi]]])-Sum[weightlist[[j]]*currents[[i-j]],{j,1,i-1}];*)
(*AppendTo[currents,currentvalue];]; *)
(**)
(*Thread[{voltagetable,currents}]]*)


(* ::Subchapter:: *)
(*modelEqEr[] - assuming the first electron transfer is sluggish and the second is "free"*)


(* ::Input:: *)
(*Options[modelEqEr]={electrons->2,estart->0.2, erev->-0.2,f->96485,r->8.314,temp->298,scanrate->0.025,steps->1200,area->0.07065*10^-4,cbulk -> 1,dreactant->3.8*10^-10,dproduct->3.8*10^-10,alpha->0.5};*)
(**)
(*modelEqEr[ezero_,k_,opts:OptionsPattern[]]:=*)
(*Module[{currentvalue,currents,weightlist=weights[OptionValue[steps]],voltagetable=makevoltagetable[FilterRules[{opts},Options[makevoltagetable]]],direction=(OptionValue[erev] - OptionValue[estart])/Abs[(OptionValue[erev] - OptionValue[estart])],f=OptionValue[f],area=OptionValue[area],cb=OptionValue[cbulk],dr=OptionValue[dreactant],dp=OptionValue[dproduct], delta=2*Abs[OptionValue[erev]-OptionValue[estart]]/(OptionValue[scanrate]*OptionValue[steps]),steps=OptionValue[steps],*)
(*electrons = OptionValue[electrons],*)
(*alpha = OptionValue[alpha]},*)
(**)
(*For[i=1;currents={};,i<=steps, i++,*)
(*currentvalue = (direction*electrons*area*f*cb*Sqrt[dr/delta]-(1+Sqrt[dr/dp]/xi[delta*i, ezero, FilterRules[{opts},Options[xi]]])Sum[weightlist[[j]]*currents[[i-j]],{j,1,i-1}])/(1+(Sqrt[dr/dp]/xi[delta*i, ezero, FilterRules[{opts},Options[xi]]])+Sqrt[dr/delta]/(k*xi[delta*i, ezero, FilterRules[{opts},Options[xi]]]^alpha));*)
(*AppendTo[currents,currentvalue];]; *)
(**)
(*Thread[{voltagetable,currents}]]*)


(* ::Subchapter:: *)
(*modelEqEq[] - two one-electron waves separated by some potential, each having its own rate constant*)


(* ::Input:: *)
(*(* Diffusivities in m^2/s*)
(*  Area in m^2*)
(*  Subscript[k, 0] in m/s*)
(*  You must specify electrons\[Rule]1 or it doesn't get passed to xi*)
(*	Also make sure you correctly specify e1 and e2. If you specify them out of order you get funky results!!*)
(**)*)


(* ::Input:: *)
(*Options[modelEqEq]={electrons->1,estart->0.2, erev->-0.2,f->96485,r->8.314,temp->298,scanrate->0.025,steps->1200,area->0.07065*10^-4,cbulk -> 1,dreactant->3.8*10^-10,dproduct->3.8*10^-10,di->3.8*10^-10,alpha1->0.5,alpha2 -> 0.5};*)
(*modelEqEq[e1_,k1_,e2_,k2_,opts:OptionsPattern[]]:=*)
(*Module[{z,y,w,v,i1,i2,i1s,i2s,itotal,wave1,wave2,fullwave,*)
(*weightlist=weights[OptionValue[steps]],voltagetable=makevoltagetable[FilterRules[{opts},Options[makevoltagetable]]],direction=(OptionValue[erev] - OptionValue[estart])/Abs[(OptionValue[erev] - OptionValue[estart])],f=OptionValue[f],area=OptionValue[area],cb=OptionValue[cbulk],*)
(*dr=OptionValue[dreactant],*)
(*dp=OptionValue[dproduct],*)
(*di=OptionValue[di],*)
(*delta=2*Abs[OptionValue[erev]-OptionValue[estart]]/(OptionValue[scanrate]*OptionValue[steps]),steps=OptionValue[steps],*)
(*electrons = OptionValue[electrons],*)
(*alpha1 = OptionValue[alpha1],*)
(*alpha2 = OptionValue[alpha2],*)
(*xiopts = FilterRules[{opts},Options[xi]]},*)
(**)
(*z[increment_]:=(xi[delta*increment,e1,xiopts]^(1-alpha1)*Sqrt[di]/(k1*Sqrt[delta]));*)
(*y[increment_]:= 1 + xi[delta*increment,e1,xiopts]*Sqrt[di/dr];*)
(*w[increment_]:=Sqrt[di]/(k2*Sqrt[delta]*xi[delta*increment,e2,xiopts]^alpha2);*)
(*v[increment_]:=1+1/xi[delta*increment,e2,xiopts]*Sqrt[di/dp];*)
(*i1[increment_]:=(direction(w[increment]+v[increment])*xi[delta*increment,e1,xiopts]f area cb Sqrt[di/delta]-(y[increment](w[increment]+v[increment])-1)*Sum[weightlist[[j]]*i1s[[increment-j]],{j,1,increment-1}]+w[increment]*Sum[weightlist[[j]]*i2s[[increment-j]],{j,1,increment-1}])/((z[increment]+y[increment])(w[increment]+v[increment])-1);*)
(*i2[increment_]:=(direction xi[delta*increment,e1,xiopts]f area cb Sqrt[di/delta]+z[increment]Sum[weightlist[[j]]*i1s[[increment-j]],{j,1,increment-1}]-(v[increment](z[increment]+y[increment])-1)*Sum[weightlist[[j]]*i2s[[increment-j]],{j,1,increment-1}])/((z[increment]+y[increment])(w[increment]+v[increment])-1);*)
(*For[i=1;i1s={};i2s={};,i<=steps, i++,*)
(*AppendTo[i1s,i1[i]];*)
(*AppendTo[i2s,i2[i]];];*)
(*itotal=i1s+i2s;*)
(*wave1=Thread[{voltagetable,i1s}];*)
(*wave2=Thread[{voltagetable,i2s}];*)
(*fullwave=Thread[{voltagetable,itotal}];*)
(*{wave1,wave2,fullwave}*)
(*]*)


(* ::Input:: *)
(*test = modelEqEq[-.87,7.*10.^-5,-.93,7.*10.^-5,{estart ->-0.7,erev->-1.15,dreactant->4.8*10^-10, dproduct->4.8*10^-10,di->4.8*10^-10,electrons->1}];*)


(* ::Input:: *)
(*ListPlot[test]*)


(* ::Input:: *)
(**)


(* ::Chapter:: *)
(*Computation work, figures, etc*)


(* ::Subchapter:: *)
(*Proving a statement in Bard and Faulkner about the appearance of 1e- separate waves separated by ~36 mV.*)


(* ::Input:: *)
(*reversiblemodel = modelEr[0.];*)
(*rev1e = modelEr[0., electrons -> 1];*)


(* ::Input:: *)
(*rev1etimes2 = currentdensityconvert[rev1e,2];*)


(* ::Input:: *)
(*35.6/2*)


(* ::Input:: *)
(*8.314*298/(2*96485)*Log[4.]*)


(* ::Input:: *)
(* *)


(* ::Input:: *)
(*entropicmodel = modelEqEq[0.0178, 0.01,-0.0178, 0.01,electrons->1];*)


(* ::Input:: *)
(*ListPlot[entropicmodel]*)


(* ::Input:: *)
(*ListPlot[{entropicmodel[[3]], reversiblemodel,rev1e ,rev1etimes2}]*)


(* ::Text:: *)
(*The above figure shows four CV models: one of a reversible, two electron process (orange), one of a reversible, one-electron process (green), the green curve stretched along the y-axis by a factor of 2 (red), and two successive one-electron transfers spaced RT/F ln 4 volts apart. This magic factor causes the curve to completely align with the red curve so you don't see it. See Bard and Faulkner ch 6, or Saveant's paper on convolution potential sweep voltammetry, part II.*)


(* ::Subchapter:: *)
(*AQDS model fitting*)
