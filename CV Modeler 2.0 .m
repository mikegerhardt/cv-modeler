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


(* ::Section::Closed:: *)
(*AQDS reversible model, pH 0*)


(* ::Input:: *)
(*aqdsmodel = modelEr[0.004];*)
(*findpeakvoltages[aqds]*)


(* ::Input:: *)
(*aqdsmodelplot=ListPlot[Table[currentdensityconvert[voltadjust[x,-0.004],1000/0.07065],{x,{aqds,aqdsmodel}}],*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*aqdsmodeleq = modelEqEr[0.004, 7.2*10^-5];*)


(* ::Input:: *)
(*aqdsmodeleqplot = ListPlot[Table[currentdensityconvert[voltadjust[x,-0.004],1000/0.07065],{x,{aqds,aqdsmodeleq}}],*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*Export["AQDS Model Reversible.png",aqdsmodelplot];Export["AQDS Model Eq measured k0.png",aqdsmodeleqplot];*)


(* ::Section::Closed:: *)
(*AQDS models at pH 1 and 14*)


(* ::Input:: *)
(*aqdsph1=Import["aqds-ph1-50mVs.csv"];*)


(* ::Input:: *)
(*findpeakvoltages[aqdsph1]*)


(* ::Input:: *)
(*aqdsph1opshift = voltadjust[aqdsph1, 0.0155];*)


(* ::Input:: *)
(*findpeakvoltages[aqdsph14]*)


(* ::Input:: *)
(*aqdsmodelEqHiScanRate=currentdensityconvert[modelEqEr[0.,7.2*10^-5,scanrate->0.05,estart->0.3,erev->-0.3],1000/0.07065];*)


(* ::Input:: *)
(*aqdsph14opshift = voltadjust[aqdsph14,(0.464+0.513)/2];*)


(* ::Input:: *)
(*aqdsphplot=ListPlot[{currentdensityconvert[aqdsph1opshift,1000/0.07065], aqdsph14opshift,aqdsmodelEqHiScanRate},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*PlotRange->{{-0.3,0.3},{-0.4,0.3}},*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"pH 1", "pH 14","Model"}]]*)


(* ::Input:: *)
(*Export["AQDS_pH_comparison.png",aqdsphplot]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(0.017-0.048)/2*)


(* ::Input:: *)
(*ListPlot[{aqdsph1,modelEqEr[-0.0155,7.2*10^-5, scanrate->.050, erev->-0.7]},PlotRange->All]*)


(* ::Subchapter:: *)
(*Proving that the EqEq model does not work for AQDS*)


(* ::Input:: *)
(*aqdsopshift = voltadjust[aqds,-0.004];*)


(* ::Input:: *)
(*aqdsmodeleq = modelEqEr[0.004, 7.2*10^-5];*)


(* ::Input:: *)
(*aqdsmodelentropic = modelEqEq[0.004 + 0.0178, 7.2*10^-5, 0.004 -0.0178, 7.2*10^-5, electrons -> 1];*)


(* ::Input:: *)
(*aqdsmodeleqeq = modelEqEq[0.004, 7.2*10^-5, 0.004, 7.2*10^-5, electrons ->1];*)
(**)


(* ::Input:: *)
(*eqeqFarMid=modelEqEq[0.05,0.01,-0.05,0.01,electrons->1];*)


(* ::Input:: *)
(*eqeqFar = modelEqEq[0.1,0.01,-0.1,0.01, electrons ->1];*)
(*eqeqFarMid=modelEqEq[0.05,0.01,-0.05,0.01,electrons->1];*)
(*eqeqMid = modelEqEq[0.0178, 0.01, -0.0178, 0.01, electrons ->1];*)
(*eqeqClose = modelEqEq[0,0.01, 0, 0.01, electrons ->1];*)
(*eqeqAutoReduce = modelEqEq[0,0.01, 0.05,0.01, electrons ->1];*)


(* ::Input:: *)
(*eqeqcurves={eqeqFar[[3]],eqeqFarMid[[3]],eqeqMid[[3]],eqeqClose[[3]],eqeqAutoReduce[[3]]};*)


(* ::Input:: *)
(*ListPlot[eqeqcurves]*)


(* ::Input:: *)
(*eqeqplots=Table[*)
(*ListPlot[{currentdensityconvert[aqdsopshift,1000/0.07065],*)
(*currentdensityconvert[x,1000/0.07065]},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]*)
(*],{x,eqeqcurves}];*)


(* ::Input:: *)
(*?Vertical*)


(* ::Input:: *)
(*?Line*)


(* ::Input:: *)
(*eqeqplots[[2]]*)


(* ::Input:: *)
(*StringReplace["abcdef","d"->ToString[1]]*)


(* ::Input:: *)
(*Table[Export[StringReplace["eqeqN.png","N"->ToString[x]],eqeqplots[[x]]],{x,1,Length[eqeqplots]}]*)


(* ::Input:: *)
(*?Export*)


(* ::Input:: *)
(*ListPlot[{aqds,eqeqFarMid[[3]]}]*)


(* ::Text:: *)
(*As the second reduction potential shifts more positive, the intermediate is consumed immediately upon its reduction. This causes the apparent reduction potential to shift positive.*)
(*To prove this model doesn't work with AQDS:*)
(*start with the waves really far apart,*)
(*gradually bring them in,*)
(*shift the second wave past the first.*)


(* ::Input:: *)
(*ListPlot[dhaq1mm]*)


(* ::Input:: *)
(**)


(* ::Subchapter:: *)
(*AQS model fitting & figures*)


(* ::Input:: *)
(*findpeakvoltages[aqs]*)


(* ::Input:: *)
(*aqsopshift = currentdensityconvert[voltadjust[aqs,0.015],1000/0.07065];*)


(* ::Input:: *)
(*aqsmodelEq =currentdensityconvert[modelEqEr[0.,2.*10^-5],1000/0.07065];*)


(* ::Input:: *)
(*aqseqplot=ListPlot[{aqsopshift,aqsmodelEq},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*aqserplot=ListPlot[{aqsopshift,currentdensityconvert[aqdsmodel,1000/0.07065]},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*Export["AQS Er model comparison.png",aqserplot]*)


(* ::Input:: *)
(*Export["AQS Eq model comparison.png",aqseqplot]*)


(* ::Subchapter:: *)
(*DHAQDS model fitting*)


(* ::Input:: *)
(*dhaqds2=Import["dhaqds_V_vs_SHE_25mVs_from_MRS_talk.csv"];*)


(* ::Input:: *)
(*ListPlot[{dhaqds,dhaqds2}]*)


(* ::Input:: *)
(*findpeakvoltages[dhaqds]*)


(* ::Input:: *)
(*findpeakvoltages[dhaqds2]*)


(* ::Input:: *)
(*(.143+.1)/2*)


(* ::Input:: *)
(*dhaqdsopshift = voltadjust[dhaqds, -0.118];findpeakvoltages[dhaqdsopshift]*)


(* ::Input:: *)
(*dhaqdsopshift2=voltadjust[dhaqds2,-0.121];*)


(* ::Input:: *)
(*dhaqdsmodelEr = currentdensityconvert[modelEr[0.],1000/0.07065];*)


(* ::Input:: *)
(*dhaqdsErplot=ListPlot[{dhaqdsopshift2,dhaqdsmodelEr},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Arial",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*dhaqdsEqplot=ListPlot[{dhaqdsopshift2,dhaqdsmodelEq},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Arial",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model"}]]*)


(* ::Input:: *)
(*Export["DHAQDS Eq model comparison.png",dhaqdsEqplot]*)


(* ::Input:: *)
(*Export["DHAQDS Er model.png",dhaqdsErplot]*)


(* ::Input:: *)
(*Export["DHAQDS Eq model higher conc.png",dhaqdsEqplothigherconc]*)


(* ::Input:: *)
(*Export["DHAQDS EQ higher conc lower k.png",dhaqdsEqplotlowerk]*)


(* ::Input:: *)
(*dhaqdsEqplothigherconc = ListPlot[{dhaqdsopshift2,dhaqdsmodelEq,dhaqdsmodelEqhigherconc},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Arial",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model", "Model, [DHAQDS] = 1.2 mM"}]]*)


(* ::Input:: *)
(*dhaqdsEqplotlowerk = ListPlot[{dhaqdsopshift2,dhaqdsmodelEq,dhaqdsmodelEqlowerkhigherconc},*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Overpotential (V)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Arial",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"Experiment","Model", "Model, [DHAQDS] = 1.2 mM, \!\(\*SubscriptBox[\(k\), \(0\)]\) = 3*\!\(\*SuperscriptBox[\(10\), \(-3\)]\) cm/s"}]]*)


(* ::Input:: *)
(*dhaqdsmodelEq=currentdensityconvert[modelEqEr[0.,1.2*10^-4],1000/0.07065];*)


(* ::Input:: *)
(*dhaqdsmodelEqhigherconc=currentdensityconvert[modelEqEr[0.,1.2*10^-4,cbulk->1.2],1000/0.07065];*)


(* ::Input:: *)
(*dhaqdsmodelEqlowerk=currentdensityconvert[modelEqEr[0.,3*10^-5],1000/0.07065];*)


(* ::Input:: *)
(*dhaqdsmodelEqlowerkhigherconc=currentdensityconvert[modelEqEr[0.,3*10^-5,cbulk->1.2],1000/0.07065];*)


(* ::Input:: *)
(*dhaqdsmodelEqEq=modelEqEq[0.,1.2*10^-4,0.,1.2*10^-4,estart->0.2,electrons->1];*)


(* ::Input:: *)
(*dhaqdsmodelEqEqcv=currentdensityconvert[dhaqdsmodelEqEq[[3]],1000/0.07068];*)


(* ::Input:: *)
(*ListPlot[{dhaqdsopshift2,dhaqdsmodelEqlowerkhigherconc}]*)


(* ::Input:: *)
(*|*)


(* ::Subchapter:: *)
(*2,6-DHAQ in base model fitting*)


(* ::Input:: *)
(*ListPlot[dhaq1mm]*)


(* ::Input:: *)
(*dhaq50mVs = Import["dhaq50.csv"];*)


(* ::Input:: *)
(*dhaqvsaqdsplot=ListPlot[Table[voltadjust[x,0.213],{x,{dhaq50mVs,aqdsph14}}],*)
(*Frame->True,*)
(*AxesStyle->Dashed,*)
(*FrameLabel->{"Potential (V vs SHE)","Current (mA/\!\(\*SuperscriptBox[\(cm\), \(2\)]\))"},*)
(*BaseStyle->{FontSize->18,FontFamily->"Calibri",FontWeight->"Plain"},*)
(*ImageSize->640,*)
(*(*PlotStyle\[Rule]{Black,Red,Red,Red,Red,Red},*)*)
(*PlotLegends->SwatchLegend[{"2,6-DHAQ","AQDS"}]]*)


(* ::Input:: *)
(*Export["DHAQ vs AQDS pH 14.png", dhaqvsaqdsplot]*)


(* ::Input:: *)
(*findpeakvoltages[dhaq1mm]*)


(* ::Input:: *)
(*(0.855+0.940)/2-.213*)


(* ::Input:: *)
(*.657+.213*)


(* ::Input:: *)
(*(.657+.717)/2*)


(* ::Input:: *)
(*dhaqmodelEqEq=modelEqEq[0.03,7.*10^-5,-0.03,7.*10^-5,electrons->1,dreactant->4.8*10^-10,dproduct->4.8*10^-10,di->4.8*10^-10,erev->-0.25];*)


(* ::Input:: *)
(*dhaqcv=currentdensityconvert[voltadjust[dhaq1mm,0.687+0.213],1000/0.07065];*)


(* ::Input:: *)
(*dhaqmodelcv=currentdensityconvert[dhaqmodelEqEq[[3]],1000/0.07065];*)


(* ::Input:: *)
(*ListPlot[{dhaqcv,dhaqmodelcv}]*)


(* ::Input:: *)
(*ListPlot[voltadjust[dhaq1mm,0.687+0.213]]*)


(* ::Subchapter:: *)
(*Data export for ACS talk*)


(* ::Input:: *)
(*aqdsexport=Table[currentdensityconvert[voltadjust[x,-0.004],1000/0.07065],{x,{aqds,aqdsmodel}}];*)


(* ::Input:: *)
(*Export["AQDS CV pH 0.txt", aqdsexport[[1]],"Table"]*)


(* ::Input:: *)
(*Export["Reversible model overpotential basis.txt", aqdsexport[[2]],"Table"]*)


(* ::Input:: *)
(*aqdsinbaseexport = {currentdensityconvert[aqdsph1opshift,1000/0.07065], aqdsph14opshift,aqdsmodelEqHiScanRate};*)


(* ::Input:: *)
(*Export["AQDS CV pH 1 50 mVs.txt", aqdsinbaseexport[[1]],"Table"];*)
(*Export["AQDS CV pH 14 50 mVs.txt", aqdsinbaseexport[[2]],"Table"];*)
(*Export["Reversible model 50 mVs.txt", aqdsinbaseexport[[3]],"Table"];*)


(* ::Input:: *)
(*aqsexport = {aqsopshift,aqsmodelEq};*)


(* ::Input:: *)
(*Export["AQS cv.txt", aqsexport[[1]],"Table"];*)
(*Export["AQS model Eq.txt",aqsexport[[2]],"Table"];*)


(* ::Input:: *)
(*dhaqdsexport={dhaqdsopshift2,dhaqdsmodelEq,dhaqdsmodelEr,dhaqdsmodelEqhigherconc,dhaqdsmodelEqlowerkhigherconc};*)


(* ::Input:: *)
(*Export["DHAQDS cv.txt",dhaqdsexport[[1]],"Table"];*)
(*Export["DHAQDS Eq model.txt",dhaqdsexport[[2]],"Table"];*)
(*Export["DHAQDS Er model.txt", dhaqdsexport[[3]],"Table"];*)
(*Export["DHAQDS Eq model higher conc.txt", dhaqdsexport[[4]],"Table"];*)
(*Export["DHAQDS Eq model higher conc lower k.txt",dhaqdsexport[[5]],"Table"];*)


(* ::Input:: *)
(*dhaqexport=Table[voltadjust[x,0.213],{x,{dhaq50mVs,aqdsph14}}];*)


(* ::Input:: *)
(*Export["DHAQ cv 50 mVs vs SHE.txt",dhaqexport[[1]],"Table"];*)
(*Export["AQDS cv 50 mVs vs SHE ph 14.txt", dhaqexport[[2]],"Table"];*)


(* ::Input:: *)
(*dhaqmodeleqer = Table[currentdensityconvert[modelEqEr[0., 1.*10^(-x)],1000/0.07065],{x,{4,5,6}}];*)


(* ::Input:: *)
(*dhaqexport2 = {dhaqcv,dhaqmodelcv,dhaqmodelEqEq[[1]],dhaqmodelEqEq[[2]],dhaqmodeleqer[[1]],dhaqmodeleqer[[2]],dhaqmodeleqer[[3]]};*)


(* ::Input:: *)
(*Export["DHAQ cv 25 mVs overpotential.txt", dhaqexport2[[1]],"Table"];*)
(*Export["DHAQ cv model total.txt", dhaqexport2[[2]],"Table"];*)
(*Export["DHAQ cv model i1.txt", dhaqexport2[[3]],"Table"];*)
(*Export["DHAQ cv model i2.txt", dhaqexport2[[4]],"Table"];*)
(*Export["DHAQ cv model eqer k 10e-2.txt", dhaqexport2[[5]],"Table"];*)
(*Export["DHAQ cv model eqer k 10e-3.txt", dhaqexport2[[6]],"Table"];*)
(*Export["DHAQ cv model eqer k 10e-4.txt", dhaqexport2[[7]],"Table"];*)
