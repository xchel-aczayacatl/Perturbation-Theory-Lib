BeginPackage["PerturbationTheory`"]

(* Declaration of the package public symbols. *)

MultiscaleExpansion

LindstedtPoincareExpansion

NayfehTrigToExp

AragonTrigToExp

BruteForceTrigToExp

FredholmAlternative

Begin["`Private`"]

(* Definition of the package symbols. *)

MultiscaleExpansion[diffEqns_,depVars_,indepVar_, param_,expansionOrder_,numberScales_] := 
Module[{cfsParam={},expansionVars={},expansionEqns={},scales,expansion,qi,i,j,m,n},
	If[Dimensions[Dimensions[depVars]][[1]]===1&&Dimensions[Dimensions[diffEqns]][[1]]===1,
		scales=Table[indepVar*param^n,{n,0,numberScales-1}];
		For[i=1,i<=Dimensions[depVars][[1]],i++,
			expansion=Table[Subscript[depVars[[i]], n]@@scales,{n,0,expansionOrder}] . Table[param^m,{m,0,expansionOrder}];
			AppendTo[expansionVars,expansion];
			Clear[expansion]];
		expansionEqns=diffEqns;
		For[i=1,i<=Dimensions[diffEqns][[1]],i++,
			For[j=1,j<=Dimensions[expansionVars][[1]],j++,
				expansionEqns[[i]]=Evaluate[expansionEqns[[i]]/.depVars[[j]]->Evaluate[expansionVars[[j]]/.indepVar-># &]]]];
		For[i=1,i<=Dimensions[expansionEqns][[1]],i++,
			For[j=(numberScales-1),j>0,j--,
				expansionEqns[[i]]=Evaluate[expansionEqns[[i]]/.param^j*indepVar->Subscript[indepVar, j]]];
			expansionEqns[[i]]=expansionEqns[[i,1]]+O[param]^(expansionOrder+1)==0];
		For[i=1,i<=Dimensions[expansionEqns][[1]],i++,
			For[j=0,j<=expansionOrder,j++,
				AppendTo[cfsParam,{param^j,SeriesCoefficient[expansionEqns[[i,1]],j]==0}]]]];
	Clear[expansionVars,expansionEqns,scales,expansion,i,j,m,n];
	Return[cfsParam]]

LindstedtPoincareExpansion[diffEqn_,depVar_,indepVar_,param_,expansionOrder_,frec_,nIndVar_] :=
Module[{shockwave={},ecExp,expF,expFrec, vectorSigma ,n},
	vectorSigma=Table[param^n,{n,0,expansionOrder}];
	expF=Table[Subscript[depVar, n][frec*indepVar],{n,0,expansionOrder}];
	expFrec=Table[Subscript[frec, n],{n,1,expansionOrder}];
	expFrec=Prepend[expFrec,1];
	expFrec=vectorSigma . expFrec;
	expF=vectorSigma . expF;ecExp=diffEqn/.depVar->(expF/.indepVar->#&);
	ecExp=(ecExp/.(frec indepVar)->nIndVar)/.frec->expFrec;ecExp=ExpandAll[ecExp[[1]]]+O[param]^(expansionOrder+1)==ecExp[[2]];
	shockwave=CoefficientList[ecExp[[1]],param];
	For[n=1,n<=Dimensions[shockwave][[1]],n++,
		shockwave[[n]]=shockwave[[n]]==0];
	Return[shockwave]]
	
NayfehTrigToExp[expr_,subscr1_,subscr2_,var1_,var2_] :=
Module[{nayfehExp},
	nayfehExp=expr/.{Cos->(E^(I #1)&),Sin->(E^(-I#1)&),C[subscr1]->(var1[#1]E^(I*var2[#1])&),C[subscr2]->(var1[#1]E^(-I*var2[#1])&)};
	Return[nayfehExp]]

AragonTrigToExp[expr_,subscr1_,subscr2_,var_] :=
Module[{aragonExp},
	aragonExp=expr/.{Cos->(E^(I #1)&),Sin->(E^(-I#1)&),C[subscr1]->(var[#1]/2&),C[subscr2]->(Conjugate[var[#1]/2]&)};
	Return[aragonExp]]

BruteForceTrigToExp[expr_] :=
Module[{fexp},
	fexp=expr/.{Cos->(E^(I #)&),Sin->(E^(-I#)&)};
	Return[fexp]]

FredholmAlternative[adjointOperatorSolution_,function_,integrationVar_,lowerLimit_,upperLimit_] :=
Module[{fredholmAlt},
	fredholmAlt=Integrate[adjointOperatorSolution*function,{integrationVar,lowerLimit,upperLimit}];
	Return[fredholmAlt]]

End[] (* End `Private` *)

EndPackage[]
