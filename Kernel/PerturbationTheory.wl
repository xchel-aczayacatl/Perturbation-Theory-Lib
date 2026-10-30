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

MultiscaleExpansion2[diffEqns_List, depVars_List, indepVar_, param_, 
  expansionOrder_Integer, numberScales_Integer] := 
 Module[{s, scalesList, finalScales, timeDerivOp, ySeriesList, rules, 
   residuals, expr, 
   result},(*1.\
 Símbolos auxiliares independientes para derivar sin efectos colateral\
es*)scalesList = Table[Unique["sc"], {numberScales}];
  finalScales = 
   Prepend[Table[Subscript[indepVar, k], {k, 1, numberScales - 1}], 
    indepVar];
  (*2. Operador:d/dt->D0+eps*D1+eps^2*D2+...*)
  timeDerivOp[e_, 1] := 
   Sum[param^k*D[e, scalesList[[k + 1]]], {k, 0, numberScales - 1}];
  timeDerivOp[e_, n_Integer /; n > 1] := 
   Nest[timeDerivOp[#, 1] &, e, n];
  (*3. Expansión asintótica de cada variable dependiente*)
  ySeriesList = 
   Table[Sum[
     param^k*Subscript[depVars[[i]], k] @@ scalesList, {k, 0, 
      expansionOrder}], {i, Length[depVars]}];
  (*4. Reglas para derivadas de cualquier orden y para la variable*)
  rules = Flatten@
    Table[{Derivative[n_Integer][depVars[[i]]][indepVar] :> 
       timeDerivOp[ySeriesList[[i]], n], 
      depVars[[i]][indepVar] :> ySeriesList[[i]]}, {i, Length[depVars]}];
  (*5. Residuo de cada ecuación*)
  residuals = (#[[1]] - #[[2]]) & /@ diffEqns;
  (*6. Coeficientes por potencia de param,
  mismo formato que el original*)
  result = 
   Flatten[Table[
     expr = (residuals[[i]] /. rules) + O[param]^(expansionOrder + 1);
     Table[{param^j, SeriesCoefficient[expr, j] == 0}, {j, 0, 
       expansionOrder}], {i, Length[residuals]}], 1];
  (*7. Ahora sí,nombres finales de las escalas*)
  result /. Thread[scalesList -> finalScales]]

LindstedtPoincareExpansion[diffEqn_, depVar_, indepVar_, param_, 
  expansionOrder_Integer, frec_, nIndVar_] :=
 Module[{ySeries, omegaSeries, rescaledEq, seriesPoly, equationsByOrder},

  (* 1. Expansiones en potencias del parametro de perturbacion *)
  ySeries = Sum[param^i*Subscript[depVar, i][nIndVar], {i, 0, expansionOrder}];
  omegaSeries = 1 + Sum[param^i*Subscript[frec, i], {i, 1, expansionOrder}];

  (* 2. Regla de la cadena generalizada: d/dt -> omega * d/dtau *)
  rescaledEq = diffEqn /. {
     Derivative[k_][depVar][indepVar] :> (omegaSeries^k)*D[ySeries, {nIndVar, k}],
     depVar[indepVar] :> ySeries,
     indepVar -> nIndVar/omegaSeries
  };

  (* 3. Reduccion de la ecuacion a la forma residual f(eps) == 0 *)
  rescaledEq = rescaledEq[[1]] - rescaledEq[[2]];

  (* 4. Truncamiento asintotico al orden deseado *)
  seriesPoly = Normal[Series[rescaledEq, {param, 0, expansionOrder}]];

  (* 5. Coeficientes por cada orden del parametro *)
  equationsByOrder = CoefficientList[seriesPoly, param];

  (* 6. Jerarquia de ecuaciones diferenciales *)
  Thread[equationsByOrder == 0]
 ]
	
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
