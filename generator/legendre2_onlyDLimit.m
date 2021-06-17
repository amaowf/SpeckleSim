LP[l_, m_, x_] := If[m >= 0, (-1)^m 2^l (1 - x^2)^(m/2) \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = m\), \(l\)]\((
\*FractionBox[\(k!\), \(\((k - m)\)!\)] 
\*SuperscriptBox[\(x\), \(k - m\)] Binomial[l, k] Binomial[
\*FractionBox[\(l + k - 1\), \(2\)], l])\)\), (l + m)!/(l - m)! 2^
    l (1 - x^2)^(-m/2) \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = \(-m\)\), \(l\)]\((
\*FractionBox[\(k!\), \(\((k + m)\)!\)] 
\*SuperscriptBox[\(x\), \(k + m\)] Binomial[l, k] Binomial[
\*FractionBox[\(l + k - 1\), \(2\)], l])\)\)];

pw = 255 ;(*pageWidth*)
maxl = 100; 
SetDirectory[NotebookDirectory[]];
fileD = "tmp/d_LPcos.f90";
fileDLimit = "tmp/d_LPcos_limit.f90"

(*generate D*)

dBuffer = 
  FullSimplify[D[LP[lReplace, mReplace, Cos[xReplace]], xReplace]];
fd[l_, m_, x_] := 
  dBuffer /. {lReplace -> l, mReplace -> m, xReplace -> x};

str = OpenWrite[fileDLimit, PageWidth -> pw]

WriteString[str, "select case(l)\n"]
For[l1 = 1, l1 <= maxl, l1++ ,
 WriteString[str, "\ncase(", l1, ") ********************* l=", l1, 
   "\n"]
  WriteString[str, "  select case(m)\n"]
  For[m1 = -l1, m1 <= l1, m1++ ,
   (*Print["l=",l1,"  m=",
   m1];*)
   (*Polynome*)
   (*generate D Limit*)
   
   erg = N[FullSimplify[Limit[fd[l1, m1, xVar], xVar -> 0]]];
   
   WriteString[str, "  case(", m1, ")\n"]
    WriteString[str, "    erg="]
    Write[str, FortranForm[erg]]
   ]
  WriteString[str, "  case default\n"]
  WriteString[str, "    call io_error(\"m out of range\")\n"]
  WriteString[str, "  endselect\n"]
 ]
WriteString[str, "case default\n"]
WriteString[str, "  call io_error(\"l out of range\")\n"]
WriteString[str, "endselect\n"]
Close[str];
(*FilePrint[%]*)