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
file1 = "tmp/poly.f90";
file2 = "tmp/limit.f90";
file3 = "tmp/d.f90";

(*generate poly*)
str = OpenWrite[file1, PageWidth -> pw]
WriteString[str, "select case(l)\n"]
For[l1 = 1, l1 <= maxl, l1++ ,
 WriteString[str, "\ncase(", l1, ") ********************* l=", l1, 
   "\n"]
  WriteString[str, "  select case(m)\n"]
  For[m1 = -l1, m1 <= l1, m1++ ,
   (*Print["l=",l1,"  m=",m1];*)
   (*Polynome*)
   
   erg = N[FullSimplify[m1 LP[l1, m1, Cos[x]]/Sin[x]]];
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
Close[str]

(*generate limit*)
str = OpenWrite[file2, PageWidth -> pw]
WriteString[str, "select case(l)\n"]
For[l1 = 1, l1 <= maxl, l1++ ,
 WriteString[str, "\ncase(", l1, ") ********************* l=", l1, 
   "\n"]
  WriteString[str, "  select case(m)\n"]
  For[m1 = -l1, m1 <= l1, m1++ ,
   (*Print["l=",l1,"  m=",m1];*)
   (*Polynome*)
   
   erg = N[Limit[FullSimplify[m1 LP[l1, m1, Cos[x]]/Sin[x]], 
      x -> 0]];
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
Close[str]
(*FilePrint[%]*)

(*generate D*)
str = OpenWrite[file3, PageWidth -> pw]
WriteString[str, "select case(l)\n"]
For[l1 = 1, l1 <= maxl, l1++ ,
 WriteString[str, "\ncase(", l1, ") ********************* l=", l1, 
   "\n"]
  WriteString[str, "  select case(m)\n"]
  For[m1 = -l1, m1 <= l1, m1++ ,
   (*Print["l=",l1,"  m=",m1];*)
   (*Polynome*)
   
   erg = N[FullSimplify[D[m1 LP[l1, m1, Cos[x]]/Sin[x], x]]];
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
Close[str]
(*FilePrint[%]*)