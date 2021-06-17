#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 17:52:55 2019

@author: itodaiber
"""

import re

file_input = "generated/PiPoly.f90"
file_output = "compilable/PiPoly.f90"

new_data = []

# add header
new_data = ["module PiPoly\n",
            "  implicit none\n",
            "\n",
            "  private\n",
            "  \n",
            "  public :: get_associated_legendre_polynomial\n",
            "\n",
            "  contains\n",
            "\n",
            "    subroutine io_error(str)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      character(len=*) :: str\n",
            "\n",
            "      print *, \"io ERROR\"\n",
            "      print *, \"  \", str\n",
            "    end subroutine\n",
            "\n",
            "    function cot(x)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      double precision :: x\n",
            "      double precision :: cot\n",
            "\n",
            "      cot = 1.0/tan(x)\n",
            "    end function\n",
            "\n",
            "    function csc(x)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      double precision :: x\n",
            "      double precision :: csc\n",
            "\n",
            "      csc = 1.0_8/sin(x)\n",
            "    end function\n",
            "\n",
            "    function get_associated_legendre_polynomial(l,m,x) result(erg)\n",
            "      implicit none\n",
            "      ! dummy\n",
            "      integer(kind=4) :: l\n",
            "      integer(kind=4) :: m\n",
            "      double precision :: x\n",
            "\n",
            "      double precision :: erg\n",
            "\n",
            "      !auxiliary\n",
            "      double precision :: sin_x\n",
            "      double precision :: cos_x\n",
            "      \n",
            "      sin_x = sin(x)\n",
            "      cos_x = cos(x)\n",
            "      \n",
            "\n"]

def line_truncation(line, length):
    
    def add_space(line, n):
        line = (' ' * n) + line
        return line
    rv = []
    
    line_in = "".join(line)
    
    count = 0
    
    for i in line_in:
        if i == ' ':
            count += 1
        else:
            break
    
    buffer = []
    buffer_line = add_space("", count)
    
    print(len(line_in))
    
    if (len(line_in) > length):
        buffer = line_in.split()
        
        for b in buffer:            
            if (len(buffer_line) + len(b) < length):
                buffer_line = buffer_line + " " + "".join(b)
            else:
                if (len(buffer_line.lstrip()) > 0):
                    rv.append(buffer_line + " &\n")
                buffer_line = add_space("", count+4)
                buffer_line = buffer_line + "" + "".join(b)
            
        rv.append(buffer_line)
    else:
        rv = line
    return rv

with open(file_input, "r") as f:
    # read file line by line without '\n'-char
    data = f.read().splitlines()
     
    flag_at_erg_line = False
    buffer = []
    
    for line in data:
        # adapt some function calls
        line = line.replace("\\", "")   # remove backslash
        line = line.replace("endselect","end select")
        line = line.replace("*********************", "!*********************")
        line = line.replace("Sin(x)", "sin_x")
        line = line.replace("Cos(x)", "cos_x")
        
        # make double
        pattern = r"(-?[\.,0-9]+e[-,0-9]+|-?[\.,0-9]+)"
        repl = r"\1_8"
        line = re.sub(pattern, repl, line)
        
        # remove kind definition at case()
        pattern = r"case\(([-,0-9]+)_[0-9]\)"
        repl = r"case(\1)"
        line = re.sub(pattern, repl, line)
        
        # add indentation
        line = "      " + line
        
        # find "erg="-lines and append all further lines up to "case"-line
        if (line.find("erg=") > 0):
            flag_at_erg_line = True
            buffer.append(line)
        elif(line.find("case") > 0):
            if (flag_at_erg_line is True):
                flag_at_erg_line = False
                b = line_truncation(buffer, 80)
                new_data.append(" ".join(b) + "\n")
                buffer=[]
                new_data.append(line + "\n")
            else:
                new_data.append(line + "\n")
        else:
            if (flag_at_erg_line is True):
                buffer.append(line.lstrip())
            else:
                new_data.append(line + "\n")

# add footer
new_data.append("".join(["  end function get_associated_legendre_polynomial\n",
                         "end module\n"]))

# write to file
with open(file_output, "w") as f:
    f.writelines(new_data)

print("Ready")