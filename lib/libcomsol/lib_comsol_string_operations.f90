!>This module is a collection of string processing subroutines.
!!The subroutines in here are used from the lib_comsol_reader
!!modlue to process the comsol file.
    module lib_comsol_string_operations
    implicit none
    
    
    contains
    
    !>This subroutine counts how often the countChar occures in string
    !!
    !!@param:    number An integer type varaible.
    !!           After the execution of the subroutine this varaible contains how often countChar occures in string
    !!
    !!@param:    string An string type varaible. This variable contains THE string.
    !!
    !!@param:    countChar An string type variable. This is the variable stores the character, that is counted in the string.
    subroutine count_character(number,string,countChar)
        character(len=*) :: string,countChar
        integer :: number,index,stringLen
        index = 1
        number = 1
        stringLen = len(string)
        do while(index < stringLen)
            if(string(index:index) == countChar) then
                number = number + 1
            endif
            index = index + 1
        end do
    end subroutine count_character
    
    
    

    !>This subroutine removes every "remChar" in a string.
    !!
    !!@param:    string This is a string type variable. From this string all "remChar"-characters will get removed.
    !!           So after executing this subroutine, string will contain all characters like before, except "remChar"
    !!
    !!@param:    remChar This is a string type variable. It contains a character. Every occurence of this character in
    !!           string will get removed.
    subroutine remove_Character(string,remChar)
		    character(len=*) :: string,remChar
		    integer :: stringLen
		    integer :: last, actual, nofoccurence

		    stringLen = len(string)
		    last = 1
		    actual = 1
            nofoccurence = 0
            call count_character(nofoccurence,string,remChar)
            nofoccurence = nofoccurence -1

		    do while (actual < stringLen)
		        if (string(last:last) == remChar) then
		            actual = actual + 1
		            string(last:last) = string(actual:actual)
		            string(actual:actual) = remChar
		        else
		            last = last + 1
		            if (actual < last) &
		                actual = last
		        endif
		    end do
		    string = trim(string(1:len(string)-nofoccurence))

    end subroutine remove_Character
    
    !>This subroutine splits the string orig_string at the first occurence of the character index
    !!This substring is stored in the string substring. And the substring gets removed from the 
    !!orig_string.
    !!
    !!Example:
    !!orig_string = "abc,def,ghi"
    !!index = ","
    !!call pop_substring_char(orig_string,substring,index)
    !!
    !!orig_string -> def,ghi
    !!substring -> abc
    !!
    !!call pop_substring_char(orig_string,substring,index)
    !!
    !!orig_string -> ghi
    !!substring -> def
    !!
    !!
    !!
    !!@param:    orig_string This is a string type variable. It stores the string, you want to pop the substring from
    !!
    !!@param:    substring This is a string type variable, it contains the poped substring, after executing the subroutine
    !!
    !!@param:    indexg This is a string type variable, it contains the character, on which the orig_string should be splitted
    
    subroutine pop_substring_char(orig_string,substring,indexg)
        implicit none
        
         character(len=30)::orig_string,substring,indexg ,h  
         integer :: indi,g

         indi = index(orig_string,indexg)
         if(indi > 0 .AND. indi<(len(orig_string)-1)) then
            substring = orig_string(1:indi-1)
            orig_string = orig_string(indi+1:len(orig_string)-1)
        elseif(indi > 0 .AND. indi == (len(orig_string)-1)) then
            substring = orig_string(1:indi-1)
            orig_string = orig_string(indi:len(orig_string)-1)
        else
            substring = ""
            orig_string = orig_string
        endif     
        
    end subroutine pop_substring_char
    
    
    !>This subroutine checks, if the given string only contains spaces and no other characters.
    !!
    !!@param:    line This is the string, that should be checked if it is blank or not
    !!
    !!@param:    log This is a variable of the type logical. After executing the subroutine
    !!           it is .TRUE. if the given string only contains blank characters or is empty, and
    !!           .FALSE. if it contains any other characters.
    subroutine is_blank_line(line, log)
        implicit none
        
        character(len = *)::line
        character(len = 200)::hl
        logical :: log
        
        hl = line
        call remove_Character(hl," ")
        call remove_Character(hl,"\t")
        call remove_Character(hl,"\n")
        if(line == "" .OR. hl == "")then
            log = .TRUE.
            return
        else
            log = .FALSE.
            return
        endif
        log = .FALSE.
    end subroutine is_blank_line
    
    !>This function converts a string to an integer.
    !!
    !!@param:    str This is the string value, consitiong of a integer number (as string)
    !!
    !!@param:    inti This has to be a integer type variable. The converted value is stored in this
    !!           variable
    !!
    !!@param:    stat This is an integer type variable. It an error variable. If this variable
    !!           is unequals 0, an error occured during the conversion.
    subroutine str2int(str,inti,stat)
        implicit none

        character(len=*),intent(in) :: str
        integer,intent(out)         :: inti
        integer,intent(out)         :: stat
    
        read(str,*,iostat=stat)  inti
    end subroutine str2int
    
    !>This function converts a string to an real number.
    !!
    !!@param:    str This is the string value, consitiong of a real number (as string)
    !!
    !!@param:    r This has to be a real type variable. The converted value is stored in this
    !!           variable
    !!
    !!@param:    stat This is an integer type variable. It an error variable. If this variable
    !!           is unequals 0, an error occured during the conversion.
    subroutine str2real(str, r,stat)
    	implicit none

    	character(len=*),intent(in) 	:: str
    	real(dp), intent(out)		:: r
    	integer,intent(out)		:: stat
    
    	read(str,*,iostat=stat) r
    end subroutine str2real
        
        
end module lib_comsol_string_operations