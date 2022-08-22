!> This module reads the json configuration file and extracts the values required.
!! To use this module properly, you first need to initialise the module by calling
!! the init_json_lib function. This function asks the user for the filepath and name
!! of the specification json file.
!! It also does some basic check, if the program can open the file.
!! If not, the program will be stopped.
!!
!! After that you can get the values stored in the json file with the
!! read_real_data
!! read_int_data
!! read_string_data
!! functions. Just choose the matching subroutine to your datatype.
!! If the requested value cannot be found in the json file, the program will also
!! be interrupted, and an error message gets displayed.
!!
!! Therefore it is required, to put all values you need in the json file.
!! You should also pay attention,that there is only one value-name pair per line in the json file.
!! This means the library actually does not support several value-name pairs in one line.
!!
!!
module json_input_reader

	implicit none
    
    !This is the value, that stores the file path and filename. It is initialised in the init_json_lib subroutine
    character(len = 300) , private :: specification_file_path_name
    !This variable stores if the library is initialised or not. This prevents from using the read functions, while
    !init_json_lib wasn't called.
    logical, private :: initialised
    
	
    contains
    
    
    
    !The subroutines below are not supposed to be used out of this module. They are used within the main subroutines of the
    !module, and provide somehow "helper-tasks"
    
    
    
    
    
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
        index = 0
        number = 0
        stringLen = len(string)
        do while(index < stringLen)
            if(string(index:index) == countChar) then
                number = number + 1
            endif
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



!>This function extracts the name of the parameter stored in the actual json-file line.
!!Before using this function you need to pre-process the string such, that it has the following format:
!!   parameter_name:value
!!There are no blanks, brackets or something similar allowed. Then this function will return a string,
!!containing the parameter_name
!!
!!@param:    line This must be a string type variable containing the pre-processed json-entry-line.
!!           After executing this subroutine the variable will only contain the parameter_name part.
subroutine get_Value_Name(line)
	character(len = *) :: line
	integer :: i
	i = 0
	!splitPos = index(line,':')
	if(index(line,':') <= 1)then
		line = ""
	else
		line = line(1:index(line,':')-1)
	endif

end subroutine get_Value_Name




!>This subroutine can be used to finally extract the value from a line of the json file.
!!Before using this function you need to pre-process the string such, that it has the following format:
!!   parameter_name:value
!!There are no blanks, brackets or something similar allowed. Then this function will return a string,
!!containing the value
!!
!!@param:    line This must be a string type variable containing the pre-processed json-entry-line.
!!           After executing the subroutine this variable will hold only the value part of the string.
subroutine get_Value(line)
	character(len=*) :: line

	call remove_Character(line,',')
	if(index(line,":") <=1)then
		line = ""
	else
		line = line(index(line,":")+1:len(line))
	endif


end subroutine get_Value





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
	real, intent(out)		:: r
	integer,intent(out)		:: stat

	read(str,*,iostat=stat) r
end subroutine str2real




!The subroutines below are supposed to be used outside this module.
!They can be used, to initialise the module and read the data from the json file.






!>This subroutine reads real type parameters from the json file.
!!If the specified parameter name cannot be found in the json file,
!!The program will be interrupted and an error message is printed.
!!
!!@param:    parameter_name This is the parameter-name used in the json file
!!           to specify the value.
!!@param:    real_variable This has to be a real variable. The value, stored
!!           in the specification file is written to this variable
subroutine read_real_data(parameter_name, real_variable)
	implicit none
	
	character (len = *) :: parameter_name
	real :: real_variable
	!Declaring the variables
	character (len=60) :: filecontent,valueName,valueString
	integer :: v,s,ioStatus
	real	:: rVal
    
    if(.NOT. initialised) then
        write(*,*) "Error, the json-library was used without previous initialisation\nPlease call init_json_lib, before trying to read a value\nExiting"
        stop
    endif


	!Opening the data.json File
	open(1,file = trim(specification_file_path_name),status = "old",iostat = ioStatus)

	if(ioStatus /= 0) then
		write(*,*) "Error while opening json-Configuration file\nExiting"
        stop
	endif
	

	do
		read(1,'(A)',iostat = ioStatus) filecontent
		if(ioStatus /= 0) exit

		call remove_Character(filecontent,' ')
		call remove_Character(filecontent,',')
		call remove_Character(filecontent,'"')

		valueName = filecontent
		valueString = filecontent
		call get_Value_Name(valueName)

		if(trim(valueName) == trim(parameter_name)) then
                call get_Value(valueString)
      			call str2real(valueString,rVal,s)
      			real_variable = rVal
      			return
      			
		endif


	end do

	close(1)
    
    write(*,*) "Error\nThe parameter name "//parameter_name//" could not be found in the specified\nconfiguration file."
    stop

end subroutine read_real_data



!>This subroutine reads int type parameters from the json file.
!!If the specified parameter name cannot be found in the json file,
!!The program will be interrupted and an error message is printed.
!!
!!@param:    parameter_name This is the parameter-name used in the json file
!!           to specify the value.
!!@param:    int_variable This has to be a integer variable. The value, stored
!!           in the specification file is written to this variable
subroutine read_int_data(parameter_name, int_variable)
	implicit none
	
	
	integer :: int_variable
	character (len = *) :: parameter_name
	!Declaring the variables
	character (len=60) :: filecontent,valueName,valueString
	integer :: v,s,ioStatus,iVal

    if(.NOT. initialised) then
        write(*,*) "Error, the json-library was used without previous initialisation\nPlease call init_json_lib, before trying to read a value\nExiting"
        stop
    endif
    
    
	!Opening the data.json File
	open(1,file = trim(specification_file_path_name),status = "old",iostat = ioStatus)

	if(ioStatus /= 0) then
        write(*,*) "Error while opening json-Configuration file\nExiting"
        stop
	endif
	

	do
		read(1,'(A)',iostat = ioStatus) filecontent
		if(ioStatus /= 0) exit

		call remove_Character(filecontent,' ')
		call remove_Character(filecontent,',')
		call remove_Character(filecontent,'"')

		valueName = filecontent
		valueString = filecontent
		call get_Value_Name(valueName)
		
		if(trim(valueName) == trim(parameter_name)) then
                call get_Value(valueString)
      			call str2int(valueString,iVal,s)
      			int_variable = iVal
      			return
      			
		endif


	end do

	close(1)
    
    write(*,*) "Error\nThe parameter name "//parameter_name//" could not be found\nin the specified configuration file.\nExiting"
    stop

end subroutine read_int_data

!>This subroutine reads string type parameters from the json file.
!!If the specified parameter name cannot be found in the json file,
!!The program will be interrupted and an error message is printed.
!!
!!@param:    parameter_name This is the parameter-name used in the json file
!!           to specify the value.
!!@param:    string_variable This has to be a string variable. The value, stored
!!           in the specification file is written to this variable
subroutine read_string_data(parameter_name, string_variable)
	implicit none
	
	
	character (len = *) :: parameter_name, string_variable
	character (len=60) :: filecontent,valueName,valueString
	integer :: v,s,ioStatus,iVal

    
    if(.NOT. initialised) then
        write(*,*) "Error, the json-library was used without previous initialisation\nPlease call init_json_lib, before trying to read a value\nExiting"
        stop
    endif
    
	!Opening the data.json File
	open(1,file = trim(specification_file_path_name),status = "old",iostat = ioStatus)

	if(ioStatus /= 0) then
        write(*,*) "Error while opening json-Configuration file\nExiting"
        stop
	endif

	do
		read(1,'(A)',iostat = ioStatus) filecontent
		if(ioStatus /= 0) exit

		call remove_Character(filecontent,' ')
		call remove_Character(filecontent,',')
		call remove_Character(filecontent,'"')

		valueName = filecontent
		valueString = filecontent
		call get_Value_Name(valueName)

		if(trim(valueName) == trim(parameter_name)) then
                call get_Value(valueString)
      			string_variable = trim(valueString)
      			return
		endif

	end do

	close(1)
    
    write(*,*) "Error\nThe parameter name "//parameter_name//" could not be found in the specified\nconfiguration file."
    stop

end subroutine read_string_data



!>This subroutine initialises this library.
!!It requests the filepath and file-name from the user and "pre-checks" the input.
!!
!!If the given file does not exist, or cannot be opened, the program will be interrupted, and an 
!!Error-Message is printed.

subroutine init_json_lib()
    implicit none
    
    integer :: ioStatus
    
    write(*,*)"Enter the configuration file-path and file-name:"
    read(*,*) specification_file_path_name
    
    open(1,file = trim(specification_file_path_name),status = "old",iostat = ioStatus)
	if(ioStatus /= 0) then
        write(*,*) "Error with specified configuration file/path.\nPlease check your input and try again.\n\nYour input was: "//specification_file_path_name//"\n\nExiting"
        close(1)
        stop
    endif
    close(1)
    
    initialised = .TRUE.
    
    
   
end subroutine init_json_lib

	



end module json_input_reader
