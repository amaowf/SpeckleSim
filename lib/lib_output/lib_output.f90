!>The lib_output module gives some basic output and logfile features.
!!The user can select between different output functions. Some of them are
!!writing only to the standard output, others are also writing to a logfile.
!!This gives the possibility to easily write messages to the output and a
!!log file without any more coding required.
module lib_output

use lib_progressbar

implicit none

character (len=255) :: filePathAndName= "SpeckleSimulator.log"
logical :: progressbar = .FALSE.



contains
  !>This subroutine set's the logfile name and the path to it.
  !!The standard logfile is called SpeckleSimulator.log and is created in
  !!the directory where the executable is located.
  !!
  !!@param: fpn The logfile path and name
  subroutine set_log_file(fpn)
    character (len=*) :: fpn
    filePathAndName = fpn
  end subroutine set_log_file


  !>This subroutine prints the given string to the standard output.
  !!That means it is only printed to the user and not written to the log file.
  !!
  !!@param: outputtouser The string which should be written to the standard output
  subroutine print_only_u(outputtouser)
    character (len=*) :: outputtouser

    if(progressbar == .TRUE.)then
        call hide_progress_bar()
        write(6,*) outputtouser
        call show_progress_bar()
    else
        write(6,*) outputtouser
    endif
  end subroutine print_only_u
  
  !>This subroutine reads a user input string from the standard input.
  !!The subroutine also doesn't store any log messages.
  !!
  !!@param: The input from the user
  subroutine read_only_u(userString)
    character (len=*) :: userString
    if(progressbar == .TRUE.)then
        call hide_progress_bar()
        read(6,*) userString
        call show_progress_bar()
    else
        write(6,*) userString
    endif
  end subroutine read_only_u


  !>This subroutine prints error messages. It additionally figures out the actual
  !!date and time and prints the actual time and date together with the
  !!given error message to the user. Error messages are also written to the logfile
  !!
  !!@param: errorMessage The message you want to be submitted to the user. This
  !!        Message is also written to the log file.
  subroutine print_error(errorMessage)
    character (len=*) :: errorMessage
    character (len=50) :: dateTimestring

    call get_time_string(dateTimestring)
    
    if(progressbar == .TRUE.)then
        call hide_progress_bar()
        write(1,*)trim(dateTimestring)//"Error Message:"
        write(1,*) errorMessage
        call show_progress_bar()
    else
        write(1,*)trim(dateTimestring)//"Error Message:"
        write(1,*) errorMessage
    endif

    call print_to_log("Error Message: "//errorMessage)

  end subroutine print_error


  !>This subroutine writes the given string to the standard output and to the log file.
  !!
  !!@param: outputtouser The string you want to be printed to the standard output
  !!        and to the logfile.
  subroutine print_u(outputtouser)
    character (len=*) :: outputtouser
    
    if(progressbar == .TRUE.)then
        call hide_progress_bar()
        write(6,*) outputtouser
        call show_progress_bar()
    else
        write(6,*) outputtouser
    endif

    call print_to_log(outputtouser)


  end subroutine print_u
  
  !>This subroutine reads the user input to the standard input.
  !!If logEN is set to .TRUE. it will also store the user input in the log file.
  !!
  !!@param userString A character type parameter, that stores the user input string afterwards
  !!
  !!@param logEN A logical type parameter, that specifies, if the user input should be stored in the log file.
  subroutine read_u(userString, logEN)
    character (len=*) :: userString
    logical ::logEN
    
    if(progressbar == .TRUE.)then
        call hide_progress_bar()
        read(6,*) userString
        call show_progress_bar()
    else
        write(6,*) userString
    endif
    
    if(logEN == .TRUE.)then
        call print_to_log("User input: <"//userString//">")
    endif
    
  end subroutine read_u


  !>This function writes strings to the logfile. It opens the file specified before
  !!and prints the given string to it. All logfile entries are appended with a
  !!timestamp.
  !!
  !!@param: text The string you want to be written to the logfile
  subroutine print_to_log(text)
    character (len=*) :: text
    character (len=50) :: dateTime
    integer :: ioStatus
    ioStatus = 0

    call get_time_string(dateTime)

    open(19,file = filePathAndName,iostat=ioStatus,position="append")
      if(ioStatus /= 0) stop "Error opening log file..."
      write(19,*) trim(dateTime)//trim(text)
    close(19)

  end subroutine print_to_log


  !>This function figures out the actual date and time. The result is converted in
  !!The following format:
  !!DD.MM.YYYY at HH:MM:SS  ->
  !!This string is then written to the timeString given by the user.
  !!The result is used in the print_to_log and in the print_error function.
  !!
  !!@param: timeString The date and time is stored in this string after executing this subroutine.
  subroutine get_time_string(timeString)
    character (len=*) :: timeString
    character (len=15) :: dfate, tfime
    character (len=20) :: readableDate, readableTime

    call date_and_time(DATE=dfate)
    call date_and_time(TIME=tfime)

    readableDate = dfate(7:8)//"."//dfate(5:6)//"."//dfate(1:4)
    readableTime = tfime(1:2)//":"//tfime(3:4)//":"//tfime(5:6)

    timeString = trim(readableDate)//" at "//trim(readableTime)//" -> "

  end subroutine get_time_string
  
  !>This subroutine creates the progress bar.
  !!If you have a progress bar running, you can still use the print functions above.
  !!This won't destroy the progress bar, everything you are printing is printed
  !!above the progress bar.
  !!If you use only the write or print funktions provided by fortran, this
  !!will mess up your output.
  !!
  !!@param: pbName This is the process-Name, that is shown in realtion with
  !!         the progress bar. This is a string variable. Only the first 15
  !!         characters will be displayed
  subroutine create_progressbar(pbName)
    character (len=*)::pbName
    if(progressbar == .TRUE.) return
    call init_progress_bar(pbName)
    progressbar = .TRUE.
  end subroutine create_progress_bar
  
  !>If the process shown by the progress bar is progressing you can
  !!submit the new progress value with this function. And the progress
  !!bar will be updated accordingly.
  !!
  !!@param: percVal This is an integer type variable. You have to submit a
  !!         value in percent. So allowed are numbers from 1,2,3...99,100
  !!         if the value is 100, the progress bar will be closed automatically
  subroutine update_progressbar(percVal)
    integer :: percVal
    
    if(progressbar == .FALSE.) return
    call update_progress_bar(percVal)
    if(percVal >= 100)then
        progressbar = .FALSE.
    endif
  end subroutine update_progressbar
  
  !>This subroutine closes the progress bar. Normally you don't need to call this
  !!subroutine. But if your process ends for any reason with less then 100% progress,
  !!you can call this subroutine and the progress bar will be closed.
  !!But in every other case, you do not need to call it.
  subroutine end_progressbar()
    if(progressbar == .FALSE.) return
    call end_progress_bar()
    progressbar = .FALSE.
  end subroutine end_progressbar
  
  !>Do not use this subroutine.
  !!This subroutine should only be used, in case, that the progress bar is doing strange things.
  !!So if the format of the progress bar on any point sucks, you can try to call this subroutine.
  !!But if everythings works like expected, don't touch this.
  subroutine reset_progressbar()
    call reset_progress_bar()
    progressbar = .FALSE.
  end subroutine reset_progressbar
    


end module lib_output
