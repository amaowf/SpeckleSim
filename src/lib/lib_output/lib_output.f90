!>The lib_output module gives some basic output and logfile features.
!!The user can select between different output functions. Some of them are
!!writing only to the standard output, others are also writing to a logfile.
!!This gives the possibility to easily write messages to the output and a
!!log file without any more coding required.
module lib_output

implicit none

character (len=255) :: filePathAndName= "SpeckleSimulator.log"



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

    write(6,*)outputtouser
  end subroutine print_only_u


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

    write(1,*)trim(dateTimestring)//"Error Message:"
    write(1,*) errorMessage

    call print_to_log("Error Message: "//errorMessage)

  end subroutine print_error


  !>This subroutine writes the given string to the standard output and to the log file.
  !!
  !!@param: outputtouser The string you want to be printed to the standard output
  !!        and to the logfile.
  subroutine print_u(outputtouser)
    character (len=*) :: outputtouser

    write(6,*) outputtouser

    call print_to_log(outputtouser)


  end subroutine print_u


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


end module lib_output
