    !>This module provides subroutines to create a progress bar.
    !!To create a progress bar, you first need to initialise the
    !!module by calling the init_progress_bar subroutine.
    !!Now your progress bar is created and you can update the
    !!progress by calling update_progress_bar. If the progress value
    !!reaches 100[%], then the progress bar will end itself automatically.
    !!Afterwards you can create a new progress bar. 
    !!You can only create one progress bar at a time.
    !!If you want to print some other text while the progress bar is running,
    !!You can use the hide_progress_bar function. Now you can print any message
    !!you want. After calling show_progress_bar, the progress bar will be visible again, 
    !!underneath everything, that was printed in between.
    !!Please use the lib_output functions to control the progress bar!!!
    
    
    module lib_progressbar
    
    implicit none
    
    character (len=15) :: process_name
    character (len=22) :: progress_string
    character (len=4) :: percentage_string
    character (len=70) :: cleaning_string
    integer :: progress_value
    logical :: initialised = .FALSE.
    logical :: hiding = .FALSE.
    
    
    
    contains
    
    !>This subroutine initialises the variables of this module.
    !!It also writes the first progress-bar to the output.
    !!
    !!@param processName This is a string type value containing the
    !!       Name of the process, that is shown by the progress bar.
    !!       This string is only allowed to have 15 characters. If
    !!       it is longer, only the first 15 characters will be shown.
    subroutine init_progress_bar(processName)
        character(len=*)::processName
        
        
        if(initialised == .TRUE.)then
            return
        endif

        if(len(processName)>15)then
            process_name = processName(1:15)
        else
            process_name = trim(processName)
        endif
        
        progress_value = 5
        
        call build_progress_string()
        call build_percentage_string()
        if(hiding == .FALSE.)then
            write(6,'(A)',advance='NO') achar(13)//process_name//" | Progress: "//percentage_string//" "//progress_string
            flush(6)
        endif
        initialised = .TRUE.
        
    end subroutine init_progress_bar
    
    !>This subroutine creates the progress_string.
    !!Therefore it uses the module intern progress_value to calculate how many
    !!Hashtags are required for the progress bar in the actual progress state.
    !!If you only want to use this library you do not have to care about this subroutine
    subroutine build_progress_string()
        integer ::numberofhashtags,help_percentage,iwl
        
        numberofhashtags=0
        progress_string = "[                    ]"
        help_percentage = progress_value - modulo(progress_value,5)
        numberofhashtags = help_percentage/5
        iwl = 2
        do while(iwl<(numberofhashtags+2))
            progress_string(iwl:iwl)='#'
            iwl = iwl + 1
        end do
    end subroutine build_progress_string
    
    !>This subroutine creates the percentage_string.
    !!Therefore it uses the module intern progress_value and turns it into a string.
    !!If you only want to use this library you do not have to care about this subroutine
    subroutine build_percentage_string()
        character(len=3) :: x1
        write (x1,'(I3.0)') progress_value
        percentage_string = trim(x1//"%")
    end subroutine build_percentage_string
    
    
    !>This subroutine will clean the line containing the progress bar.
    !!All the variables like process name and progress value will not be deleted.
    !!After executing this subroutine you can print normal output messages.
    !!After that you can then execute show_progress_bar to unhide the
    !!progress bar and make it visible again. Executing any other function,
    !!like update_progress_bar will not have any visible effect to the output
    !!while the progress bar is hidden. If the progress bar is hidden and the process
    !!will finish (with a percentage of 100%) then there is a message printed, and this
    !!module get's reseted.
    subroutine hide_progress_bar()
        if(hiding==.TRUE. .OR. initialised == .FALSE.)then
            return
        endif
        
        cleaning_string=" "
        write(6,'(A)',advance='NO') achar(13)//cleaning_string
        flush(6)
        write(6,'(A)',advance='NO') achar(13)
        flush(6)
        hiding = .TRUE.
    end subroutine hide_progress_bar
    
    !>This subroutine recreates the progress bar after it has been hidden.
    !!It is neccesarry to call this subroutine to unhide the progress bar.
    !!Only calling update_progress_bar does not unhide the progress bar.
    subroutine show_progress_bar()
        if(initialised == .TRUE. .AND. hiding==.TRUE.)then
            call build_progress_string()
            call build_percentage_string()
            write(6,'(A)',advance='NO') achar(13)//process_name//" | Progress: "//percentage_string//" "//progress_string
            flush(6)
        else
            return
        endif
        hiding = .FALSE.
    end subroutine show_progress_bar
    
    
    !>This subroutine updates the progress bar with a new value.
    !!
    !!@param percentage This parameter is an integer value.
    !!       It contains the new percentage, that should be displayed by the
    !!       progress bar.
    subroutine update_progress_bar(percentage)
        integer :: percentage
        
        if(initialised == .FALSE.)then
            return
        endif
        progress_value = percentage
        if(hiding == .FALSE.)then
            call build_progress_string()
            call build_percentage_string()
            write(6,'(A)',advance='NO') achar(13)//process_name//" | Progress: "//percentage_string//" "//progress_string
            flush(6)
        endif
        if(percentage >= 100)then
            call end_progress_bar()
        endif
    end subroutine update_progress_bar
    
    
    !>This subroutine ends the progress bar and resets this module afterwards.
    !!When using this library you don't need to call these subroutine. It is called
    !!by the update_progress_bar subroutine automatically, when a percentage value
    !!bigger or equals to 100 is submitted.
    !!However, if you want to end the process of the progress bar before it reaches 100%
    !!You can call this function and it will interrupt the progress bar.
    !!This subroutine also automatically reset's the module so you can use it afterwards
    !!for a new progress bar.
    
    subroutine end_progress_bar()
        if(initialised == .FALSE.)then
            return
        endif
        
        call build_progress_string()
        call build_percentage_string()
        write(6,'(A)',advance='YES') achar(13)//process_name//" | Progress: "//percentage_string//" "//progress_string
        flush(6)
        if(progress_value < 100)then
            write(6,*)"Process '"//trim(process_name)//"' finished unexpected early with less than 100%"
        else
            write(6,*)"Process '"//trim(process_name)//"' finished."
        endif
        call reset_progress_bar()
    end subroutine end_progress_bar
    
    
    !>This subroutine cleans the module, after a progress bar was shown.
    !!As user you don't have to take care of these function. You never need to call it.
    !!This function is called by the end_progress_bar function, after it ended the progress bar
    !!properly.
    subroutine reset_progress_bar()
        percentage_string="    "
        progress_string = "[                    ]"
        progress_value = 0
        process_name = " "
        initialised = .FALSE.
        hiding = .FALSE.
    end subroutine reset_progress_bar
    
    
end module lib_progressbar