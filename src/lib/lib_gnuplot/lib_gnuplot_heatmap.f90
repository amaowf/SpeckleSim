!>This module can be used to print heatmap diagrams. To be able to use this module
!!you need to have gnuplot installed on your system. Gnuplot must also be available in
!!your system path.
    
module lib_gnuplot_heatmap
    
    implicit none
    character (len=200)::DataFilePath
    character (len=50)::diagram_title,diagram_x_text,diagram_y_text,output_file_name
    character (len=80)::outputType,heatmap_colour,fontType,pltfilename
    character (len = 1000)::customParameter
    integer :: sizeX,sizeY,textSize
    logical :: grid,colour_set
    logical :: initialised = .FALSE.
    
    
    contains
    
    !>This subroutine is used to initialise the module.
    !!To print the most basic lineplot, you need to call this subroutine
    !!Then you need to call create_plot_file_lineplot and afterwards
    !!plot_diagram_lineplot
    !!
    !!
    !!All parameters get initialised by standard values.
    !!The standard parameters are:
    !!Diagram-Title: 'Heatmap'
    !!X-Achsis Text: 'x'
    !!Y-Achsis Text: 'y'
    !!Output file name: lineplot
    !!Colour of the plot: gnuplot's default colour scheme
    !!Output file type: png
    !!font: arial
    !!font size: 15
    !!Size in X-Direction: 700px
    !!Size in Y-Direction: 700px
    !!A grid is not displayed
    !!
    !!
    !!To be able to plot a lineplot from data stored in a file, the file must have
    !!The following format:
    !!
    !!  #X-Val          Y-Val       Z-Val(This value defines the colour at the point)     
    !!  #(Every line starting with '#' is a comment)
    !!  -0.9800000     -0.9411921    0
    !!  -0.9600000     -0.8847359    1
    !!  -0.9400000     -0.8305840    2
    !!  -0.9200000     -0.7786880    3
    !!  -0.9000000     -0.7289999    0
    !!  -0.8800000     -0.6814720    1
    !!  -0.8600000     -0.6360560    2
    !!  -0.8400000     -0.5927041    3
    !!  -0.8200000     -0.5513680    0
    !!  -0.8000000     -0.5120000    1
    !!^^-0.7800000     -0.4745519    2
    !!    
    !!@param filenameandpath This is a string type parameter with the filename (and it's path) to the data file
    subroutine set_data_heatmap(filenameandpath)
        character(len=*)::filenameandpath
        DataFilePath = filenameandpath
        grid = .FALSE.
        colour_set = .FALSE.
        diagram_title = "Heatmap"
        diagram_x_text = "x"
        diagram_y_text = "y"
        sizeX = 700
        sizeY = 700
        outputType = "png"
        fontType = "arial"
        textSize = 15
        output_file_name = "heatmap_plot"
        initialised = .TRUE.
        
    end subroutine set_data_heatmap
    
    !>By default no grid net is shown on the diagram. You
    !!can controll that with this subroutine
    !!
    !!@param booleaniforifnot A logical type paramter.
    !!                       setting this to .TRUE. enables the grid. .FALSE. disables it.
    subroutine set_grid_heatmap(booleaniforifnot)
        logical:: booleaniforifnot
        grid = booleaniforifnot
    end subroutine set_grid_heatmap
    
    !>With this subroutine you can specify the colour scheme of the heatmap.
    !!You need to submit a string formated in the "gnuplot format" for colour
    !!strings.
    !!This means you need to use the following format:
    !!   ( 0 'red',1 'black')"
    !!You can add as many colours as you want. Do do this you can just add another
    !!number to the list and specify the colour behind. To specify the colour you can either use
    !!one of gnuplots default colour or specify a color by it's hex value.
    !!This means you can generate a heatmap with the colours green, blue and yellow with the following
    !!colour string:
    !!   (0 'green', 1 'blue', 2 'yellow')
    !!
    !!@param colorVal A string type parameter with the colour information, in the above described format.
    subroutine colour_heatmap(colorVal)
        character(len=*)::colorVal
        colour_set = .TRUE.
        heatmap_colour = colorVal
    end subroutine colour_heatmap

    !>This subroutine makes it possible to customise the title of the plot.
    !!
    !!@param titletext A string type parameter. This should be the title you want to
    !!               be displayed above the diagram.
    subroutine title_text_heatmap(titletext)
        character(len=*)::titletext
        diagram_title = titletext
    end subroutine title_text_heatmap
    
    
    !>This subroutine makes it possible to customise the X-Achsis text of the plot.
    !!
    !!@param titletext A string type parameter. This should be the x-Achsis text you want to
    !!               be displayed above the diagram.
    subroutine xAchsis_text_heatmap(xtext)
        character(len=*)::xtext
        diagram_x_text = xtext
    end subroutine xAchsis_text_heatmap
    
    
    !>This subroutine makes it possible to customise the Y-Achsis text of the plot.
    !!
    !!@param titletext A string type parameter. This should be the y-Achsis text you want to
    !!               be displayed above the diagram.
    subroutine yAchsis_text_heatmap(ytext)
        character(len=*)::ytext
        diagram_y_text = ytext
    end subroutine yAchsis_text_heatmap
    
    !>With this subroutine you can specify the resolution of the output.
    !!The values are in pixels.
    !!
    !!@param x_resolution An integer type parameter, which contains the resolution in
    !!                   x-direction
    !!@param y_resolution An integer type parameter, which contains the resolution in
    !!                   y-direction
    subroutine output_resolution_heatmap(x_resolution,y_resolution)
        integer :: x_resolution,y_resolution
        sizeX = x_resolution
        sizeY = y_resolution
    end subroutine output_resolution_heatmap
    
    
    !>With this subroutine you can specify, which format the output file should have.
    !!
    !!@param otype A string type parameter with the type of the output file. (png,svg,pdf,ps...)
    !!           See the possible gnuplot outputs for further details.
    subroutine output_type_heatmap(otype)
        character(len=*)::otype
        outputType = otype
    end subroutine output_type_heatmap
    
    !>With this subroutine you can specify a custom font and font-size.
    !!Please pay attention, that GNUPLOT relies to external fonts.
    !!Therefore you can only choose fonts, that are available on that system.
    !!If the specified font is not available, GNUPLOT will use it's default font.
    !!The size parameter works anyway.
    !!
    !!@param font A string type parameter for the name of the font
    !!
    !!@param size An Integer type parameter for the text size.
    subroutine set_font_and_size_heatmap(font,size)
        character(len=*)::font
        integer :: size
        fontType = font
        textSize = size
    end subroutine set_font_and_size_heatmap
    
    !>This subroutine creates a .plt file with all the specifications from above.
    !!This file is stored in the same directory in which the output file will be generated to.
    subroutine create_plot_file_heatmap()
        if(.NOT. initialised)then
            return
        endif
        
        pltfilename = trim(output_file_name)//"hm-plotting.plt"
        OPEN(UNIT=48,FILE=trim(pltfilename))
            WRITE(48,*) "set view map"
            WRITE(48,*) "set dgrid3d 100,100,2"
            WRITE(48,*) "set pm3d"
            !set terminal png size 1200,900 font "arial,20"
            WRITE(48,'(A)',advance = 'no') "set terminal "
            WRITE(48,'(A)',advance = 'no') trim(outputType)
            WRITE(48,'(A)',advance = 'no') " size "
            WRITE(48,'(A)',advance = 'no') trim(intToStr(sizeX))//","//trim(intToStr(sizeY))
            WRITE(48,'(A)',advance = 'no') ' font "'
            WRITE(48,*) trim(fontType)//","//trim(intToStr(textSize))//'"'
            !set output "output2.png"
            WRITE(48,'(A)',advance = 'no') 'set output "'
            WRITE(48,*) trim(output_file_name)//"."//trim(outputType)//'"'
            !set title 'Affensattel'
            WRITE(48,'(A)',advance = 'no') 'set title "'
            WRITE(48,*) trim(diagram_title)//'"'
            !set xlabel "x"
            WRITE(48,'(A)',advance = 'no') 'set xlabel "'
            WRITE(48,*) trim(diagram_x_text)//'"'
            !set ylabel "y"
            WRITE(48,'(A)',advance = 'no') 'set ylabel "'
            WRITE(48,*) trim(diagram_y_text)//'"'
            !set grid
            if(grid)then
                WRITE(48,*) "set grid"
            else
                WRITE(48,*) "unset grid"
            endif
            !set palette defined ( 0 'red',1 'black')
            if(colour_set)then
                WRITE(48,'(A)',advance = 'no') "set palette defined "
                WRITE(48,*) trim(heatmap_colour)
            endif
            !unset key #Hierdurch wird der Name der Datei im Diagram angezeigt
            !unset surface#Hierdurch wird an jedem Datenpunkt ein Kreuz angezeigt.
            WRITE(48,*) "unset key"
            WRITE(48,*) "unset surface"
            !splot "data2"
            WRITE(48,'(A)',advance = 'no') 'splot "'
            WRITE(48,'(A)',advance = 'no') trim(DataFilePath)//'"'
        CLOSE(48)
    end subroutine create_plot_file_heatmap
    
    
    !>After generating the .plt file, you can use this subroutine to call gnuplot and create
    !!the plot from the .plt file.
    subroutine plot_diagram_heatmap()
        if(.NOT. initialised)then
            return
        endif
        CALL SYSTEM('gnuplot -p '//trim(pltfilename))
        initialised = .FALSE.
    end subroutine plot_diagram_heatmap
    
    !This is an module internal function, that is used to convert integer values to strings
    !Don't use ore care about when you only want to use this module.
    !
    !@param k An integer type parameter. This integer is converted into a string number
    !
    !@returns An string containing the number of the integer given to the function
    character(len=20) function intToStr(k)
        integer, intent(in) :: k
        write (intToStr, *) k
        intToStr = adjustl(intToStr)
    end function intToStr
    
    

end module