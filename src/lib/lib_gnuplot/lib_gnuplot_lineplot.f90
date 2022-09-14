!>This module provides functions to plot a 2D lineplot from data, stored in a file.
!!To be able to use this module you need to have gnuplot installed. It should also be available on your system path.
!!otherwise gnuplot cannot be found by fortran.
    module lib_gnuplot_lineplot
    
    implicit none
    
    
    character (len=200)::DataFilePath
    character (len=50)::diagram_title,diagram_x_text,diagram_y_text,output_file_name,linecolour
    character (len=80)::outputType,lineplot_linetype,fontType,typeLine,pltfilename
    integer :: sizeX,sizeY,textSize
    logical :: grid
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
    !!Diagram-Title: 'Lineplot'
    !!X-Achsis Text: 'x'
    !!Y-Achsis Text: 'y'
    !!Output file name: lineplot
    !!Colour of the plot: Blue
    !!Output file type: png
    !!font: arial
    !!font size: 15
    !!Size in X-Direction: 700px
    !!Size in Y-Direction: 700px
    !!A grid is displayed
    !!
    !!
    !!To be able to plot a lineplot from data stored in a file, the file must have
    !!The following format:
    !!
    !!   #X-Val          Y-Val      (Every line starting with '#' is a comment)
    !!   -0.9800000     -0.9411921    
    !!   -0.9600000     -0.8847359    
    !!   -0.9400000     -0.8305840    
    !!   -0.9200000     -0.7786880    
    !!   -0.9000000     -0.7289999    
    !!   -0.8800000     -0.6814720    
    !!   -0.8600000     -0.6360560    
    !!   -0.8400000     -0.5927041    
    !!   -0.8200000     -0.5513680    
    !!   -0.8000000     -0.5120000    
    !!   -0.7800000     -0.4745519    
    !!    
    !!@param filenameandpath This is a string type parameter with the filename (and it's path) to the data file
    subroutine set_data_lineplot(filenameandpath)
        character(len=*)::filenameandpath
        
        DataFilePath = filenameandpath
        diagram_title = "Lineplot"
        diagram_x_text = "x"
        diagram_y_text = "y"
        output_file_name = "lineplot"
        linecolour = "blue"
        outputType = "png"
        fontType = "arial"
        typeLine = "lines"
        sizeX = 700
        sizeY = 700
        textSize = 15
        grid = .TRUE.
        initialised = .TRUE.
        
    end subroutine set_data_lineplot
    
    !>By default a grid net is shown on the diagram. You
    !!can controll that with this subroutine
    !!
    !!@param booleaniforifnot A logical type paramter.
    !!                       setting this to .TRUE. enables the grid. .FALSE. disables it.
    subroutine set_grid_lineplot(booleaniforifnot)
        logical:: booleaniforifnot
        grid = booleaniforifnot
    end subroutine set_grid_lineplot
    
    !>With this subroutine you can specify the colour of the line which is
    !!plotted between the data points. By default this line will be blue.
    !!
    !!@param: colorVal A string type parameter. You can either use hex-colour values
    !!               or use one of gnuplots default colour values.(red,green,blue,black,grey,gold,.......)
    subroutine colour_lineplot(colorVal)
        character(len=*)::colorVal
        linecolour = colorVal
    end subroutine colour_lineplot

    !>This subroutine makes it possible to customise the title of the plot.
    !!
    !!@param titletext A string type parameter. This should be the title you want to
    !!               be displayed above the diagram.
    subroutine title_text_lineplot(titletext)
        character(len=*)::titletext
        diagram_title = titletext
    end subroutine title_text_lineplot
    !>This subroutine makes it possible to customise the X-Achsis text of the plot.
    !!
    !!@param titletext A string type parameter. This should be the x-Achsis text you want to
    !!               be displayed above the diagram.
    subroutine xAchsis_text_lineplot(xtext)
        character(len=*)::xtext
        diagram_x_text = xtext
    end subroutine xAchsis_text_lineplot
    
    !>This subroutine makes it possible to customise the Y-Achsis text of the plot.
    !!
    !!@param titletext A string type parameter. This should be the y-Achsis text you want to
    !!               be displayed above the diagram.
    subroutine yAchsis_text_lineplot(ytext)
        character(len=*)::ytext
        diagram_y_text = ytext
    end subroutine yAchsis_text_lineplot
    
    !>With this subroutine you can specify the resolution of the output.
    !!The values are in pixels.
    !!
    !!@param x_resolution An integer type parameter, which contains the resolution in
    !!                   x-direction
    !!@param y_resolution An integer type parameter, which contains the resolution in
    !!                   y-direction
    subroutine output_resolution_lineplot(x_resolution,y_resolution)
        integer :: x_resolution,y_resolution
        sizeX = x_resolution
        sizeY = y_resolution
    end subroutine output_resolution_lineplot
    
    !>With this subroutine you can specify, which format the output file should have.
    !!
    !!@param otype A string type parameter with the type of the output file. (png,svg,pdf,ps...)
    !!           See the possible gnuplot outputs for further details.
    subroutine output_type_lineplot(otype)
        character(len=*)::otype
        outputType = otype
    end subroutine output_type_lineplot
    
    subroutine linetype_lineplot(linetype)
        character(len=*)::linetype
        typeLine = linetype
    end subroutine linetype_lineplot
    

    !>With this subroutine you can specify a custom font and font-size.
    !!Please pay attention, that GNUPLOT relies to external fonts.
    !!Therefore you can only choose fonts, that are available on that system.
    !!If the specified font is not available, GNUPLOT will use it's default font.
    !!The size parameter works anyway.
    !!
    !!@param font A string type parameter for the name of the font
    !!
    !!@param size An Integer type parameter for the text size.
    subroutine set_font_and_size_lineplot(font,size)
        character(len=*)::font
        integer :: size
        textSize = size
        fontType = font
    end subroutine set_font_and_size_lineplot
    
    
    !>This subroutine creates a .plt file with all the specifications from above.
    !!This file is stored in the same directory in which the output file will be generated to.
    subroutine create_plot_file_lineplot()
        pltfilename = trim(output_file_name)//"lp-plotting.plt"
        OPEN(48,FILE=trim(pltfilename))
        !unset key #Hierdurch wird der Name der Datei im Diagram angezeigt
        WRITE(48,*)"unset key"
        !set title 'X-Achse Affensattel'
        WRITE(48,'(A)',advance = 'no') 'set title "'
        WRITE(48,*) trim(diagram_title)//'"'
        !set output "output1.png"
        WRITE(48,'(A)',advance = 'no') 'set output "'
        WRITE(48,*) trim(output_file_name)//"."//trim(outputType)//'"'
        !set grid
        if(grid)then
            WRITE(48,*) "set grid"
        else
            WRITE(48,*) "unset grid"
        endif
        !set terminal png size 1000,1000 font "arial,20"
        WRITE(48,'(A)',advance = 'no') "set terminal "
        WRITE(48,'(A)',advance = 'no') trim(outputType)
        WRITE(48,'(A)',advance = 'no') " size "
        WRITE(48,'(A)',advance = 'no') trim(intToStr_lp(sizeX))//","//trim(intToStr_lp(sizeY))
        WRITE(48,'(A)',advance = 'no') ' font "'
        WRITE(48,*) trim(fontType)//","//trim(intToStr_lp(textSize))//'"'
        !set xlabel "x"
        WRITE(48,'(A)',advance = 'no') 'set xlabel "'
        WRITE(48,*) trim(diagram_x_text)//'"'
        !set ylabel "y"
        WRITE(48,'(A)',advance = 'no') 'set ylabel "'
        WRITE(48,*) trim(diagram_y_text)//'"'
        !plot "data1" with lines lt rgb "green"
        WRITE(48,'(A)',advance = 'no') 'plot "'
        WRITE(48,'(A)',advance = 'no') trim(DataFilePath)//'" with '//trim(typeLine)//' lt rgb "'//linecolour//'"'
        CLOSE(48)
    end subroutine create_plot_file_lineplot
    
    
    !>After generating the .plt file, you can use this subroutine to call gnuplot and create
    !!the plot from the .plt file.
    subroutine plot_diagram_lineplot()
        if(.NOT. initialised)then
            return
        endif
        CALL SYSTEM('gnuplot -p '//trim(pltfilename))
        initialised = .FALSE.
    end subroutine plot_diagram_lineplot
    
    !This is an module internal function, that is used to convert integer values to strings
    !Don't use ore care about when you only want to use this module.
    !
    !@param k An integer type parameter. This integer is converted into a string number
    !
    !@returns An string containing the number of the integer given to the function
    character(len=20) function intToStr_lp(k)
        integer, intent(in) :: k
        write (intToStr_lp, *) k
        intToStr_lp = adjustl(intToStr_lp)
    end function intToStr_lp
    
    
    
end module lib_gnuplot_lineplot