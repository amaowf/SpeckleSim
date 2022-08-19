
module lib_comsol_reader
    use lib_comsol_string_operations
    
    implicit none
    
    
    
    contains
    
    !>This subroutine figures out, how many points, and how many triangles are specified by
    !!the comsol file.
    !!To do that the subroutine opens the file and iterates over the lines.
    !!The basic structure of this subroutine is pretty similar to the structure of the
    !!read_comsol_file subroutine. So if you are interested in how this subroutine works
    !!it is recomended to read the comments on how the read_comsol_file subroutine works first.
    !!
    !!@param:    filename_path This is a string type variable, containing the filename and path to that file
    !!
    !!@param:    numb_coord This is an integer type variable. It holds the information on how many points are
    !!           specified, by it's coordinates in the file(This value is 1 based)
    !!
    !!@param:    numb_tri This is an integer type variable. It holds the information on how many triangles are
    !!
    
    
    specified in the file(This value is 1 based)
    subroutine number_of_coordinates_triangels(filename_path, numb_coord,numb_tri)
        implicit none
        !Initialisation, Opening the comsol file [...]
        character(len = *) :: filename_path
        character(len = 200) :: linecontent
        character(len = 100) :: numberSubstring
        integer :: ioStatus,numberofHashtags
        integer :: numb_coord,numb_tri
        logical :: isBlankLine
        
        numb_coord = 0
        numb_tri = 0
        
        
        open(2,file = trim(filename_path),status = "old",iostat = ioStatus)
        if(ioStatus /= 0) then
            write(*,*) "Error while opening comsol file.\nThe specified comsol file was: "//filename_path//"\n\nExiting"
            close(1)
            stop
        endif
    
        
        !Searching for the "# Mesh point coordinates" line, that indicates the beginning of the coordinates-point-list
        read(2,'(A)',iostat = ioStatus) linecontent
		if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call remove_Character(linecontent," ")
        do while(trim(linecontent) /= "#Meshpointcoordinates")
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call remove_Character(linecontent," ")
        end do
        
        !Handling and saving the coordinates for each mesh point
        read(2,'(A)',iostat = ioStatus) linecontent
        if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call count_character(numberofHashtags,linecontent,"#")
        do while(numberofHashtags == 0)
            call is_blank_line(linecontent, isBlankLine)
            if(.NOT. isBlankLine) then
                numb_coord = numb_coord + 1
            endif
            
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call count_character(numberofHashtags,linecontent,"#")
        end do
        !Removing the lines, between the coordinate list, and the triangle list.
        read(2,'(A)',iostat = ioStatus) linecontent
		if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call remove_Character(linecontent," ")
        do while(trim(linecontent) /= "#Elements")
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call remove_Character(linecontent," ")
        end do
        !Handling the triangle list. And increasing every point index by 1 (to make the data compatible with the 1 based fortran index system)
        read(2,'(A)',iostat = ioStatus) linecontent
        if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call count_character(numberofHashtags,linecontent,"#")
        
        do while(numberofHashtags == 0)
            call is_blank_line(linecontent, isBlankLine)
            
            if(.NOT. isBlankLine) then
                numb_tri = numb_tri + 1
            endif
            
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call count_character(numberofHashtags,linecontent,"#")
        end do
        
        !Closing the file and cleaning up[...]
        close(2)
    end subroutine number_of_coordinates_triangels
    
    
    !This subroutine reads a comsol file and writes all the data stored in there to a
    !structure_tri structure. This subroutine is a little bit longer, and probably not
    !that intuitive. Therefore I divided the subroutine with code-comments in several
    !parts. Hopefully it got better understandable through that
    !
    !How a comsol file is structured:
    !A comsol file is kind of a text file, that stores Point coordinates and the indexes
    !of 3 Points, that define a triangle. All the points and triangles together building up the mesh.
    !At the beginning of the text file is a header. The inforamtions stored in this header
    !aren't required for our purpose. Therefore we will ignore them.
    !Before the list of coordinates for each point starts, there is one line containing
    !the following: "# Mesh point coordinates"
    !This subroutine searches for this line and starts handling the coordinate list, that
    !comes below. This means, that in every line are 3 real numbers defining the position
    !of a point. The numbers are seperated by a blank.
    !There are no indices, they are defined implicitely by the order of the coordintes.
    !This means, that the first coordinate triple defines point 1(0 in Comsol, but we add 1 to 
    !every point indice to be fortran compatible)
    !After that block there are once again some lines. The information stored by this lines is
    !also irrelevant for us, so we will ignore them too. We can detect these lines, by the presence
    !of an #.
    !The beginning of the list, storing wich points define a triangle is merked by the following line:
    !"# Elements". The subroutine also searches for that line.
    !Then it stores all the triangle data. This data is an integer triple, separated by a blank.
    !
    !Comsol File example for clarification:
    !
    !# Created by COMSOL Multiphysics Fri Jul 29 11:17:54 2022
    !
    !
    !# Major & minor version
    !0 1 
    !1 # number of tags
    !# Tags
    !5 mesh1 
    !1 # number of types
    !# Types
    !3 obj 
    !
    !# --------- Object 0 ----------
    !
    !0 0 1 
    !4 Mesh # class
    !4 # version
    !3 # sdim
    !7360 # number of mesh points
    !0 # lowest mesh point index
    !
    !# Mesh point coordinates <-This is where the header ends
    !9.7545161008064488e-008 2.4999999999999998e-006 4.9039264020161532e-007 
    !4.900868465586856e-008 2.5000000000000002e-006 4.9759235687141682e-007 
    !1.4514233756944688e-007 2.4999999999999998e-006 4.7847016490641793e-007 
    !1.2149008987738374e-007 2.4573349182922882e-006 4.8501562661587051e-007
    ![... Many more node coordinates ...]
    !-1.6844492662404178e-007 -2.4574051918749036e-006 -4.7077203261729685e-007 
    !-1.9134171618254481e-007 -2.5000000000000006e-006 -4.6193976625564335e-007 
    !-9.754516100806413e-008 -2.5000000000000006e-006 -4.9039264020161521e-007 
    !-1.4514233756944653e-007 -2.5000000000000006e-006 -4.7847016490641772e-007 
    !
    !1 # number of element types    <- This # is detected
    !
    !# Type #0
    !
    !3 tri # type name
    !
    !
    !3 # number of nodes per element
    !14592 # number of elements
    !# Elements <- This is where the triangle data starts
    !2 0 3 
    !0 1 7 
    !3 0 7 
    !3 7 6 
    !7 1 8 
    !1 9 8 
    ![... Many more triangles ...]
    !7353 7358 7355 
    !7358 7359 7355 
    !7359 7357 7356 
    !7355 7359 7356 
    !
    !0 # number of geometric entity indices <-This # markes the end
    !
    !All these numbers are stored in matrixes. These matrixes are then stored in the struct given
    !
    !@param:    filename_path A string type variable, that stores the filename and path to the comsol file
    !
    !@param:    struct A struct of the type structure_tri. In this struct all the data is stored in the end
    !           and you can use it from there
    
    

    subroutine read_comsol_file(filename_path, struct)
        implicit none
        !Initialisation, Opening the comsol file [...]
        type(structure_tri) :: struct		
		integer, dimension(:, :), allocatable :: Element_NodeIndices	
        real(dp), dimension(:, :), allocatable :: Node_coordinates
        character(len = *) :: filename_path
        character(len = 200) :: linecontent,hs
        character(len = 100) :: numberSubstring
        integer :: ioStatus,numberofHashtags,n_p,n_t
        integer :: indexOfPoint,point1,point2,point3
        real(dp) :: coordinate1,coordinate2,coordinate3
        logical :: isBlankLine
        
        call number_of_coordinates_triangels(filename_path, n_p,n_t)
        
        allocate(Node_coordinates(3, n_p))
		allocate(struct%points(n_p))
        allocate(Element_NodeIndices(3, n_t))
		allocate(struct%elements(n_t))
        
        open(2,file = trim(filename_path),status = "old",iostat = ioStatus)
        if(ioStatus /= 0) then
            write(*,*) "Error while opening comsol file.\nThe specified comsol file was: "//filename_path//"\n\nExiting"
            close(1)
           stop
        endif
    
        
        !Searching for the "# Mesh point coordinates" line, that indicates the beginning of the coordinates-point-list
        read(2,'(A)',iostat = ioStatus) linecontent
		if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call remove_Character(linecontent," ")
        do while(trim(linecontent) /= "#Meshpointcoordinates")
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call remove_Character(linecontent," ")
        end do
        
        !Handling and saving the coordinates for each mesh point
        indexOfPoint = 1
        read(2,'(A)',iostat = ioStatus) linecontent
        if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call count_character(numberofHashtags,linecontent,"#")
        
        do while(numberofHashtags == 0)
            call is_blank_line(linecontent, isBlankLine)
            
            if(.NOT. isBlankLine) then
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2real(numberSubstring, coordinate1,ioStatus)
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2real(numberSubstring, coordinate2,ioStatus)
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2real(numberSubstring, coordinate3,ioStatus)
                Node_coordinates(1,indexOfPoint) = coordinate1
                Node_coordinates(2,indexOfPoint) = coordinate2
                Node_coordinates(3,indexOfPoint) = coordinate3
            endif
            
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call count_character(numberofHashtags,linecontent,"#")
            indexOfPoint = indexOfPoint + 1
        end do
        !Removing the lines, between the coordinate list, and the triangle list.
        read(2,'(A)',iostat = ioStatus) linecontent
		if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call remove_Character(linecontent," ")
        do while(trim(linecontent) /= "#Elements")
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call remove_Character(linecontent," ")
        end do
        !Handling the triangle list. And increasing every point index by 1 (to make the data compatible with the 1 based fortran index system)
        indexOfPoint = 1
        read(2,'(A)',iostat = ioStatus) linecontent
        if(ioStatus /= 0)then
            write(*,*) "Unexpected EOF while reading comsol file "//filename_path
            stop
        endif
        call count_character(numberofHashtags,linecontent,"#")
        
        do while(numberofHashtags == 0)
            call is_blank_line(linecontent, isBlankLine)
            if(.NOT. isBlankLine) then
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2int(numberSubstring, point1,ioStatus)
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2int(numberSubstring, point2,ioStatus)
                call pop_substring_char(linecontent,numberSubstring," ")
                call str2int(numberSubstring, point3,ioStatus)
                point1 = point1 + 1
                point2 = point2 + 1
                point3 = point3 + 1
                Element_NodeIndices(1, indexOfPoint) = point1
                Element_NodeIndices(1, indexOfPoint) = point2
                Element_NodeIndices(1, indexOfPoint) = point3
            endif
            
            read(2,'(A)',iostat = ioStatus) linecontent
		    if(ioStatus /= 0) exit
            call count_character(numberofHashtags,linecontent,"#")
            indexOfPoint = indexOfPoint + 1
        end do
        
        !Closing the file and cleaning up[...]
        do j = 1, n_p
			struct%points(j)%point = Node_coordinates(:, j)
        end do
        
        do j = 1, n_t
			allocate(struct%elements(j)%corners(3))
			struct%elements(j)%corners = Element_NodeIndices(:, j)
        end do
        
        deallocate(Element_NodeIndices)
		deallocate(Node_coordinates)
        close(2)
        
    end subroutine read_comsol_file
    
end module lib_comsol_reader