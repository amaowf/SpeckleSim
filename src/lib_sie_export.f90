!>This module contains subroutines to save the simultion data in hdf5 files.
!!Therefore the hdf5 library is required.
!!When trying to compile this please make shure, that your include paths are setup properly
!!And that your compiler creates an 64 bit executable. 32 bit won't work with a 64 bit hdf5 library
module lib_sie_export
    use hdf5
    
    implicit none
    
    contains
    
    
    !>This subroutine saves the simulation result in a HDF5 file.
    !!The simultion data file must be an e-field of the type vector_c.
    !!This datatype consists of a 3 dimensional Array of complex numbers.
    !!The e-field is represented by an 2 dimensional array of vector_c.
    !!To store the e-field the subroutine extracts the x,y and z data from
    !!The given e-field and stores these information seperately to the hdf5 file.
    !!The parameter wavelength gives the possiblity to store the wavelength as well.
    !!
    !!@param f_name The name of the resulting hdf5 file, including the ending (e.g.: .h)
    !!               As well as the path to the file
    !!@param e_field The e-field that should be stored. This has to be a 2 dimensional array
    !!               of the type vector_c
    !!@param wavelength This is the wavelength used in the simulation. (real(dp))
    subroutine hdf5_EXPORT(f_name,e_field,wavelength)
    
        character(len=*)::f_name
        type(vector_c),dimension(:,:)::e_field
        real(dp)::wavelength
        real(dp),dimension(1),target::wl
        integer, parameter :: r_k8 = KIND(0.0d0)
        integer:: dimensions_of_field_x,dimensions_of_field_y,i,j,error
        integer(HID_T)   :: sample_type_id, X_set_id,Y_set_id,Z_set_id, X_space_id,Y_space_id,Z_space_id,WL_space_id,WL_set_id, file_id
        integer(HSIZE_T),dimension(2):: dims
        integer(HSIZE_T),dimension(1):: dimsWL = (1) 
        complex(dp), dimension(:,:),allocatable,target::dat_x,dat_y,dat_z
        type(C_PTR) :: cpointer
        integer(8) :: real_size, real_complex_size
        
        wl(1) = wavelength
        
        dimensions_of_field_x =size(e_field,1)
        dimensions_of_field_y =size(e_field,2)
        dims(1) = dimensions_of_field_x
        dims(2) = dimensions_of_field_y
        
        if(dimensions_of_field_x <= 0 .or. dimensions_of_field_y <= 0)then
            write(*,*)"Error: Output E-Field empty"
            return
        endif
        
        real_size = storage_size(1_dp)/8
        real_complex_size = 2* real_size
        
        !Splitting up the e-field in the x,y and z components
        allocate(dat_x(dimensions_of_field_x,dimensions_of_field_y))
        
        do, i = 1,dimensions_of_field_x
            do, j = 1,dimensions_of_field_y
                dat_x(i,j) = e_field(i,j)%vector(1)
            end do
        end do
        
        allocate(dat_y(dimensions_of_field_x,dimensions_of_field_y))
        
        do, i = 1,dimensions_of_field_x
            do, j = 1,dimensions_of_field_y
                dat_y(i,j) = e_field(i,j)%vector(2)
            end do
        end do
        
        allocate(dat_z(dimensions_of_field_x,dimensions_of_field_y))
        
        do, i = 1,dimensions_of_field_x
            do, j = 1,dimensions_of_field_y
                dat_z(i,j) = e_field(i,j)%vector(3)
            end do
        end do
        
        !Here starts the actual hdf5 Library stuff...
        CALL h5open_f(error)
        CALL h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, error)
        if(error /= 0) then
            write(*,*)"Error while creating hdf5 output File..."
            write(*,*)"...Exiting. Nothing was stored!"
            return
        endif
        
        !Create a complex data type
        CALL H5Tcreate_f(H5T_COMPOUND_F, real_complex_size, sample_type_id, error)
        CALL H5Tinsert_f( sample_type_id, "r", 0, h5kind_to_type(dp,H5_REAL_KIND), error)
        CALL H5Tinsert_f( sample_type_id, "i", real_size, h5kind_to_type(dp,H5_REAL_KIND), error)
        if(error /= 0) then
            write(*,*)"Error while creating hdf5 complex data type..."
            write(*,*)"...Exiting. Nothing was stored!"
            return
        endif
        !create dataspaces
        CALL h5screate_simple_f(2, dims, X_space_id, error)
        CALL h5screate_simple_f(2, dims, Y_space_id, error)
        CALL h5screate_simple_f(2, dims, Z_space_id, error)
        CALL h5screate_simple_f(1, dimsWL, WL_space_id, error)
        if(error /= 0) then
            write(*,*)"Error while creating hdf5 dataspaces"
            write(*,*)"...Exiting. Nothing was stored!"
            return
        endif
        ! create datasets
        CALL H5Dcreate_f(file_id, "data_x",  sample_type_id,X_space_id, X_set_id, error)
        CALL H5Dcreate_f(file_id, "data_y",  sample_type_id,Y_space_id,Y_set_id, error)
        CALL H5Dcreate_f(file_id, "data_z",  sample_type_id,Z_space_id,Z_set_id, error)
        CALL H5Dcreate_f(file_id, "wavelength",  H5T_NATIVE_DOUBLE,WL_space_id,WL_set_id, error)  
        if(error /= 0) then
            write(*,*)"Error while creating hdf5 datasets"
            write(*,*)"...Exiting. Nothing was stored!"
            return
        endif
        !Write data to the datasets
        cpointer = C_LOC(dat_x(1,1))
        CALL H5Dwrite_f(X_set_id, sample_type_id, cpointer ,error)
        cpointer = C_LOC(dat_y(1,1))
        CALL H5Dwrite_f(Y_set_id, sample_type_id, cpointer ,error)
        cpointer = C_LOC(dat_z(1,1))
        CALL H5Dwrite_f(Z_set_id, sample_type_id, cpointer ,error)
        cpointer = C_LOC(wl)
        CALL H5Dwrite_f(WL_set_id, H5T_NATIVE_DOUBLE, cpointer ,error)
        if(error /= 0) then
            write(*,*)"Error while writing to hdf5 file"
        endif
        !close everything...
        CALL h5dclose_f(X_set_id, error)
        CALL h5sclose_f(X_space_id, error)
        CALL h5dclose_f(Y_set_id, error)
        CALL h5sclose_f(Y_space_id, error)
        CALL h5dclose_f(Z_set_id, error)
        CALL h5sclose_f(Z_space_id, error)
        CALL h5dclose_f(WL_set_id, error)
        CALL h5sclose_f(WL_space_id, error)
        CALL H5Tclose_f(sample_type_id, error)
        CALL h5fclose_f(file_id, error)
        CALL h5close_f(error)
    end subroutine hdf5_EXPORT
    
    
    
end module lib_sie_export