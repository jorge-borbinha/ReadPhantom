!*********************************************************************
!**********************************************************************
!*                     READ PHANTOM PROGRAM                           *
!*                                                                    *
!*     This program reads and ASCII/binary voxel phantom file, where  *
!*   the organ numbers are listed countinously in x, y and z, and     *
!*   generates the files necessary for visualization and simulation   *
!*   with PENELOPE/penEasy.                                           *
!*                                                                    *
!*  Created by: Jorge Cebola Borbinha, December 2017                  *
!*    @website: jorge-borbinha.github.io                              *
!*                                                                    *
!*  Preferred citations:                                              *
!*     (1) Borbinha J. Organ Dose Estimates in Thorax CT: Voxel       *
!*   Phantom Organ Matching With Individual Patient Anatomy. MSc      *
!*   Thesis, 2017. Available from: http://hdl.handle.net/10362/       *
!*     (2) Borbinha J, Di Maria S, Madeira P, Belchior A, Baptista M, *
!*   Vaz P. Increasing organ dose accuracy through voxel phantom      *
!*   organ matching with individual patient anatomy. Radiat Phys      *
!*   Chem.2019 Jun;159:35–46. DOI: 10.1016/j.radphyschem.2019.02.014  *
!*                                                                    *
!* Copyright (C) 2017 Jorge Cebola Borbinha,                          *
!* Faculty of Sciences and Technology, NOVA University of Lisbon      *
!* This program is free software: you can redistribute it and/or      *
!* modify it under the terms of the MIT License.                      *
!* For more information, see:                                         *
!* https://github.com/jorge-borbinha/ReadPhantom?tab=MIT-1-ov-file    *
!*                                                                    *
!**********************************************************************
!**********************************************************************

        program readPhantom
        implicit none
        integer::x,y,z,counter,file_type,nvoxx,nvoxy,nvoxz,nvoxtot,material,max_material, n_materials, n_organs,organ,max_organ
        integer ::i,j,unit1,k
        real :: voxResx,voxResy,voxResz,lenx,leny,lenz,density,max_density
        integer,allocatable, dimension(:) :: mat,organ_vector,organlist_organ,organlist_mat
        real, allocatable, dimension(:) :: den, organlist_den
        character(20) :: phantom_file
        character(1) :: choice1, choice2, B
        max_density=0
        counter=0

        write(*,*)' Please insert the name of the phantom file (max 20 characters).'
        read (*,*) phantom_file
        write(*,*)' Is your file in binary or in ASCCI? (type 0 for bin or 1 for ASCII)'
        read (*,*) file_type
        do while (file_type /= 0 .and. file_type /= 1)
                write(*,*)' The number you typed is not valid. Type 0 for bin or 1 for ASCII.'
                read(*,*) file_type
        end do
        write(*,*) ' Type in the number of voxels of the phantom in x,y,z.'
        read (*,*) nvoxx,nvoxy,nvoxz
        write(*,*) ' Type in the voxel resolution in x,y,z /cm.'
        read (*,*) voxResx,voxResy,voxResz
        write(*,*) ' Type in the number of different materials in the phantom.'
        read (*,*) n_materials
        write(*,*) ' Type in the number of organ IDs - number of lines in organlist.dat.'
        read (*,*) n_organs
        nvoxtot = nvoxx*nvoxy*nvoxz
        allocate(mat(nvoxtot))
        allocate(den(nvoxtot))
        allocate(organ_vector(nvoxtot))
        allocate(organlist_organ(n_organs))
        allocate(organlist_mat(n_organs))
        allocate(organlist_den(n_organs))
        lenx=voxResx*nvoxx
        leny=voxResy*nvoxy
        lenz=voxResz*nvoxz
        write(*,*)'  '
        write(*,*)'>>> Characteristics of your phantom file:'
        30 format(' > Number of voxels in x,y,z:',1x,i4,1x,i4,1x,i4)
        write(*,30) nvoxx, nvoxy, nvoxz
        31 format(' > Voxel resolution in x,y,z /cm:' ,1x,f10.7,2x,f10.7,2x,f10.7)
        write(*,31) voxResx,voxResy,voxResz
        32 format(' > Phantom size (approximate value) in x,y,z /cm: ',1x,f12.7,2x,f12.7,2x,f12.7)
        write(*,32) lenx,leny,lenz
        33 format(' > Total number of voxels:',1x,i9)
        write(*,33) nvoxtot
        write(*,*)
        write(*,*) ' Please check if these values are correct. Do you wish to continue? (y/n)'
        do while (choice1/='n'.and.choice1/='N'.and.choice1/='y'.and.choice1/='Y')
                write(*,*)' Please type y or n.'
                read (*,*)choice1
        end do
        if ((choice1.eq.'n').or.(choice1.eq.'N')) then
                write(*,*) ' The script will stop at your request.'
                stop              
        end if
        write(*,*)'  ' 
        write(*,*)' Reading the phantom file...'

        open(unit=11, file='organlist.dat', status='old')
        open(unit=12, file='organlistAsRead.dat')
        34 format(2X,I4,1X,I4,2X,E13.5)
        do i=1,n_organs
                read (11,*) organ,material,density
                organlist_organ(i) = organ
                organlist_mat(i) = material
                organlist_den(i) = density
                write(12,34) organ,material,density
        end do
        close(11)
        close(12)
        40 format('  Finished reading the phantom file.',1x,i9,1x,'voxels were read.')

!       for binary file
        if (file_type == 0) then
                open (1, file=phantom_file, status='old', form='unformatted', access='stream')
                write(*,*)' Opened phantom file in binary...'
                do z=1,nvoxz
                        do y=1,nvoxy
                                do x=1,nvoxx
                                        counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                                        read(1) B
                                        organ_vector(counter)= ICHAR(B)
                                        do j=1, n_organs
                                        if (organ_vector(counter)==organlist_organ(j)) then
                                                mat(counter)=organlist_mat(j)
                                                den(counter)=organlist_den(j)
                                        end if
                                        end do
                                end do      
                        end do
                        call read_progress(z)
                end do
                close(1)
        write(*,40) nvoxtot
        end if

!       for ASCII file
        if (file_type == 1) then
                open (2, file=phantom_file, status='old', action='read')
                write(*,*)' Opened phantom file in ASCII...'
                read (2,*) organ_vector
                do z=1,nvoxz
                        do y=1,nvoxy
                                do x=1,nvoxx
                                        counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                                        do j=1, n_organs
                                        if (organ_vector(counter)==organlist_organ(j)) then
                                                mat(counter)=organlist_mat(j)
                                                den(counter)=organlist_den(j)
                                        end if
                                        end do
                                end do
                        end do
                        call read_progress(z)
                end do
                close(2)
                write(*,40) nvoxtot
                write(*,*)'  ' 

        end if

!       Check if materials and densities are okay, max values. 
!       (just for verification)
!
        write(*,*)' Reading the phantom data...'
        if (n_materials==1) then
                write(*,*)' The number of materials must be larger than 1.'
                write(*,*)' The script will stop.'
                stop
        end if
        do k=1,nvoxtot
                if (den(k) < 0) then
                        write(*,*)' You have negative densities in your organlist.dat file.'
                        write(*,*)' This cannot happen.'
                        write(*,*)' Do you wish to change the densities in organlist.dat to their absolute value?'
                        do while (choice2/='n'.and.choice2/='N'.and.choice2/='y'.and.choice2/='Y')
                                write(*,*)' Please type y or n.'
                                read (*,*)choice2
                        end do
                        if ((choice2.eq.'n').or.(choice2.eq.'N')) then
                                write(*,*) ' The script will stop.'
                                stop
                        end if
                        if ((choice2.eq.'y').or.(choice2.eq.'Y')) then
                                den=abs(den)
                                write(*,*) 'All densities were changed to their absolute value.'
                                write(*,*) 'The script will continue.'
                                write(*,*) '  '
                        end if
                end if
        end do
        max_material=maxval(mat)
        max_density=maxval(den)
        max_organ=maxval(organ_vector)
        41  format('  > Number of materials:',1x,i4)
        write(*,41) n_materials
        42  format('  > Maximum value of material ID:',1x,i4)
        write(*,42) max_material
        43  format('  > Maximum value of organ (tag) ID:',1x,i4)
        write(*,43) max_organ
        44  format('  > Maximum density: ',2x,e13.5)
        write(*,44) max_density
        write(*,*)' The phantom data was read.'
        write(*,*)'  '
!
!       Write VOX file
!
        open(4, file='phantom.vox', status='replace')
        write(*,*)' Creating the VOX file...'
!       VOX file header
        write(4,*)'[SECTION VOXELS HEADER v.2008-04-13]'
        52 format(1x,3I4)
        53 format(1x,f10.7,1X,f10.7,1X,f10.7)
        write(4,52) nvoxx,nvoxy,nvoxz
        write(4,53) voxResx,voxResy,voxResz
        write(4,*)' 1'
        write(4,*)' 2'
        write(4,*)' 0'
        write(4,*)'[END OF VXH SECTION]'
        do z=1,nvoxz
                do y=1,nvoxy
                        do x=1,nvoxx
                                counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                                write(4,'(1X,I3,2X,F7.4)')mat(counter), den(counter)
                        end do      
                end do
                call read_progress(z)
        end do
        write(*,*)' The phantom.vox file was created.'
        write(*,*)'  '
        close(4)

!
!
!       Create ct-den-mat.dat files
!
        write(*,*)' Creating the ct-den-mat.dat files...'        
        open(7, file = 'ct-den-matXY.dat', status='replace')
        open(8, file = 'ct-den-matXZ.dat', status='replace')
        open(9, file = 'ct-den-matYZ.dat', status='replace')
!       File Headers
        60  format(1x,'#  CT enclosure limits:  XL = ',1p,e13.6,' cm,  XU = ',e13.6,' cm')
        61  format(1x,'#',20x,'YL = ',1p,e13.6,' cm,  YU = ',e13.6,' cm')
        62  format(1x,'#',20x,'ZL = ',1p,e13.6,' cm,  ZU = ',e13.6,' cm')
        64  format (1x,'# Numbers of voxels:   NVX = ',i3,', NVY = ',i3,', NVZ = ',i3)
        do unit1=7,9
                write(unit1,*)'#  CT structure (GNUPLOT format).'
                write(unit1,60) 0.0D0,lenx
                write(unit1,61) 0.0D0,leny
                write(unit1,62) 0.0D0,lenz
                write(unit1,64) nvoxx,nvoxy,nvoxz
                write(unit1,*)'#'
                write(unit1,*) 
                write(unit1,*)'#  columns 1 to 3: bin indices IX, IY and IZ'
                write(unit1,*)'#  4th column: density (g/cm**3).'
                write(unit1,*)'#  5th column: material.'
                write(unit1,*)'#  CT structure (GNUPLOT format).'
        end do
!       Create ct-den-matXY.dat (should have 14 303 159 lines (11 +14 255 124 +137*348 +348))
        65 format (1x,3i4,e13.5,1x,i4,3x,i4)
        do z=1,nvoxz
                do y=1,nvoxy
                        do x=1,nvoxx
                                counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                                write(7,65)x,y,z,den(counter),mat(counter),organ_vector(counter)
                        end do
                        write (7,'(''  '')')
                end do
                write(7,'(''  '')')
                call read_progress(z)
        end do  
        write (*,*) ' File ct-den-matXY.dat created.'
        write (*,*) ' '		


            
!       Create ct-den-matXZ.dat
        do y=1,nvoxy
                do x=1,nvoxx
                       do z=1,nvoxz
                               counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                               write(8,65)x,y,z,den(counter),mat(counter),organ_vector(counter)   
                        end do
                        write (8,'(''  '')')
                end do
                write(8,'(''  '')')
                call read_progress(y)
        end do   
        write (*,*) ' File ct-den-matXZ.dat created.'
        write (*,*) ' '			
!       Create ct-den-matYZ.dat
        do x=1,nvoxx
                do z=1,nvoxz
                        do y=1,nvoxy
                                counter=x+nvoxx*((y-1)+nvoxy*(z-1))
                                write(9,65)x,y,z,den(counter),mat(counter),organ_vector(counter)   
                        end do
                        write (9,'(''  '')')
                end do
                write(9,'(''  '')')
                call read_progress(x)
        end do   
        write (*,*) ' File ct-den-matYZ.dat created.'
        write (*,*) ' '			
        close(7)
        close(8)
        close(9)

        write (*,*) ' All Files for visualization and simulation with this phantom have been created.' 

        end program readPhantom

        subroutine read_progress(z)
                integer z
                if (z==50) then
                        write(*,*)' > Number of slices read:',z
                else if (z==100) then
                        write(*,*)' > Number of slices read:',z
                else if (z==150) then
                        write(*,*)' > Number of slices read:',z
                else if (z==200) then
                        write(*,*)' > Number of slices read:',z
                else if (z==250) then
                        write(*,*)' > Number of slices read:',z
                else if (z==300) then
                        write(*,*)' > Number of slices read:',z
                else if (z==350) then
                        write(*,*)' > Number of slices read:',z
				else if (z==400) then
                        write(*,*)' > Number of slices read:',z
                end if
        return
        end


