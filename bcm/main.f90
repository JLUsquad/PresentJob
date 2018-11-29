program main  

  implicit none
  integer(kind=4)             :: k, a, b, c, nele, natom, model, nconf,ls, nn, i, j
  integer(kind=4),allocatable :: ele_num(:)
  parameter(nconf=5755, ls=1,nn=56)
  real(kind=8)                :: qall(11,1,1,nconf),E(nconf), Distance(nconf),cfrd(nconf),fpd(nconf), xxd(nconf), E0, temp
  real(kind=8),allocatable    :: q0(:,:,:), fp0(:,:),xx0(:,:,:),q1(:,:,:), fp1(:,:), xx1(:,:,:),cfr0(:),cfr1(:)
  logical                     :: alive
  character(100)              :: filename,tempstr
  character(50),dimension(2857) ::filenameall
  character(3)                :: c1, c2, c3
  parameter(model=1)

!  nele=0
!  natom=0
!  filename='1-ICSD-431.vasp'
!  filename=trim(filename)
!  print* ,len(filename)
!  call getpara(filename,nele,natom)
!  write(*,*)'nele=',nele,'natom=',natom
!  allocate(q0(11,nele,nele))
!  allocate(ele_num(nele))
!  call fingerprint(filename,nele,natom,q0,model)
!  deallocate(ele_num)

!do i = 1, 2
!do j = 1, 2
!  write(*,*)q0(:,i,j),'(',i,j,')'
!enddo
!enddo

open(50,file = 'qresult.txt')
open(30,file = 'Dfrom0andE.txt')
k = 0

do a = 1, 20
  write(c1,'(i3)')a
  c1 = adjustl(c1)
  do b = 1, 10
    write(c2,'(i3)')b
    c2 = adjustl(c2)
    do c = 1, 100
      write(c3,'(i3)')c
      c3 = adjustl(c3)
      filename = trim('Calypso_opt/calypso'//trim(c1)//'/a_'//trim(c2)//'/'//trim(c3)//'/POSCAR')
      inquire(file=filename,exist=alive)
      if(alive)then

        k = k + 1
        print* ,k
        nele=0
        natom=0

        call getpara(filename,nele,natom)

        write(*,*)"nele,natom",nele,natom

        allocate(q1(11,nele,nele))

        q1=0
    
        allocate(q0(11,nele,nele))
        q0=0
  

        call fingerprint(filename,nele,natom,E(k),q1,model)

!        print* ,E(k),q1

  !      do i=1,11
 !       qall(i,1,1,k)=q1(i,1,1)
   !     enddo
        write(50,*)k
        write(50,*)E(k)
        
        
        write(50,*)q1
      
        call dist_BCM(Distance(k),q0,q1,nele)      
        write(30,*)Distance(k),'  ',E(k)
        deallocate(q1)
        deallocate(q0)
    else
    endif 

enddo
enddo
enddo

close(50)
close(30)

!open(10,file = 'filenamelist.txt')
!open(30,file = 'result01.txt')
!do a = 1, 2857
!        read(10,*)filename
!        k = k + 1
!        filenameall(k)=trim(filename)
!        print* ,'k=',k
!
!        nele=0
!        natom=0
!        call getpara(filename,nele,natom)
!        allocate(q1(11,nele,nele))
!!        allocate(ele_num(nele))
!        q1=0
!        call fingerprint(filename,nele,natom,q1,model)
!        write(30,*)filename
!        do i = 1,2
!        do j = 1,2
!        write(30,*)q1(:,i,j)
!        enddo
!        enddo
!        call dist_BCM(Distance(k),q0,q1,nele)
 !       print* ,'Distance = ',Distance(k)
!        write(30,*)'from 344 D=',Distance(k)
!        deallocate(q1)
!        deallocate(ele_num)
!end do
!close(30)
!close(10)

!else
!endif

!do i= 1,nconf
!    do j= 1,nconf-1
!      if(Distance(j) .gt. Distance(j+1))then
!      temp=Distance(j+1)
!      Distance(j+1)=Distance(j)
!      Distance(j)=temp
!      do a=1,2
!        do b=1,2
!          do c=1,11
!            temp2=qall(c,a,b,j+1)
!            qall(c,a,b,j+1)=qall(c,a,b,j)
!            qall(c,a,b,j)=temp2
!          enddo
!        enddo
!      enddo
!      tempstr=filenameall(j+1)
!      filenameall(j+1)=filenameall(j)
!      filenameall(j)=tempstr
!      else
!      end if
!    end do
!  end do

!open(50,file = 'result.txt')

!do i = 1, nconf

!write(50,*)Distance(i),filenameall(i)

!enddo

!close(50)


end program main




