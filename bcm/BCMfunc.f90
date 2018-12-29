
subroutine BCM(nele,natom,tau0,lat_matrix,ele_num,q0)

  implicit none
  integer(kind=4)             :: i, j, k, ix, iy, iz, i1, i2, i3, k1
  integer(kind=4)             :: l, m, n, j1, j2, j3
  integer(kind=4)             :: nele, natom
  integer(kind=4)             :: ele_num(nele), super_ele_num(nele), counter(nele,nele)
  real(kind=8)                :: lat_matrix(3,3), pi
  real(kind=8)                :: eps,eps1,tmp,temp(3)
  real(kind=8)                :: rms,x,y,z,r,theta,phi,maxr,rp(400*natom),rpmin(natom,nele)
  real(kind=8)                :: pos (3,400*natom,nele), mo(nele,nele), tau1(3,natom), tau(3,natom), tau0(natom,3)
  real(kind=8)                :: dis(nele,nele), q(11,nele,nele), q0(6,nele,nele)
  complex                     :: Yx(nele,nele),sumY(nele,nele)

!f2py intent(in) tau0
!f2py intent(in) lat_matrix
!f2py intent(in) ele_num
!f2py intent(out) q0
!f2py depend(ele_num) nele
!f2py depend(tau0) natom

  do i = 1,3
        do j=1,natom
                tau(i,j)=tau0(j,i)
        enddo
  enddo


  pi = 3.1415926

  eps = 1e-1
  eps1 = 1e-4 
  Yx  = (0.0,0.0)
  sumY = (0.0,0.0)
  mo = 0.0

  ix = 1
  iy = 1
  iz = 1
  super_ele_num = 0
  k=0
 
  do i=1,nele
     do j=1,ele_num(i)
        k=k+1
        do i1 = -ix, ix
           do i2 = -iy, iy
              do i3 = -iz, iz
                 super_ele_num(i) = super_ele_num(i) + 1
 
                 temp(1) = tau(1,k) + i1
                 temp(2) = tau(2,k) + i2
                 temp(3) = tau(3,k) + i3
                 pos(1,super_ele_num(i),i)=temp(1)*lat_matrix(1,1) + temp(2)*lat_matrix(2,1) +temp(3)*lat_matrix(3,1) 
                 pos(2,super_ele_num(i),i)=temp(1)*lat_matrix(1,2) + temp(2)*lat_matrix(2,2) +temp(3)*lat_matrix(3,2)
                 pos(3,super_ele_num(i),i)=temp(1)*lat_matrix(1,3) + temp(2)*lat_matrix(2,3) +temp(3)*lat_matrix(3,3)
              end do
           end do
        end do
     end do
  end do
  
  k=0
    do i=1,nele
       do j=1,ele_num(i)
       k=k+1
        tau1(1,k)=tau(1,k)*lat_matrix(1,1) + tau(2,k)*lat_matrix(2,1) +tau(3,k)*lat_matrix(3,1)
        tau1(2,k)=tau(1,k)*lat_matrix(1,2) + tau(2,k)*lat_matrix(2,2) +tau(3,k)*lat_matrix(3,2)
        tau1(3,k)=tau(1,k)*lat_matrix(1,3) + tau(2,k)*lat_matrix(2,3) +tau(3,k)*lat_matrix(3,3)
       end do
    end do
  dis=100.0

k=0
!-----------------
rpmin(:,:)=10000.0
rp = 0
       do i1=1, nele
           do j1 = 1,ele_num(i1)
              k = k + 1
              k1 = 0
              rp = 0
              do i2=1,nele
                 do j2=1,super_ele_num(i2)
                    k1 = k1 + 1
                    x = tau1(1,k) - pos(1,j2,i2)
                    y = tau1(2,k) - pos(2,j2,i2)
                    z = tau1(3,k) - pos(3,j2,i2)
                    rp(k1)=sqrt(x*x+y*y+z*z)
                    if (rp(k1) .gt. eps)then
                   
                    if (rpmin(k,i2) .gt. rp(k1))then
                    rpmin(k,i2) = rp(k1) 
                    else
                    endif
                    else
                    endif
                 enddo
              enddo
           enddo
       enddo
do i=1,natom
rpmin(i,:)=rpmin(i,:)+1e-2
enddo

  do l = 0, 10
     mo = 0.0
     do m = -l, l
        sumY = (0.0, 0.0)
        counter = 0
       k = 0
        do i1=1, nele
           do j1 = 1,ele_num(i1)
              k = k + 1
              do i2=1,nele
                 do j2=1,super_ele_num(i2)
                    x = tau1(1,k) - pos(1,j2,i2)
                    y = tau1(2,k) - pos(2,j2,i2)
                    z = tau1(3,k) - pos(3,j2,i2)
                    r = sqrt(x*x+y*y+z*z)
                    if (r>eps .and. r<rpmin(k,i2) ) then
                       call car2sphe(x,y,z,r,theta,phi)
                       counter(i1,i2) = counter(i1,i2) + 1
                       call spherical(l,m,theta,phi,Yx(i1,i2))
                       sumY(i1,i2) = sumY(i1,i2) + Yx(i1,i2)
                    end if
                 end do
              end do
           end do
        end do
        
        do i1=1,nele
           do i2=1,nele
              if (counter(i1,i2)==0) then
                  sumY(i1,i2) = 0.0
              else
                  sumY(i1,i2) = sumY(i1,i2)/real(counter(i1,i2))
              end if              
              mo(i1,i2)=mo(i1,i2)+cabs(sumY(i1,i2))**2
           end do
        end do
     end do
     do i1=1,nele
        do i2=1,nele
           q(l+1,i1,i2) = sqrt ( mo(i1,i2) * 4 * pi / ( 2 * l + 1 ) )
        end do
     end do
  end do


  do i = 1,6
    do j = 1,nele
      do k = 1,nele
         q0(i,j,k) = q(2*i-1,j,k)
      enddo
    enddo
  enddo




  end subroutine BCM

  subroutine dist_BCM(Distance,q1,q2,nele)
  implicit none
  integer(kind=4) ::l,i,j,nele,Ntype
  real(kind=8)    ::q1(11,nele,nele),q2(11,nele,nele)
  real(kind=8)    ::sumdif,Distance
  
  
  
  sumdif = 0
  do l = 1, 11
     do i = 1, nele
        do j = i, nele
           sumdif = sumdif + (q1(l,i,j)-q2(l,i,j))**2
        end do
     end do
  end do

  Ntype = nele*(nele-1)/2 + nele
  Distance = sqrt(sumdif/Ntype)
end subroutine dist_BCM

