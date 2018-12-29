subroutine car2sphe(x,y,z,r,theta,psi)

  implicit none
  real(kind=8)      :: x,y,z,r,theta,psi,pi,tol

  pi=3.1415926
  tol = 0.01

  r=sqrt(x*x+y*y+z*z)
  
  if (abs(x)<tol .and. abs(y)<tol ) then 
      if (abs(z)<tol) then
          theta = 0.0
      else if (z<0.0) then
         theta=pi
      else
         theta=0.0
      endif
  else
     theta=acos(z/r)
  endif
  if (abs(x)<tol .and. abs(y)<tol ) then
      psi = 0.0
  else if (abs(x)<tol .and. y > 0.0 ) then
      psi = pi/2.0
  else if (abs(x)<tol .and. y < 0.0 ) then
      psi = 3*pi/2.0
  else if (abs(y)<tol .and. x > 0.0 ) then
      psi = 0.0
  else if (abs(y)<tol .and. x < 0.0 ) then
      psi = pi
  else if (x>0 .and. y>0) then
      psi = atan(y/x)
  else if (x<0 .and. y>0) then
      psi = atan(y/x) + pi
  else if (x<0 .and. y<0) then
      psi = atan(y/x) + pi
  else if (x>0 .and. y<0) then
      psi = atan(y/x) + 2*pi
  end if

end subroutine car2sphe
