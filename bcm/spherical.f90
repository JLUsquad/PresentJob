subroutine spherical(l,m,angle1,angle2,Y)
  ! angle1 is theta  0-pi
  ! angle2 is phi    0-2pi
implicit none

integer(kind=4) l,m
real(kind=8) angle1,angle2,sina,cosa,sinb,cosb,pi,c,d,topi
complex i,Y

i=(0,1)
sina=sin(angle1)
cosa=cos(angle1)
pi=3.1415926
topi=2.0*pi
if (l==0) then
  if(m==0)  Y=0.5*sqrt(1.0/pi)
elseif (l==1) then
  if( m==-1) then
      Y=0.5*sqrt(3.0/(topi))*(exp(-i*angle2))*sin(angle1)
  elseif (m==0)  then
      Y=0.5*sqrt(3.0/pi)*cos(angle1)
  elseif (m==1) then
      Y=-0.5*sqrt(3.0/(topi))*exp(i*angle2)*sin(angle1)
  end if
else if (l==2) then
  if (m==-2) then
     Y= 0.25*sqrt(15.0/(topi))*exp(-i*2.0*angle2)*sin(angle1)**2
  elseif (m==-1) then
     Y=0.5*sqrt(15.0/(topi))*exp(-i*angle2)*sin(angle1)*cos(angle1) 
  elseif (m==0) then
     Y=0.25*sqrt(5.0/pi)*(3.0*(cos(angle1))**2-1.0)
  elseif (m==1) then
      Y=-0.5*sqrt(15.0/(topi))*exp(i*angle2)*sin(angle1)*cos(angle1) 
  elseif (m==2) then
      Y=0.25*sqrt(15.0/topi)*exp(i*2.0*angle2)*sina**2
  end if
else if (l==3) then
  if (m==-3) then
     Y=1.0/8.0*sqrt(35.0/pi)*exp(-i*3.0*angle2)*(sin(angle1))**3
  elseif (m==-2) then
     Y=1.0/4.0*sqrt(105.0/topi)*exp(-i*2.0*angle2)*(sin(angle1))**2*cos(angle1)
  elseif (m==-1) then
     Y=1.0/8.0*sqrt(21.0/pi)*exp(-i*angle2)*sin(angle1)*(5*(cos(angle1))**2-1)
  elseif (m==0) then
     Y=1.0/4.0*sqrt(7.0/pi)*(5.0*(cos(angle1))**3-3.0*cos(angle1))
  elseif (m==1) then
     Y=-1.0/8.0*sqrt(21.0/pi)*exp(i*angle2)*sin(angle1)*(5.0*(cos(angle1))**2-1)
  elseif (m==2) then
     Y=1.0/4.0*sqrt(105.0/(topi))*exp(2.0*i*angle2)*(sin(angle1))**2*cos(angle1)
  elseif (m==3) then
     Y=-1.0/8.0*sqrt(35.0/pi)*exp(3*i*angle2)*(sin(angle1))**3
  end if
else if(l==4) then
  if (m==-4) then
      Y=3.0/16.0*sqrt(35.0/topi)*exp(-4.0*i*angle2)*(sin(angle1))**4
  elseif (m==-3) then
      Y=3.0/8.0*sqrt(35.0/pi)*exp(-3*i*angle2)*(sin(angle1))**3*cos(angle1)
  elseif (m==-2) then
      Y=3.0/8.0*sqrt(5.0/topi)*exp(-2.0*i*angle2)*(sin(angle1))**2*(7.0*(cos(angle1))**2-1)
  elseif (m==-1) then
      Y=3.0/8.0*sqrt(5.0/pi)*exp(-i*angle2)*sin(angle1)*(7.0*(cos(angle1))**3-3.0*cos(angle1))
  elseif (m==0) then
      Y=3.0/16.0*sqrt(1.0/pi)*(35.0*cosa**4-30.0*cosa**2+3.0)
  elseif (m==1)then
      Y=-3.0/8.0*sqrt(5.0/pi)*exp(i*angle2)*sina*(7.0*cosa**3-3.0*cosa)
  elseif (m==2) then
      Y=3.0/8.0*sqrt(5.0/topi)*exp(2.0*i*angle2)*sina**2*(7.0*cosa**2-1.0)
  elseif (m==3) then
      Y=-3.0/8.0*sqrt(35.0/pi)*exp(3.0*i*angle2)*sina**3*cosa
  elseif (m==4) then
      Y=3.0/16.0*sqrt(35.0/topi)*exp(4.0*i*angle2)*sina**4
  endif
elseif (l==5) then
  if (m==-5) then
      Y=3.0/32.0*sqrt(77.0/pi)*exp(-5.0*i*angle2)*sina**5
  elseif (m==-4) then
     Y=3.0/16.0*sqrt(385.0/topi)*exp(-4.0*i*angle2)*sina**4*cosa
 elseif (m==-3) then
     Y=1.0/32.0*sqrt(385.0/pi)*exp(-3.0*i*angle2)*sina**3*(9.0*cosa**2-1.0)
 elseif (m==-2)then
     Y=1.0/8.0*sqrt(1155.0/topi)*exp(-2.0*i*angle2)*sina**2*(3.0*cosa**3-cosa)
 elseif (m==-1) then
     Y=1.0/16.0*sqrt(165.0/topi)*exp(-i*angle2)*sina*(21.0*cosa**4-14.0*cosa**2+1.0)
 elseif (m==0) then
     Y=1.0/16.0*sqrt(11.0/pi)*(63.0*cosa**5-70.0*cosa**3+15.0*cosa)
 elseif (m==1) then
     Y=-1.0/16.0*sqrt(165.0/topi)*exp(i*angle2)*sina*(21.0*cosa**4-14.0*cosa**2+1.0)
 elseif (m==2) then
     Y=1.0/8.0*sqrt(1155.0/topi)*exp(i*angle2)*sina**2*(3.0*cosa**3-cosa)
 elseif (m==3) then
     Y=-1.0/32.0*sqrt(385.0/pi)*exp(i*3.0*angle2)*sina**3*(9.0*cosa**2-1.0)
 elseif (m==4) then
     Y=3.0/16.0*sqrt(385.0/topi)*exp(i*4.0*angle2)*sina**4*cosa
 elseif (m==5) then
     Y=-3.0/32.0*sqrt(77.0/pi)*exp(5.0*i*angle2)*sina**5
 end if
 elseif (l==6)then
   if (m==-6)then
       Y=1.0/64.0*sqrt(3003.0/pi)*exp(-6.0*i*angle2)*sina**6
   elseif (m==-5) then
       Y=3.0/32.0*sqrt(1001.0/pi)*exp(-5.0*i*angle2)*sina**5*cosa
   elseif (m==-4) then
       Y=3.0/32.0*sqrt(91.0/topi)*exp(-4.0*i*angle2)*sina**4*(11.0*cosa**2-1.0)
   elseif (m==-3) then
       Y=1.0/32.0*sqrt(1365.0/pi)*exp(-3.0*i*angle2)*sina**3*(11.0*cosa**3-3.0*cosa)
   elseif (m==-2) then
       Y=1.0/64.0*sqrt(1365.0/pi)*exp(-2.0*i*angle2)*sina**2*(33.0*cosa**4-18.0*cosa**2+1.0)
   elseif (m==-1) then
       Y=1.0/16.0*sqrt(273.0/topi)*exp(-i*angle2)*sina*(33.0*cosa**5-30.0*cosa**3+5.0*cosa) 
   elseif (m==0) then
       Y=1.0/32.0*sqrt(13.0/pi)*(231.0*cosa**6-315.0*cosa**4+105.0*cosa**2-5.0)
   elseif (m==1) then
       Y=-1.0/16.0*sqrt(273.0/topi)*exp(i*angle2)*sina*(33.0*cosa**5-30.0*cosa**3+5.0*cosa)
   elseif (m==2) then
       Y=1.0/64.0*sqrt(1365.0/pi)*exp(i*2.0*angle2)*sina**2*(33.0*cosa**4-18.0*cosa**2+1.0)
   elseif (m==3) then
       Y=-1.0/32.0*sqrt(1365.0/pi)*exp(3.0*i*angle2)*sina**3*(11.0*cosa**3-3.0*cosa)
   elseif (m==4)then
       Y=3.0/32.0*sqrt(91.0/topi)*exp(4.0*i*angle2)*sina**4*(11.0*cosa**2-1.0)
   elseif (m==5) then
       Y=-3.0/32.0*sqrt(1001.0/pi)*exp(5.0*i*angle2)*sina**5*cosa
   elseif (m==6) then
       Y=1.0/64.0*sqrt(3003.0/pi)*exp(6.0*i*angle2)*sina**6
   endif
 elseif (l==7) then
    if (m==-7) then
      Y=3.0/64.0*sqrt(715.0/topi)*exp(-7.0*i*angle2)*sina**7
    elseif (m==-6) then
      Y=3.0/64.0*sqrt(5005.0/pi)*exp(-6.0*i*angle2)*sina**6*cosa
    elseif (m==-5) then
      Y=3.0/64.0*sqrt(385.0/topi)*exp(-5.0*i*angle2)*sina**5*(13.0*cosa**2-1.0)
    elseif (m==-4) then
      Y=3.0/32.0*sqrt(385.0/topi)*exp(-4.0*i*angle2)*sina**4*(13.0*cosa**3-3.0*cosa)
    elseif (m==-3) then
      Y=3.0/64.0*sqrt(35.0/topi)*exp(-3.0*i*angle2)*sina**3*(143.0*cosa**4-66.0*cosa**2+3.0)
    elseif (m==-2) then
      Y=3.0/64.0*sqrt(35.0/pi)*exp(-2.0*i*angle2)*sina**2*(143.0*cosa**5-110.0*cosa**3+15.0*cosa)
    elseif (m==-1) then
      Y=1.0/64.0*sqrt(105.0/topi)*exp(-i*angle2)*sina*(429.0*cosa**6-495.0*cosa**4+135.0*cosa**2-5.0)
    elseif (m==0) then
      Y=1.0/32.0*sqrt(15.0/pi)*(429.0*cosa**7-693.0*cosa**5+315.0*cosa**3-35.0*cosa)
    elseif (m==1) then
      Y=-1.0/64.0*sqrt(105.0/topi)*exp(i*angle2)*sina*(429.0*cosa**6-495.0*cosa**4+135.0*cosa**2-5.0)
    elseif (m==2) then
      Y=3.0/64.0*sqrt(35.0/pi)*exp(2.0*i*angle2)*sina**2*(143.0*cosa**5-110.0*cosa**3+15.0*cosa)
    elseif (m==3) then
      Y=-3.0/64.0*sqrt(35.0/topi)*exp(3.0*i*angle2)*sina**3*(143.0*cosa**4-66.0*cosa**2+3.0)
    elseif (m==4) then
      Y=3.0/32.0*sqrt(385.0/topi)*exp(4.0*i*angle2)*sina**4*(13.0*cosa**3-3.0*cosa)
    elseif (m==5) then
      Y=-3.0/64.0*sqrt(385.0/topi)*exp(5.0*i*angle2)*sina**5*(13.0*cosa**2-1.0)
    elseif (m==6) then
      Y=3.0/64.0*sqrt(5005.0/pi)*exp(6.0*i*angle2)*sina**6*cosa
    elseif (m==7) then
      Y=-3.0/64.0*sqrt(715.0/topi)*exp(7.0*i*angle2)*sina**7
    endif
elseif (l==8) then
   if(m==-8) then
     Y=3.0/256.0*sqrt(12155.0/topi)*exp(-8.0*i*angle2)*sina**8
   elseif (m==-7) then
     Y=3.0/64.0*sqrt(12155.0/topi)*exp(-7.0*i*angle2)*sina**7*cosa
   elseif (m==-6) then
     Y=1.0/128.0*sqrt(7293.0/pi)*exp(-6.0*i*angle2)*sina**6*(15.0*cosa**2-1.0)
   elseif (m==-5) then
     Y=3.0/64.0*sqrt(17017.0/topi)*exp(-5.0*i*angle2)*sina**5*(5.0*cosa**3-cosa)
   elseif (m==-4) then
     Y=3.0/128.0*sqrt(1309.0/topi)*exp(-4.0*i*angle2)*sina**4*(65.0*cosa**4-26.0*cosa**2+1.0)
   elseif (m==-3) then
     Y=1.0/64.0*sqrt(19635.0/topi)*exp(-3.0*i*angle2)*sina**3*(39.0*cosa**5-26.0*cosa**3+3.0*cosa)
   elseif (m==-2) then
     Y=3.0/128.0*sqrt(595.0/pi)*exp(-2.0*i*angle2)*sina**2*(143.0*cosa**6-143.0*cosa**4+33.0*cosa**2-1.0)
   elseif (m==-1) then
     Y=3.0/64.0*sqrt(17.0/topi)*exp(-i*angle2)*sina*(715.0*cosa**7-1001.0*cosa**5+385.0*cosa**3-35.0*cosa)
   elseif (m==0) then
     Y=1.0/256.0*sqrt(17.0/pi)*(6435.0*cosa**8-12012.0*cosa**6+6930.0*cosa**4-1260.0*cosa**2+35.0)
   elseif (m==1) then
     Y=-3.0/64.0*sqrt(17.0/topi)*exp(i*angle2)*sina*(715.0*cosa**7-1001.0*cosa**5+385.0*cosa**3-35.0*cosa)
   elseif (m==2) then
     Y=3.0/128.0*sqrt(595.0/pi)*exp(2.0*i*angle2)*sina**2*(143.0*cosa**6-143.0*cosa**4+33.0*cosa**2-1.0)
   elseif (m==3) then
     Y=-1.0/64.0*sqrt(19635.0/topi)*exp(3.0*i*angle2)*sina**3*(39.0*cosa**5-26.0*cosa**3+3.0*cosa)
   elseif (m==4) then
     Y=3.0/128.0*sqrt(1309.0/topi)*exp(4.0*i*angle2)*sina**4*(65.0*cosa**4-26.0*cosa**2+1.0)
   elseif (m==5) then
     Y=-3.0/64.0*sqrt(17017.0/topi)*exp(5.0*i*angle2)*sina**5*(5.0*cosa**3-cosa)
   elseif (m==6) then
     Y=1.0/128.0*sqrt(7293.0/pi)*exp(6.0*i*angle2)*sina**6*(15.0*cosa**2-1.0)
   elseif (m==7) then
     Y=-3.0/64.0*sqrt(12155.0/topi)*exp(7.0*i*angle2)*sina**7*cosa
   elseif (m==8) then
     Y=3.0/256.0*sqrt(12155.0/topi)*exp(8.0*i*angle2)*sina**8
   endif
elseif (l==9) then
   if (m==-9) then
     Y=1.0/512.0*sqrt(230945.0/pi)*exp(-9.0*i*angle2)*sina**9
   elseif (m==-4) then
     Y=3.0/128.0*sqrt(95095.0/topi)*exp(-4.0*i*angle2)*sina**4*(17.0*cosa**5-10.0*cosa**3+cosa)
   elseif (m==-5) then
     Y=3.0/256.0*sqrt(2717.0/pi)*exp(-5.0*i*angle2)*sina**5*(85.0*cosa**4-30.0*cosa**2+1.0)
   elseif (m==-6) then
     Y=1.0/128.0*sqrt(40755.0/pi)*exp(-6.0*i*angle2)*sina**6*(17.0*cosa**3-3.0*cosa)
   elseif (m==-7)  then
     Y=3.0/512.0*sqrt(13585.0/pi)*exp(-7.0*i*angle2)*sina**7*(17.0*cosa**2-1.0)
   elseif (m==-8) then
     Y=3.0/256.0*sqrt(230945.0/topi)*exp(-8.0*i*angle2)*sina**8*cosa
   elseif (m==-3) then
      Y=1.0/256.0*sqrt(21945.0/pi)*exp(-3.0*i*angle2)*sina**3*(221.0*cosa**6-195.0*cosa**4+39.0*cosa**2-1.0)
   elseif (m==-2) then
      Y=3.0/128.0*sqrt(1045.0/pi)*exp(-2.0*i*angle2)*sina**2*(221.0*cosa**7-273.0*cosa**5+91.0*cosa**3-7.0*cosa)
   elseif (m==-1) then
      Y=3.0/256.0*sqrt(95.0/topi)*exp(-i*angle2)*sina*(2431.0*cosa**8-4004.0*cosa**6+2002.0*cosa**4-308.0*cosa**2+7.0)
   elseif (m==0) then
     Y=1.0/256.0*sqrt(19.0/pi)*(12155.0*cosa**9-25740.0*cosa**7+18018.0*cosa**5-4620.0*cosa**3+315.0*cosa)  
   elseif (m==1) then
    Y=-3.0/256.0*sqrt(95.0/topi)*exp(i*angle2)*sina*(2431.0*cosa**8-4004.0*cosa**6+2002.0*cosa**4-308.0*cosa**2+7.0)
   elseif (m==2) then
     Y=3.0/128.0*sqrt(1045.0/pi)*exp(2.0*i*angle2)*sina**2*(221.0*cosa**7-273.0*cosa**5+91.0*cosa**3-7.0*cosa)
   elseif (m==3) then
     Y=-1.0/256.0*sqrt(21945.0/pi)*exp(3.0*i*angle2)*sina**3*(221.0*cosa**6-195.0*cosa**4+39.0*cosa**2-1.0)
   elseif (m==4) then
     Y=3.0/128.0*sqrt(95095.0/topi)*exp(4.0*i*angle2)*sina**4*(17.0*cosa**5-10*cosa**3+cosa)
   elseif (m==5) then
     Y=-3.0/256.0*sqrt(2717.0/pi)*exp(5.0*i*angle2)*sina**5*(85.0*cosa**4-30.0*cosa**2+1.0)
   elseif (m==6) then
     Y=1.0/128.0*sqrt(40755.0/pi)*exp(6.0*i*angle2)*sina**6*(17.0*cosa**3-3.0*cosa)
   elseif (m==7)  then
     Y=-3.0/512.0*sqrt(13585.0/pi)*exp(7.0*i*angle2)*sina**7*(17.0*cosa**2-1.0)
   elseif (m==8) then
     Y=3.0/256.0*sqrt(230945.0/topi)*exp(8.0*i*angle2)*sina**8.0*cosa
   elseif (m==9) then
     Y=-1.0/512.0*sqrt(230945.0/pi)*exp(9.0*i*angle2)*sina**9  
   endif
elseif (l==10) then
  if (m==-10) then
      Y=1.0/1024.0*sqrt(969969.0/pi)*exp(-10.0*i*angle2)*sina**10
  elseif (m==-9) then
      Y=1.0/512.0*sqrt(4849845.0/pi)*exp(-9.0*i*angle2)*sina**9*cosa
  elseif (m==-8) then
     Y=1.0/512.0*sqrt(255255.0/topi)*exp(-8.0*i*angle2)*sina**8*(19.0*cosa**2-1.0)
  elseif (m==-7) then
     Y=3.0/512.0*sqrt(85085.0/pi)*exp(-7.0*i*angle2)*sina**7*(19.0*cosa**3-3.0*cosa)
  elseif  (m==-6) then
     Y=3.0/1024.0*sqrt(5005.0/pi)*exp(-6.0*i*angle2)*sina**6*(323.0*cosa**4-102.0*cosa**2+3.0)
  elseif (m==-5) then
     Y=3.0/256.0*sqrt(1001.0/pi)*exp(-5.0*i*angle2)*sina**5*(323.0*cosa**5-170.0*cosa**3+15*cosa)
  elseif (m==-4) then
     Y=3.0/256.0*sqrt(5005.0/topi)*exp(-4.0*i*angle2)*sina**4*(323.0*cosa**6-255.0*cosa**4+45.0*cosa**2-1.0)
  elseif (m==-3) then
     Y=3.0/256.0*sqrt(5005.0/pi)*exp(-3.0*i*angle2)*sina**3*(323.0*cosa**7-357.0*cosa**5+105.0*cosa**3-7.0*cosa)
  elseif (m==-2) then
     Y=3.0/512.0*sqrt(385.0/topi)*exp(-2.0*i*angle2)*sina**2*(4199.0*cosa**8-6188.0*cosa**6+2730.0*cosa**4-364.0*cosa**2+7.0)
  elseif (m==-1) then
     Y=1.0/256.0*sqrt(1155.0/topi)*exp(-i*angle2)*sina*(4199.0*cosa**9-7956.0*cosa**7+4914.0*cosa**5-1092.0*cosa**3+63.0*cosa)
  elseif (m==0) then
     Y=1.0/512.0*sqrt(21.0/pi)*(46189.0*cosa**10-109395.0*cosa**8+90090.0*cosa**6-30030.0*cosa**4+3465.0*cosa**2-63.0)
  elseif (m==1) then
   Y=-1.0/256.0*sqrt(1155.0/topi)*exp(i*angle2)*sina*(4199.0*cosa**9-7956.0*cosa**7+4914.0*cosa**5-1092.0*cosa**3+63.0*cosa) 
  elseif (m==2) then
    Y=3.0/512.0*sqrt(385.0/topi)*exp(2.0*i*angle2)*sina**2*(4199.0*cosa**8-6188.0*cosa**6+2730.0*cosa**4-364.0*cosa**2+7.0)
  elseif (m==3)then
    Y=-3.0/256.0*sqrt(5005.0/pi)*exp(3.0*i*angle2)*sina**3*(323.0*cosa**7-357.0*cosa**5+105.0*cosa**3-7.0*cosa)
  elseif (m==4) then
    Y=3.0/256.0*sqrt(5005.0/topi)*exp(4.0*i*angle2)*sina**4*(323.0*cosa**6-255.0*cosa**4+45.0*cosa**2-1.0)
  elseif (m==5)then
    Y=-3.0/256.0*sqrt(1001.0/pi)*exp(5.0*i*angle2)*sina**5*(323.0*cosa**5-170.0*cosa**3+15.0*cosa)
  elseif (m==6) then
    Y=3.0/1024.0*sqrt(5005.0/pi)*exp(6.0*i*angle2)*sina**6*(323.0*cosa**4-102.0*cosa**2+3.0)
  elseif (m==7)then
    Y=-3.0/512.0*sqrt(85085.0/pi)*exp(7.0*i*angle2)*sina**7*(19.0*cosa**3-3.0*cosa)
  elseif (m==8)then
    Y=1.0/512.0*sqrt(255255.0/topi)*exp(8.0*i*angle2)*sina**8*(19.0*cosa**2-1.0)
  elseif(m==9) then
    Y=-1.0/512.0*sqrt(4849845.0/pi)*exp(9.0*i*angle2)*sina**9*cosa
  elseif (m==10) then
    Y=1.0/1024.0*sqrt(969969.0/pi)*exp(10.0*i*angle2)*sina**10 
  endif
  end if
  end
