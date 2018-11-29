subroutine MagFieldRuf(xTor,yTor,zTor)
implicit double precision (a-h,o-z)

common /dat1/ bx,by,bz,wx,wy,wz
common /ParameterTor/ rc,a0,wTor,bTor0,rTorMax,alfa,beta,s1,s2
common /const/ pi,zk1,zk2,re,omega,w,b0,e,ep,ek,p0,bre,rlar,tt,c

xyTor=sqrt(xTor**2.+yTor**2.)-rc   !расстояние от оси тора до проекции точки на плоскость ХY
rTor=sqrt(xyTor**2.+zTor**2.)      !расстояние от оси тора до точки

tetTor=acos(xyTor/rTor)            !угол teta в квазитороидальной системе координат
if(zTor.lt.0.)tetTor=2.*pi-tetTor

feTor=acos(xTor/(rc+rTor*cos(tetTor)))  !угол fee в квазитороидальной системе координат
if(yTor.lt.0.)feTor=-abs(feTor)

t0=wTor*pi/2.                                   !вспомогательные коэффициенты
t2=((rTor/a0)**2.+(cos(feTor/2.))**2.)**(-1.5)
t3=rTor/(rc+rTor*cos(tetTor))

!Расчет компонент магнитного поля в квазитороидальной системе координат согласно модели Руффоло
br=-s1*0.5*sin(feTor/2.)*t2*t3*bTor0                                   !радиальная компонента
bfe=s1*cos(feTor/2.)*t2*bTor0                                          !fee-товская компонента
btet=s2*t0*t2*t3*bTor0*exp(-rTor/(a0*cos(feTor/2.)))*(cos(feTor/2.))**2.  !teta- компонента          !модель Руффоло в которой спиральность уменьшается к поверхности тора, также как и  напряженность.
!btet=s2*t0*t2*t3*bTor0*exp(rTor/(a0*cos(feTor/2.)))*(cos(feTor/2.))**2.  !teta- компонента          !Видоизмененная модель Руффоло в которой спиральность увеличивается к поверхности тора, также как и  напряженность.

b=sqrt(br**2.+bfe**2.+btet**2.)                                   !расчет напряженности магнитного поля

bx=br*cos(tetTor)*cos(feTor)-bfe*sin(feTor)-btet*sin(tetTor)*cos(feTor)   !расчет компонент магнитного поля в декартовой системе координат
by=br*cos(tetTor)*sin(feTor)+bfe*cos(feTor)-btet*sin(tetTor)*sin(feTor)
bz=br*sin(tetTor)+btet*cos(tetTor)

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Определение находится ли точка внутри квазитора.
subroutine INorOUT(xTor,yTor,zTor,location)   !система координат с центром в центре тора, ось Х от Солнца к земле, Z вертикально вверх к плоскости тора масштаб 1=RCT
implicit double precision (a-h,o-z)

common /ParameterTor/ rc,a0,wTor,bTor0,rTorMax,alfa,beta,s1,s2
common /const/ pi,zk1,zk2,re,omega,w,b0,e,ep,ek,p0,bre,rlar,tt,c

xyTor=sqrt(xTor**2.+yTor**2.)-rc   !расстояние от оси тора до проекции точки на плоскость ХY
rTor=sqrt(xyTor**2.+zTor**2.)      !расстояние от оси тора до точки

feTor=acos(xTor/(rc+rTor*cos(tetTor)))  !угол fee в квазитороидальной системе координат
if(yTor.lt.0.)feTor=2.*pi-feTor

rTorGran=rTorMax*cos(feTor/2.)      !расстояние от оси тора до поверхности для угла feTor выбранной точки в пространстве

if(abs(rTorGran).ge.rTor)then   !если точка находится внутри квазитора то location=1 если снаружи location=2
location=1                      
else
location=2
endif

return
end subroutine
