subroutine sech   !Подпрограмма для получения контура поверхности квазитора.
implicit double precision (a-h,o-z)
parameter(N=100000)
dimension x(N),y(N),x1(N),y1(N),gran(2,12),x2(N),y2(N)

common /kontur/ zkoord(2,N),Xcentr,Ycentr,mmm
 common /const/ pi

Xcentr=6.4 !2.3     !координаты центра вспомогательной системы координат. Должны лежать внутри контура. Для того что бы угол менялся от 0 до 2 пи. Каждому углу соответствовала одна точка на контуре. Нужную координату надо преобразовывать также.
Ycentr=0.           !Можно взять координаты оси образующего тора, если же квазитор лежит вне образующего тора, то можно взять координаты произвольно

x=0.  !Обнуление массивов
y=0.
zkoord=0.

nnn=1
open(unit=11,file='dataINT-1.dat')   !считываем ранее насчитанную силовую линию, которую принимаем за лежащую на поверхности
do

read(11,*,end=14)znach1,znach2,znach3,x2(nnn),y2(nnn),znach4,znach5 !Для определения контура необходима только проекция силовой линии на сечение.
x2(nnn)=x2(nnn)-Xcentr     !Переводим координаты проекции сечения силовой линии в систему координат центр которой лежит внутри контура
y2(nnn)=y2(nnn)-Ycentr
nnn=nnn+1
enddo

14 continue
 close(11)

if(y2(1).gt.y2(2))then   !точки сечения силовой линии переводятся в новые массивы которые заполнены против часовой стрелки не смотря на накрутку силовой линии
do i=1,nnn
x(i)=x2(nnn+1-i)
y(i)=y2(nnn+1-i)
enddo
else
x=x2
y=y2
endif

do i=2,1000000    !находится точка с максимальным x, ее номер nachalo, контур будет начинаться с нее
if(x(i-1).lt.x(i).and.x(i+1).lt.x(i))then
nachalo=i
if(y(i).lt.0.)nachalo=nachalo+1
goto 1
endif
enddo

print*,'AHTUNG ne hvataet kontura'   !если контур не замкнется
pause

1 continue

print*,'nachalo',nachalo

ugolSTAR=0.   !начальный угол отсчитывается от точки с максимальным x
do i=1,nnn-nachalo
ugol=atan(abs(y(nachalo+i-1))/abs(x(nachalo+i-1)))   !расчитываются углы каждой точки сечения 

if(y(nachalo+i-1).gt.0.and.x(nachalo+i-1).gt.0.)ugol=ugol    
if(y(nachalo+i-1).gt.0.and.x(nachalo+i-1).lt.0.)ugol=pi-ugol
if(y(nachalo+i-1).lt.0.and.x(nachalo+i-1).lt.0.)ugol=pi+ugol
if(y(nachalo+i-1).lt.0.and.x(nachalo+i-1).gt.0.)ugol=2*pi-ugol

if(ugolSTAR.gt.ugol)then
goto 2    !когда новый угол станет меньше предыдущего контур замкнут
else
ugolSTAR=ugol
endif

rr=sqrt(y(nachalo+i-1)**2.+x(nachalo+i-1)**2.)    !расстояние от центра вспомогательной системы координат до контура
zkoord(1,i)=ugol
zkoord(2,i)=rr
enddo

2 continue

do i=1,nnn        !определяется сколько точек в контуре
if(zkoord(2,i).eq.0.)then
goto 3
endif
enddo

3 continue
mmm=i-1


open(unit=211,file='kontur3.dat')   !записываются для каждой точки контура угол и расстояние до центра вспомогательной системы координат
do i=1,mmm
write(211,*)zkoord(1,i),zkoord(2,i)
enddo
 close(211)

return
end subroutine




!!! Подпрограмма определяющая находится ли точка внутри квазитора или снаружи
subroutine INorOUT(x,y,nObl)    !nObl=1     !внутри области; nObl=0     !снаружи области
implicit double precision (a-h,o-z)
parameter(N=100000)

 common /const/ pi
common /kontur/ zkoord(2,N),Xcentr,Ycentr,mmm

xTek=x-Xcentr    !Координаты текущей точки переводятся во вспомогательную систему координат
yTek=y-Ycentr

ugolTek=atan(abs(yTek)/abs(xTek))   !определяется угол между радиус вектором в эту точку и осью Х

if(yTek.gt.0.and.xTek.gt.0.)ugolTek=ugolTek
if(yTek.gt.0.and.xTek.lt.0.)ugolTek=pi-ugolTek
if(yTek.lt.0.and.xTek.lt.0.)ugolTek=pi+ugolTek
if(yTek.lt.0.and.xTek.gt.0.)ugolTek=2*pi-ugolTek

rTek=sqrt(yTek**2.+xTek**2.)       !Определяется расстояние от центра вспомогательной системы координат до текущей точки


!Так как для точек контура угол от оси Х монотонно увеличивается, то можно брать точку по середине и сравнивать ее угол с искомым, если углы не равны, берем подходящую часть и делим пополам, повторяем итерацию.
kmax=mmm
kmin=1
ksr=0
do k=1,1000
ksr=kmin+(kmax-kmin)/2

if(zkoord(1,ksr).eq.ugolTek)then  
 if(zkoord(2,ksr).ge.rTek)then   !если текущий угол равняется одному из углов точек контура, сравниваются расстояния от центра контура до контура для этого угла и расстояние до текущей точки, если точка ближе, то она находится внутри квазитора
 nObl=1     !внутри области
 else
 nObl=0     !снаружи области
 endif
 goto 5
else
 if(zkoord(1,ksr).gt.ugolTek)then  
 kmax=ksr
 else
 kmin=ksr
 endif
endif

if(kmin+1.eq.kmax)then   !итерации повторяются пока соседние точки контура не будут располагаться по обе стороны текущей точки
 if(zkoord(2,kmax).ge.rTek)then !Для точки с большим углом сравнивается расстояние от центра контура до границы контура и до текущей точки, и определяется лежит ли точка внутри контура или вне.
 nObl=1     !внутри области
 else
 nObl=0     !снаружи области
 endif
 goto 5
endif

!print*,k,kmin,kmax,zkoord(1,kmin),zkoord(1,kmax),ugolTek
enddo

5 continue

return
end subroutine




subroutine MagFieldInt(x,y,z)

implicit double precision (a-h,o-z)
common /torparam/ R0,s1,s2,b0
common /const/ pi
common /BInt/ bx,by,bz

parameter(nf=180)

dimension f(nf)

hf=2.*pi/float(nf)

f(1)=0.

do i=1,nf-1
f(i+1)=f(i)+hf
enddo

sum1=0.
sum2=0.
sum3=0.

do i=1,nf
r=sqrt(z**2.+(x*cos(f(i))+y*sin(f(i))-R0)**2.)
rr=s2*2.41*r
bes0=fun0(rr)               !Функции Бесселя
bes1=fun1(rr)

sum1=sum1+bes1*z*cos(f(i))/r-bes0*sin(f(i))    !Интегрирование
sum2=sum2+bes1*z*sin(f(i))/r+bes0*cos(f(i))
sum3=sum3-bes1*(x*cos(f(i))+y*sin(f(i))-R0)/r
enddo


bx=sum1*hf*s1*b0    !Расчет компонент магнитного поля Интегральным методом
by=sum2*hf*s1*b0
bz=sum3*hf*s1*b0

return
end subroutine


!!! Функция Бесселя первого рода, нулевого порядка аппроксимируется многочленами приведенными в "Справочнике по спецфункциям" под редакцией М. Абрамовица и И. Стиган. Москва. Наука. 1979.
function fun0(x)    !Функция Бесселя
implicit double precision (a-h,o-z)
if(x.le.3.)then
a=x/3.
a1=1.-2.2499997*a**2.
a2=1.2656208*a**4.-0.3163866*a**6.
a3=0.0444479*a**8.-0.0039444*a**10.+0.00021*a**12.
fun0=a1+a2+a3
else
a=3./x
f1=0.79788456-0.00000077*a
f2=-0.0055274*a**2.-0.00009512*a**3.
f3=0.00137237*a**4.-0.00072805*a**5.+0.00014476*a**6.
f0=f1+f2+f3

t1=x-0.78539816-0.04166397*a-0.00003954*a**2.
t2=0.00262573*a**3.-0.00054125*a**4.
t3=-0.00029333*a**5.+0.00013558*a**6.
t0=t1+t2+t3
fun0=x**(-0.5)*f0*cos(t0)
endif

return
end function

!!! Функция Бесселя первого рода, первого порядка аппроксимируется многочленами приведенными в "Справочнике по спецфункциям" под редакцией М. Абрамовица и И. Стиган. Москва. Наука. 1979.
function fun1(x)     !Функция Бесселя
implicit double precision (a-h,o-z)
if(x.le.3.)then
a=x/3.
a1=0.5-0.56249985*a**2.+0.21093573*a**4.
a2=-0.03954289*a**6+0.00443319*a**8.
a3=-0.00031761*a**10+0.00001109*a**12.
fun1=x*(a1+a2+a3)
else
a=3./x
f2=0.79788456+0.00000156*a+0.01659667*a**2.
f3=0.00017105*a**3.-0.00249511*a**4.
f4=0.00113653*a**5.-0.00020033*a**6.
f1=f2+f3+f4

t2=x-2.35619449+0.12499612*a
t3=0.0000565*a**2.-0.00637879*a**3.
t4=0.00074348*a**4.+0.00079824*a**5.-0.00029166*a**6.
t1=t2+t3+t4
fun1=x**(-0.5)*f1*cos(t1)
endif

return
end
