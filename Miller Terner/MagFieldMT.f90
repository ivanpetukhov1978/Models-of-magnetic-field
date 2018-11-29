subroutine MagFieldMT(x,y,z) !Внешняя подпрограмма по расчету компонент магнитного поля в торе методо Миллера-Тернера

implicit double precision (a-h,o-z)

common /MFMT/ bx,by,bz           !Компоненты магнитного поля
common /initial/ rct,rt,s1,s2,b0 !Параметры тора и накрутки магнитного поля
common /const/ pi                !Константы

r=sqrt(x**2.+y**2.+z**2.)   !координаты точки в сферической системе координат
teta=acos(z/r)
fe=acos(x/sqrt(x**2.+y**2.))
if(y.lt.0.)fe=2.*pi-fe

rpt=sqrt((sqrt(x**2.+y**2.)-rct)**2.+z**2.) !координаты в тороидальной системе координат
ft=acos((sqrt(x**2.+y**2.)-rct)/(rpt))
if(z.lt.0.)ft=2.*pi-ft

a1=0.5*rpt*cos(ft)/rct !коэффициенты для функции Бесселя
a2=s2*2.41*rpt/rt
a3=s2*rt/(4.82*rct)

zJ0=fun0(a2)           !значение функции Бесселя
zJ1=fun1(a2)

br=a3*sin(ft)*zJ0   !отнормированные компоненты магнитного поля в тороидальной системе координат рассчитанные согласно методу Миллера-Тернера
bft=-((1.-a1)*zJ1-a3*cos(ft)*zJ0)
blt=(1.-a1)*zJ0

br=s1*b0*br               !компоненты магнитного поля с учетом значения поля на оси тора
bft=s1*b0*bft
blt=s1*b0*blt

!компоненты магнитного поля в сферической системе координат
brSpher=br*sin(teta+ft)+bft*cos(teta+ft)
btetSpher=br*cos(teta+ft)-bft*sin(teta+ft)
bfeSpher=blt

!компоненты магнитного поля в декартовой системе координат, передаются в основную программу
bx=brSpher*cos(fe)*sin(teta)+btetSpher*cos(fe)*cos(teta)-bfeSpher*sin(fe)
by=brSpher*sin(fe)*sin(teta)+btetSpher*sin(fe)*cos(teta)+bfeSpher*cos(fe)
bz=brSpher*cos(teta)-btetSpher*sin(teta)

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
function fun1(x)   !Функция Бесселя
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
