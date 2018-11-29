subroutine MagFieldRV(x,y,z)  !Расчет магнитного поля в торе методом Ромашец-Вандас
implicit double precision (a-h,o-z)

common /MFRV2/ bx,by,bz
common /initial/ R0,s1,s2,b0
common /const/ pi

a=sqrt(r0**2.-1.)   !Вспомогательная величина

 call fun1(r0,eps)  !

r2=x**2.+y**2.+z**2.   !квадрат радиуса точки
f1=(r2+a**2.)**2./(4.*a**2.*(x**2.+y**2.))  !Вспомогательная функция
zmu=asinh(1./sqrt(f1-1.))                   !
f2=(r2-a**2.)/(2.*a*z)                      !Вспомогательная функция
psi=-1./(sinh(zmu)**2.)                     !

if(z.ge.0.)then
eta=acos(f2/sqrt(1.+f2**2.))                !
else
eta=-pi+acos(f2/sqrt(1.+f2**2.))
endif

argcos=x*(cosh(zmu)-cos(eta))/(a*sinh(zmu)) !

if(argcos.ge.1.)argcos=0.99999999
if(argcos.le.-1.)argcos=-0.99999999

if(y.ge.0.)then
fe=acos(argcos)                             !
else
fe=-acos(argcos)
endif

 call fun2(psi,fgip0,fgip1,eps)            !

beta=s2*eps*cosh(zmu)*(cosh(zmu)-cos(eta))*fgip1/(2.*(sinh(zmu))**3.)
bfe=s1*(cosh(zmu)-cos(eta))*fgip0/sinh(zmu)

dd=1./(cosh(zmu)-cos(eta))                             !Вспомогательная функция

bx=-dd*beta*sinh(zmu)*sin(eta)*cos(fe)-bfe*sin(fe)    !Расчет компонент магнитного поля согласно модели Ромашец-Вандас
by=-dd*beta*sinh(zmu)*sin(eta)*sin(fe)+bfe*cos(fe)
bz=dd*beta*(cosh(zmu)*cos(eta)-1.)

return
end subroutine



!!! Гипергеометрическая функция аппроксимируется гипергеометрическим рядом Гаусса приведенным в "Справочнике по спецфункциям" под редакцией М. Абрамовица и И. Стиган. Москва. Наука. 1979.
subroutine fun1(r0,eps)  !Подпрограмма для определения eps, корня гипергеометрической функции. Методом подбора определяется значение eps где гипергеометрическая функция меняет знак.
implicit double precision (a-h,o-z)

eps=0.1   !Начальное значение
deps=0.01  !Шаг с которым подбирается корень функции

psi=1./(1.-(r0)**2.)   !Вспомогательная функция определяется заданными параметрами тора (отношением расстояния от центра тора до его оси к радиусу поперечного сечения тора)

do keps=1,5000
d1=1.-4.*eps**2.   !Вспомогательная функция

if(d1.ge.0.)then   !
delta=sqrt(d1)/4.
zk=-1.
else
delta=sqrt(-d1)/4.
zk=1.
endif

fgip0=1.        !Гипергеометрическая функция заданная только первым членом ряда 
gold0=1.        !Первый член ряда

do k=1,500
zk1=float(k)

gnew0=gold0*((zk1-1.+0.25)**2.+zk*delta**2.)/zk1*(psi/zk1)   !член ряда, которым аппроксимируется гипергеометрическая функция
fgip0=fgip0+gnew0                                            !гипергеометрическая функция
dd0=abs(gnew0)/abs(fgip0)                                    !Приращение функции очередным членом ряда

if(dd0.lt.1.e-4)goto 10                                      !Если приращение меньше 0.01% расчет заканчивается, ряд считается сошедшимся
gold0=gnew0
enddo
10 continue

if(keps.eq.1)then                                             !Определяется смена знака гипергеометрической функции
g1=fgip0
else
g2=g1*fgip0
 if(g2.lt.0.)goto 15                                          !Если знак изменился запоминается eps и подпрограмма заканчивает работу
endif

eps=eps+deps
enddo !keps

print*,'SOS fun1'      !Если не обнаружен корень
pause
15 continue

return
end subroutine


!!! Гипергеометрическая функция аппроксимируется гипергеометрическим рядом Гаусса приведенным в "Справочнике по спецфункциям" под редакцией М. Абрамовица и И. Стиган. Москва. Наука. 1979.
subroutine fun2(psi,fgip0,fgip1,eps) !Расчет гипергеометрической функции
implicit double precision (a-h,o-z)

d1=1.-4.*eps**2.

if(d1.ge.0.)then
delta=sqrt(d1)/4.
zk=-1.
else
delta=sqrt(-d1)/4.
zk=1.
endif

fgip0=1.     !Гипергеометрическая функция заданная только первым членом ряда 
fgip1=1.     !Гипергеометрическая функция заданная только первым членом ряда, все аргументы больше на 1

gold0=1.     !Первый член ряда
gold1=1.     !Первый член ряда

do k=1,500
zk1=float(k)

gnew0=gold0*((zk1-1.+0.25)**2.+zk*delta**2.)/zk1*(psi/zk1)       !Член ряда, которым аппроксимируется гипергеометрическая функция
gnew1=gold1*((zk1+0.25)**2.+zk*delta**2.)/(zk1+1.)*(psi/zk1)     !Член ряда, которым аппроксимируется гипергеометрическая функция, все аргументы больше на 1

fgip0=fgip0+gnew0    !Гипергеометрическая функция
fgip1=fgip1+gnew1    !Гипергеометрическая функция

dd0=abs(gnew0)/abs(fgip0)    !Приращение функции очередным членом ряда
dd1=abs(gnew1)/abs(fgip1)    !Приращение функции очередным членом ряда

if(dd0.lt.1.e-4.and.dd1.lt.1.e-4)goto 10   !Если приращение обоих функций меньше 0.01% расчет заканчивается, ряды считаются сошедшимися
gold0=gnew0
gold1=gnew1
enddo
print*,'SOS fun2'      !Если хотя бы один ряд не сошелся
pause

10 continue



return
end subroutine











