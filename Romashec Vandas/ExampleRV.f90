program ExampleRV !Программа - пример использования подпрограммы по расчету компонент магнитного поля в торе методом Ромашец-Вандас, строится силовая линия магнитного поля в торе и распределение напряженности по сечению тора.
implicit double precision (a-h,o-z)

common /MFRV2/ bx,by,bz           !Компоненты магнитного поля
common /initial/ R0,s1,s2,b0 !Параметры тора и накрутки магнитного поля
common /const/ pi                !Константы
pi=3.1415926535897932384626433832795

!Параметры определяющие тор и магнитное поле в нем.
rct=10. !Расстояние от центра тора до оси
rt=1.   !Расстояние от оси до поверхности тора
s1=1.   !коэффициенты определяющие тип накрутки магнитного поля
s2=1.
b0=1.   !величина напряженности магнитного поля на оси тора

R0=rct/rt !Так как в формулах Ромашец-Вандас принято что расстояние от оси до поверхности тора является масштабом переходим к величине R0 расстояние от центра тора до оси выраженное в этих расстояниях


!!! Постороение силовой линии магнитного поля
open(unit=10,file='dataRV-1.dat')

dR=0.00001 !шаг с которым рассчитывается магнитное поле

x=r0+0.1   !координаты начальной точки силовой линии
y=0.
z=0.0001

do iR=1,100000000 !00

 call MagFieldRV(x,y,z) !вызов подпрограммы для расчета магнитного поля методом Ромашец-Вандас

bb=sqrt(bx**2.+by**2.+bz**2.) !расчет модуля магнитного поля

sinAlpha1=y/sqrt(x**2.+y**2.) !синус угла между радиус вектором текущей точки на силовой линии и осью Х

x=x+dr*bx/bb !координаты новой точки на силовой линии
y=y+dr*by/bb
z=z+dr*bz/bb

sinAlpha2=y/sqrt(x**2.+y**2.) !синус угла между радиус вектором новой точки на силовой линии и осью Х

if(sinAlpha1.lt.0.and.sinAlpha2.gt.0.)goto 1 !Конец построения силовой линии при пересечении оси Х. Для параметров s1=1 s2=1 это один оборот силовой линии в торе.


if(ir/10000.eq.float(ir)/10000.)write(10,*)sqrt(x**2.+y**2.)-r0,z,x,y,z,bx,by,bz !запись рассчитанной проекции силовой линии на вертикальное сечение прореженая через 10000 точек, силовой линии в единицах rt,  и компоненты магнитного поля вдоль силовой линии.

enddo

1 continue
!!!конец блока по построению силовой линии

!!!Блок по получению распределения магнитного поля в сечении тора плоскостью XZ
open(unit=10,file='dataRV-2.dat')

y=0.   !в точках с этими координатами будет определено магнитное поле.
z=1.5
do i1=1,501
z=z-3./501. !Z будет меняться в пределах от -1.5 до 1.5, с разбиением на 500 шагов
x=r0+1.5
do j1=1,501
x=x-3./501.  !Х будет меняться в пределах от 8.5 до 11.5, с разбиением на 500 шагов


rTek=sqrt((x-r0)**2.+z**2.)  !Определяется расстояние от оси тора до текущей точки
if(rTek.le.1.)then 
 call MagFieldRV(x,y,z)   !Если точка лежит в пределах тора, вызывается подпрограмма по расчету напряженности магнитного поля
bb=sqrt(bx**2.+by**2.+bz**2.)   !Рассчитывается модуль магнитного поля
else
bb=0.  !Если точка лежит вне тора, то напряженность магнитного поля берется равна 0
endif

if(j1.eq.501)then    !Для удобства построения результаты записываются в матрицу 501х501 точек, которая в дальнейшем будет обработана в программе QtiPlot или аналогичной.
write(10,'(f15.10)')bb
else
write(10,'(f15.10$)')bb
endif

enddo
enddo

!!!! Конец блока по построению распределения магнитного поля в сечении тора.


stop
end

include "MagFieldRV.f90"  !Внешняя подпрограмма по расчету компонент магнитного поля в торе методом Ромашец-Вандас
