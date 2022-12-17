!> solve double mach reflection problem
      program main
        !> パラメータ等設定　includeを用いるのは本来好ましくない
        include "parameters.f"
        !> 最大反復数
        parameter(ilmt = 1000)
        !> 終了時間
        parameter(tfin = 0.2d0)
        !> 格子点座標 0とmax+1は便宜上　計算領域は0-max
        common /mesh/ x(0:jmax+1,0:kmax+1),y(0:jmax+1,0:kmax+1)
        !>  \(\Delta x, \Delta y\)
        !!@todo 逆数を定義しておく　intgがちょっと速くなるはず
        common /dlxy/ dltx,dlty
        !> 保存量　\(\rho,\rho u,\rho v,e \)
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        !> 基本量　\(\rho,u,v,p \)
        common /pris/ pris(4,0:jmax+1,0:kmax+1)
        !> t:現在時間, tt:サブステップの時間
        common /time/ tn,tt
        !> 反復数
        common /loop/ im
        !> before shock state
        common /stt0/ r0,p0,u0,v0,ru0,rv0,e0
        !> after shock state
        common /stt1/ r1,p1,u1,v1,ru1,rv1,e1
        !> 非粘性流束
        common /flux/ xflx(4,jmax,kmax),yflx(4,jmax,kmax)
        !> 流束評価面に補間された値 l:left, r:right
        common /sval/ pxl(4,0:jmax+1,0:kmax+1),pxr(4,0:jmax+1,0:kmax+1),
     &                pyl(4,0:jmax+1,0:kmax+1),pyr(4,0:jmax+1,0:kmax+1)
        !> 時間増分
        common /dflx/ df(4,jmax,kmax)
        !> 非粘性流束評価で用いる値
        common /sflr/ conL(4), conR(4),dnflx(4)
        !> 流束評価面の法線ベクトル
        common /nvec/ dnx,dny
        !> n stepの保存量の一時的置き場
        double precision conn(4,0:jmax+1,0:kmax+1)

! pre-process
        ! 格子生成
        call grid
        ! 初期条件代入
        call init
        im = 0
        tn = 0.0d0
        df = 0.0d0
        conn = cons
        ! 書き出し
        call outf
! main iteration
        do im=1,ilmt
!         runge kutta
!         (1)
          tt = tn+dt
          ! 境界条件代入
          call stbc
          ! 保存量から基本量へ変換
          call c2p
          ! 界面まで補間
          call intp
          ! 非粘性流束評価
          call iflx
          ! 時間増分の計算
          call intg
          do k=1,kmax
            do j=1,jmax
              cons(1,j,k) = conn(1,j,k)-dt*df(1,j,k)
              cons(2,j,k) = conn(2,j,k)-dt*df(2,j,k)
              cons(3,j,k) = conn(3,j,k)-dt*df(3,j,k)
              cons(4,j,k) = conn(4,j,k)-dt*df(4,j,k)
            end do
          end do
!         (2)
          tt = tn+0.5d0*dt
          call stbc
          call c2p
          call intp
          call iflx
          call intg
          do k=1,kmax
            do j=1,jmax
              cons(1,j,k) = 0.75* conn(1,j,k)
     &                     +0.25*(cons(1,j,k)-dt*df(1,j,k))
              cons(2,j,k) = 0.75* conn(2,j,k)
     &                     +0.25*(cons(2,j,k)-dt*df(2,j,k))
              cons(3,j,k) = 0.75* conn(3,j,k)
     &                     +0.25*(cons(3,j,k)-dt*df(3,j,k))
              cons(4,j,k) = 0.75* conn(4,j,k)
     &                     +0.25*(cons(4,j,k)-dt*df(4,j,k))
            end do
          end do
 
!         (3)
          tt = tn+dt
          call stbc
          call c2p
          call intp
          call iflx
          call intg
          do k=1,kmax
            do j=1,jmax
              cons(1,j,k) = onethird* conn(1,j,k)
     &               +2.0d0*onethird*(cons(1,j,k)-dt*df(1,j,k))
              cons(2,j,k) = onethird* conn(2,j,k)
     &               +2.0d0*onethird*(cons(2,j,k)-dt*df(2,j,k))
              cons(3,j,k) = onethird* conn(3,j,k)
     &               +2.0d0*onethird*(cons(3,j,k)-dt*df(3,j,k))
              cons(4,j,k) = onethird* conn(4,j,k)
     &               +2.0d0*onethird*(cons(4,j,k)-dt*df(4,j,k))
            end do
          end do
          conn = cons
          tn = dt*im
          !call stbc
          if(mod(im,100)==0) then
            call stbc
            call outf
          endif
          if(tn>=tfin) then
            call stbc
            call outf
            exit
          end if
        end do
! post-process

      end program main

      subroutine grid
        include "parameters.f"
        parameter (xmin=0.00d0,xmax=4.00d0)
        parameter (ymin=0.00d0,ymax=1.00d0)
        common /mesh/ x(0:jmax+1,0:kmax+1),y(0:jmax+1,0:kmax+1)
        common /dlxy/ dltx,dlty

        write(*,*) jmax,kmax
        dltx=(xmax-xmin)/(jmax-1)
        dlty=(ymax-ymin)/(kmax-1)
        write(*,*) 'dx,dy: ', dltx,dlty
        cfl = 10.0d0*dt/dmax1(dltx,dlty)
        write(*,*) 'cfl: ', cfl

        do k=0,kmax+1
          do j=0,jmax+1
            x(j,k) = (j-1)*dltx
            y(j,k) = (k-1)*dlty
          end do
        end do

        open(7,form='unformatted',file='grid.xyz')
        write(7) jmax,kmax
        write(7) ((x(j,k),j=1,jmax),k=1,kmax),
     &             ((y(j,k),j=1,jmax),k=1,kmax) 
        close(7)

      end subroutine grid

      subroutine init
        include "parameters.f"
        common /stt0/ r0,p0,u0,v0,ru0,rv0,e0
        common /stt1/ r1,p1,u1,v1,ru1,rv1,e1
        common /mesh/ x(0:jmax+1,0:kmax+1),y(0:jmax+1,0:kmax+1)
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        common /pris/ pris(4,0:jmax+1,0:kmax+1)

        r0 = 1.40d0
        p0 = 1.00d0
        u0 = 0.00d0
        v0 = 0.00d0
        ru0 = r0*u0
        rv0 = r0*v0
        e0 = p0*gm1i+0.5d0*r0*(u0*u0+v0*v0)

        r1 = 8.0d0
        p1 = 116.5d0
        u1 = 4.125d0*dsqrt(3.0d0)
        v1 = -4.125d0
        ru1 = r1*u1
        rv1 = r1*v1
        e1 = p1*gm1i+0.5d0*r1*(u1*u1+v1*v1)

        sqr3 = dsqrt(3.0d0)
        sqr3x2i = 1.0d0/(sqr3*2.0d0)
        do k=0,kmax+1
          do j=0,jmax+1
            xn = x(j,k)
            yn = y(j,k)
            yl = sqr3*xn-sqr3x2i
            if(yn>=yl) then
              cons(1,j,k) = r1
              cons(2,j,k) = ru1
              cons(3,j,k) = rv1
              cons(4,j,k) = e1
              pris(1,j,k) = r1
              pris(2,j,k) = u1
              pris(3,j,k) = v1
              pris(4,j,k) = p1
            else
              cons(1,j,k) = r0
              cons(2,j,k) = ru0
              cons(3,j,k) = rv0
              cons(4,j,k) = e0
              pris(1,j,k) = r0
              pris(2,j,k) = u0
              pris(3,j,k) = v0
              pris(4,j,k) = p0
            end if
          end do
        end do

!       b.c.
! j=1 (inflow)
        j=1
        do k=1,kmax
          cons(1,j-1,k) = cons(1,j,k)
          cons(2,j-1,k) = cons(2,j,k)
          cons(3,j-1,k) = cons(3,j,k)
          cons(4,j-1,k) = cons(4,j,k)
        end do
!
        k=1
        onesix=onethird*0.5d0
        do j=1,jmax
          if(x(j,k)>onesix) then
            vl1 = dsqrt(pris(2,j,k+1)*pris(2,j,k+1)
     &                 +pris(3,j,k+1)*pris(3,j,k+1))
            vl0 = dsign(vl1,cons(2,j,k))
            rho = cons(1,j,k+1)
            cons(1,j,k) = rho
            cons(2,j,k) = rho*vl0
            cons(3,j,k) = 0.0d0
            cons(4,j,k) = cons(4,j,k+1)
            cons(1,j,k-1) = rho
            cons(2,j,k-1) = cons(2,j,k+1)
            cons(3,j,k-1) =-cons(3,j,k+1)
            cons(4,j,k-1) = cons(4,j,k+1)
          end if
        end do

      end subroutine init

      subroutine outf
        include "parameters.f"
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        common /time/ tn,tt
        common /loop/ im
        character(len=50) fnam

        write(fnam,'("flowfield_",i5.5,".q")'),im
        open(8,form='unformatted',file=fnam)
        write(8) jmax,kmax
        write(8) 10.0d0,0.0d0,1.0d15,tn
        write(8) (((cons(n,j,k),j=1,jmax),k=1,kmax),n=1,4)
        close(8)
        write(*,*) 'output',tn,im

      end subroutine outf

      subroutine stbc
        include "parameters.f"
        common /stt0/ r0,p0,u0,v0,ru0,rv0,e0
        common /stt1/ r1,p1,u1,v1,ru1,rv1,e1
        common /mesh/ x(0:jmax+1,0:kmax+1),y(0:jmax+1,0:kmax+1)
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        common /time/ tn,tt
        common /loop/ im
! j=jmax (subsonic outflow)
        j=jmax
        do k=1,kmax
          cons(1,j,k) = cons(1,j-1,k)
          cons(2,j,k) = cons(2,j-1,k)
          cons(3,j,k) = cons(3,j-1,k)
          dke = 0.5d0*(cons(2,j,k)**2+cons(3,j,k)**2)/cons(1,j,k)
          pout =gmm1*( cons(4,j-1,k)-dke)
          cons(4,j,k) = gm1i*pout+dke
          cons(1,j+1,k) = cons(1,j,k)
          cons(2,j+1,k) = cons(2,j,k)
          cons(3,j+1,k) = cons(3,j,k)
          cons(4,j+1,k) = cons(4,j,k)
        end do
! k=kmax (pre/post shock state)
        sqrt3 = dsqrt(3.0d0)
        xshock = tt*20.0d0/sqrt3+0.5d0*onethird+1.0d0/sqrt3
        k=kmax
        do j=1,jmax
          if(x(j,k)<=xshock) then
              cons(1,j,k) = r1
              cons(2,j,k) = ru1
              cons(3,j,k) = rv1
              cons(4,j,k) = e1
              cons(1,j,k+1) = cons(1,j,k)
              cons(2,j,k+1) = cons(2,j,k)
              cons(3,j,k+1) = cons(3,j,k)
              cons(4,j,k+1) = cons(4,j,k)
            else
              cons(1,j,k) = r0
              cons(2,j,k) = ru0
              cons(3,j,k) = rv0
              cons(4,j,k) = e0
              cons(1,j,k+1) = cons(1,j,k)
              cons(2,j,k+1) = cons(2,j,k)
              cons(3,j,k+1) = cons(3,j,k)
              cons(4,j,k+1) = cons(4,j,k)
            end if
        end do

        k=1
        do j=1,jmax
          if(x(j,k)>0.5d0*onethird) then
            vl1 = dsqrt(cons(2,j,k+1)*cons(2,j,k+1)
     &                 +cons(3,j,k+1)*cons(3,j,k+1))/cons(1,j,k+1)
            vl0 = dsign(vl1,cons(2,j,k))
            rho = cons(1,j,k+1)
            cons(1,j,k) = rho
            cons(2,j,k) = rho*vl0
            cons(3,j,k) = 0.0d0
            cons(4,j,k) = cons(4,j,k+1)
            cons(1,j,k-1) = rho
            cons(2,j,k-1) = cons(2,j,k+1)
            cons(3,j,k-1) =-cons(3,j,k+1)
            cons(4,j,k-1) = cons(4,j,k+1)
          end if
        end do


      end subroutine stbc

      subroutine c2p
        include "parameters.f"
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        common /pris/ pris(4,0:jmax+1,0:kmax+1)
 
        do k=0,kmax+1
          do j=0,jmax+1
            rho = cons(1,j,k)
            rhoi = 1.0d0/rho
            u = cons(2,j,k)*rhoi
            v = cons(3,j,k)*rhoi
            p = gmm1*(cons(4,j,k)-0.5d0*rho*(u*u+v*v))
            pris(1,j,k) = rho
            pris(2,j,k) = u
            pris(3,j,k) = v
            pris(4,j,k) = p
          end do
        end do

      end subroutine c2p

      subroutine intp
        include "parameters.f"
        common /pris/ pris(4,0:jmax+1,0:kmax+1)
        common /sval/ pxl(4,0:jmax+1,0:kmax+1),pxr(4,0:jmax+1,0:kmax+1),
     &                pyl(4,0:jmax+1,0:kmax+1),pyr(4,0:jmax+1,0:kmax+1)
        parameter(eps = 1.0d0,dkpa=onethird)
        parameter(b=(3.0d0-dkpa)/(1.0d0-dkpa))
! 1st order
!       x
!        do k=2,kmax-1
!          do j=1,jmax-1
!            pxl(1,j,k) = pris(1,j,k)
!            pxl(2,j,k) = pris(2,j,k)
!            pxl(3,j,k) = pris(3,j,k)
!            pxl(4,j,k) = pris(4,j,k)
!            pxr(1,j,k) = pris(1,j+1,k)
!            pxr(2,j,k) = pris(2,j+1,k)
!            pxr(3,j,k) = pris(3,j+1,k)
!            pxr(4,j,k) = pris(4,j+1,k)
!          end do
!        end do         
!!       y
!        do j=2,jmax-1
!          do k=1,kmax-1
!            pyl(1,j,k) = pris(1,j,k)
!            pyl(2,j,k) = pris(2,j,k)
!            pyl(3,j,k) = pris(3,j,k)
!            pyl(4,j,k) = pris(4,j,k)
!            pyr(1,j,k) = pris(1,j,k+1)
!            pyr(2,j,k) = pris(2,j,k+1)
!            pyr(3,j,k) = pris(3,j,k+1)
!            pyr(4,j,k) = pris(4,j,k+1)
!          end do
!        end do
! high(3rd) order
!       x
        do k=2,kmax-1
          do j=1,jmax
            dm1 = -pris(1,j-1,k)+pris(1,j  ,k)
            dp1 = -pris(1,j  ,k)+pris(1,j+1,k)
            dm2 = -pris(2,j-1,k)+pris(2,j  ,k)
            dp2 = -pris(2,j  ,k)+pris(2,j+1,k)
            dm3 = -pris(3,j-1,k)+pris(3,j  ,k)
            dp3 = -pris(3,j  ,k)+pris(3,j+1,k)
            dm4 = -pris(4,j-1,k)+pris(4,j  ,k)
            dp4 = -pris(4,j  ,k)+pris(4,j+1,k)

            bdm1 = dsign(1.0d0,dm1)*dmax1(0.0d0,
     &              dmin1(dabs(dm1),dsign(b,dm1)*dp1))
            bdm2 = dsign(1.0d0,dm2)*dmax1(0.0d0,
     &              dmin1(dabs(dm2),dsign(b,dm2)*dp2))
            bdm3 = dsign(1.0d0,dm3)*dmax1(0.0d0,
     &              dmin1(dabs(dm3),dsign(b,dm3)*dp3))
            bdm4 = dsign(1.0d0,dm4)*dmax1(0.0d0,
     &              dmin1(dabs(dm4),dsign(b,dm4)*dp4))

            bdp1 = dsign(1.0d0,dp1)*dmax1(0.0d0,
     &              dmin1(dabs(dp1),dsign(b,dp1)*dm1))
            bdp2 = dsign(1.0d0,dp2)*dmax1(0.0d0,
     &              dmin1(dabs(dp2),dsign(b,dp2)*dm2))
            bdp3 = dsign(1.0d0,dp3)*dmax1(0.0d0,
     &              dmin1(dabs(dp3),dsign(b,dp3)*dm3))
            bdp4 = dsign(1.0d0,dp4)*dmax1(0.0d0,
     &              dmin1(dabs(dp4),dsign(b,dp4)*dm4))

            pxl(1,j,k) = pris(1,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm1+(1.0d0+dkpa)*bdp1)
            pxl(2,j,k) = pris(2,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm2+(1.0d0+dkpa)*bdp2)
            pxl(3,j,k) = pris(3,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm3+(1.0d0+dkpa)*bdp3)
            pxl(4,j,k) = pris(4,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm4+(1.0d0+dkpa)*bdp4)

            pxr(1,j-1,k) = pris(1,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp1+(1.0d0+dkpa)*bdm1)
            pxr(2,j-1,k) = pris(2,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp2+(1.0d0+dkpa)*bdm2)
            pxr(3,j-1,k) = pris(3,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp3+(1.0d0+dkpa)*bdm3)
            pxr(4,j-1,k) = pris(4,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp4+(1.0d0+dkpa)*bdm4)
          end do
        end do         
!       y
        do j=2,jmax-1
          do k=1,kmax
            dm1 = -pris(1,j,k-1)+pris(1,j,k  )
            dp1 = -pris(1,j,k  )+pris(1,j,k+1)
            dm2 = -pris(2,j,k-1)+pris(2,j,k  )
            dp2 = -pris(2,j,k  )+pris(2,j,k+1)
            dm3 = -pris(3,j,k-1)+pris(3,j,k  )
            dp3 = -pris(3,j,k  )+pris(3,j,k+1)
            dm4 = -pris(4,j,k-1)+pris(4,j,k  )
            dp4 = -pris(4,j,k  )+pris(4,j,k+1)

            bdm1 = dsign(1.0d0,dm1)*dmax1(0.0d0,
     &              dmin1(dabs(dm1),dsign(b,dm1)*dp1))
            bdm2 = dsign(1.0d0,dm2)*dmax1(0.0d0,
     &              dmin1(dabs(dm2),dsign(b,dm2)*dp2))
            bdm3 = dsign(1.0d0,dm3)*dmax1(0.0d0,
     &              dmin1(dabs(dm3),dsign(b,dm3)*dp3))
            bdm4 = dsign(1.0d0,dm4)*dmax1(0.0d0,
     &              dmin1(dabs(dm4),dsign(b,dm4)*dp4))
            bdp1 = dsign(1.0d0,dp1)*dmax1(0.0d0,
     &              dmin1(dabs(dp1),dsign(b,dp1)*dm1))
            bdp2 = dsign(1.0d0,dp2)*dmax1(0.0d0,
     &              dmin1(dabs(dp2),dsign(b,dp2)*dm2))
            bdp3 = dsign(1.0d0,dp3)*dmax1(0.0d0,
     &              dmin1(dabs(dp3),dsign(b,dp3)*dm3))
            bdp4 = dsign(1.0d0,dp4)*dmax1(0.0d0,
     &              dmin1(dabs(dp4),dsign(b,dp4)*dm4))

            pyl(1,j,k) = pris(1,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm1+(1.0d0+dkpa)*bdp1)
            pyl(2,j,k) = pris(2,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm2+(1.0d0+dkpa)*bdp2)
            pyl(3,j,k) = pris(3,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm3+(1.0d0+dkpa)*bdp3)
            pyl(4,j,k) = pris(4,j,k)+0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdm4+(1.0d0+dkpa)*bdp4)

            pyr(1,j,k-1) = pris(1,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp1+(1.0d0+dkpa)*bdm1)
            pyr(2,j,k-1) = pris(2,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp2+(1.0d0+dkpa)*bdm2)
            pyr(3,j,k-1) = pris(3,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp3+(1.0d0+dkpa)*bdm3)
            pyr(4,j,k-1) = pris(4,j,k)-0.25d0*eps*
     &                    ((1.0d0-dkpa)*bdp4+(1.0d0+dkpa)*bdm4)
 
          end do
        end do



      end subroutine intp

      subroutine slau
        implicit double precision(a-h,o-z)
        common /sflr/ conL(4), conR(4),dnflx(4)
        common /nvec/ dnx,dny
        parameter(gmma = 1.4d0, gmm1 = gmma-1.0d0, gm1i = 1.0d0/gmm1)

!  Left state
        rl = conL(1)
        ul = conL(2)/conL(1)
        vl = conL(3)/conL(1)
        pl = gmm1*( conL(4) - 0.5d0*rL*(uL*uL+vL*vL) )
        al = dsqrt(gmma*pL/rL)
        dhl = ( conL(4) + pL ) / rL
        vnp = ul*dnx+vl*dny
!  Right state
        rr = conR(1)
        ur = conR(2)/conR(1)
        vr = conR(3)/conR(1)
        pr = gmm1*( conR(4) - 0.5d0*rR*(uR*uR+vR*vR) )
        ar = dsqrt(gmma*pR/rR)
        dhr = ( conR(4) + pR ) / rR
        vnm = ur*dnx+vr*dny
        
        cbar = 0.5d0*(al+ar)
        cbi = 1.0d0/cbar
        dmp = vnp*cbi
        dmm = vnm*cbi
        g = -dmax1(dmin1(dmp,0.0d0),-1.0d0)*
     &       dmin1(dmax1(dmm,0.0d0), 1.0d0)
        dmht = dmin1(1.0d0,cbi*dsqrt(0.5d0*(ul*ul+vl*vl+ur*ur+vr*vr)))
        xi = (1.0d0-dmht)*(1.0d0-dmht)
        absbv = (rl*dabs(vnp)+rr*dabs(vnm))/(rl+rr)
        dmdt = 0.5d0*(rl*vnp+rr*vnm-absbv*(rr-rl))*(1.0d0-g)
     &          -xi*0.5d0*cbi*(pr-pl)
        admd = dabs(dmdt)

        if(dabs(dmp)<1.0d0) then
          bp = 0.25d0*(2.0d0-dmp)*(dmp+1.0d0)*(dmp+1.0d0)
        else
          bp = 0.5d0*(1.0d0+dsign(1.0d0, dmp))
        end if
        if(dabs(dmm)<1.0d0) then
          bm = 0.25d0*(2.0d0+dmm)*(dmm-1.0d0)*(dmm-1.0d0)
        else
          bm = 0.5d0*(1.0d0+dsign(1.0d0,-dmm))
        end if

        ptld = 0.5d0*( (pl+pr)+(bp-bm)*(pl-pr)
     &                +(1.0d0-xi)*(bp+bm-1.0d0)*(pl+pr))

        dm1 = dmdt+admd
        dm2 = dmdt-admd
        dnflx(1) = dmdt
        dnflx(2) = 0.5d0*(dm1* ul+dm2* ur)+ptld*dnx
        dnflx(3) = 0.5d0*(dm1* vl+dm2* vr)+ptld*dny
        dnflx(4) = 0.5d0*(dm1*dhl+dm2*dhr)


      end subroutine slau

      subroutine iflx
        include "parameters.f"
        common /sval/ pxl(4,0:jmax+1,0:kmax+1),pxr(4,0:jmax+1,0:kmax+1),
     &                pyl(4,0:jmax+1,0:kmax+1),pyr(4,0:jmax+1,0:kmax+1)
        common /flux/ xflx(4,jmax,kmax),yflx(4,jmax,kmax)
!        double precision conl(4),conr(4),dnflx(4)
        common /sflr/ conL(4), conR(4),dnflx(4)
        common /nvec/ dnx,dny
 
        one = 1.0d0
        zero = 0.0d0
!       x flux
        dnx = 1.0d0
        dny = 0.0d0
        do k=2,kmax-1
          do j=1,jmax-1
!           left state            
            rl = pxl(1,j,k)
            ul = pxl(2,j,k)
            vl = pxl(3,j,k)
            pl = pxl(4,j,k)
            conl(1) = rl 
            conl(2) = rl*ul 
            conl(3) = rl*vl 
            conl(4) = gm1i*pl+0.5d0*rl*(ul*ul+vl*vl)
!           right state
            rr = pxr(1,j,k)
            ur = pxr(2,j,k)
            vr = pxr(3,j,k)
            pr = pxr(4,j,k) 
            conr(1) = rr
            conr(2) = rr*ur 
            conr(3) = rr*vr 
            conr(4) = gm1i*pr+0.5d0*rr*(ur*ur+vr*vr)
            call slau
            xflx(1,j,k) = dnflx(1)
            xflx(2,j,k) = dnflx(2)
            xflx(3,j,k) = dnflx(3)
            xflx(4,j,k) = dnflx(4)
          end do
        end do        
!       y flux
        dnx = 0.0d0
        dny = 1.0d0
        do j=2,jmax-1
          do k=1,kmax-1
!           left state            
            rl = pyl(1,j,k)
            ul = pyl(2,j,k)
            vl = pyl(3,j,k)
            pl = pyl(4,j,k)
            conl(1) = rl 
            conl(2) = rl*ul 
            conl(3) = rl*vl 
            conl(4) = gm1i*pl+0.5d0*rl*(ul*ul+vl*vl)
!           right state
            rr = pyr(1,j,k)
            ur = pyr(2,j,k)
            vr = pyr(3,j,k)
            pr = pyr(4,j,k) 
            conr(1) = rr
            conr(2) = rr*ur 
            conr(3) = rr*vr 
            conr(4) = gm1i*pr+0.5d0*rr*(ur*ur+vr*vr)
            call slau
            yflx(1,j,k) = dnflx(1)
            yflx(2,j,k) = dnflx(2)
            yflx(3,j,k) = dnflx(3)
            yflx(4,j,k) = dnflx(4)
          end do
        end do
        !stop
      end subroutine iflx

      subroutine intg
        include "parameters.f"
        common /flux/ xflx(4,jmax,kmax),yflx(4,jmax,kmax)
        common /cons/ cons(4,0:jmax+1,0:kmax+1)
        common /dflx/ df(4,jmax,kmax)
        common /dlxy/ dltx,dlty


        do k=2,kmax-1
          do j=2,jmax-1
 
            df(1,j,k) = (-xflx(1,j-1,k)+xflx(1,j,k))/dltx
            df(2,j,k) = (-xflx(2,j-1,k)+xflx(2,j,k))/dltx
            df(3,j,k) = (-xflx(3,j-1,k)+xflx(3,j,k))/dltx
            df(4,j,k) = (-xflx(4,j-1,k)+xflx(4,j,k))/dltx
            df(1,j,k) = df(1,j,k)+(-yflx(1,j,k-1)+yflx(1,j,k))/dlty
            df(2,j,k) = df(2,j,k)+(-yflx(2,j,k-1)+yflx(2,j,k))/dlty
            df(3,j,k) = df(3,j,k)+(-yflx(3,j,k-1)+yflx(3,j,k))/dlty
            df(4,j,k) = df(4,j,k)+(-yflx(4,j,k-1)+yflx(4,j,k))/dlty

          end do
        end do

      end subroutine intg