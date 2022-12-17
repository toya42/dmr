      implicit double precision (a-h,o-z)
      !> x方向格子点数
      parameter(jmax=1001)
      !> y方向格子点数
      parameter(kmax= 251)
      !> 時間刻み幅
      parameter(dt=5.0d-5)
      !> 比熱比と関連定数
      parameter(gmma = 1.4d0, gmm1 = gmma-1.0d0, gm1i = 1.0d0/gmm1)
      !> 定数（計算速度性能に寄与するかはよくわからない）
      parameter(onethird=1.0d0/3.0d0)