!USE buffersize 12M

program state_dummy

      use kinds_mod
      use constants
      use omp_lib

      implicit none

      integer (int_kind), parameter :: &
      nx_block = 324, &
      ny_block = 388, &
      km = 60

      

      real (r8), dimension(nx_block,ny_block,km) :: TMIX,TMIX1,TMIX2,TMIX3
 
      real (r8), dimension(nx_block,ny_block,km) :: DRDT,DRDS

      real (r8), dimension(2880,2619) :: TTMIX,TEMPDRDT,TEMPDRDS,TEMPTMIX1,TEMPTMIX2,KARRAY

      integer (int_kind) kk, this_block,counter,k

      integer (int_kind),dimension(nx_block,ny_block,km) :: K_3darray

      integer (int_kind),dimension(2880,2619) :: K_2darray

      this_block = 1
      kk = 1

      do k=1,km
      K_3darray(:,:,k)=k
      enddo 
      
      call random_number(TMIX)  
      call random_number(TMIX1)
      call random_number(TMIX2)
      call random_number(TMIX3)

      K_2darray = RESHAPE(K_3darray , (/2880, 2619/))
      TTMIX = RESHAPE(TMIX, (/2880, 2619/))
      TEMPDRDT = RESHAPE(DRDT,(/2880,2619/))
      TEMPDRDS = RESHAPE(DRDS,(/2880,2619/))
      TEMPTMIX1 = RESHAPE(TMIX1,(/2880,2619/)) 
      TEMPTMIX2 = RESHAPE(TMIX2,(/2880,2619/))

      do counter=1,10 
      call state (kk, kk, TTMIX, TTMIX, this_block ,k_2Darray,RHOOUT=TEMPDRDT,RHOFULL=TEMPDRDS,DRHODT=TEMPTMIX1,DRHODS=TEMPTMIX2)  
      enddo     

      TMIX = RESHAPE( TTMIX, (/324, 388,60/) )
      DRDT = RESHAPE( TEMPDRDT, (/324,388,60/) )
      DRDS = RESHAPE( TEMPDRDS, (/324,388,60/) )
      TMIX1 = RESHAPE( TEMPTMIX1, (/324,388,60/) )
      TMIX2 = RESHAPE( TEMPTMIX2, (/324,388,60/) )

      print *,"sum of DRDT is",sum(DRDT)
                

contains

subroutine state(k, kk, TEMPK, SALTK, this_block, K_2Darray, RHOOUT, RHOFULL, DRHODT, DRHODS)




      integer (int_kind), intent(in) :: this_block 

      !integer (int_kind), parameter :: &
      !nx_block = 2880, & 
      !ny_block = 2619, &
      !km = 60

      integer (int_kind) i,j  

      real (r8) start_time, end_time,first_time

      real (r8), parameter ::                  &
      mwjfnp0s0t0 =   9.99843699e+2_r8 * p001, &
      mwjfnp0s0t1 =   7.35212840e+0_r8 * p001, &
      mwjfnp0s0t2 =  -5.45928211e-2_r8 * p001, &
      mwjfnp0s0t3 =   3.98476704e-4_r8 * p001, &
      mwjfnp0s1t0 =   2.96938239e+0_r8 * p001, &
      mwjfnp0s1t1 =  -7.23268813e-3_r8 * p001, &
      mwjfnp0s2t0 =   2.12382341e-3_r8 * p001, &
      mwjfnp1s0t0 =   1.04004591e-2_r8 * p001, &
      mwjfnp1s0t2 =   1.03970529e-7_r8 * p001, &
      mwjfnp1s1t0 =   5.18761880e-6_r8 * p001, &
      mwjfnp2s0t0 =  -3.24041825e-8_r8 * p001, &
      mwjfnp2s0t2 =  -1.23869360e-11_r8* p001

   !*** these constants will be used to construct the denominator

      real (kind=r8), parameter ::       &
      mwjfdp0s0t0 =   1.0e+0_r8,         &
      mwjfdp0s0t1 =   7.28606739e-3_r8,  &
      mwjfdp0s0t2 =  -4.60835542e-5_r8,  &
      mwjfdp0s0t3 =   3.68390573e-7_r8,  &
      mwjfdp0s0t4 =   1.80809186e-10_r8, &
      mwjfdp0s1t0 =   2.14691708e-3_r8,  &
      mwjfdp0s1t1 =  -9.27062484e-6_r8,  &
      mwjfdp0s1t3 =  -1.78343643e-10_r8, &
      mwjfdp0sqt0 =   4.76534122e-6_r8,  &
      mwjfdp0sqt2 =   1.63410736e-9_r8,  &
      mwjfdp1s0t0 =   5.30848875e-6_r8,  &
      mwjfdp2s0t3 =  -3.03175128e-16_r8, &
      mwjfdp3s0t1 =  -1.27934137e-17_r8
 

      integer (int_kind) ,intent(in) :: &
      k,                    &! depth level index
      kk                     ! level to which water is adiabatically displaced

      integer (int_kind),dimension(2880,2619),intent(in) :: K_2darray 

      real (r8), dimension(2880,2619) ,optional, intent(out) :: & 
      RHOOUT,  &! perturbation density of water
      RHOFULL, &! full density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity


      real (r8), dimension(2880,2619) ,intent(in) :: & 
      TEMPK,             &! temperature at level k/
      SALTK               ! salinity    at level k
      
      real (r8), dimension(2880,2619) :: & 
      WORK1, WORK2, WORK3, WORK4

      real(r8) twork1,twork2,twork3,twork4,TDENOMK,TQ,SQ,SQR 
        

      real (r8), dimension(km) :: & 
      tmin, tmax,        &! valid temperature range for level k
      smin, smax,        &! valid salinity    range for level k
      pressz              ! ref pressure (bars) at each level

     real (r8) :: p,p2, &! temporary pressure scalars 
      mwjfnums0t1, mwjfnums0t3,              &
      mwjfnums1t1, mwjfnums2t0,              &
      mwjfdens1t0, mwjfdens0t2, mwjfdens0t4, &
      mwjfdens1t1, mwjfdens1t3,                           &
      mwjfdensqt0, mwjfdensqt2
 
      real (r8), dimension(2880,2619) :: mwjfnums0t0,mwjfnums0t2,mwjfnums1t0
      real (r8), dimension(2880,2619) :: mwjfdens0t0,mwjfdens0t1,mwjfdens0t3
 
      tmin =  -2.0_r8  ! limited   on the low  end
      tmax = 999.0_r8  ! unlimited on the high end
      smin =   0.0_r8  ! limited   on the low  end
      smax = 0.999_r8  ! unlimited on the high end

      call random_number(pressz)

      do j=1,2619
      do i=1,2880 
      p = c10 * pressz( k_2darray(i,j) )
      mwjfnums0t0(i,j) = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0)
      mwjfnums0t2(i,j) = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2)
      mwjfnums1t0(i,j) = mwjfnp0s1t0 + p*mwjfnp1s1t0
      mwjfdens0t0(i,j) = mwjfdp0s0t0 + p*mwjfdp1s0t0
      mwjfdens0t1(i,j) = mwjfdp0s0t1 + p**3 * mwjfdp3s0t1
      mwjfdens0t3(i,j) = mwjfdp0s0t3 + p**2 * mwjfdp2s0t3
      enddo
      enddo 

      !***
      !*** first calculate numerator of MWJF density [P_1(S,T,p)]
      !***

      mwjfnums0t1 = mwjfnp0s0t1 
      mwjfnums0t3 = mwjfnp0s0t3
      mwjfnums1t1 = mwjfnp0s1t1
      mwjfnums2t0 = mwjfnp0s2t0
       
       
      mwjfdens0t2 = mwjfdp0s0t2
      mwjfdens0t4 = mwjfdp0s0t4
      mwjfdens1t0 = mwjfdp0s1t0
      mwjfdens1t1 = mwjfdp0s1t1
      mwjfdens1t3 = mwjfdp0s1t3
      mwjfdensqt0 = mwjfdp0sqt0
      mwjfdensqt2 = mwjfdp0sqt2 

      !dir$ offload begin target(mic:0)

        
      first_time = omp_get_wtime() 

      !dir$ assume_aligned SALTK: 64
      !dir$ assume_aligned RHOOUT: 64
      !dir$ assume_aligned RHOFULL: 64
      !dir$ assume_aligned DRHODT: 64 
      !dir$ assume_aligned RHOFULL: 64
      !dir$ assume_aligned DRHODS: 64
      !dir$ assume_aligned TEMPK: 64

      !$omp parallel default(none)shared(TEMPK,SALTK,WORK1,WORK2,WORK3,WORK4) &
      !$omp shared(RHOOUT,RHOFULL,DRHODS,DRHODT,mwjfnums0t0,mwjfnums0t2)&
      !$omp shared(mwjfnums1t0,mwjfdens0t0,mwjfdens0t1,mwjfdens0t3)&  
      !$omp firstprivate(tmax,tmin,smax,smin,kk) &
      !$omp firstprivate(mwjfnums0t1,mwjfnums0t3) &
      !$omp firstprivate(mwjfnums1t1,mwjfnums2t0) &
      !$omp firstprivate(mwjfdens0t2) &
      !$omp firstprivate(mwjfdens0t4,mwjfdens1t0,mwjfdens1t1,mwjfdens1t3) &
      !$omp firstprivate(mwjfdensqt0,mwjfdensqt2) &
      !$omp private(TWORK1,TWORK2,TWORK3,TWORK4,TDENOMK,TQ,SQ,SQR)

      !$omp do schedule(static)  
      do j=1,2619
      !dir$ simd
      !dir$ ivdep
      !dir$ vector nontemporal 
      do i=1,2880

      TQ = min(TEMPK(i,j),tmax(kk))
      TQ = max(TQ,tmin(kk))
      SQ = min(SALTK(i,j),smax(kk))
      SQ = max(SQ,smin(kk))
      SQ  = c1000*SQ
      SQR = sqrt(SQ)

      TWORK1 = mwjfnums0t0(i,j) + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2(i,j) + &
              mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0(i,j) +              &
              mwjfnums1t1 * TQ + mwjfnums2t0 * SQ )
      
      TWORK2 = mwjfdens0t0(i,j) + TQ * (mwjfdens0t1(i,j) + TQ * (mwjfdens0t2 +    &
           TQ * (mwjfdens0t3(i,j) + mwjfdens0t4 * TQ ))) +                   &
           SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ * TQ * mwjfdens1t3)+ &
           SQR * (mwjfdensqt0 + TQ * TQ * mwjfdensqt2))

      TDENOMK = c1/TWORK2

      if (present(RHOOUT)) then
          RHOOUT(i,j)  = TWORK1 * TDENOMK
      endif

      if (present(RHOFULL)) then
          RHOFULL(i,j) = TWORK1 * TDENOMK
      endif

      if (present(DRHODT)) then
         TWORK3 = &! dP_1/dT
                 mwjfnums0t1 + TQ * (c2*mwjfnums0t2(i,j) +    &
                 c3*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ

        TWORK4 = &! dP_2/dT
                 mwjfdens0t1(i,j) + SQ * mwjfdens1t1 +               &
                 TQ * (c2*(mwjfdens0t2 + SQ * SQR * mwjfdensqt2) +  &
                 TQ * (c3*(mwjfdens0t3(i,j) + SQ * mwjfdens1t3) +    &
                 TQ *  c4*mwjfdens0t4))
        
          DRHODT(i,j) = (TWORK3 - TWORK1 * TDENOMK * TWORK4)* TDENOMK
      endif

      if (present(DRHODS)) then
         TWORK3 = &! dP_1/dS
                 mwjfnums1t0(i,j) + mwjfnums1t1 * TQ + c2*mwjfnums2t0 * SQ

         TWORK4 = mwjfdens1t0 +   &! dP_2/dS
                 TQ * (mwjfdens1t1 + TQ * TQ * mwjfdens1t3) +   &
                 c1p5 * SQR *(mwjfdensqt0 + TQ * TQ * mwjfdensqt2)

         DRHODS(i,j) = ( TWORK3 - TWORK1 * TDENOMK * TWORK4 ) * TDENOMK * c1000
     endif  
       enddo
       enddo
       !$omp end do

      !$omp end parallel         
      end_time = omp_get_wtime()

     !dir$ end offload   
 
      print *,"total time is",end_time - first_time

      print *, sum(WORK1)
      print *, sum(WORK2)
      print *, sum(WORK3)
      print *, sum(WORK4)
      print *, sum(DRHODS)
      print *, sum(DRHODT)
      print *, sum(RHOOUT)
      print *, sum(RHOFULL)
       
end subroutine state

end program state_dummy
