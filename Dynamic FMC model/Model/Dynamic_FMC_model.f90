! --------------------------------------------------
! Fuel moisture model, Matthews model / Koba model
! --------------------------------------------------

PROGRAM MoistureModel
IMPLICIT NONE

INTEGER :: i
INTEGER :: j
INTEGER :: k
INTEGER :: local_index
INTEGER, parameter :: Inp_N_layers = 5;
INTEGER, parameter :: n_data = 3106;
doubleprecision, allocatable, dimension(:,:) :: q, T_m, mo,	T_a, liq, T_l

doubleprecision :: temp_a_all(n_data,1)
doubleprecision :: DNI_all(n_data,1)
doubleprecision :: RelH_all(n_data,1)
doubleprecision :: pressure_all(n_data,1)
doubleprecision :: wind_all(n_data,1)
doubleprecision :: rain_all(n_data,1)
doubleprecision :: T_s(n_data,1)

doubleprecision :: temp_a
doubleprecision :: DNI
doubleprecision :: RelH
doubleprecision :: pressure
doubleprecision :: wind
doubleprecision :: rain

doubleprecision :: f1, f2, f3, f4, f5, f6
doubleprecision, allocatable, dimension(:,:) :: m_old, q_old, T_a_old, T_m_old, liq_old, T_liq_old
doubleprecision, allocatable, dimension(:,:) :: F, F1_new
doubleprecision :: f_m1, f_q1, f_Ta1, f_Tm1, f_li1, f_li2, f_tliq1, f_tliq2
doubleprecision :: f_m2, f_q2, f_Ta2, f_Tm2, f_m3, f_q3, f_Ta3, f_Tm3
doubleprecision :: f_m4, f_q4, f_Ta4, f_Tm4, f_m5, f_q5, f_Ta5, f_Tm5
doubleprecision, allocatable, dimension(:,:):: f_old, f_old1
doubleprecision, allocatable, dimension(:,:):: f_new, f_new1
doubleprecision :: dt, dh1, dh2, dh3, dh4
integer :: data_num, it
doubleprecision, allocatable, dimension(:,:) :: itera_m, itera_q, itera_ta, &
                        itera_tm, itera_l, itera_tl
doubleprecision, allocatable, dimension(:,:) :: al, Y_old, Y_guess, Y_new
doubleprecision :: tol, error1
doubleprecision, allocatable, dimension(:,:) :: err_mat
doubleprecision, allocatable, dimension(:,:) :: m_fin, q_fin, T_a_fin, T_m_fin, liq_fin, T_l_fin
doubleprecision :: dfm1, dq1, dTa1, dTm1, dfm2, dq2, dTa2, dliq1, dTl1, dTm2, dliq2, &
                        dTl2, dfm3, dq3, dTa3, dTm3, dliq3, dTl3
doubleprecision :: dfm4, dq4, dTa4, dTm4, dfm5, dq5, dTa5, dTm5
doubleprecision, allocatable, dimension(:,:) :: err_m, err_q, err_ta, err_tm, err_l, err_tl
doubleprecision :: max_itera
doubleprecision, allocatable, DIMENSION(:,:) :: q_sat_check

! Ambient weather condition inputs
open (unit=1, file='Air_temperature.csv')
read(1,*) temp_a_all
close (1)

open (unit=2, file='Solar_radiation.csv')
read(2,*) DNI_all
close (2)

open (unit=3, file='Relative_humidity.csv')
read(3,*) RelH_all
close (3)

open (unit=4, file='Pressure.csv')
read(4,*) pressure_all
close (4)

open (unit=5, file='Wind_speed.csv')
read(5,*) wind_all
close (5)

open (unit=7, file='Rainfall.csv')
read(7,*) rain_all
close (7)

dt = 112.5;
data_num = 1;

dh1 = 0.000001;
dh2 = 0.000001;
dh3 = 0.00001;
dh4 = 0.00001;

allocate(q(Inp_N_layers,1))
allocate(T_m(Inp_N_layers,1))
allocate(mo(Inp_N_layers,1))
allocate(T_a(Inp_N_layers,1))
allocate(liq(Inp_N_layers,1))
allocate(T_l(Inp_N_layers,1))

allocate(m_old(Inp_N_layers,1))
allocate(q_old(Inp_N_layers,1))
allocate(T_a_old(Inp_N_layers,1))
allocate(T_m_old(Inp_N_layers,1))
allocate(liq_old(Inp_N_layers,1))
allocate(T_liq_old(Inp_N_layers,1))

allocate(F(6,Inp_N_layers))
allocate(F1_new(6,Inp_N_layers))

allocate(f_old(6,Inp_N_layers))
allocate(f_old1(6,Inp_N_layers))
allocate(f_new(6,Inp_N_layers))
allocate(f_new1(6,Inp_N_layers))

allocate(itera_m(Inp_N_layers,1))
allocate(itera_q(Inp_N_layers,1))
allocate(itera_ta(Inp_N_layers,1))
allocate(itera_tm(Inp_N_layers,1))
allocate(itera_l(Inp_N_layers,1))
allocate(itera_tl(Inp_N_layers,1))

allocate(al(6,Inp_N_layers))
allocate(Y_old(6,Inp_N_layers))
allocate(Y_guess(6,Inp_N_layers))
allocate(Y_new(6,Inp_N_layers))

allocate(err_mat(6,Inp_N_layers))

allocate(m_fin(n_data,Inp_N_layers))
allocate(q_fin(n_data,Inp_N_layers))
allocate(T_a_fin(n_data,Inp_N_layers))
allocate(T_m_fin(n_data,Inp_N_layers))
allocate(liq_fin(n_data,Inp_N_layers))
allocate(T_l_fin(n_data,Inp_N_layers))

allocate(err_m(Inp_N_layers,1))
allocate(err_q(Inp_N_layers,1))
allocate(err_ta(Inp_N_layers,1))
allocate(err_tm(Inp_N_layers,1))
allocate(err_l(Inp_N_layers,1))
allocate(err_tl(Inp_N_layers,1))

allocate(q_sat_check(Inp_N_layers,1))

temp_a = temp_a_all(data_num,1);
DNI = DNI_all(data_num,1);
RelH = RelH_all(data_num,1);
wind = abs(wind_all(data_num,1));
pressure = pressure_all(data_num,1);
rain = rain_all(data_num,1);

do i = 1,Inp_N_layers
mo(i,1) = 0.03;
q(i,1) = 0.03;
T_a(i,1) = 290;
T_m(i,1) = 285 ;
liq(i,1) = 0;
T_l(i,1) = 290;
end do

mo(1,1) = 0.02;
T_m = 287;
m_old = mo;
q_old = q;
T_a_old = T_a;
T_m_old = T_m;
liq_old = liq;
T_liq_old = T_l;

max_itera =1000;

it = 1;
!tol = 0.5*(10**(-1));
tol = (10**(-3))

! Newton Raphson Loop for Equations 11-16
do while (it <= n_data)

    error1 = 1;

    do i = 1,Inp_N_layers
    Y_old(1,i) = m_old(i,1);
    Y_old(2,i) = q_old(i,1);
    Y_old(3,i) = T_a_old(i,1);
    Y_old(4,i) = T_m_old(i,1);
    Y_old(5,i) = liq_old(i,1);
    Y_old(6,i) = T_m_old(i,1);

    Y_guess(1,i) = m_old(i,1);
    Y_guess(2,i) = q_old(i,1);
    Y_guess(3,i) = T_a_old(i,1);
    Y_guess(4,i) = T_m_old(i,1);
    Y_guess(5,i) = liq_old(i,1);
    Y_guess(6,i) = T_m_old(i,1);
    end do

tol = (10**(-2));
data_num = it;

do i = 1,Inp_N_layers
err_m(i,1) = 1;
err_q(i,1) = 1;
err_ta(i,1) = 1;
err_tm(i,1) = 1;
err_l(i,1) = 1;
err_tl(i,1) = 1;

itera_m(i,1) = 1;
itera_q(i,1) = 1;
itera_ta(i,1) = 1;
itera_tm(i,1) = 1;
itera_l(i,1) = 1;
itera_tl(i,1) = 1;
end do

do i = 1,6
    do j = 1,Inp_N_layers
		al(i,j) = 1;
    end do
end do

temp_a = temp_a_all(data_num,1);
DNI = DNI_all(data_num,1);
RelH = RelH_all(data_num,1);
wind = abs(wind_all(data_num,1));
pressure = pressure_all(data_num,1);
rain = rain_all(data_num,1);
T_s = temp_a_all(data_num,1);

		tol = 0.001;
		data_num = it;

			error1 = 1;

			do i = 1,Inp_N_layers

				Y_old(1,i) = m_old(i,1);
				Y_old(2,i) = q_old(i,1);
				Y_old(3,i) = T_a_old(i,1);
				Y_old(4,i) = T_m_old(i,1);
				Y_old(5,i) = liq_old(i,1);
			!    Y_old(6,1) = T_liq_old;
				Y_old(6,i) = T_m_old(i,1);

				Y_guess(1,i) = m_old(i,1);
				Y_guess(2,i) = q_old(i,1);
				Y_guess(3,i) = T_a_old(i,1);
				Y_guess(4,i) = T_m_old(i,1);
				Y_guess(5,i) = liq_old(i,1);
			!    Y_guess(6,1) = T_liq_old;
				Y_guess(6,i) = T_m_old(i,1);
			end do

			do i = 1,Inp_N_layers
				err_m(i,1) = 1;
				err_q(i,1) = 1;
				err_ta(i,1) = 1;
				err_tm(i,1) = 1;
				err_l(i,1) = 1;
				err_tl(i,1) = 1;

				itera_m(i,1) = 1;
				itera_q(i,1) = 1;
				itera_ta(i,1) = 1;
				itera_tm(i,1) = 1;
				itera_l(i,1) = 1;
				itera_tl(i,1) = 1;
			end do

			do i = 1,6
				do j = 1,Inp_N_layers
					al(i,j) = 1;
				end do
			end do

			do local_index = 1,Inp_N_layers

				do while (err_q(local_index,1) > tol)

					do i = 1,Inp_N_layers
						mo(i,1) = Y_guess(1,i);
						q(i,1) = Y_guess(2,i);
						T_a(i,1) = Y_guess(3,i);
						T_m(i,1) = Y_guess(4,i);
						liq(i,1) = Y_guess(5,i);
					!    T_l = Y_guess(6,1);
						T_l(i,1) = Y_guess(4,i);
					end do

					call vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,wind,DNI, &
					temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,F);

					q(local_index,1) = q(local_index,1) + dh2;
					call f_q_calc_1a(mo,q,T_a,T_m,q_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_q1);
					q(local_index,1) = q(local_index,1) - dh2;

					call f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,wind,DNI, &
					temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,f_old);

					dq1 = (f_q1-f_old(2,local_index))/dh2;

					if (dq1 == 0) then

						Y_new(2,local_index) = Y_guess(2,local_index);

					else

						Y_new(2,local_index) = Y_guess(2,local_index) - al(2,local_index)*F(2,local_index)/dq1;

					endif

					q_sat_check(local_index,1) = (0.622*0.611*(exp((17.3*(T_a(local_index,1)-273))/ &
                            (T_a(local_index,1)-273+237.3))))/(pressure*0.001);

					if (Y_new(2,local_index) .GT. q_sat_check(local_index,1)) then

					Y_new(2,local_index) = q_sat_check(local_index,1);

					elseif (Y_new(2,local_index) .LT. 1.E-4) then

					Y_new(2,local_index) = 0.0001

					else

					Y_new(2,local_index) = Y_new(2,local_index);

					endif

					q(local_index,1) = Y_new(2,local_index);

					err_q(local_index,1) = abs(Y_new(2,local_index)-Y_guess(2,local_index));

					Y_guess(2,local_index) = Y_new(2,local_index);

					if (itera_q(local_index,1) .GE. max_itera) exit

					itera_q(local_index,1) = itera_q(local_index,1) + 1;

				end do

			end do

			do local_index = 1,Inp_N_layers

				do while (err_ta(local_index,1) > tol)

					do i = 1,Inp_N_layers
						mo(i,1) = Y_guess(1,i);
						q(i,1) = Y_guess(2,i);
						T_a(i,1) = Y_guess(3,i);
						T_m(i,1) = Y_guess(4,i);
						liq(i,1) = Y_guess(5,i);
					!    T_l = Y_guess(6,1);
						T_l(i,1) = Y_guess(4,i);
					end do

					call vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt, &
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,F);

					T_a(local_index,1) = T_a(local_index,1) + dh3;
					call f_Ta_calc_1a(mo,q,T_a,T_m,T_a_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_Ta1);
					T_a(local_index,1) = T_a(local_index,1) - dh3;

					call f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt, &
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,f_old);

					dTa1 = (f_Ta1-f_old(3,local_index))/dh3;

					if (dTa1 == 0) then

						Y_new(3,local_index) = Y_guess(3,local_index);

					else

						Y_new(3,local_index) = Y_guess(3,local_index) - al(3,local_index)*F(3,local_index)/dTa1;

					endif

					T_a(local_index,1) = Y_new(3,local_index);

					err_Ta(local_index,1) = abs(Y_new(3,local_index)-Y_guess(3,local_index));

					Y_guess(3,local_index) = Y_new(3,local_index);

					if (itera_ta(local_index,1) .GE. max_itera) exit

					itera_ta(local_index,1) = itera_ta(local_index,1) + 1;

				end do

			!print*, 't_a', local_index, err_ta(local_index,1), itera_ta(local_index,1)
			!print*, 't_a', T_a(local_index,1), dTa1

				! if (Y_new(3,local_index) .LE. 283) then
					! Y_new(3,local_index) = 283
					! Y_guess(3,local_index) = Y_new(3,local_index)
					! T_a(local_index,1) = Y_new(3,local_index);
				! end if

			!	print*, 'ta', local_index, Y_new(3, local_index), itera_ta(local_index,1)

			end do

			do local_index = 1,Inp_N_layers

				do while (err_tm(local_index,1) > tol)

					do i = 1,Inp_N_layers
						mo(i,1) = Y_guess(1,i);
						q(i,1) = Y_guess(2,i);
						T_a(i,1) = Y_guess(3,i);
						T_m(i,1) = Y_guess(4,i);
						liq(i,1) = Y_guess(5,i);
					!    T_l = Y_guess(6,1);
						T_l(i,1) = Y_guess(4,i);
					end do

					call vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,F); ! Additionally output E_ma(T_m)

					T_m(local_index,1) = T_m(local_index,1) + dh4; 				! T_m + delta_Tm
					call f_Tm_calc_1a(mo,q,T_a,T_m,T_m_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_Tm1); ! f(T_m + delta_Tm)
					T_m(local_index,1) = T_m(local_index,1) - dh4;				! output E_ma(T_m + delta_Tm)

					call f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,f_old);

					dTm1 = (f_Tm1-f_old(4,local_index))/dh4;

					if (dTm1 == 0) then

						Y_new(4,local_index) = Y_guess(4,local_index);

					else

						Y_new(4,local_index) = Y_guess(4,local_index) - al(4,local_index)*F(4,local_index)/dTm1;
						!print*, Y_new(4,1), Y_guess(4,1), al(4,1), F(4,1), dTm1
					endif

					T_m(local_index,1) = Y_new(4,local_index);

					err_Tm(local_index,1) = abs(Y_new(4,local_index)-Y_guess(4,local_index));

					Y_guess(4,local_index) = Y_new(4,local_index);

					if (Y_new(4,local_index) .LE. 275.0) then
						Y_new(4,local_index) = 275.0
						Y_guess(4,local_index) = Y_new(4,local_index);
						T_m(local_index,1) = Y_new(4,local_index);
					endif

					if (itera_tm(local_index,1) .GE. max_itera) exit

					itera_tm(local_index,1) = itera_tm(local_index,1) + 1;


				end do

			end do

			do local_index = 1,Inp_N_layers

				do while (err_m(local_index,1) > tol)

					do i = 1,Inp_N_layers
						mo(i,1) = Y_guess(1,i);
						q(i,1) = Y_guess(2,i);
						T_a(i,1) = Y_guess(3,i);
						T_m(i,1) = Y_guess(4,i);
						liq(i,1) = Y_guess(5,i);
					!    T_l = Y_guess(6,1);
						T_l(i,1) = Y_guess(4,i);
					end do

					call vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,F);

					mo(local_index,1) = mo(local_index,1) + dh1;
					call f_m_calc_1a(mo,q,T_a,T_m,m_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_m1);
					mo(local_index,1) = mo(local_index,1) - dh1;

					call f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
					wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,f_old);

					dfm1 = (f_m1-f_old(1,local_index))/dh1;

					if (dfm1 == 0) then

						Y_new(1,local_index) = Y_guess(1,local_index);

					else

						Y_new(1,local_index) = Y_guess(1,local_index) - al(1,1)*F(1,local_index)/dfm1;

					endif

					mo(local_index,1) = Y_new(1,local_index);

					err_m(local_index,1) = abs(Y_new(1,local_index)-Y_guess(1,local_index));

					Y_guess(1,local_index) = Y_new(1,local_index);

					if (Y_new(1,local_index) .LE. 1.E-4) then
						Y_new(1,local_index) = 0.0001
						Y_guess(1,local_index) = Y_new(1,local_index)
						mo(local_index,1) = Y_new(1,local_index);
					end if

					if (itera_m(local_index,1) .GE. max_itera) exit

					itera_m(local_index,1) = itera_m(local_index,1) + 1;

				end do

			end do

			call vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,wind,DNI,temp_a,pressure,RelH,rain,&
							T_s,INP_N_layers,local_index,F);

			call f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
						wind,DNI,temp_a,pressure,RelH,rain,T_s,INP_N_layers,local_index,f_old);

			do local_index = 1,Inp_N_layers

				mo(local_index,1) = mo(local_index,1) + dh1;
				call f_m_calc_1a(mo,q,T_a,T_m,m_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_m1);
				mo(local_index,1) = mo(local_index,1) - dh1;

				dfm1 = (f_m1-f_old(1,local_index))/dh1;

				Y_new(1,local_index) = Y_guess(1,local_index) - al(1,local_index)*F(1,local_index)/dfm1;

				if (Y_new(1,local_index) .LE. 1.E-4) then
					Y_new(1,local_index) = 0.0001
					Y_guess(1,local_index) = Y_new(1,local_index)
					mo(local_index,1) = Y_new(1,local_index);
				end if
			end do

			do local_index = 1,Inp_N_layers

				q(local_index,1) = q(local_index,1) + dh2;
				call f_q_calc_1a(mo,q,T_a,T_m,q_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_q1);
				q(local_index,1) = q(local_index,1) - dh2;

				dq1 = (f_q1-f_old(2,local_index))/dh2;

				Y_new(2,local_index) = Y_guess(2,local_index) - al(2,local_index)*F(2,local_index)/dq1;
				q_sat_check(local_index,1) = (0.622*0.611*(exp((17.3*(T_a(local_index,1)-273)) &
				/(T_a(local_index,1)-273+237.3))))/(pressure*0.001);

				if (Y_new(2,local_index) .GT. q_sat_check(local_index,1)) then

					Y_new(2,local_index) = q_sat_check(local_index,1);

				elseif (Y_new(2,local_index) .LT. 0) then

					Y_new(2,local_index) = 0.001

				else

					Y_new(2,local_index) = Y_new(2,local_index);

				endif
			end do

			do local_index = 1,Inp_N_layers

				T_a(local_index,1) = T_a(local_index,1) + dh3;
				call f_Ta_calc_1a(mo,q,T_a,T_m,T_a_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_Ta1);
				T_a(local_index,1) = T_a(local_index,1) - dh3;

				dTa1 = (f_Ta1-f_old(3,local_index))/dh3;

				Y_new(3,local_index) = Y_guess(3,local_index) - al(3,local_index)*F(3,local_index)/dTa1;

				! if (Y_new(3,local_index) .LE. 283) then
					! Y_new(3,local_index) = 283
					! Y_guess(3,local_index) = Y_new(3,local_index)
					! T_a(local_index,1) = Y_new(3,local_index);
				! end if
			end do

			do local_index = 1,Inp_N_layers

				T_m(local_index,1) = T_m(local_index,1) + dh4;
				call f_Tm_calc_1a(mo,q,T_a,T_m,T_m_old,liq,T_l,dt,wind,DNI,temp_a,pressure,RelH,rain,T_s,local_index,INP_N_layers,f_Tm1);
				T_m(local_index,1) = T_m(local_index,1) - dh4;

				dTm1 = (f_Tm1-f_old(4,local_index))/dh4;

				Y_new(4,local_index) = Y_guess(4,local_index) - al(4,local_index)*F(4,local_index)/dTm1;

				if (Y_new(4,local_index) .LE. 275.0) then
					Y_new(4,local_index) = 275.0
					Y_guess(4,local_index) = Y_new(4,local_index);
					T_m(local_index,1) = Y_new(4,local_index);
				endif

			end do

			do i = 1,Inp_N_layers

				m_old(i,1) = Y_new(1,i);
				q_old(i,1) = Y_new(2,i);
				T_a_old(i,1) = Y_new(3,i);
				T_m_old(i,1) = Y_new(4,i);
				liq_old(i,1) = Y_new(5,i);
			!    T_liq_old = Y_new(6,1);
				T_liq_old(i,1) = Y_new(4,i);

				m_fin(it,i) = Y_new(1,i);
				q_fin(it,i) = Y_new(2,i);
				T_a_fin(it,i) = Y_new(3,i);
				T_m_fin(it,i) = Y_new(4,i);
				liq_fin(it,i) = Y_new(5,i);
				T_l_fin(it,i) = Y_new(6,i);

			end do

    it = it + 1;

end do


!! Write Data

! Fuel Moisture Content
open (unit=10,file="m.csv",action="write",status="replace")
do i=1,n_data
write (10,"(32(f0.6,',',:))") m_fin(i,:)
end do

! Specific Humidity
open (unit=15,file="q.csv",action="write",status="replace")
do i=1,n_data
write (15,"(32(f0.6,',',:))") q_fin(i,:)
end do

! Air Temperature
open (unit=12,file="Ta.csv",action="write",status="replace")
do i=1,n_data
write (12,"(32(f0.6,',',:))") T_a_fin(i,:)
end do

! Fuel Temperature
open (unit=13,file="Tm.csv",action="write",status="replace")
do i=1,n_data
write (13,"(32(f0.6,',',:))") T_m_fin(i,:)
end do

! Liquid water
open (unit=16,file="liq.csv",action="write",status="replace")
do i=1,n_data
write (16,"(32(f0.6,',',:))") liq_fin(i,:)
end do

! Liquid water temperature
open (unit=17,file="Tl.csv",action="write",status="replace")
do i=1,n_data
write (17,"(32(f0.6,',',:))") T_l_fin(i,:)
end do

END PROGRAM MoistureModel

subroutine f_m_calc_1a(mo,q,T_a,T_m,m_old,liq,T_l,dt,wind,DNI,temp_a,P_INF,RelH,rain,T_s,local_index,INP_N_layers,f_m1)
IMPLICIT NONE

doubleprecision, intent(in) :: mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
									m_old(Inp_N_layers,1), liq(Inp_N_layers,1), T_l(Inp_N_layers,1)
doubleprecision, INTENT(IN) :: DT
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH, rain, T_s
INTEGER, INTENT(IN) :: local_index,INP_N_layers
doubleprecision, intent(out) :: f_m1

doubleprecision :: a_con, b_con, q_sat, RH_fuel, calc, V_a, V_m
doubleprecision :: D_T0, zeta, RiB, zama_h, gamma_h, q0, dz

INTEGER :: i, j

doubleprecision :: E_ma, E_T0, E_T, H_T0, H_T, H_ma, R_net0, mu_fmc_ma1, K_T0, temp_aK

doubleprecision :: x(Inp_N_layers)

doubleprecision :: f1

doubleprecision :: H_ml, H_la, E_la, E_ml, K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1, &
					S_lma, S_ml, S_ma, S_la, Drain, D_p, C_hl, V_l
doubleprecision :: check(2,1)
doubleprecision :: mo_local, T_m_local, q_local, T_a_local, T_l_local, l_local
!!!! Moisture model parameter
doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
doubleprecision, PARAMETER :: T_ref = 273.16

doubleprecision, PARAMETER :: rho_air = 1.2

doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann

!< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;

dz = hei/(Inp_N_layers)

zama_h = log10(z_screen/z_0)

check(1,1) = (liq(1,1)/(D_s*rho_bulk));
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;

mu_fmc_ma = S_ma/V_m;

mu_fmc_ma1 = S_ma/V_a;

mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

calc = exp((mo(local_index,1)*B_fmc) + A_fmc);

RH_fuel = exp((-4.19*Mass_fmc*calc)/((R0/1000)*T_m(local_index,1)));
q_sat = (0.622*0.611*(exp((17.3*(T_m(local_index,1)-273))/(T_m(local_index,1)-273+237.3))))/(P_INF*0.001);
q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);
E_ml = 0;

E_ma = rho_air*K_maE*(RH_fuel*q_sat-q(local_index,1));

f1 = - m_old(local_index,1) + (dt/rho_litter)*(mu_fmc_ma*E_ma) + mo(local_index,1);

f_m1 = f1;

end subroutine f_m_calc_1a

subroutine f_q_calc_1a(mo,q,T_a,T_m,q_old,liq,T_l,dt,wind,DNI,temp_a,P_INF,RelH,rain,T_s,local_index,INP_N_layers,f_q1)
IMPLICIT NONE

doubleprecision, intent(in) :: mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
									q_old(Inp_N_layers,1), liq(Inp_N_layers,1), T_l(Inp_N_layers,1)
!integer, intent(in) :: data_num
doubleprecision, INTENT(IN) :: DT
INTEGER, INTENT(IN) :: local_index,INP_N_layers
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH, rain, T_s

doubleprecision, intent(out) :: f_q1

doubleprecision :: a_con, b_con, q_sat, RH_fuel, calc, V_a


doubleprecision :: D_T0, zeta, V_m, RiB, zama_h, gamma_h, q0, dz

INTEGER :: i, j

doubleprecision :: S_dn, L_up0, S_up0, tau_i, S_up, L_dn, L_up, R_net, K_maH, H_c

doubleprecision :: E_ma, E_T0, E_T, H_T0, H_T, H_ma, R_net0, mu_fmc_ma1, &
					K_T0, temp_aK, RH_soil, q_soil, theta

doubleprecision :: x(Inp_N_layers)
doubleprecision :: D_T(Inp_N_layers,1)

doubleprecision :: f2
doubleprecision :: H_ml, H_la, E_la, E_ml, K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1,S_lma, &
					S_ml, S_ma, S_la, Drain, D_p, C_hl, V_l
doubleprecision :: check(2,1)

doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
							! universal gas constant
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
doubleprecision, PARAMETER :: T_ref = 273.16
doubleprecision, PARAMETER :: rho_air = 1.2
doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann
!REAL(EB), PARAMETER :: LAI = 1._EB				                 !< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;

!tau_i = EXP(gama*h)

dz = hei/(INP_N_layers)
zama_h = log10(z_screen/z_0)

D_T0 = D_T0a*exp(wind*D_T0b);
zeta = zeta_a + (wind*zeta_b);
K_maH = 3.23*(10**(-3))*(wind/d_fmc)**(0.5);

do i = 1,Inp_N_layers
	if (i == 1) then
		x(i) = dz/2.0;
	else
		x(i) = x(1) + (i-1)*dz;
	endif
		D_T(i,1) = D_T0*exp(zeta*((x(i)/hei)-1));
end do

!check(1,1) = (liq/(D_s*rho_bulk));
check(1,1) = 0;
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;

mu_fmc_ma = S_ma/V_m;
mu_fmc_ma1 = S_ma/V_a;
mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

K_laE = (kappa_h/kappa_v)**(0.67);

RiB = (z_screen*9.81*(temp_aK-T_a(1,1)))/((wind**(2))*temp_aK);
gamma_h = (zama_h-10*RiB*zama_h-zama_h)/(2-10*RiB);
K_T0 = (((K_von)**(2))*wind)/((zama_h-gamma_h)*(zama_h-gamma_h));

calc = exp((mo(local_index,1)*B_fmc) + A_fmc);
RH_fuel = exp((-4.19*Mass_fmc*calc)/((R0/1000)*T_m(local_index,1)));
q_sat = (0.622*0.611*(exp((17.3*(T_m(local_index,1)-273))/(T_m(local_index,1)-273+237.3))))/(P_INF*0.001);
q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

E_la = 0;

E_ma = rho_air*K_maE*(RH_fuel*q_sat-q(local_index,1));

if (local_index == 1) then
	E_T0 = rho_air*K_T0*(q(local_index,1)-q0);
	E_T = rho_air*D_T0*exp(zeta*(((hei/((Inp_N_layers)*hei))-1)))*(q(local_index+1,1)-q(local_index,1))/(dz);

	f2 = - q_old(local_index,1) - (dt/rho_air)*(((E_T-E_T0)/(V_a*0.5*dz)) + mu_fmc_ma1*E_ma) + q(local_index,1);

elseif (local_index == Inp_N_layers) then

	E_T0 = rho_air*D_T0*exp(zeta*((((local_index-1)*hei)/((Inp_N_layers)*hei))-1))*(q(local_index,1)-q(local_index-1,1))/dz;
	!E_T = -(10**(-6));
	E_T = rho_air*D_T0*exp(zeta*((((local_index-1)*hei)/((Inp_N_layers)*hei))-1))*((0.0001)-q(local_index,1))/(dz/2);
	f2 = - q_old(local_index,1) - (dt/rho_air)*(((E_T-E_T0)/(V_a*dz)) + mu_fmc_ma1*E_ma) + q(local_index,1);

else
	E_T0 = rho_air*D_T0*exp(zeta*((((local_index-1)*hei)/((Inp_N_layers)*hei))-1))*(q(local_index,1)-q(local_index-1,1))/dz;
	E_T = rho_air*D_T0*exp(zeta*((((local_index)*hei)/((Inp_N_layers)*hei))-1))*(q(local_index+1,1)-q(local_index,1))/dz;

	f2 = - q_old(local_index,1) - (dt/rho_air)*(((E_T-E_T0)/(V_a*dz)) + mu_fmc_ma1*E_ma) + q(local_index,1);

end if

!print*, 'local_index, E_T0, E_T ', local_index, E_T0, E_T

f_q1 = f2;

END subroutine f_q_calc_1a

subroutine f_Ta_calc_1a(mo,q,T_a,T_m,T_a_old,liq,T_l,dt,wind,DNI,temp_a,P_INF,RelH,rain,T_s,local_index,INP_N_layers,f_Ta1)
IMPLICIT NONE

doubleprecision, intent(in) ::  mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
									T_a_old(Inp_N_layers,1), liq(Inp_N_layers,1), T_l(Inp_N_layers,1)
!integer, intent(in) :: data_num
doubleprecision, INTENT(IN) :: DT
INTEGER, INTENT(IN) :: local_index,INP_N_layers
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH, rain, T_s
doubleprecision, intent(out) :: f_Ta1

doubleprecision :: a_con, b_con, q_sat, RH_fuel, calc, V_a
doubleprecision :: D_T0
doubleprecision :: zeta

INTEGER :: i
INTEGER :: j

doubleprecision :: V_m
doubleprecision :: RiB
doubleprecision :: zama_h
doubleprecision :: gamma_h
doubleprecision :: q0
doubleprecision :: dz
doubleprecision :: S_dn, L_up0, S_up0, tau_i, S_up, L_dn, L_up, R_net, K_maH, H_c
doubleprecision :: E_ma, E_T0, E_T, H_T0, H_T, H_ma, R_net0, mu_fmc_ma1, K_T0, temp_aK
doubleprecision :: x(Inp_N_layers), D_T(Inp_N_layers,1)
doubleprecision :: f3
doubleprecision :: H_ml, H_la, E_la, E_ml, K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1,S_lma, S_ml, S_ma, S_la, &
					Drain, D_p, C_hl, V_l
doubleprecision :: check(2,1)
!!!! Moisture model parameter
doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
!REAL(EB), PARAMETER :: R = 8.314								! universal gas constant
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
!REAL(EB), PARAMETER :: P = 100000.0
doubleprecision, PARAMETER :: T_ref = 273.16

doubleprecision, PARAMETER :: rho_air = 1.2

doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
!REAL(EB), PARAMETER :: sigma_st = 5.67/(10**(8))
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann
!REAL(EB), PARAMETER :: LAI = 1._EB				                 !< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;

!tau_i = EXP(gama*h)

dz = hei/(INP_N_layers)
zama_h = log10(z_screen/z_0)

D_T0 = D_T0a*exp(wind*D_T0b);
zeta = zeta_a + (wind*zeta_b);
K_maH = 3.23*(10**(-3))*(wind/d_fmc)**(0.5);

do i = 1,Inp_N_layers
	if (i == 1) then
		x(i) = dz/2.0;
	else
		x(i) = x(1) + (i-1)*dz;
	endif
		D_T(i,1) = D_T0*exp(zeta*((x(i)/hei)-1));
end do

check(1,1) = 0;
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;

mu_fmc_ma = S_ma/V_m;
mu_fmc_ma1 = S_ma/V_a;
mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

RiB = (z_screen*9.81*(temp_aK-T_a(1,1)))/((wind**(2))*temp_aK);
gamma_h = (zama_h-10.0*RiB*zama_h-zama_h)/(2-10.0*RiB);
K_T0 = (((K_von)**(2))*wind)/((zama_h-gamma_h)*(zama_h-gamma_h));

K_laH = K_maH;
K_laE = (kappa_h/kappa_v)**(0.67);

calc = exp((mo(1,1)*B_fmc) + A_fmc);
RH_fuel = exp((-4.19*Mass_fmc*calc)/((R0/1000)*T_m(1,1)));
q_sat = (0.622*0.611*(exp((17.3*(T_m(1,1)-273))/(T_m(1,1)-273+237.3))))/(P_INF*0.001);

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

H_la = 0;

if (local_index == 1) then
	H_T0 = rho_air*C_p*K_T0*(T_a(local_index,1)-temp_aK);
	H_T = rho_air*C_p*D_T0*exp(zeta*(((hei)/((Inp_N_layers)*hei))-1))*(T_a(local_index+1,1)-T_a(local_index,1))/dz;
	H_ma = rho_air*C_p*K_maH*(T_m(local_index,1)-T_a(local_index,1));

	f3 = - T_a_old(local_index,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T-H_T0)/(0.5*dz) + &
        mu_fmc_ma1*H_ma + mu_fmc_la1*H_la) + T_a(local_index,1);
else if (local_index == Inp_N_layers) then
	H_T0 = rho_air*C_p*D_T0*exp(zeta*((((local_index-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(local_index,1)-T_a(local_index-1,1))/dz;
	H_T = rho_air*C_p*D_T0*exp(zeta*((((local_index)*hei)/((Inp_N_layers)*hei))-1))*(T_s+273-T_a(local_index,1))/(0.5*dz);	! Should be dz/2
	!H_T = 0;
	H_ma = rho_air*C_p*K_maH*(T_m(local_index,1)-T_a(local_index,1));

	f3 = - T_a_old(local_index,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T-H_T0)/dz + &
        mu_fmc_ma1*H_ma + mu_fmc_la1*H_la) + T_a(local_index,1);
else
	H_T0 = rho_air*C_p*D_T0*exp(zeta*((((local_index-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(local_index,1)-T_a(local_index-1,1))/dz;
	H_T = rho_air*C_p*D_T0*exp(zeta*((((local_index)*hei)/((Inp_N_layers)*hei))-1))*(T_a(local_index+1,1)-T_a(local_index,1))/dz;
	H_ma = rho_air*C_p*K_maH*(T_m(local_index,1)-T_a(local_index,1));

	f3 = - T_a_old(local_index,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T-H_T0)/dz + &
        mu_fmc_ma1*H_ma + mu_fmc_la1*H_la) + T_a(local_index,1);
end if

f_Ta1 = f3;

END subroutine f_Ta_calc_1a

subroutine f_Tm_calc_1a(mo,q,T_a,T_m,T_m_old,liq,T_l,dt,wind,DNI,temp_a,P_INF,RelH,rain,T_s,local_index,INP_N_layers,f_Tm1)
IMPLICIT NONE

doubleprecision, intent(in) ::  mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
									T_m_old(Inp_N_layers,1), liq(Inp_N_layers,1), T_l(Inp_N_layers,1)
!integer, intent(in) :: data_num
doubleprecision, INTENT(IN) :: DT
INTEGER, INTENT(IN) :: local_index,INP_N_layers
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH, rain, T_s
doubleprecision, intent(out) :: f_Tm1

doubleprecision :: a_con, b_con, q_sat(Inp_N_layers,1), RH_fuel(Inp_N_layers,1), calc(Inp_N_layers,1), V_a

doubleprecision :: D_T0, zeta, V_m, RiB, zama_h, gamma_h, q0, dz

INTEGER :: i, j

doubleprecision :: S_dn(Inp_N_layers+1,1), L_up0, S_up0, tau_i(Inp_N_layers+1,1), S_up(Inp_N_layers+1,1), &
							L_dn(Inp_N_layers+1,1), L_up(Inp_N_layers+1,1), R_net(Inp_N_layers+1,1), &
							K_maH, H_c, D_T(Inp_N_layers,1), H_C0

doubleprecision :: E_ma, E_T0, E_T(Inp_N_layers,1), H_T0, H_T(Inp_N_layers,1), H_ma, R_net0, mu_fmc_ma1, K_T0, temp_aK, &
						S_dn_st, L_up_st, L_dn_st

doubleprecision :: x(Inp_N_layers)

doubleprecision :: f4
doubleprecision :: H_ml, H_la, E_la, E_ml, K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1,S_lma, S_ml, S_ma, S_la, Drain, D_p, C_hl, V_l
doubleprecision :: check(2,1)
!!!! Moisture model parameter
doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
!REAL(EB), PARAMETER :: R = 8.314								! universal gas constant
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
!REAL(EB), PARAMETER :: P = 100000.0
doubleprecision, PARAMETER :: T_ref = 273.16

doubleprecision, PARAMETER :: rho_air = 1.2

doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
!REAL(EB), PARAMETER :: sigma_st = 5.67/(10**(8))
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann
!REAL(EB), PARAMETER :: LAI = 1._EB				                 !< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;

!tau_i = EXP(gama*h)

dz = hei/(INP_N_layers)
zama_h = log10(z_screen/z_0)

!check(1,1) = (liq/(D_s*rho_bulk));
check(1,1) = 0;
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;


mu_fmc_ma = S_ma/V_m;
mu_fmc_ma1 = S_ma/V_a;
mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

! correct index x(i)
! forward +1 index tau_i(i) S_up

D_T0 = D_T0a*exp(wind*D_T0b);
zeta = zeta_a + (wind*zeta_b);
K_maH = 3.23*(10**(-3))*(wind/d_fmc)**(0.5);

do i = 1,Inp_N_layers
	if (i == 1) then
		x(i) = dz/2.0;
	else
		x(i) = x(1) + (i-1)*dz;
	endif
		D_T(i,1) = D_T0*exp(zeta*((x(i)/hei)-1));
end do


do i = 1,Inp_N_layers+1
	if (i == 1) then
		tau_i(i,1) = 1;
	else
		!tau_i(i,1) = exp(-gama*1.5*x(i-1));
		tau_i(i,1) = exp(-gama*LAI*x(i-1));
	end if
end do

do i = 1,Inp_N_layers+1
S_dn(i,1) = DNI*tau_i(i,1);
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		S_up(i,1) = alpha_s*S_dn(i,1)
	else
		S_dn_st = 0
		do j = i,Inp_N_layers
			S_dn_st = S_dn_st + alpha_l*(S_dn(j,1)-S_dn(j+1,1))*tau_i(j-i+1,1)
		end do
		S_up(i,1) = (alpha_s*S_dn(Inp_N_layers+1,1)*tau_i(Inp_N_layers-i+2,1)) + S_dn_st
	end if

end do

do i = 1,Inp_N_layers+1
	L_dn_st = 0;
	if (i == 1) then
		L_dn(i,1) = SIGMA*(temp_aK**4);
	else
		do j = 2,i
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_dn_st = L_dn_st + (tau_i(i-j+1,1)-tau_i(i-j+2,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_dn_st = L_dn_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_dn(i,1)  = tau_i(i,1)*SIGMA*(temp_aK**4) + L_dn_st
	end if
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		L_up(i,1) = SIGMA*((T_s+273)**4);
	else
		L_up_st = 0;
		do j = i+1,Inp_N_layers+1
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_up_st = L_up_st + (tau_i(i-j+2,1)-tau_i(i-j+3,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_up_st = L_up_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_up(i,1) = tau_i(Inp_N_layers-i+2,1)*SIGMA*((T_s+273)**4)+L_up_st
	end if

end do

do i = 1,Inp_N_layers+1
R_net(i,1) = S_dn(i,1) - S_up(i,1) + L_dn(i,1) - L_up(i,1);
end do
!print*, R_net(Inp_N_layers+1,1);
RiB = (z_screen*9.81*(temp_aK-T_a(1,1)))/((wind**(2))*temp_aK);
gamma_h = (zama_h-10.0*RiB*zama_h-zama_h)/(2-10.0*RiB);
K_T0 = (((K_von)**(2))*wind)/((zama_h-gamma_h)*(zama_h-gamma_h));
K_laH = K_maH;
K_laE = (kappa_h/kappa_v)**(0.67);

do i = 1,Inp_N_layers
calc(i,1) = exp((mo(i,1)*B_fmc) + A_fmc);
RH_fuel(i,1) = exp((-4.19*Mass_fmc*calc(i,1))/((R0/1000)*T_m(i,1)));
q_sat(i,1) = (0.622*0.611*(exp((17.3*(T_m(i,1)-273))/(T_m(i,1)-273+237.3))))/(P_INF*0.001);
end do

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

H_ml = 0;
E_la = 0;

E_ma = rho_air*K_maE*(RH_fuel(local_index,1)*q_sat(local_index,1)-q(local_index,1));
H_ma = rho_air*C_p*K_maH*(T_m(local_index,1)-T_a(local_index,1));

if (local_index == 1) then
	H_c = (0.14+0.22*mo(local_index,1))*(T_m(local_index+1,1)-T_m(local_index,1))/dz;
	H_C0 = (0.14+0.22*mo(local_index,1))*(T_m(local_index,1)-temp_aK)/(0.5*dz);	! Should be dz/2

	f4 =  - T_m_old(local_index,1) - (dt/C_hm)*((1/V_m)*((R_net(local_index+1,1)-R_net(local_index,1))/(0.5*dz)) + &
				(1/V_m)*((H_c-H_C0)/(0.5*dz)) - mu_fmc_ma*H_ma - mu_fmc_ma*lambda*E_ma - mu_fmc_ml1*H_ml) + T_m(local_index,1);

elseif (local_index == Inp_N_layers) then
	H_c = (0.14+0.22*mo(local_index,1))*(T_s+273-T_m(local_index,1))/(0.5*dz);	! Should be dz/2
	H_C0 = (0.14+0.22*mo(local_index,1))*(T_m(local_index,1)-T_m(local_index-1,1))/dz;

	f4 =  - T_m_old(local_index,1) - (dt/C_hm)*((1/V_m)*((R_net(local_index+1,1)-R_net(local_index,1))/dz) + (1/V_m)*((H_c-H_C0)/dz) &
			- mu_fmc_ma*H_ma - mu_fmc_ma*lambda*E_ma - mu_fmc_ml1*H_ml) + T_m(local_index,1);

else

	H_c = (0.14+0.22*mo(local_index,1))*(T_m(local_index+1,1)-T_m(local_index,1))/dz;
	H_C0 = (0.14+0.22*mo(local_index,1))*(T_m(local_index,1)-T_m(local_index-1,1))/dz;

	f4 =  - T_m_old(local_index,1) - (dt/C_hm)*((1/V_m)*((R_net(local_index+1,1)-R_net(local_index,1))/dz) + (1/V_m)*((H_c-H_C0)/dz) &
        - mu_fmc_ma*H_ma - mu_fmc_ma*lambda*E_ma - mu_fmc_ml1*H_ml) + T_m(local_index,1);

end if

! IF (INT(T) .GT. 65) THEN
! print*, 'H_c, H_c0, K_maH', local_index, H_c, H_c0, K_maH
! print*, 'T_m(local_index+1), T_m(local_index), T_s, temp_aK', T_m(local_index+1,1), T_m(local_index,1), T_s, temp_aK
! print*, 'R_net(local_index+1,1), R_net(local_index,1)', local_index, R_net(local_index+1,1), R_net(local_index,1)
! ENDIF

f_Tm1 = f4;

END subroutine f_Tm_calc_1a

subroutine f_olda(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
            wind,DNI,temp_a,P_INF,RelH,rain,T_s,INP_N_layers,local_index,f_old)
IMPLICIT NONE

doubleprecision, intent(in) :: mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
							m_old(Inp_N_layers,1), q_old(Inp_N_layers,1), T_a_old(Inp_N_layers,1), T_m_old(Inp_N_layers,1), &
                                liq(Inp_N_layers,1), T_l(Inp_N_layers,1), liq_old(Inp_N_layers,1), T_liq_old(Inp_N_layers,1)
!integer, intent(in) :: data_num
doubleprecision, INTENT(IN) :: DT
INTEGER, intent(in) :: INP_N_layers,local_index
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH,rain,T_s
doubleprecision, intent(out) :: f_old(6,Inp_N_layers)

doubleprecision :: a_con, b_con, q_sat(Inp_N_layers,1), RH_fuel(Inp_N_layers,1), calc(Inp_N_layers,1), V_a


doubleprecision :: D_T0, zeta, V_m, RiB(Inp_N_layers,1), zama_h, gamma_h(Inp_N_layers,1), q0, dz

INTEGER :: i, j

doubleprecision :: S_dn(Inp_N_layers+1,1), L_up0, S_up0, tau_i(Inp_N_layers+1,1), S_up(Inp_N_layers+1,1), &
							L_dn(Inp_N_layers+1,1), L_up(Inp_N_layers+1,1), R_net(Inp_N_layers+1,1), &
							K_maH, H_c(Inp_N_layers,1), D_T(Inp_N_layers,1)

doubleprecision :: E_ma(Inp_N_layers,1), E_T0(Inp_N_layers,1), E_T(Inp_N_layers,1), H_T0(Inp_N_layers,1), &
						H_T(Inp_N_layers,1), H_ma(Inp_N_layers,1), &
						R_net0, mu_fmc_ma1, K_T0(Inp_N_layers,1), temp_aK, H_C0(Inp_N_layers,1)

doubleprecision :: x(Inp_N_layers)

doubleprecision :: f1(Inp_N_layers,1), f2(Inp_N_layers,1), f3(Inp_N_layers,1), f4(Inp_N_layers,1), &
						f5(Inp_N_layers,1), f6(Inp_N_layers,1)
doubleprecision :: H_ml(Inp_N_layers,1), H_la(Inp_N_layers,1), E_la(Inp_N_layers,1), E_ml(Inp_N_layers,1), K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1,S_lma, S_ml, S_ma, S_la, Drain, D_p, C_hl, V_l, &
						S_dn_st, L_up_st, L_dn_st
doubleprecision :: check(2,1)
!!!! Moisture model parameter
doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
!REAL(EB), PARAMETER :: R = 8.314								! universal gas constant
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
!REAL(EB), PARAMETER :: P = 100000.0
doubleprecision, PARAMETER :: T_ref = 273.16

doubleprecision, PARAMETER :: rho_air = 1.2

doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
!REAL(EB), PARAMETER :: sigma_st = 5.67/(10**(8))
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann
!REAL(EB), PARAMETER :: LAI = 1._EB				                 !< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;
!tau_i = EXP(gama*h)

dz = hei/(Inp_N_layers)
zama_h = log10(z_screen/z_0)

!check(1,1) = (liq/(D_s*rho_bulk));
check(1,1) = 0;
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;

mu_fmc_ma = S_ma/V_m;
mu_fmc_ma1 = S_ma/V_a;
mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

D_T0 = D_T0a*exp(wind*D_T0b);
zeta = zeta_a + (wind*zeta_b);
K_maH = 3.23*(10**(-3))*(wind/d_fmc)**(0.5);


do i = 1,Inp_N_layers
	if (i == 1) then
		x(i) = dz/2.0;
	else
		x(i) = x(1) + (i-1)*dz;
	endif
		D_T(i,1) = D_T0*exp(zeta*((x(i)/hei)-1));
end do


do i = 1,Inp_N_layers+1
	if (i == 1) then
		tau_i(i,1) = 1;
	else
		!tau_i(i,1) = exp(-gama*1.5*x(i-1));
		tau_i(i,1) = exp(-gama*LAI*x(i-1));
	end if
end do

do i = 1,Inp_N_layers+1
S_dn(i,1) = DNI*tau_i(i,1);
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		S_up(i,1) = alpha_s*S_dn(i,1)
	else
		S_dn_st = 0
		do j = i,Inp_N_layers
			S_dn_st = S_dn_st + alpha_l*(S_dn(j,1)-S_dn(j+1,1))*tau_i(j-i+1,1)
		end do
		S_up(i,1) = (alpha_s*S_dn(Inp_N_layers+1,1)*tau_i(Inp_N_layers-i+2,1)) + S_dn_st
	end if

end do

do i = 1,Inp_N_layers+1
	L_dn_st = 0;
	if (i == 1) then
		L_dn(i,1) = SIGMA*(temp_aK**4);
	else
		do j = 2,i
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_dn_st = L_dn_st + (tau_i(i-j+1,1)-tau_i(i-j+2,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_dn_st = L_dn_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_dn(i,1)  = tau_i(i,1)*SIGMA*(temp_aK**4) + L_dn_st
	end if
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		L_up(i,1) = SIGMA*((T_s+273)**4);
	else
		L_up_st = 0;
		do j = i+1,Inp_N_layers+1
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_up_st = L_up_st + (tau_i(i-j+2,1)-tau_i(i-j+3,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_up_st = L_up_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_up(i,1) = tau_i(Inp_N_layers-i+2,1)*SIGMA*((T_s+273)**4)+L_up_st
	end if

end do

do i = 1,Inp_N_layers+1
R_net(i,1) = S_dn(i,1) - S_up(i,1) + L_dn(i,1) - L_up(i,1);
end do

do i = 1,Inp_N_layers
RiB(i,1) = (z_screen*9.81*(temp_aK-T_a(i,1)))/((wind**(2))*temp_aK);
gamma_h(i,1) = (zama_h-10.0*RiB(i,1)*zama_h-zama_h)/(2-10.0*RiB(i,1));
K_T0(i,1) = (((K_von)**(2))*wind)/((zama_h-gamma_h(i,1))*(zama_h-gamma_h(i,1)));
end do

do i = 1,Inp_N_layers
calc(i,1) = exp((mo(i,1)*B_fmc) + A_fmc);
RH_fuel(i,1) = exp((-4.19*Mass_fmc*calc(i,1))/((R0/1000)*T_m(i,1)));
q_sat(i,1) = (0.622*0.611*(exp((17.3*(T_m(i,1)-273))/(T_m(i,1)-273+237.3))))/(P_INF*0.001);
end do

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

K_laH = K_maH;
K_laE = (kappa_h/kappa_v)**(0.67);

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

do i = 1,INP_N_layers
    H_ml(i,1) = 0;
    H_la(i,1) = 0;
    E_la(i,1) = 0;
    E_ml(i,1) = 0;
end do


do i = 1,Inp_N_layers
E_ma(i,1) = rho_air*K_maE*(RH_fuel(i,1)*q_sat(i,1)-q(i,1));
f1(i,1) = - m_old(i,1) + (dt/rho_litter)*(mu_fmc_ma*E_ma(i,1) + mu_fmc_ml*E_ml(i,1)) + mo(i,1);
end do

do i = 1,Inp_N_layers

	if (i == 1) then
		E_T0(i,1) = rho_air*K_T0(i,1)*(q(i,1)-q0);
		E_T(i,1) = rho_air*D_T0*exp(zeta*((((hei/(Inp_N_layers)))/hei)-1))*(q(i+1,1)-q(i,1))/dz;
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*0.5*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	else if (i == Inp_N_layers) then
		E_T0(i,1) = rho_air*D_T0*exp(zeta*(((Inp_N_layers*hei/(Inp_N_layers))/hei)-1))*(q(i,1)-q(i-1,1))/dz;
		!E_T(i,1) = 0;
		E_T(i,1) = rho_air*D_T0*exp(zeta*(((Inp_N_layers*hei/(Inp_N_layers))/hei)-1))*((0.0001)-q(i,1))/(0.5*dz);
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	else
		E_T0(i,1) = rho_air*D_T0*exp(zeta*((((i-1)*hei)/(Inp_N_layers))/hei-1))*(q(i,1)-q(i-1,1))/dz;
		E_T(i,1) = rho_air*D_T0*exp(zeta*(((i*hei)/(Inp_N_layers))/hei-1))*(q(i+1,1)-q(i,1))/dz;
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	end if

end do

do i = 1,Inp_N_layers

	if (i == 1) then
		H_T0(i,1) = rho_air*C_p*K_T0(i,1)*(T_a(i,1)-temp_aK);
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*(((hei)/((Inp_N_layers)*hei))-1))*(T_a(i+1,1)-T_a(i,1))/dz;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/(0.5*dz) + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);

	else if (i == Inp_N_layers) then
		H_T0(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i,1)-T_a(i-1,1))/dz;
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i)*hei)/((Inp_N_layers)*hei))-1))*(T_s+273-T_a(i,1))/(0.5*dz);	! Should be dz/2
		!H_T(i,1) = 0;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/dz + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);

	else
		H_T0(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i,1)-T_a(i-1,1))/dz;
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i+1,1)-T_a(i,1))/dz;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/dz + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);
	end if

end do

do i = 1,Inp_N_layers

	if (i == 1) then
		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_m(i+1,1)-T_m(i,1))/dz;
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-temp_aK)/(0.5*dz);	! Should be dz/2

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/(0.5*dz)) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/(0.5*dz)) &
				- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	elseif (i == Inp_N_layers) then
		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_s+273-T_m(i,1))/(0.5*dz);	! Should be dz/2
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-T_m(i-1,1))/dz;

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/dz) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/dz) &
				- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	else

		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_m(i+1,1)-T_m(i,1))/dz;
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-T_m(i-1,1))/dz;

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/dz) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/dz) &
			- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	end if

end do

do i = 1,Inp_N_layers
f_old(1,i) = f1(i,1); f_old(2,i) = f2(i,1); f_old(3,i) = f3(i,1); f_old(4,i) = f4(i,1); f_old(5,i) = f5(i,1); f_old(6,i) = f6(i,1);
end do

END subroutine f_olda

subroutine vector_f(mo,q,T_a,T_m,m_old,q_old,T_a_old,T_m_old,liq,T_l,liq_old,T_liq_old,dt,&
            wind,DNI,temp_a,P_INF,RelH,rain,T_s,INP_N_layers,local_index,F)
IMPLICIT NONE

doubleprecision, intent(in) :: mo(Inp_N_layers,1), q(Inp_N_layers,1), T_a(Inp_N_layers,1), T_m(Inp_N_layers,1), &
							m_old(Inp_N_layers,1), q_old(Inp_N_layers,1), T_a_old(Inp_N_layers,1), T_m_old(Inp_N_layers,1), &
                                liq(Inp_N_layers,1), T_l(Inp_N_layers,1), liq_old(Inp_N_layers,1), T_liq_old(Inp_N_layers,1)
!integer, intent(in) :: data_num
doubleprecision, INTENT(IN) :: DT
doubleprecision, intent(in) :: wind, DNI, temp_a, P_INF, RelH, rain, T_s
INTEGER, intent(in) :: INP_N_layers,local_index
doubleprecision, intent(out) :: F(6,Inp_N_layers)

doubleprecision :: a_con, b_con, q_sat(Inp_N_layers,1), RH_fuel(Inp_N_layers,1), calc(Inp_N_layers,1), V_a


doubleprecision :: D_T0, zeta, V_m, RiB(Inp_N_layers,1), zama_h, gamma_h(Inp_N_layers,1), q0, dz

INTEGER :: i, j

doubleprecision :: S_dn(Inp_N_layers+1,1), L_up0, S_up0, tau_i(Inp_N_layers+1,1), S_up(Inp_N_layers+1,1), &
							L_dn(Inp_N_layers+1,1), L_up(Inp_N_layers+1,1), R_net(Inp_N_layers+1,1), K_maH, &
							H_c(Inp_N_layers,1), D_T(Inp_N_layers,1)

doubleprecision :: E_ma(Inp_N_layers,1), E_T0(Inp_N_layers,1), E_T(Inp_N_layers,1), H_T0(Inp_N_layers,1), &
							H_T(Inp_N_layers,1), H_ma(Inp_N_layers,1), &
							R_net0, mu_fmc_ma1, K_T0(Inp_N_layers,1), temp_aK,H_C0(Inp_N_layers,1)

doubleprecision :: x(Inp_N_layers)

doubleprecision :: f1(Inp_N_layers,1), f2(Inp_N_layers,1), f3(Inp_N_layers,1), &
    f4(Inp_N_layers,1), f5(Inp_N_layers,1), f6(Inp_N_layers,1)
doubleprecision :: H_ml(Inp_N_layers,1), H_la(Inp_N_layers,1), E_la(Inp_N_layers,1), E_ml(Inp_N_layers,1), K_laH, K_laE
doubleprecision :: mu_fmc_ma,mu_fmc_ml, mu_fmc_la, mu_fmc_ml1, mu_fmc_la1,S_lma, S_ml, S_ma, S_la, Drain, D_p, C_hl, V_l,&
 S_dn_st, L_up_st, L_dn_st
doubleprecision :: check(2,1)
!!!! Moisture model parameter
doubleprecision, PARAMETER :: rho_litter = 512.0						! litter density
!REAL(EB), PARAMETER :: R = 8.314								! universal gas constant
doubleprecision, PARAMETER :: Mass_fmc = 18.0153							! molecular mass of water
doubleprecision, PARAMETER :: A_fmc = 5.2								! Nelson model constant
doubleprecision, PARAMETER :: B_fmc = -19.0								! Nelson model constant
!REAL(EB), PARAMETER :: P = 100000.0
doubleprecision, PARAMETER :: T_ref = 273.16

doubleprecision, PARAMETER :: rho_air = 1.2

doubleprecision, PARAMETER :: gama = 1.363
doubleprecision, PARAMETER :: hei = 0.02
doubleprecision, PARAMETER :: alpha_s = 0.2
!REAL(EB), PARAMETER :: sigma_st = 5.67/(10**(8))
doubleprecision, PARAMETER :: alpha_l = 0.27

doubleprecision, PARAMETER :: mu_fmc = 9770.0
doubleprecision, PARAMETER :: K_maE = 0.0006
doubleprecision, PARAMETER :: rho_bulk = 201.91

doubleprecision, PARAMETER :: C_p = 1004.5
doubleprecision, PARAMETER :: lambda = 2.45*1000000
doubleprecision, PARAMETER :: K_c = 0.2
doubleprecision, PARAMETER :: tau_i0 = 1.0

doubleprecision, PARAMETER :: D_T0a = 0.00002
doubleprecision, PARAMETER :: D_T0b = 2.6
doubleprecision, PARAMETER :: zeta_a = 2.08
doubleprecision, PARAMETER :: zeta_b = 2.38
doubleprecision, PARAMETER :: d_fmc = 0.03

doubleprecision, PARAMETER :: z_screen = 2.0
doubleprecision, PARAMETER :: z_0 = 0.01
doubleprecision, PARAMETER :: K_von = 0.4

doubleprecision, PARAMETER :: C_hm = 1.0*1000000
doubleprecision, PARAMETER :: kappa_h = 2.08
doubleprecision, PARAMETER :: kappa_v = 2.34

doubleprecision, PARAMETER :: rho_water = 1000.0

doubleprecision, PARAMETER :: D_s = 1.153
doubleprecision, PARAMETER :: D_a = 0.00003
doubleprecision, PARAMETER :: L_a = 0.23;
doubleprecision, PARAMETER :: L_b = -1.63;
doubleprecision, PARAMETER :: K_mlH = 700.0;
doubleprecision, PARAMETER :: R0=8314.472                       !< Gas constant (J/K/kmol)
doubleprecision, PARAMETER :: SIGMA=5.670373/100000000                 !< Stefan-Boltzmann
!REAL(EB), PARAMETER :: LAI = 1._EB				                 !< Leaf area index
doubleprecision, PARAMETER :: LAI = 0.5				                 ! CC = 94.00*[1-exp(-0.43*LAI)]^(0.52)

temp_aK = temp_a + 273.0;
!tau_i = EXP(gama*h)

dz = hei/(INP_N_layers)
zama_h = log10(z_screen/z_0)

check(1,1) = 0;
check(2,1) = 1;

S_lma = mu_fmc*(rho_bulk/rho_litter);
S_la = minval(check);
S_ml = S_la;
S_ma = S_lma - S_ml;
V_m = rho_bulk/rho_litter;
V_l = 0;
V_a = 1-V_m-V_l;

mu_fmc_ma = S_ma/V_m;
mu_fmc_ma1 = S_ma/V_a;
mu_fmc_la = 0;
mu_fmc_la1 = 0;
mu_fmc_ml = 0;
mu_fmc_ml1 = 0;

D_T0 = D_T0a*exp(wind*D_T0b);
zeta = zeta_a + (wind*zeta_b);
K_maH = 3.23*(10**(-3))*(wind/d_fmc)**(0.5);

do i = 1,Inp_N_layers
	if (i == 1) then
		x(i) = dz/2.0;
	else
		x(i) = x(1) + (i-1)*dz;
	endif
	D_T(i,1) = D_T0*exp(zeta*((x(i)/hei)-1));
end do

do i = 1,Inp_N_layers+1
	if (i == 1) then
		tau_i(i,1) = 1;
	else
		!tau_i(i,1) = exp(-gama*1.5*x(i-1));
		tau_i(i,1) = exp(-gama*LAI*x(i-1));
	end if
end do

do i = 1,Inp_N_layers+1
S_dn(i,1) = DNI*tau_i(i,1);
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		S_up(i,1) = alpha_s*S_dn(i,1)
	else
		S_dn_st = 0
		do j = i,Inp_N_layers
			S_dn_st = S_dn_st + alpha_l*(S_dn(j,1)-S_dn(j+1,1))*tau_i(j-i+1,1)
		end do
		S_up(i,1) = (alpha_s*S_dn(Inp_N_layers+1,1)*tau_i(Inp_N_layers-i+2,1)) + S_dn_st
	end if

end do

do i = 1,Inp_N_layers+1
	L_dn_st = 0;
	if (i == 1) then
		L_dn(i,1) = SIGMA*(temp_aK**4);
	else
		do j = 2,i
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_dn_st = L_dn_st + (tau_i(i-j+1,1)-tau_i(i-j+2,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_dn_st = L_dn_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_dn(i,1)  = tau_i(i,1)*SIGMA*(temp_aK**4) + L_dn_st
	end if
end do

do i = 1,Inp_N_layers+1

	if (i == Inp_N_layers+1) then
		L_up(i,1) = SIGMA*((T_s+273)**4);
	else
		L_up_st = 0;
		do j = i+1,Inp_N_layers+1
			if (j .GE. 2 .AND. j .LE. i+1) then
				L_up_st = L_up_st + (tau_i(i-j+2,1)-tau_i(i-j+3,1))*SIGMA*(T_m(j-1,1)**4);
			else if (j .GT. i+1 .AND. j .LT. (Inp_N_layers+2)) then
				L_up_st = L_up_st + (tau_i(j-i,1)-tau_i(j-i+1,1))*SIGMA*(T_m(j-1,1)**4);
			end if
		end do
		L_up(i,1) = tau_i(Inp_N_layers-i+2,1)*SIGMA*((T_s+273)**4)+L_up_st
	end if

end do

do i = 1,Inp_N_layers+1
R_net(i,1) = S_dn(i,1) - S_up(i,1) + L_dn(i,1) - L_up(i,1);
end do

do i = 1,Inp_N_layers
RiB(i,1) = (z_screen*9.81*(temp_aK-T_a(i,1)))/((wind**(2))*temp_aK);
gamma_h(i,1) = (zama_h-10.0*RiB(i,1)*zama_h-zama_h)/(2-10.0*RiB(i,1));
K_T0(i,1) = (((K_von)**(2))*wind)/((zama_h-gamma_h(i,1))*(zama_h-gamma_h(i,1)));
end do

do i = 1,Inp_N_layers
calc(i,1) = exp((mo(i,1)*B_fmc) + A_fmc);
RH_fuel(i,1) = exp((-4.19*Mass_fmc*calc(i,1))/((R0/1000)*T_m(i,1)));
q_sat(i,1) = (0.622*0.611*(exp((17.3*(T_m(i,1)-273))/(T_m(i,1)-273+237.3))))/(P_INF*0.001);
end do

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

K_laH = K_maH;
K_laE = (kappa_h/kappa_v)**(0.67);

q0 = (RelH*exp((17.67*((temp_a+273)-T_ref))/((temp_a+273)-29.65)))/(0.263*P_INF);

do i = 1,INP_N_layers
    H_ml(i,1) = 0;
    H_la(i,1) = 0;
    E_la(i,1) = 0;
    E_ml(i,1) = 0;
end do

do i = 1,Inp_N_layers
E_ma(i,1) = rho_air*K_maE*(RH_fuel(i,1)*q_sat(i,1)-q(i,1));
f1(i,1) = - m_old(i,1) + (dt/rho_litter)*(mu_fmc_ma*E_ma(i,1) + mu_fmc_ml*E_ml(i,1)) + mo(i,1);
end do

do i = 1,Inp_N_layers

	if (i == 1) then
		E_T0(i,1) = rho_air*K_T0(i,1)*(q(i,1)-q0);
		E_T(i,1) = rho_air*D_T0*exp(zeta*((((hei/(Inp_N_layers)))/hei)-1))*(q(i+1,1)-q(i,1))/dz;
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*0.5*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	else if (i == Inp_N_layers) then
		E_T0(i,1) = rho_air*D_T0*exp(zeta*(((Inp_N_layers*hei/(Inp_N_layers))/hei)-1))*(q(i,1)-q(i-1,1))/dz;
		!E_T(i,1) = 0;
		E_T(i,1) = rho_air*D_T0*exp(zeta*(((Inp_N_layers*hei/(Inp_N_layers))/hei)-1))*((0.0001)-q(i,1))/(0.5*dz);
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	else
		E_T0(i,1) = rho_air*D_T0*exp(zeta*((((i-1)*hei)/(Inp_N_layers))/hei-1))*(q(i,1)-q(i-1,1))/dz;
		E_T(i,1) = rho_air*D_T0*exp(zeta*(((i*hei)/(Inp_N_layers))/hei-1))*(q(i+1,1)-q(i,1))/dz;
		f2(i,1) = - q_old(i,1) - (dt/rho_air)*(((E_T(i,1)-E_T0(i,1))/(V_a*dz)) + mu_fmc_ma1*E_ma(i,1)) + q(i,1);

	end if

end do

do i = 1,Inp_N_layers

	if (i == 1) then
		H_T0(i,1) = rho_air*C_p*K_T0(i,1)*(T_a(i,1)-temp_aK);
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*(((hei)/((Inp_N_layers)*hei))-1))*(T_a(i+1,1)-T_a(i,1))/dz;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/(0.5*dz) + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);

	else if (i == Inp_N_layers) then
		H_T0(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i,1)-T_a(i-1,1))/dz;
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i)*hei)/((Inp_N_layers)*hei))-1))*(T_s+273-T_a(i,1))/(0.5*dz);	! Should be dz/2
		!H_T(i,1) = 0;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/dz + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);

	else
		H_T0(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i-1)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i,1)-T_a(i-1,1))/dz;
		H_T(i,1) = rho_air*C_p*D_T0*exp(zeta*((((i)*hei)/((Inp_N_layers)*hei))-1))*(T_a(i+1,1)-T_a(i,1))/dz;
		H_ma(i,1) = rho_air*C_p*K_maH*(T_m(i,1)-T_a(i,1));

		f3(i,1) = - T_a_old(i,1) - (dt/(C_p*rho_air))*((1/V_a)*(H_T(i,1)-H_T0(i,1))/dz + mu_fmc_ma1*H_ma(i,1) + &
					mu_fmc_la1*H_la(i,1)) + T_a(i,1);
	end if

end do

do i = 1,Inp_N_layers

	if (i == 1) then
		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_m(i+1,1)-T_m(i,1))/dz;
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-temp_aK)/(0.5*dz);	! Should be dz/2

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/(0.5*dz)) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/(0.5*dz)) &
				- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	elseif (i == Inp_N_layers) then
		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_s+273-T_m(i,1))/(0.5*dz);	! Should be dz/2
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-T_m(i-1,1))/dz;

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/dz) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/dz) &
				- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	else

		H_c(i,1) = (0.14+0.22*mo(i,1))*(T_m(i+1,1)-T_m(i,1))/dz;
		H_C0(i,1) = (0.14+0.22*mo(i,1))*(T_m(i,1)-T_m(i-1,1))/dz;

		f4(i,1) =  - T_m_old(i,1) - (dt/C_hm)*((1/V_m)*((R_net(i+1,1)-R_net(i,1))/dz) + (1/V_m)*((H_c(i,1)-H_C0(i,1))/dz) &
			- mu_fmc_ma*H_ma(i,1) - mu_fmc_ma*lambda*E_ma(i,1) - mu_fmc_ml1*H_ml(i,1)) + T_m(i,1);

	end if

end do

do i = 1,Inp_N_layers
F(1,i) = f1(i,1); F(2,i) = f2(i,1); F(3,i) = f3(i,1); F(4,i) = f4(i,1); F(5,i) = f5(i,1); F(6,i) = f6(i,1);
end do

END subroutine vector_f
