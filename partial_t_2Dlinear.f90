    subroutine partial_t_2Dlinear()
    use module_solvers
    implicit doubleprecision (a-h,o-z)
    tvtt_p=0d0
    tvtt_z=0d0
    tpse_e=0d0
    
    call partial_initial_linear()
        !$omp parallel num_threads(max_thread)
        !$omp do collapse(2) private(i,j)
        do j=1,nj
            do i=1,ni
                tps1_a(i,j)=-vtt_y0(i,j)*cpt_v1(i,j) &
                            +vtt_x1(i,j)*vlt_y0_dx(i,j) &
                            +vtt_y0(i,j)*vlt_y1_dy(i,j) &
                            -vlt_x1(i,j)*vtt_y0_dx(i,j) &
                            -vlt_y0(i,j)*vtt_y1_dy(i,j) &
                            -eta*mgt_y1_d2(i,j)

                tvtt_z1(i,j)=-vtt_z0(i,j)*cpt_v1(i,j) &
                            +vtt_x1(i,j)*vlt_z0_dx(i,j) &
                            +vtt_y0(i,j)*vlt_z1_dy(i,j) &
                            -vlt_x1(i,j)*vtt_z0_dx(i,j) &
                            -vlt_y0(i,j)*vtt_z1_dy(i,j) &
                            -pse_e1_dy(i,j)*dst_e0_dx(i,j)/dst_e0(i,j)**2d0 &
                            -eta*mgt_z1_d2(i,j)

                !tpse_e1(i,j)=-vlt_y0(i,j)*pse_e1_dy(i,j) &
                !            -gamma*pse_e0(i,j)*cpt_v1(i,j)
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    
    return
    end