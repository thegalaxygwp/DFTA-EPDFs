    module dimens_analysis
    doubleprecision,allocatable :: &
        vot_x0(:,:),vot_y0(:,:),vot_z0(:,:),    &
        vot_x0_x1(:,:),vot_y0_x1(:,:),vot_z0_x1(:,:),    &
        
        vot_x1(:,:),vot_y1(:,:),vot_z1(:,:),    &
        vot_x1_x1(:,:),vot_y1_x1(:,:),vot_z1_x1(:,:),    &
        vot_x1_y1(:,:),vot_y1_y1(:,:),vot_z1_y1(:,:),    &
        
        tvot_x1(:,:),tvot_y1(:,:),tvot_z1(:,:),    &
        tvot_x1_x1(:,:),tvot_y1_x1(:,:),tvot_z1_x1(:,:),    &
        tvot_x1_y1(:,:),tvot_y1_y1(:,:),tvot_z1_y1(:,:)
    
    doublecomplex,allocatable :: &        
        fvot_x1(:,:),fvot_y1(:,:),fvot_z1(:,:),    &
        fvot_x1_x1(:,:),fvot_y1_x1(:,:),fvot_z1_x1(:,:),    &
        fvot_x1_y1(:,:),fvot_y1_y1(:,:),fvot_z1_y1(:,:),    &
        
        ftvot_x1(:,:),ftvot_y1(:,:),ftvot_z1(:,:),    &
        ftvot_x1_x1(:,:),ftvot_y1_x1(:,:),ftvot_z1_x1(:,:),    &
        ftvot_x1_y1(:,:),ftvot_y1_y1(:,:),ftvot_z1_y1(:,:)
        
    doubleprecision,allocatable :: &
        tpsi_z1_fr1(:,:),tpsi_z1_fr2(:,:),    &
        tpsi_z1_de1(:,:),tpsi_z1_de2(:,:),tpsi_z1_de3(:,:),    &
        
        
        tmgt_z1_fr1(:,:),tmgt_z1_fr2(:,:),tmgt_z1_fr3(:,:),    &
        tmgt_z1_bir(:,:),    &
        tmgt_z1_de1(:,:),tmgt_z1_de2(:,:),tmgt_z1_de3(:,:),    &
        tmgt_z1_de4(:,:),tmgt_z1_de0(:,:),   &
        
        tvtt_x1(:,:),tvtt_y1(:,:),tdst_e1(:,:),    &
        tpsi_z1(:,:),    &
        tmgt_x1(:,:),tmgt_y1(:,:),tmgt_z1(:,:),    &
        tcrt_x1(:,:),tcrt_y1(:,:),tcrt_z1(:,:),    &
        tvlt_x1(:,:),tvlt_y1(:,:),tvlt_z1(:,:),    &
        
        tvtt_x1_x1(:,:),tvtt_y1_x1(:,:),tvtt_z1_x1(:,:),    &
        tvtt_p1_x1(:,:),tpse_e1_x1(:,:),tdst_e1_x1(:,:),    &
        tpsi_z1_x1(:,:),    &
        tmgt_x1_x1(:,:),tmgt_y1_x1(:,:),tmgt_z1_x1(:,:),    &
        tcrt_x1_x1(:,:),tcrt_y1_x1(:,:),tcrt_z1_x1(:,:),    &
        tvlt_x1_x1(:,:),tvlt_y1_x1(:,:),tvlt_z1_x1(:,:),    &
        
        tvtt_x1_y1(:,:),tvtt_y1_y1(:,:),tvtt_z1_y1(:,:),    &
        tvtt_p1_y1(:,:),tpse_e1_y1(:,:),tdst_e1_y1(:,:),    &
        tpsi_z1_y1(:,:),    &
        tmgt_x1_y1(:,:),tmgt_y1_y1(:,:),tmgt_z1_y1(:,:),    &
        tcrt_x1_y1(:,:),tcrt_y1_y1(:,:),tcrt_z1_y1(:,:),    &
        tvlt_x1_y1(:,:),tvlt_y1_y1(:,:),tvlt_z1_y1(:,:),    &    
        
        tvtt_x1_d2(:,:),tvtt_y1_d2(:,:),tvtt_z1_d2(:,:),    &
        tvtt_p1_d2(:,:),tpse_e1_d2(:,:),tdst_e1_d2(:,:),    &
        tpsi_z1_d2(:,:),    &
        tmgt_x1_d2(:,:),tmgt_y1_d2(:,:),tmgt_z1_d2(:,:),    &
        tcrt_x1_d2(:,:),tcrt_y1_d2(:,:),tcrt_z1_d2(:,:),    &
        tvlt_x1_d2(:,:),tvlt_y1_d2(:,:),tvlt_z1_d2(:,:)
        
    doublecomplex,allocatable ::    &
        ftvtt_x1(:,:),ftvtt_y1(:,:),ftvtt_z1(:,:),    &
        ftvtt_p1(:,:),ftpse_e1(:,:),ftdst_e1(:,:),    &
        ftpsi_z1(:,:),    &
        ftmgt_x1(:,:),ftmgt_y1(:,:),ftmgt_z1(:,:),    &
        ftcrt_x1(:,:),ftcrt_y1(:,:),ftcrt_z1(:,:),    &
        ftvlt_x1(:,:),ftvlt_y1(:,:),ftvlt_z1(:,:),    &
        
        ftvtt_x1_x1(:,:),ftvtt_y1_x1(:,:),ftvtt_z1_x1(:,:),    &
        ftvtt_p1_x1(:,:),ftpse_e1_x1(:,:),ftdst_e1_x1(:,:),    &
        ftpsi_z1_x1(:,:),    &
        ftmgt_x1_x1(:,:),ftmgt_y1_x1(:,:),ftmgt_z1_x1(:,:),    &
        ftcrt_x1_x1(:,:),ftcrt_y1_x1(:,:),ftcrt_z1_x1(:,:),    &
        ftvlt_x1_x1(:,:),ftvlt_y1_x1(:,:),ftvlt_z1_x1(:,:),    &
        
        ftvtt_x1_y1(:,:),ftvtt_y1_y1(:,:),ftvtt_z1_y1(:,:),    &
        ftvtt_p1_y1(:,:),ftpse_e1_y1(:,:),ftdst_e1_y1(:,:),    &
        ftpsi_z1_y1(:,:),    &
        ftmgt_x1_y1(:,:),ftmgt_y1_y1(:,:),ftmgt_z1_y1(:,:),    &
        ftcrt_x1_y1(:,:),ftcrt_y1_y1(:,:),ftcrt_z1_y1(:,:),    &
        ftvlt_x1_y1(:,:),ftvlt_y1_y1(:,:),ftvlt_z1_y1(:,:),    &
        
        ftvtt_x1_d2(:,:),ftvtt_y1_d2(:,:),ftvtt_z1_d2(:,:),    &
        ftvtt_p1_d2(:,:),ftpse_e1_d2(:,:),ftdst_e1_d2(:,:),    &
        ftpsi_z1_d2(:,:),    &
        ftmgt_x1_d2(:,:),ftmgt_y1_d2(:,:),ftmgt_z1_d2(:,:),    &
        ftcrt_x1_d2(:,:),ftcrt_y1_d2(:,:),ftcrt_z1_d2(:,:),    &
        ftvlt_x1_d2(:,:),ftvlt_y1_d2(:,:),ftvlt_z1_d2(:,:)
    
    doublecomplex,allocatable ::  &
        fans01(:),fans02(:),fans03(:),fans04(:),fans05(:),fans06(:),    &
        fans07(:),fans08(:),fans09(:),fans10(:),fans11(:),fans12(:),    &
        fans13(:),fans14(:),fans15(:),fans16(:),fans17(:),fans18(:),    &
        fans19(:),fans20(:),fans21(:),fans22(:),fans23(:),fans24(:),fans25(:),fans26(:)
    end module