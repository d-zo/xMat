! ------------------------------------------------------------------ ! ----------------------------------------------- !
module PlaxisInformationPool
   use General_Settings, only: dp, setting_num_statevariables, setting_max_mat_params, setting_len_id, &
                               setting_id_elasticity, setting_id_hypo_wu92, setting_id_hypo_vw96, &
                               setting_id_viscohypo_ni03, setting_id_barodesy_ko15, setting_id_barodesy_sc18, &
                               setting_id_barodesy_ko21
   implicit none

   ! ATTENTION: The order of the constitutive models has to be the same for all variables. Currently, if any changes
   !            are made to the constitutive models (or which constitutive models are included) adjustments may be
   !            necessary here
   character(len=setting_len_id), parameter, dimension(7) :: ModelList = [ &
      setting_id_elasticity, setting_id_hypo_wu92, setting_id_hypo_vw96, setting_id_viscohypo_ni03, &
      setting_id_barodesy_ko15, setting_id_barodesy_sc18, setting_id_barodesy_ko21]
   integer, parameter, dimension(7) :: numParameters =     [2, 4, 16, 15, 7, 9, 10]
   integer, parameter, dimension(7) :: numStateVariables = [0, 1, 11, 11, 2, 2, 13]

   character(len=16), parameter, dimension(16, 7) :: ModelParameterNamesUnits = reshape([ &
      ! Elasticity
      'E        :F/L^2#', '@nu#     :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-Wu92
      'C_1#     :-     ', 'C_2#     :-     ', 'C_3#     :-     ', 'C_4#     :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Hypo-VW96
      '@j#_c#   :rad   ', '@nu#     :-     ', 'h_s#     :F/L^2#', 'n_H      :-     ', &
      'e_d0#    :-     ', 'e_c0#    :-     ', 'e_i0#    :-     ', '@a#_H#   :-     ', &
      '@b#_H#   :-     ', 'm_T#     :-     ', 'm_R#     :-     ', 'R_max#   :-     ', &
      '@b#_r#   :-     ', '@c#      :-     ', '-        :-     ', 'e_0#     :-     ', &
      !
      ! ViHy-Ni03
      'e_100#   :-     ', '@nu#     :-     ', '@l#      :-     ', '@k#      :-     ', &
      '@b#_b#   :-     ', 'I_v#     :-     ', 'D_r#     :-     ', '@j#_c#   :rad   ', &
      '-        :-     ', 'm_T#     :-     ', 'm_R#     :-     ', 'R_max#   :-     ', &
      '@b#_r#   :-     ', '@c#      :-     ', 'OCR      :-     ', '         :-     ', &
      !
      ! Baro-Ko15
      '@j#_c#   :rad   ', 'c_2#     :-     ', 'c_3#     :-     ', 'c_4#     :-     ', &
      'c_5#     :-     ', 'e_c0#    :-     ', 'e_min#   :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Sc18
      '@j#_c#   :rad   ', 'k_r#     :-     ', '@x#      :-     ', '@k#      :-     ', &
      'e_c0#    :-     ', 'c_e#     :-     ', 'c_5#     :-     ', 'c_6#     :-     ', &
      'c_7#     :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      !
      ! Baro-Ko21
      '@j#_c#   :rad   ', 'c_2#     :-     ', 'c_3#     :-     ', 'c_4#     :-     ', &
      'c_5#     :-     ', '@k_1#    :-     ', '@k_2#    :-     ', 'k_1#     :-     ', &
      'k_2#     :-     ', 'e_c0#    :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ' &
      ], [16, 7])

   character(len=16), parameter, dimension(13, 7) :: ModelStatevarNamesUnits = reshape([ &
      ! Elasticity
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', &
      !
      ! Hypo-Wu92
      'voidratio:-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', &
      !
      ! Hypo-VW96
      'voidratio:-     ', 'igran11  :-     ', 'igran22  :-     ', 'igran33  :-     ', &
      'igran12  :-     ', 'igran21  :-     ', 'igran23  :-     ', 'igran32  :-     ', &
      'igran13  :-     ', 'igran31  :-     ', 'young_rep:F/L^2#', '         :-     ', &
      '         :-     ', &
      !
      ! ViHy-Ni03
      'voidratio:-     ', 'igran11  :-     ', 'igran22  :-     ', 'igran33  :-     ', &
      'igran12  :-     ', 'igran21  :-     ', 'igran23  :-     ', 'igran32  :-     ', &
      'igran13  :-     ', 'igran31  :-     ', 'young_rep:F/L^2#', '         :-     ', &
      '         :-     ', &
      !
      ! Baro-Ko15
      'voidratio:-     ', 'young_rep:F/L^2#', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', &
      !
      ! Baro-Sc18
      'voidratio:-     ', 'young_rep:F/L^2#', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', '         :-     ', '         :-     ', '         :-     ', &
      '         :-     ', &
      !
      ! Baro-Ko21
      'voidratio:-     ', 'e_c1#    :-     ', 'e_c2#    :-     ', 'D_old#11 :-     ', &
      'D_old#22 :-     ', 'D_old#33 :-     ', 'D_old#12 :-     ', 'D_old#21 :-     ', &
      'D_old#23 :-     ', 'D_old#32 :-     ', 'D_old#13 :-     ', 'D_old#31 :-     ', &
      'young_rep:F/L^2#' &
      ], [13, 7])
end module PlaxisInformationPool


#addfile subroutine USER_MOD


#addfile subroutine GetModelCount


#addfile subroutine GetModelName


#addfile subroutine GetParamCount


#addfile subroutine GetParamName


#addfile subroutine GetParamUnit


#addfile subroutine GetStateVarCount


#addfile subroutine GetStateVarName


#addfile subroutine GetStateVarNameInternal


#addfile subroutine GetStateVarUnit


#addfile subroutine GetStateVarUnitInternal
