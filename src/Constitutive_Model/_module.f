#adddirectory module Baseclass


#adddirectory module Elasticity


#adddirectory module Hypoplasticity_Wu92


#adddirectory module Hypoplasticity_VW96


#adddirectory module Viscohypoplasticity_Ni03


#adddirectory module Barodesy_Ko15


#adddirectory module Barodesy_Sc18


#adddirectory module Barodesy_Ko21


#adddirectory module Test_DGL


! ------------------------------------------------------------------ ! ----------------------------------------------- !
module Constitutive_Model_Class
   use General_Settings, only: dp, Write_Error_And_Exit, setting_len_id, setting_id_elasticity, setting_id_hypo_wu92, &
                               setting_id_hypo_vw96, setting_id_viscohypo_ni03, setting_id_barodesy_ko15, &
                               setting_id_barodesy_sc18, setting_id_barodesy_ko21, setting_id_test_dgl
   use Constitutive_Model_Baseclass
   implicit none

   private
   public :: Select_Constitutive_Model, Constitutive_Model


   contains


#addfile subroutine Select_Constitutive_Model
end module Constitutive_Model_Class
