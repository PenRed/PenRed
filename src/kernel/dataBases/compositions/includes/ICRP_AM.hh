 
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//
//    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
//
//    PenRed is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PenRed is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//
//

#ifndef _PENRED_ICRP_AM_TISSUE_DATABASE_
#define _PENRED_ICRP_AM_TISSUE_DATABASE_

#include "dataBasesCommon.hh"

namespace penred{

  namespace dataBases{

    namespace compositions{

      namespace ICRP{

	struct AM : public materials{

	  static constexpr const std::array<massFraction, 9> Adrenal_left = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.228),
	    massFraction(7, 0.028),
	    massFraction(8, 0.63),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Adrenal_right = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.228),
	    massFraction(7, 0.028),
	    massFraction(8, 0.63),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET1_0_8_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET1_8_40_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET1_40_50_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET1_50_Surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 2> ET2_15_0_ = {
	    massFraction(1, 0.112),
	    massFraction(8, 0.888)
	  };

	  static constexpr const std::array<massFraction, 9> ET2_0_40_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET2_40_50_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET2_50_55_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET2_55_65_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> ET2_65_Surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Oral_mucosa_tongue = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Oral_mucosa_mouth_floor = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Oral_mucosa_lips_and_cheeks = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Trachea = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 2> BB_11_6_ = {
	    massFraction(1, 0.112),
	    massFraction(8, 0.888)
	  };

	  static constexpr const std::array<massFraction, 9> BB_6_0_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_0_10_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_10_35_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_35_40_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_40_50_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_50_60_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_60_70_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> BB_70_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_arteries_head = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_veins_head = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_arteries_trunk = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_veins_trunk = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_arteries_arms = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_veins_arms = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_arteries_legs = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_large_veins_legs = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Humeri_upper_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Humeri_upper_spogiosa = {
	    massFraction(1, 0.081),
	    massFraction(6, 0.354),
	    massFraction(7, 0.028),
	    massFraction(8, 0.41),
	    massFraction(11, 0.002),
	    massFraction(12, 0.001),
	    massFraction(15, 0.037),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.083)
	  };

	  static constexpr const std::array<massFraction, 7> Humeri_upper_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Humeri_lower_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Humeri_lower_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 7> Humeri_lower_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Ulnae_and_radii_cotical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Ulnae_and_radii_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 7> Ulnae_and_radii_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Wrists_and_hand_bones_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Wrists_and_hand_bones_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 9> Clavicles_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 10> Clavicles_spongiosa = {
	    massFraction(1, 0.089),
	    massFraction(6, 0.409),
	    massFraction(7, 0.025),
	    massFraction(8, 0.385),
	    massFraction(11, 0.001),
	    massFraction(15, 0.027),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.06)
	  };

	  static constexpr const std::array<massFraction, 9> Cranium_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Cranium_spongiosa = {
	    massFraction(1, 0.088),
	    massFraction(6, 0.395),
	    massFraction(7, 0.026),
	    massFraction(8, 0.395),
	    massFraction(11, 0.001),
	    massFraction(12, 0.001),
	    massFraction(15, 0.028),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.062)
	  };

	  static constexpr const std::array<massFraction, 9> Femora_upper_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Femora_upper_spongiosa = {
	    massFraction(1, 0.093),
	    massFraction(6, 0.441),
	    massFraction(7, 0.023),
	    massFraction(8, 0.365),
	    massFraction(11, 0.001),
	    massFraction(12, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.05)
	  };

	  static constexpr const std::array<massFraction, 7> Femora_upper_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Femora_lower_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Femora_lower_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 7> Femora_lower_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Tibiae_fibulae_and_patellae_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Tibiae_fibulae_and_patellae_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 7> Tibiae_fibulae_and_patellae_medullary_cavity = {
	    massFraction(1, 0.115),
	    massFraction(6, 0.636),
	    massFraction(7, 0.007),
	    massFraction(8, 0.239),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Ankles_and_foot_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 9> Ankles_and_foot_spongiosa = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.504),
	    massFraction(7, 0.017),
	    massFraction(8, 0.308),
	    massFraction(11, 0.001),
	    massFraction(15, 0.022),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(20, 0.049)
	  };

	  static constexpr const std::array<massFraction, 9> Mandible_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Mandible_spongiosa = {
	    massFraction(1, 0.077),
	    massFraction(6, 0.332),
	    massFraction(7, 0.03),
	    massFraction(8, 0.42),
	    massFraction(11, 0.002),
	    massFraction(12, 0.001),
	    massFraction(15, 0.041),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.093)
	  };

	  static constexpr const std::array<massFraction, 9> Pelvis_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Pelvis_spongiosa = {
	    massFraction(1, 0.094),
	    massFraction(6, 0.409),
	    massFraction(7, 0.026),
	    massFraction(8, 0.4),
	    massFraction(11, 0.001),
	    massFraction(12, 0.001),
	    massFraction(15, 0.02),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.045)
	  };

	  static constexpr const std::array<massFraction, 9> Ribs_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 12> Ribs_spongiosa = {
	    massFraction(1, 0.088),
	    massFraction(6, 0.346),
	    massFraction(7, 0.031),
	    massFraction(8, 0.444),
	    massFraction(11, 0.001),
	    massFraction(12, 0.001),
	    massFraction(15, 0.026),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.058),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Scapulae_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Scapulae_spongiosa = {
	    massFraction(1, 0.084),
	    massFraction(6, 0.373),
	    massFraction(7, 0.027),
	    massFraction(8, 0.404),
	    massFraction(11, 0.001),
	    massFraction(12, 0.001),
	    massFraction(15, 0.033),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001),
	    massFraction(19, 0.001),
	    massFraction(20, 0.073)
	  };

	  static constexpr const std::array<massFraction, 9> Cervical_spine_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Cervical_spine_spongiosa = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.416),
	    massFraction(7, 0.028),
	    massFraction(8, 0.428),
	    massFraction(11, 0.001),
	    massFraction(15, 0.006),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(20, 0.012),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Thoracic_spine_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Thoracic_spine_spongiosa = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.403),
	    massFraction(7, 0.028),
	    massFraction(8, 0.431),
	    massFraction(11, 0.001),
	    massFraction(15, 0.01),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(20, 0.021),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Lumbar_spine_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Lumbar_spine_spongiosa = {
	    massFraction(1, 0.095),
	    massFraction(6, 0.38),
	    massFraction(7, 0.03),
	    massFraction(8, 0.436),
	    massFraction(11, 0.001),
	    massFraction(15, 0.016),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(20, 0.036),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Sacrum_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 11> Sacrum_spongiosa = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.426),
	    massFraction(7, 0.027),
	    massFraction(8, 0.426),
	    massFraction(11, 0.001),
	    massFraction(15, 0.003),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(20, 0.006),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Sternum_cortical = {
	    massFraction(1, 0.036),
	    massFraction(6, 0.159),
	    massFraction(7, 0.042),
	    massFraction(8, 0.448),
	    massFraction(11, 0.003),
	    massFraction(12, 0.002),
	    massFraction(15, 0.094),
	    massFraction(16, 0.003),
	    massFraction(20, 0.213)
	  };

	  static constexpr const std::array<massFraction, 10> Sternum_spongiosa = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.421),
	    massFraction(7, 0.028),
	    massFraction(8, 0.427),
	    massFraction(15, 0.005),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(20, 0.009),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Cartilage_costal = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.099),
	    massFraction(7, 0.022),
	    massFraction(8, 0.744),
	    massFraction(11, 0.005),
	    massFraction(15, 0.022),
	    massFraction(16, 0.009),
	    massFraction(17, 0.003)
	  };

	  static constexpr const std::array<massFraction, 8> Cartilage_discs = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.099),
	    massFraction(7, 0.022),
	    massFraction(8, 0.744),
	    massFraction(11, 0.005),
	    massFraction(15, 0.022),
	    massFraction(16, 0.009),
	    massFraction(17, 0.003)
	  };

	  static constexpr const std::array<massFraction, 9> Brain = {
	    massFraction(1, 0.107),
	    massFraction(6, 0.143),
	    massFraction(7, 0.023),
	    massFraction(8, 0.713),
	    massFraction(11, 0.002),
	    massFraction(15, 0.004),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 7> Breast_left_adipose_tissue = {
	    massFraction(1, 0.114),
	    massFraction(6, 0.581),
	    massFraction(7, 0.008),
	    massFraction(8, 0.294),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Breast_left_glandular_tissue = {
	    massFraction(1, 0.106),
	    massFraction(6, 0.324),
	    massFraction(7, 0.03),
	    massFraction(8, 0.535),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 7> Breast_right_adipose_tissue = {
	    massFraction(1, 0.114),
	    massFraction(6, 0.581),
	    massFraction(7, 0.008),
	    massFraction(8, 0.294),
	    massFraction(11, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Breast_right_glandular_tissue = {
	    massFraction(1, 0.106),
	    massFraction(6, 0.324),
	    massFraction(7, 0.03),
	    massFraction(8, 0.535),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Eye_lens_sensitive_left = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.195),
	    massFraction(7, 0.057),
	    massFraction(8, 0.646),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Eye_lens_insensitive_left = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.195),
	    massFraction(7, 0.057),
	    massFraction(8, 0.646),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Cornea_left = {
	    massFraction(1, 0.101),
	    massFraction(6, 0.125),
	    massFraction(7, 0.037),
	    massFraction(8, 0.732),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 4> Aqueous_left = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.004),
	    massFraction(7, 0.001),
	    massFraction(8, 0.883)
	  };

	  static constexpr const std::array<massFraction, 4> Vitreous_left = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.004),
	    massFraction(7, 0.001),
	    massFraction(8, 0.883)
	  };

	  static constexpr const std::array<massFraction, 8> Eye_lens_sensitive_right = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.195),
	    massFraction(7, 0.057),
	    massFraction(8, 0.646),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Eye_lens_insensitive_right = {
	    massFraction(1, 0.096),
	    massFraction(6, 0.195),
	    massFraction(7, 0.057),
	    massFraction(8, 0.646),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> Cornea_right = {
	    massFraction(1, 0.101),
	    massFraction(6, 0.125),
	    massFraction(7, 0.037),
	    massFraction(8, 0.732),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 4> Aqueous_right = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.004),
	    massFraction(7, 0.001),
	    massFraction(8, 0.883)
	  };

	  static constexpr const std::array<massFraction, 4> Vitreous_right = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.004),
	    massFraction(7, 0.001),
	    massFraction(8, 0.883)
	  };

	  static constexpr const std::array<massFraction, 9> Gall_bladder_wall = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Gall_bladder_contents = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.256),
	    massFraction(7, 0.027),
	    massFraction(8, 0.602),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Stomach_wall_0_60_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Stomach_wall_60_100_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Stomach_wall_100_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Stomach_wall_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Stomach_contents = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Small_intestine_wall_0_130_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Small_intestine_wall_130_150_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Small_intestine_wall_150_200_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Small_intestine_wall_200_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Small_intestine_contents_500_0_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Small_intestine_contents_centre_500_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Ascending_colon_wall_0_280_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Ascending_colon_wall_280_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Ascending_colon_wall_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Ascending_colon_content = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_right_0_280_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_right_280_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_right_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Transverse_colon_contents_right = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_left_0_280_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_left_280_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Transverse_colon_wall_left_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Transverse_colon_content_left = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Descending_colon_wall_0_280_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Descending_colon_wall_280_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Descending_colon_wall_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Descending_colon_content = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Sigmoid_colon_wall_0_280_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Sigmoid_colon_wall_280_300_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Sigmoid_colon_wall_300_surface_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Sigmoid_colon_contents = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Rectum_wall = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.114),
	    massFraction(7, 0.025),
	    massFraction(8, 0.75),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Heart_wall = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.135),
	    massFraction(7, 0.029),
	    massFraction(8, 0.722),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 10> Blood_in_heart_chamber = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.11),
	    massFraction(7, 0.033),
	    massFraction(8, 0.745),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_left_cortex = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_left_medulla = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_left_pelvis = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_right_cortex = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_right_medulla = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Kidney_right_pelvis = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.126),
	    massFraction(7, 0.031),
	    massFraction(8, 0.729),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Liver = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.132),
	    massFraction(7, 0.031),
	    massFraction(8, 0.723),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 10> Lung_AI__left = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.108),
	    massFraction(7, 0.032),
	    massFraction(8, 0.748),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 10> Lung_AI__right = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.108),
	    massFraction(7, 0.032),
	    massFraction(8, 0.748),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.002),
	    massFraction(26, 0.001)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_ET = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_thoracic = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_head = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_trunk = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_arms = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 7> Lymphatic_nodes_legs = {
	    massFraction(1, 0.108),
	    massFraction(6, 0.045),
	    massFraction(7, 0.012),
	    massFraction(8, 0.827),
	    massFraction(11, 0.003),
	    massFraction(16, 0.001),
	    massFraction(17, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Muscle_head = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Muscle_trunk = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Muscle_arms = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Muscle_legs = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Oesophagus_wall_0_190_ = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.223),
	    massFraction(7, 0.028),
	    massFraction(8, 0.635),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Oesophagus_wall_190_200_ = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.223),
	    massFraction(7, 0.028),
	    massFraction(8, 0.635),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Oesophagus_wall_200_surface_ = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.223),
	    massFraction(7, 0.028),
	    massFraction(8, 0.635),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 10> Oesophagus_contents = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Pancreas = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.158),
	    massFraction(7, 0.024),
	    massFraction(8, 0.704),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Pituitary_gland = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Prostate = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 8> RST_head = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.517),
	    massFraction(7, 0.011),
	    massFraction(8, 0.355),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> RST_trunk = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.517),
	    massFraction(7, 0.011),
	    massFraction(8, 0.355),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> RST_arms = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.517),
	    massFraction(7, 0.011),
	    massFraction(8, 0.355),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 8> RST_legs = {
	    massFraction(1, 0.112),
	    massFraction(6, 0.517),
	    massFraction(7, 0.011),
	    massFraction(8, 0.355),
	    massFraction(11, 0.001),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Salivary_glands_left = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Salivary_glandss_right = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_head_insensitive = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_head_sensitive_50_100_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_trunk_insensitive = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_trunk_sensitive_50_100_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_arms_insensitive = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_arms_sensitive_50_100_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_legs_insensitive = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Skin_legs_sensitive_50_100_ = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.199),
	    massFraction(7, 0.042),
	    massFraction(8, 0.65),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Spinal_cord = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Spleen = {
	    massFraction(1, 0.103),
	    massFraction(6, 0.112),
	    massFraction(7, 0.032),
	    massFraction(8, 0.743),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 7> Teeth = {
	    massFraction(1, 0.023),
	    massFraction(6, 0.095),
	    massFraction(7, 0.029),
	    massFraction(8, 0.426),
	    massFraction(12, 0.007),
	    massFraction(15, 0.135),
	    massFraction(20, 0.285)
	  };

	  static constexpr const std::array<massFraction, 10> Teeth_retention_region = {
	    massFraction(1, 0.1),
	    massFraction(6, 0.222),
	    massFraction(7, 0.022),
	    massFraction(8, 0.644),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004),
	    massFraction(20, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Testis_left = {
	    massFraction(1, 0.106),
	    massFraction(6, 0.099),
	    massFraction(7, 0.021),
	    massFraction(8, 0.765),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Testis_right = {
	    massFraction(1, 0.106),
	    massFraction(6, 0.099),
	    massFraction(7, 0.021),
	    massFraction(8, 0.765),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.002),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Thymus = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 10> Thyroid = {
	    massFraction(1, 0.104),
	    massFraction(6, 0.118),
	    massFraction(7, 0.025),
	    massFraction(8, 0.745),
	    massFraction(11, 0.002),
	    massFraction(15, 0.001),
	    massFraction(16, 0.001),
	    massFraction(17, 0.002),
	    massFraction(19, 0.001),
	    massFraction(53, 0.001)
	  };

	  static constexpr const std::array<massFraction, 9> Tongue_upper_food_ = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Tongue_lower = {
	    massFraction(1, 0.102),
	    massFraction(6, 0.142),
	    massFraction(7, 0.034),
	    massFraction(8, 0.711),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.001),
	    massFraction(19, 0.004)
	  };

	  static constexpr const std::array<massFraction, 9> Tonsils = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Ureter_left = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Ureter_right = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.251),
	    massFraction(7, 0.027),
	    massFraction(8, 0.607),
	    massFraction(11, 0.001),
	    massFraction(15, 0.002),
	    massFraction(16, 0.003),
	    massFraction(17, 0.002),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 9> Urinary_bladder_wall_insensitive = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.096),
	    massFraction(7, 0.026),
	    massFraction(8, 0.761),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 9> Urinary_bladder_wall_sensitive_75_118_ = {
	    massFraction(1, 0.105),
	    massFraction(6, 0.096),
	    massFraction(7, 0.026),
	    massFraction(8, 0.761),
	    massFraction(11, 0.002),
	    massFraction(15, 0.002),
	    massFraction(16, 0.002),
	    massFraction(17, 0.003),
	    massFraction(19, 0.003)
	  };

	  static constexpr const std::array<massFraction, 7> Urinary_bladder_content = {
	    massFraction(1, 0.107),
	    massFraction(6, 0.003),
	    massFraction(7, 0.01),
	    massFraction(8, 0.873),
	    massFraction(11, 0.004),
	    massFraction(15, 0.001),
	    massFraction(19, 0.002)
	  };

	  static constexpr const std::array<massFraction, 2> Air_inside_body = {
	    massFraction(7, 0.8),
	    massFraction(8, 0.2)
	  };


	  static constexpr std::array<const char*, 187> names = {
	    "Adrenal_left",
	    "Adrenal_right",
	    "ET1(0-8)",
	    "ET1(8-40)",
	    "ET1(40-50)",
	    "ET1(50-Surface)",
	    "ET2(-15-0)",
	    "ET2(0-40)",
	    "ET2(40-50)",
	    "ET2(50-55)",
	    "ET2(55-65)",
	    "ET2(65-Surface)",
	    "Oral_mucosa_tongue",
	    "Oral_mucosa_mouth_floor",
	    "Oral_mucosa_lips_and_cheeks",
	    "Trachea",
	    "BB(-11--6)",
	    "BB(-6-0)",
	    "BB(0-10)",
	    "BB(10-35)",
	    "BB(35-40)",
	    "BB(40-50)",
	    "BB(50-60)",
	    "BB(60-70)",
	    "BB(70-surface)",
	    "Blood_in_large_arteries_head",
	    "Blood_in_large_veins_head",
	    "Blood_in_large_arteries_trunk",
	    "Blood_in_large_veins_trunk",
	    "Blood_in_large_arteries_arms",
	    "Blood_in_large_veins_arms",
	    "Blood_in_large_arteries_legs",
	    "Blood_in_large_veins_legs",
	    "Humeri_upper_cortical",
	    "Humeri_upper_spogiosa",
	    "Humeri_upper_medullary_cavity",
	    "Humeri_lower_cortical",
	    "Humeri_lower_spongiosa",
	    "Humeri_lower_medullary_cavity",
	    "Ulnae_and_radii_cotical",
	    "Ulnae_and_radii_spongiosa",
	    "Ulnae_and_radii_medullary_cavity",
	    "Wrists_and_hand_bones_cortical",
	    "Wrists_and_hand_bones_spongiosa",
	    "Clavicles_cortical",
	    "Clavicles_spongiosa",
	    "Cranium_cortical",
	    "Cranium_spongiosa",
	    "Femora_upper_cortical",
	    "Femora_upper_spongiosa",
	    "Femora_upper_medullary_cavity",
	    "Femora_lower_cortical",
	    "Femora_lower_spongiosa",
	    "Femora_lower_medullary_cavity",
	    "Tibiae_fibulae_and_patellae_cortical",
	    "Tibiae_fibulae_and_patellae_spongiosa",
	    "Tibiae_fibulae_and_patellae_medullary_cavity",
	    "Ankles_and_foot_cortical",
	    "Ankles_and_foot_spongiosa",
	    "Mandible_cortical",
	    "Mandible_spongiosa",
	    "Pelvis_cortical",
	    "Pelvis_spongiosa",
	    "Ribs_cortical",
	    "Ribs_spongiosa",
	    "Scapulae_cortical",
	    "Scapulae_spongiosa",
	    "Cervical_spine_cortical",
	    "Cervical_spine_spongiosa",
	    "Thoracic_spine_cortical",
	    "Thoracic_spine_spongiosa",
	    "Lumbar_spine_cortical",
	    "Lumbar_spine_spongiosa",
	    "Sacrum_cortical",
	    "Sacrum_spongiosa",
	    "Sternum_cortical",
	    "Sternum_spongiosa",
	    "Cartilage_costal",
	    "Cartilage_discs",
	    "Brain",
	    "Breast_left_adipose_tissue",
	    "Breast_left_glandular_tissue",
	    "Breast_right_adipose_tissue",
	    "Breast_right_glandular_tissue",
	    "Eye_lens_sensitive_left",
	    "Eye_lens_insensitive_left",
	    "Cornea_left",
	    "Aqueous_left",
	    "Vitreous_left",
	    "Eye_lens_sensitive_right",
	    "Eye_lens_insensitive_right",
	    "Cornea_right",
	    "Aqueous_right",
	    "Vitreous_right",
	    "Gall_bladder_wall",
	    "Gall_bladder_contents",
	    "Stomach_wall(0-60)",
	    "Stomach_wall(60-100)",
	    "Stomach_wall(100-300)",
	    "Stomach_wall(300-surface)",
	    "Stomach_contents",
	    "Small_intestine_wall(0-130)",
	    "Small_intestine_wall(130-150)",
	    "Small_intestine_wall(150-200)",
	    "Small_intestine_wall(200-surface)",
	    "Small_intestine_contents(-500-0)",
	    "Small_intestine_contents(centre--500)",
	    "Ascending_colon_wall(0-280)",
	    "Ascending_colon_wall(280-300)",
	    "Ascending_colon_wall(300-surface)",
	    "Ascending_colon_content",
	    "Transverse_colon_wall_right(0-280)",
	    "Transverse_colon_wall_right(280-300)",
	    "Transverse_colon_wall_right(300-surface)",
	    "Transverse_colon_contents_right",
	    "Transverse_colon_wall_left(0-280)",
	    "Transverse_colon_wall_left(280-300)",
	    "Transverse_colon_wall_left(300-surface)",
	    "Transverse_colon_content_left",
	    "Descending_colon_wall(0-280)",
	    "Descending_colon_wall(280-300)",
	    "Descending_colon_wall(300-surface)",
	    "Descending_colon_content",
	    "Sigmoid_colon_wall(0-280)",
	    "Sigmoid_colon_wall(280-300)",
	    "Sigmoid_colon_wall(300-surface)",
	    "Sigmoid_colon_contents",
	    "Rectum_wall",
	    "Heart_wall",
	    "Blood_in_heart_chamber",
	    "Kidney_left_cortex",
	    "Kidney_left_medulla",
	    "Kidney_left_pelvis",
	    "Kidney_right_cortex",
	    "Kidney_right_medulla",
	    "Kidney_right_pelvis",
	    "Liver",
	    "Lung(AI)_left",
	    "Lung(AI)_right",
	    "Lymphatic_nodes_ET",
	    "Lymphatic_nodes_thoracic",
	    "Lymphatic_nodes_head",
	    "Lymphatic_nodes_trunk",
	    "Lymphatic_nodes_arms",
	    "Lymphatic_nodes_legs",
	    "Muscle_head",
	    "Muscle_trunk",
	    "Muscle_arms",
	    "Muscle_legs",
	    "Oesophagus_wall(0-190)",
	    "Oesophagus_wall(190-200)",
	    "Oesophagus_wall(200-surface)",
	    "Oesophagus_contents",
	    "Pancreas",
	    "Pituitary_gland",
	    "Prostate",
	    "RST_head",
	    "RST_trunk",
	    "RST_arms",
	    "RST_legs",
	    "Salivary_glands_left",
	    "Salivary_glandss_right",
	    "Skin_head_insensitive",
	    "Skin_head_sensitive(50-100)",
	    "Skin_trunk_insensitive",
	    "Skin_trunk_sensitive(50-100)",
	    "Skin_arms_insensitive",
	    "Skin_arms_sensitive(50-100)",
	    "Skin_legs_insensitive",
	    "Skin_legs_sensitive(50-100)",
	    "Spinal_cord",
	    "Spleen",
	    "Teeth",
	    "Teeth_retention_region",
	    "Testis_left",
	    "Testis_right",
	    "Thymus",
	    "Thyroid",
	    "Tongue_upper(food)",
	    "Tongue_lower",
	    "Tonsils",
	    "Ureter_left",
	    "Ureter_right",
	    "Urinary_bladder_wall_insensitive",
	    "Urinary_bladder_wall_sensitive(75-118)",
	    "Urinary_bladder_content",
	    "Air_inside_body",
	  };

	  inline std::vector<std::string> matList() const override{
	    return std::vector<std::string>(names.cbegin(), names.cend());
	  }
	  inline unsigned getIndex(const std::string& name) const override{
	    for(size_t i = 0; i < names.size(); ++i){
	      if(name.compare(names[i]) == 0){
		return i;
	      }
	    }
	    return name.size();
	  }
	  std::vector<massFraction> getElements(const unsigned index) const override;
	  double getDensity(const unsigned index) const override;
	  
	};

      } //namespace tissues

    } //namespace compositions
    
  } //namespace dataBases

} //namespace penred
#endif
