# Generated results

Produced by `python scripts/run_analysis.py`. Do not edit by hand.

## Drive requirements
| Quantity | Value |
|---|---|
| Assist cut-off speed | 25 km/h |
| Gear ratio (motor : wheel) | 8 : 1 |
| Motor speed at cut-off (base) | 1608 rpm |
| Corner speed (end of constant torque) | 805 rpm |
| Rated torque (250 W @ base) | 1.48 N m |
| Cruise torque (flat, 25 km/h) | 0.88 N m (149 W) |
| Hill torque (6 %, 15 km/h) | 3.03 N m |
| Launch torque (8 %, 0.5 m/s^2) | 5.93 N m |
| Design peak torque | 5.93 N m |

## Pole-configuration comparison (at max. speed)
| config | phases | stroke_deg | strokes_per_rev | f_phase_at_nmax_Hz | switches_asym_bridge | self_start_both_dir |
|---|---|---|---|---|---|---|
| 6/4 | 3 | 30.0 | 12 | 128.6 | 6 | True |
| 8/6 | 4 | 15.0 | 24 | 192.9 | 8 | True |
| 12/8 | 3 | 15.0 | 24 | 257.2 | 6 | True |
| 16/12 | 4 | 7.5 | 48 | 385.8 | 8 | True |
| 24/16 | 3 | 7.5 | 48 | 514.4 | 6 | True |

## Machine design
| Quantity | Value |
|---|---|
| configuration | 12/8, 3-phase |
| outer_diameter_mm | 135.3 |
| bore_diameter_mm | 71.7 |
| stack_length_mm | 35.9 |
| air_gap_mm | 0.3 |
| beta_s_deg | 15.0 |
| beta_r_deg | 17.0 |
| stator_pole_width_mm | 9.36 |
| rotor_pole_width_mm | 10.51 |
| stator_pole_height_mm | 24.78 |
| rotor_pole_height_mm | 9.0 |
| stator_yoke_mm | 7.02 |
| rotor_yoke_mm | 19.06 |
| turns_per_phase | 164 |
| turns_per_pole | 41 |
| wire_diameter_mm | 1.659 |
| slot_area_mm2 | 394.0 |
| mean_turn_length_mm | 115.4 |
| phase_resistance_mohm | 197.9 |
| L_aligned_mH | 9.055 |
| L_unaligned_mH | 1.371 |
| L_sat_mH | 1.508 |
| inductance_ratio | 6.6 |
| psi_knee_mWb | 85.33 |
| mass_iron_kg | 2.307 |
| mass_copper_kg | 1.093 |
| mass_active_kg | 3.4 |
| rotor_inertia_kgm2 | 0.0004133 |

## Operating points
| point | speed_rpm | theta_on_deg | theta_off_deg | i_ref_A | torque_avg_Nm | torque_ripple_pct | i_peak_A | i_rms_A | psi_peak_mWb | b_pole_peak_T | f_switch_kHz | p_out_W | eta_motor_pct | eta_system_pct | continuous_conduction |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| launch (peak) | 80.5 | 4.52 | 21.68 | 30.67 | 6.813 | 113.3 | 31.43 | 18.92 | 126.51 | 2.0 | 3.62 | 55.9 | 20.7 | 18.8 | False |
| hill 6 % | 964.6 | -0.54 | 17.56 | 20.45 | 3.025 | 144.7 | 20.73 | 12.07 | 81.39 | 1.478 | 2.44 | 287.9 | 73.4 | 70.8 | False |
| rated | 1607.6 | -4.91 | 17.56 | 15.44 | 1.481 | 147.7 | 15.66 | 9.48 | 56.57 | 1.028 | 3.64 | 229.5 | 75.8 | 72.8 | False |
| cruise | 1607.6 | -4.91 | 17.56 | 10.2 | 0.882 | 142.2 | 10.35 | 6.83 | 48.63 | 0.883 | 5.79 | 133.1 | 75.6 | 71.6 | False |
| overspeed | 1929.2 | -6.35 | 17.56 | 9.82 | 0.742 | 143.9 | 9.97 | 6.59 | 43.45 | 0.789 | 5.92 | 133.0 | 75.7 | 71.8 | False |

## Converter ratings
| Quantity | Value |
|---|---|
| topology | asymmetric half bridge, 3 legs, 6 MOSFETs + 6 diodes |
| device_voltage_rating_V | 63.0 |
| selected_device_voltage_V | 100 |
| phase_current_limit_A | 30.7 |
| device_rms_current_peak_A | 18.92 |
| battery_current_peak_A | 8.17 |
| battery_current_rated_A | 8.74 |
| dc_link_ripple_current_rms_peak_A | 15.42 |
| dc_link_ripple_current_rms_rated_A | 10.89 |
| dc_link_capacitance_min_uF | 3782.0 |
| switching_frequency_rated_kHz | 3.64 |

## Thermal estimate (continuous rated point)
| Quantity | Value |
|---|---|
| motor_loss_W | 73.2 |
| cooling_area_m2 | 0.0749 |
| temperature_rise_K | 39.1 |
| winding_temp_estimate_C | 76.9 |
| current_density_rated_A_mm2 | 4.39 |

## Route simulation
| Quantity | Value |
|---|---|
| Distance | 6.40 km |
| Duration | 18.9 min |
| Rider energy | 25.5 Wh |
| Motor shaft energy | 23.4 Wh |
| Battery energy | 32.3 Wh |
| Cycle-average drive efficiency | 72.3 % |
| Consumption | 5.05 Wh/km |
| Estimated range (90 % of 360 Wh) | 64 km |
| Time with torque deficit | 0.0 s |

## Trade study: gear ratio and electric loading
| gear_ratio | A_peak_kA_per_m | D_out_mm | L_stk_mm | active_mass_kg | I_peak_A | eta_motor_rated_pct | eta_system_rated_pct |
|---|---|---|---|---|---|---|---|
| 5 | 60.0 | 143.8 | 38.1 | 4.102 | 27.7 | 70.6 | 68.5 |
| 5 | 45.0 | 158.3 | 41.9 | 5.514 | 29.2 | 76.3 | 73.6 |
| 5 | 35.0 | 172.1 | 45.6 | 7.136 | 31.2 | 81.0 | 78.0 |
| 8 | 60.0 | 122.9 | 32.6 | 2.526 | 29.3 | 70.6 | 68.1 |
| 8 | 45.0 | 135.3 | 35.9 | 3.4 | 30.7 | 75.8 | 72.8 |
| 8 | 35.0 | 147.1 | 39.0 | 4.404 | 32.0 | 79.9 | 76.6 |
| 11 | 60.0 | 110.6 | 29.3 | 1.816 | 30.3 | 69.9 | 67.3 |
| 11 | 45.0 | 121.7 | 32.2 | 2.447 | 31.6 | 74.8 | 71.6 |
| 11 | 35.0 | 132.3 | 35.1 | 3.172 | 33.7 | 80.2 | 77.0 |
