# Deployability of a ZVS-Constrained Optimal Modulation LUT for PI Current Control of a Wide-Voltage-Range DAB

Optimal triple-phase-shift (TPS) modulation of a dual active bridge (DAB) is usually computed offline, stored in a look-up table (LUT) and interpolated between its stored points. This project asks whether a plain PI current controller can run such a LUT along a 200–1000 V battery-charging profile. The LUT keeps zero-voltage switching (ZVS) on the commutations and, among those duty triplets, uses the one with the lowest RMS current.

## Contributions

1. **Command-to-current gain of an optimal-modulation LUT.** An ideal inverse modulation delivers exactly the commanded current, so its gain from command to delivered current is unity. A stored, bilinearly interpolated LUT does not. The controller samples the current once per switching period, so the gain it actually sees is the sampled gain g_s, evaluated on the exact five-state model of the converter. The gain decomposition shows that g_s can leave unity only in LUT cells that span a jump of the modulation between operating modes.

2. **g_s as a deployability metric for a fixed PI controller.** An interpolation error is only a bias, which the integral action removes. A negative g_s inverts the sign of the feedback, and the controller then amplifies the error. The widest reversal along the charging ramp limits the usable crossover frequency, and the spread of g_s sets the margins the PI is tuned for. With one PI tuned per voltage range, the metric shows that the ZVS-constrained LUT is controllable at light and heavy load and not at medium load.

## Converter

| Parameter | Value |
|---|---|
| Input voltage V_i | 800 V |
| Output voltage V_2 | 200–1000 V |
| Turns ratio n | 1.7 |
| AC-link inductance L_s / resistance R_s | 15 µH / 25 mΩ |
| Input filter L_i / R_i / C_i | 10 µH / 0.1 Ω / 30 µF |
| Output filter L_o / R_o / C_o | 10 µH / 0.1 Ω / 44 µF |
| Switching frequency f_s | 50 kHz |
| Charging profile | I_ref = min(122 A, 50 kW / V_2), ramp 8 V/ms |

## How the LUT is built

**Grid.** 40 current × 17 voltage operating points, I_cmd = 5–122.5 A and V_2 = 200–1000 V in 50 V steps. At each point, D1 and D2 are scanned on a 70 × 70 grid. For each (D1, D2) pair, every outer shift D3 that delivers the target power P = I_cmd · V_2 is found by bisection.

**Analytical model.** Each candidate is evaluated on the ideal TPS waveform of the AC-link current, with these assumptions:
- Lossless converter, constant DC voltages, filters not included.
- Secondary referred to the primary (n · V_2).
- The link current is piecewise linear, with slope (v_p − v_s)/L_s on each interval between the eight sorted switching edges, and zero mean over the period.
- Power and RMS current follow in closed form from this waveform.

**ZVS criterion.** It is checked separately at each of the eight edges:
- The link inductance and the output capacitances of the commutating leg form a resonant circuit, with characteristic impedance Z = √(L_s / 2C_oss) on the primary side and n·Z on the secondary side.
- C_oss = 1 nF is used as a constant, charge-equivalent value.
- The transition completes if the link current at the edge flows in the right direction and its magnitude reaches √((V_b − v_o)² − (V_a − v_o)²) / Z. Here V_a and V_b are the bridge voltages before and after the edge, and v_o is the voltage of the opposite bridge during the transition.
- Dead time is not modelled.

**Selection.** Among the triplets with the highest number of ZVS edges (eight wherever possible), the one with the lowest RMS current is stored. The LUT is then evaluated on the exact five-state model of the converter, including the filters and resistances, to compute the gain the controller sees.

## Sampled gain along the charging profile

![gain profile](figures/gain_profile.png)

*Figure 1. Command-to-current gain of the interpolated LUT along the charging profile, I_ref = min(122 A, 50 kW / V_2). The dashed line marks zero gain, and the dotted lines mark 500 and 700 V. Below zero, the delivered current moves opposite to the command.*

The gain stays close to unity at light and heavy load. At medium load, between 500 and 700 V, it reverses in three stretches, down to about −23. The widest reversal is 25.4 V, which the ramp crosses in 3.2 ms, so the crossover frequency of a PI controller in that range is limited to 50 Hz.

## PI controller per voltage range

Each range gets its own PI controller, tuned only on the gain inside that range. The plant is scaled to the median positive g_s of the range. The controller is accepted when the loop scaled to the 10th and 90th percentiles of g_s keeps PM ≥ 45° and GM ≥ 6 dB, and when its crossover respects the reversal limit.

| Range | Load | g_s median | g_s p10 / p90 | Reversals | f_c | Kp | Ki | PM min | GM min | Ramp worst error | Verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 200–500 V | light | 1.000 | 0.988 / 1.026 | 0 | 300 Hz | 0.01882 | 1882.2 | 86.1° | 18.5 dB | – | controllable |
| 500–700 V | medium | 0.997 | 0.972 / 1.014 | 3 | 50 Hz | 0.00315 | 315.0 | 89.4° | 34.0 dB | 69.0 A | PI fails |
| 700–1000 V | heavy | 0.998 | 0.961 / 1.013 | 0 | 300 Hz | 0.01887 | 1887.1 | 86.2° | 18.7 dB | 0.40 A | controllable |

The margins are met in every range. At medium load the PI controller still fails: on the ramp, the correction saturates at ±20 A in 163 of 1250 switching periods and the tracking error reaches 69.0 A.

![PI ramp 500-700 V](figures/pi_ramp_500_700.png)

*Figure 2. Medium load: V_2 ramped from 500 to 700 V at 8 V/ms on the exact five-state model. Red: LUT in open loop. Blue: LUT with the PI controller of this range (f_c = 50 Hz). Dashed: reference I_ref. Inside the reversal stretches the PI controller drives the current away from the reference and oscillates.*

![PI ramp 700-1000 V](figures/pi_ramp_700_1000.png)

*Figure 3. Heavy load: V_2 ramped from 700 to 1000 V at 8 V/ms on the exact five-state model. Red: LUT in open loop. Blue: LUT with the PI controller of this range (f_c = 300 Hz). Dashed: reference I_ref. The current is tracked with a worst error of 0.40 A.*

## Efficiency against SPS

The losses are estimated from the Wolfspeed CAB320M17XM3 SiC module datasheet at T_vj = 125 °C. The module is used only in this efficiency estimate.
- Conduction losses use R_DS(on) plus the package resistance.
- A ZVS commutation costs E_off.
- A hard commutation costs E_on (not below E_oss) plus E_rr.

Single phase shift (SPS) modulation is evaluated with the same model at the same power.

![efficiency](figures/efficiency_estimate.png)

*Figure 4. Estimated efficiency along the charging profile. Blue: ZVS-constrained optimal modulation LUT. Dashed orange: SPS at the same power.*

At light and heavy load the ZVS-constrained LUT is far more efficient than SPS. SPS is efficient only near the matched voltage, around 470 V, and drops steeply towards 1000 V. At medium load, around 570–650 V, the ZVS-constrained LUT falls below SPS. Holding ZVS on all eight commutations costs more there than it saves.

## Simulink validation

The switched circuit model uses ideal switches (R_on = 1 µΩ, no dead time) and the ode23t solver with a maximum step of 0.5 µs. It runs the full ramp of V_2 from 200 to 1000 V in 0.1 s, with the LUT once in closed loop and once in open loop.

![Simulink closed loop](figures/closedloopsimulink.png)

*Figure 5. Simulink, closed loop with the PI controller tuned for 500–700 V (f_c = 50 Hz), time in s. Blue: delivered current i_2. Yellow: reference I_ref. The current is tracked at light and heavy load. In the medium-load band (about 45–56 ms, 560–650 V) the loop is lost and oscillates.*

![Simulink open loop](figures/openloopsimulink.png)

*Figure 6. Simulink, LUT in open loop (no controller), time in s. Yellow: delivered current i_2. Blue: reference I_ref. The two peaks between 45 and 56 ms are the misdelivery of the interpolated LUT in the medium-load band. Elsewhere the LUT delivers the reference.*

The switched model reproduces the offline result: the current is tracked at light and heavy load, and the loop is lost in the same medium-load band.



## Conclusion

This project shows that the ZVS constraint is cheap at low and high output voltage, and expensive in the medium voltage range. Three independent results point to the same picture: the sampled gain g_s, the PI controller on the charging ramp, and the efficiency against SPS.

**Low (200–500 V) and high (700–1000 V) voltage range, light and heavy load.** Here the minimum-RMS triplet with ZVS on all eight commutations varies smoothly with the command, so the LUT behaves like an ideal inverse modulation. The sampled gain stays close to unity, between 0.96 and 1.03, with no reversal along the charging profile. A PI controller tuned on this gain reaches a crossover of 300 Hz with PM ≥ 86° and GM ≥ 18.5 dB. On the 700–1000 V ramp it tracks the reference with a worst error of 0.40 A, and the switched Simulink model tracks the current over the whole low and high range. Because every commutation is soft, the estimated efficiency stays between about 97 and 98.5 %. SPS reaches that level only near the matched voltage, around 470 V, and falls to about 83 % at 1000 V.

**Medium voltage range (500–700 V, medium load).** Here ZVS on all eight commutations can be held only by triplets with a high circulating current, and the minimum-RMS selection jumps between them from one operating point to the next. The efficiency drops to about 92.5 %, below SPS, so holding full ZVS costs more than it saves. The interpolated LUT spans these jumps and the sampled gain reverses in three stretches, down to about −23. A command correction then moves the current the wrong way, and on the 500–700 V ramp the PI controller saturates in 163 of 1250 periods with a tracking error of 69.0 A. The switched Simulink model loses the current in the same band.

**The g_s metric.** g_s identified this failure from the stored LUT alone, before any closed-loop test. It located the reversal stretches along the charging profile, and the widest one (25.4 V, crossed in 3.2 ms) limited the crossover of any PI controller in that range to 50 Hz. Its spread gave the margins the controllers were tuned for in the other two ranges. The ramps and the Simulink model then confirmed exactly what g_s predicted: tracking at light and heavy load, loss of the current at medium load. This makes g_s a practical check of whether an optimal-modulation LUT can be run by a fixed PI controller, and of where it cannot.

The result points to the next step: in the medium voltage range, the ZVS constraint should be relaxed, so that the modulation keeps a low RMS current and a sampled gain close to unity.

## How to run

MATLAB scripts, in this order:

1. `zvsconstrainedoptimalmodulation.m` builds the LUT and saves `pointwise_fullTPS_NG70.mat`. This step is slow, and the LUT is already included.
2. `gain_reversal.m` computes the gain along the charging profile and saves `figures/gain_profile.png`.
3. `controller_ramp.m` tunes one PI controller per voltage range, runs the 500–700 V and 700–1000 V ramps, and saves `controller_ramp.mat` and the `figures/pi_ramp_*.png` plots.
4. `efficiency_estimate.m` compares the LUT with SPS and saves `figures/efficiency_estimate.png`.

`dab_params.m` holds the converter parameters. `load_cab320.m` reads the digitized datasheet curves (`.csv` files) for the efficiency estimate.

Simulink:

1. Run `controller_ramp.m` once, so that `controller_ramp.mat` exists.
2. Run `init_simulink.m`. `BAND` selects which range's PI controller is loaded (1, 2 or 3; 2 was used for Figure 5).
3. Open and run `closed_TPScontrol_varduty.slx` or `openloopTPScontrol_varduty.slx`.

## Credits

Switch data digitized from the Wolfspeed CAB320M17XM3 datasheet.
