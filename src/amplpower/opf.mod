########## SETS ##########

set N;
set G;
set L;

########## PARAMETERS ##########

param OPF_TYPE symbolic;
param CONNECTIVITY symbolic;
param BASEMVA;
param MAXANGLE;
param MINANGLE;
param COST_2 {G};
param COST_1 {G};
param COST_0 {G};
param PG {G};
param QG {G};
param PMAX {G};
param PMIN {G};
param QMAX {G};
param QMIN {G};
param GEN_BUS {G};
param PD {N};
param QD {N};
param GS {N};
param BS {N};
param VMAX {N};
param VMIN {N};
param AMAX {N};
param AMIN {N};
param VOL0 {N};
param ANG0 {N};
param VOLR0 {N};
param VOLI0 {N};
param CF {L,N};
param CT {L,N};
param CG {G,N};
param GFF {L};
param BFF {L};
param GFT {L};
param BFT {L};
param GTF {L};
param BTF {L};
param GTT {L};
param BTT {L};
param BUS_I {N};
param BUS_TYPE {N};
param BUS_AREA {N};
param VM {N};
param VA {N};
param BASE_KV {N};
param ZONE {N};
param F_BUS {L};
param T_BUS {L};
param BR_R {L};
param BR_X {L};
param BR_B {L};
param RATE_A {L};
param RATE_B {L};
param RATE_C {L};
param PFMAX {L};
param PFMIN {L};
param QFMAX {L};
param QFMIN {L};
param TAP {L};
param SHIFT {L};
param BR_STATUS {L};
param BR_SWITCH {L};
param ANGMIN {L};
param ANGMAX {L};
param VG {G};
param MBASE {G};
param GEN_STATUS {G};
param PC1 {G};
param PC2 {G};
param QC1MIN {G};
param QC1MAX {G};
param QC2MIN {G};
param QC2MAX {G};
param RAMP_AGC {G};
param RAMP_10 {G};
param RAMP_30 {G};
param RAMP_Q {G};
param APF {G};
param MODEL {G};
param STARTUP {G};
param SHUTDOWN {G};
param NCOST {G};
param PG0 {G};
param QG0 {G};
param PF0 {L};
param QF0 {L};
param PT0 {L};
param QT0 {L};
param PFLODC {L};
param PFUPDC {L};
param PFLOAC {L};
param PFUPAC {L};
param PTLOAC {L};
param PTUPAC {L};
param QFLOAC {L};
param QFUPAC {L};
param QTLOAC {L};
param QTUPAC {L};
param COSFTMAX {L};
param COSFTMIN {L};
param SINFTMAX {L};
param SINFTMIN {L};
param VRMAX {N};
param VRMIN {N};
param VIMAX {N};
param VIMIN {N};

########## VARIABLES ##########
var Pg {g in G} >= PMIN[g], <= PMAX[g]:= PG0[g];
var Qg {g in G} >= QMIN[g], <= QMAX[g]:= QG0[g];
var Pf {l in L} := PF0[l];
var Pt {l in L} := PT0[l];
var Qf {l in L} := QF0[l];
var Qt {l in L} := QT0[l];
var Vm {n in N} >= VMIN[n], <= VMAX[n] := VOL0[n];
var Vr {n in N} >= VRMIN[n], <=VRMAX[n] := VOLR0[n];
var Vi {n in N} >= VIMIN[n], <=VIMAX[n] := VOLI0[n];
var V2 {n in N} >= VMIN[n]^2, <= VMAX[n]^2 := VOL0[n]^2;
var Va {n in N} >= AMIN[n], <= AMAX[n] := ANG0[n];
var cosft {l in L} >= COSFTMIN[l], <= COSFTMAX[l] := VOL0[F_BUS[l]]*VOL0[T_BUS[l]]*cos(ANG0[F_BUS[l]]-ANG0[T_BUS[l]]);
var sinft {l in L} >= SINFTMIN[l], <= SINFTMAX[l] := VOL0[F_BUS[l]]*VOL0[T_BUS[l]]*sin(ANG0[F_BUS[l]]-ANG0[T_BUS[l]]);
var Pfa {l in L} := PF0[l];
var Pta {l in L} := PT0[l];
var Qfa {l in L} := QF0[l];
var Qta {l in L} := QT0[l];
var u {n in N} >= 1, <= card(N);
var status {l in L} binary;
var statusf {l in L} binary;
var statust {l in L} binary;

########## POWER BALANCE ##########

subject to active_power_balance_1 {n in N:OPF_TYPE='dc'}:
    sum {g in G} CG[g,n] * Pg[g] - PD[n] = sum {l in L} (CF[l,n] * Pf[l] + CT[l,n] * Pt[l]);

subject to active_power_balance_2 {n in N:OPF_TYPE='acpolar'}:
    sum {g in G} CG[g,n] * Pg[g] - PD[n] = GS[n]*Vm[n]*Vm[n] + sum {l in L} (CF[l,n] * Pf[l] + CT[l,n] * Pt[l]);

subject to reactive_power_balance_2 {n in N:OPF_TYPE='acpolar'}:
    sum {g in G} CG[g,n] * Qg[g] - QD[n] = -BS[n]*Vm[n]*Vm[n] + sum {l in L} (CF[l,n] * Qf[l] + CT[l,n] * Qt[l]);

subject to active_power_balance_3_4 {n in N: OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    sum {g in G} CG[g,n] * Pg[g] - PD[n] = GS[n]*V2[n] + sum {l in L} (CF[l,n] * Pf[l] + CT[l,n] * Pt[l]);

subject to reactive_power_balance_3_4 {n in N: OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    sum {g in G} CG[g,n] * Qg[g] - QD[n] = -BS[n]*V2[n] + sum {l in L} (CF[l,n] * Qf[l] + CT[l,n] * Qt[l]);

########## POWER FLOW DEFINITIONS (BR_SWITCH = 0) ##########

subject to active_flow_from_0 {l in L:BR_SWITCH[l] == 0}:
    Pf[l] = 0;

subject to active_flow_to_0 {l in L:BR_SWITCH[l] == 0}:
    Pt[l] = 0;

subject to reactive_flow_from_0 {l in L:BR_SWITCH[l] == 0}:
    Qf[l] = 0;

subject to reactive_flow_to_0 {l in L:BR_SWITCH[l] == 0}:
    Qt[l] = 0;

########## POWER FLOW DEFINITIONS (BR_SWITCH = 1) ##########

subject to active_flow_from_1 {l in L:OPF_TYPE='dc' and BR_SWITCH[l] == 1}:
    Pf[l] = (1 / BR_X[l]) * (Va[F_BUS[l]] - Va[T_BUS[l]]);

subject to active_flow_to_1 {l in L:OPF_TYPE='dc' and BR_SWITCH[l] == 1}:
    Pt[l] = (1 / BR_X[l]) * (Va[T_BUS[l]] - Va[F_BUS[l]]);

subject to active_flow_from_2 {l in L:OPF_TYPE='acpolar' and BR_SWITCH[l] == 1}:
    Pf[l] = GFF[l]*Vm[F_BUS[l]]*Vm[F_BUS[l]] + Vm[F_BUS[l]]*Vm[T_BUS[l]]*(GFT[l]*cos(Va[F_BUS[l]]-Va[T_BUS[l]])+BFT[l]*sin(Va[F_BUS[l]]-Va[T_BUS[l]]));

subject to active_flow_to_2 {l in L:OPF_TYPE='acpolar' and BR_SWITCH[l] == 1}:
    Pt[l] = GTT[l]*Vm[T_BUS[l]]*Vm[T_BUS[l]] + Vm[F_BUS[l]]*Vm[T_BUS[l]]*(GTF[l]*cos(Va[T_BUS[l]]-Va[F_BUS[l]])+BTF[l]*sin(Va[T_BUS[l]]-Va[F_BUS[l]]));

subject to reactive_flow_from_2 {l in L:OPF_TYPE='acpolar' and BR_SWITCH[l] == 1}:
    Qf[l] = -BFF[l]*Vm[F_BUS[l]]*Vm[F_BUS[l]] - Vm[F_BUS[l]]*Vm[T_BUS[l]]*(BFT[l]*cos(Va[F_BUS[l]]-Va[T_BUS[l]])-GFT[l]*sin(Va[F_BUS[l]]-Va[T_BUS[l]]));

subject to reactive_flow_to_2 {l in L:OPF_TYPE='acpolar' and BR_SWITCH[l] == 1}:
    Qt[l] = -BTT[l]*Vm[T_BUS[l]]*Vm[T_BUS[l]] - Vm[F_BUS[l]]*Vm[T_BUS[l]]*(BTF[l]*cos(Va[T_BUS[l]]-Va[F_BUS[l]])-GTF[l]*sin(Va[T_BUS[l]]-Va[F_BUS[l]]));

subject to active_flow_from_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_SWITCH[l] == 1}:
    Pf[l] = GFF[l]*V2[F_BUS[l]] + GFT[l]*cosft[l] + BFT[l]*sinft[l];

subject to active_flow_to_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_SWITCH[l] == 1}:
    Pt[l] = GTT[l]*V2[T_BUS[l]] + GTF[l]*cosft[l] - BTF[l]*sinft[l];

subject to reactive_flow_from_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_SWITCH[l] == 1}:
    Qf[l] = -BFF[l]*V2[F_BUS[l]] - BFT[l]*cosft[l] + GFT[l]*sinft[l];

subject to reactive_flow_to_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_SWITCH[l] == 1}:
    Qt[l] = -BTT[l]*V2[T_BUS[l]] - BTF[l]*cosft[l] - GTF[l]*sinft[l];

########## POWER FLOW DEFINITIONS (BR_SWITCH = 2) ##########

subject to active_flow_from_1_switch {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 2}:
    Pf[l] = status[l] * (1 / BR_X[l]) * (Va[F_BUS[l]] - Va[T_BUS[l]]);

subject to active_flow_to_1_switch {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 2}:
    Pt[l] = status[l] * (1 / BR_X[l]) * (Va[T_BUS[l]] - Va[F_BUS[l]]);

subject to active_flow_from_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_SWITCH[l] == 2}:
    Pf[l] = status[l] * (GFF[l]*Vm[F_BUS[l]]*Vm[F_BUS[l]] + Vm[F_BUS[l]]*Vm[T_BUS[l]]*(GFT[l]*cos(Va[F_BUS[l]]-Va[T_BUS[l]])+BFT[l]*sin(Va[F_BUS[l]]-Va[T_BUS[l]])));

subject to active_flow_to_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_SWITCH[l] == 2}:
    Pt[l] = status[l] * (GTT[l]*Vm[T_BUS[l]]*Vm[T_BUS[l]] + Vm[F_BUS[l]]*Vm[T_BUS[l]]*(GTF[l]*cos(Va[T_BUS[l]]-Va[F_BUS[l]])+BTF[l]*sin(Va[T_BUS[l]]-Va[F_BUS[l]])));

subject to reactive_flow_from_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_SWITCH[l] == 2}:
    Qf[l] = status[l] * (-BFF[l]*Vm[F_BUS[l]]*Vm[F_BUS[l]] - Vm[F_BUS[l]]*Vm[T_BUS[l]]*(BFT[l]*cos(Va[F_BUS[l]]-Va[T_BUS[l]])-GFT[l]*sin(Va[F_BUS[l]]-Va[T_BUS[l]])));

subject to reactive_flow_to_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_SWITCH[l] == 2}:
    Qt[l] = status[l] * (-BTT[l]*Vm[T_BUS[l]]*Vm[T_BUS[l]] - Vm[F_BUS[l]]*Vm[T_BUS[l]]*(BTF[l]*cos(Va[T_BUS[l]]-Va[F_BUS[l]])-GTF[l]*sin(Va[T_BUS[l]]-Va[F_BUS[l]])));

subject to active_flow_from_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 2}:
    Pf[l] = status[l] * (GFF[l]*V2[F_BUS[l]] + GFT[l]*cosft[l] + BFT[l]*sinft[l]);

subject to active_flow_to_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 2}:
    Pt[l] = status[l] * (GTT[l]*V2[T_BUS[l]] + GTF[l]*cosft[l] - BTF[l]*sinft[l]);

subject to reactive_flow_from_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 2}:
    Qf[l] = status[l] * (-BFF[l]*V2[F_BUS[l]] - BFT[l]*cosft[l] + GFT[l]*sinft[l]);

subject to reactive_flow_to_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 2}:
    Qt[l] = status[l] * (-BTT[l]*V2[T_BUS[l]] - BTF[l]*cosft[l] - GTF[l]*sinft[l]);

########## POWER FLOW DEFINITIONS (BR_SWITCH = 3) ##########

subject to active_flow_from_1_bigm {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pfa[l] = (1 / BR_X[l]) * (Va[F_BUS[l]] - Va[T_BUS[l]]);

subject to active_flow_to_1_bigm {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pta[l] = (1 / BR_X[l]) * (Va[T_BUS[l]] - Va[F_BUS[l]]);

subject to Pfa_lower_1 {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    PFLODC[l] * (1 - status[l]) <= -Pf[l] + Pfa[l];

subject to Pfa_upper_1 {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    -Pf[l] + Pfa[l] <= PFUPDC[l] * (1 - status[l]);

subject to Pta_lower_1 {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    -PFUPDC[l] * (1 - status[l]) <= -Pt[l] + Pta[l];

subject to Pta_upper_1 {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    -Pt[l] + Pta[l] <= -PFLODC[l] * (1 - status[l]);

subject to active_flow_from_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Pfa[l] = GFF[l] * V2[F_BUS[l]] + GFT[l] * cosft[l] + BFT[l] * sinft[l];

subject to active_flow_to_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Pta[l] = GTT[l] * V2[T_BUS[l]] + GTF[l] * cosft[l] - BTF[l] * sinft[l];

subject to reactive_flow_from_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Qfa[l] = -BFF[l] * V2[F_BUS[l]] - BFT[l] * cosft[l] + GFT[l] * sinft[l];

subject to reactive_flow_to_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Qta[l] = -BTT[l] * V2[T_BUS[l]] - BTF[l] * cosft[l] - GTF[l] * sinft[l];

subject to Pfa_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    PFLOAC[l] * (1 - status[l]) <= -Pf[l] + Pfa[l];

subject to Pfa_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    -Pf[l] + Pfa[l] <= PFUPAC[l] * (1 - status[l]);

subject to Pta_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    PTLOAC[l] * (1 - status[l]) <= -Pt[l] + Pta[l];

subject to Pta_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    -Pt[l] + Pta[l] <= PTUPAC[l] * (1 - status[l]);

subject to Qfa_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    QFLOAC[l] * (1 - status[l]) <= -Qf[l] + Qfa[l];

subject to Qfa_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    -Qf[l] + Qfa[l] <= QFUPAC[l] * (1 - status[l]);

subject to Qta_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    QTLOAC[l] * (1 - status[l]) <= -Qt[l] + Qta[l];

subject to Qta_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    -Qt[l] + Qta[l] <= QTUPAC[l] * (1 - status[l]);

########## POWER FLOW LIMITS ##########

subject to flowf_limits_1 {l in L:OPF_TYPE='dc' and (BR_SWITCH[l] == 1 or BR_SWITCH[l] == 2)}:
    PFMIN[l] <= Pf[l] <= PFMAX[l];

subject to flowt_limits_1 {l in L:OPF_TYPE='dc' and (BR_SWITCH[l] == 1 or BR_SWITCH[l] == 2)}:
    -PFMAX[l] <= Pt[l] <= -PFMIN[l];

subject to flowf_limits_dc_lower {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pf[l] >= PFMIN[l] * status[l];

subject to flowf_limits_dc_upper {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pf[l] <= PFMAX[l] * status[l];

subject to flowt_limits_dc_lower {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pt[l] >= -PFMAX[l] * status[l];

subject to flowt_limits_dc_upper {l in L: OPF_TYPE == 'dc' and BR_SWITCH[l] == 3}:
    Pt[l] <= -PFMIN[l] * status[l];

subject to flow_limits_from_23_4 {l in L:(OPF_TYPE='acpolar' or OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr') and (BR_SWITCH[l] == 1 or BR_SWITCH[l] == 2)}:
    Pf[l]^2 + Qf[l]^2 <= RATE_A[l]^2;

subject to flow_limits_to_23_4 {l in L:(OPF_TYPE='acpolar' or OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr') and (BR_SWITCH[l] == 1 or BR_SWITCH[l] == 2)}:
    Pt[l]^2 + Qt[l]^2 <= RATE_A[l]^2;

subject to flow_limits_from_acrect {l in L: (OPF_TYPE == 'acpolar' or OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Pf[l]^2 + Qf[l]^2 <= RATE_A[l]^2 * status[l];

subject to flow_limits_to_acrect {l in L: (OPF_TYPE == 'acpolar' or OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_SWITCH[l] == 3}:
    Pt[l]^2 + Qt[l]^2 <= RATE_A[l]^2 * status[l];

########## RECTANGULAR DEFINITIONS ##########

subject to eq_vol_squared {n in N:OPF_TYPE = 'acrect'}:
    V2[n] == Vr[n]*Vr[n] + Vi[n]*Vi[n];

subject to eq_cosft {l in L:OPF_TYPE = 'acrect'}:
    cosft[l] == Vr[F_BUS[l]]*Vr[T_BUS[l]] + Vi[F_BUS[l]]*Vi[T_BUS[l]];

subject to eq_sinft {l in L:OPF_TYPE = 'acrect'}:
    sinft[l] == Vi[F_BUS[l]]*Vr[T_BUS[l]] - Vr[F_BUS[l]]*Vi[T_BUS[l]];

########## JABR RELAXATION ##########

subject to jabr_relaxation_ft {l in L:OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    cosft[l]^2 + sinft[l]^2 <= V2[F_BUS[l]]*V2[T_BUS[l]];

# TODO: Add the relaxation by Muñoz where losses are positive?

########## SLACK BUS ##########

subject to eq_slack:
    Va[0] == 0;

subject to eq_slack_imag:
    Vi[0] == 0;

########## CONNECTIVITY CONSTRAINTS ##########

subject to status_split {l in L: CONNECTIVITY = 'on' and (BR_SWITCH[l] == 2 or BR_SWITCH[l] == 3)}:
    statusf[l] + statust[l] == status[l];

subject to connectivity {n in N: CONNECTIVITY = 'on'}:
    sum {l in L: F_BUS[l] == n} statusf[l] + sum {l in L: T_BUS[l] == n} statust[l] >= 1;

subject to mtz_connectivity_f {l in L: CONNECTIVITY = 'on' and F_BUS[l] != 0 and T_BUS[l] != 0}:
    u[F_BUS[l]] - u[T_BUS[l]] + card(N) * statusf[l] <= card(N) - 1;

subject to mtz_connectivity_t {l in L: CONNECTIVITY = 'on' and F_BUS[l] != 0 and T_BUS[l] != 0}:
    u[T_BUS[l]] - u[F_BUS[l]] + card(N) * statust[l] <= card(N) - 1;

subject to eq_slack_mtz:
    u[0] == 1;

########## STATUS FIX ##########

subject to fix_status_0 {l in L: BR_SWITCH[l] == 0}:
    status[l] == 0;

subject to fix_status_1 {l in L: BR_SWITCH[l] == 1}:
    status[l] == 1;

# TODO: Split AMPL models?
