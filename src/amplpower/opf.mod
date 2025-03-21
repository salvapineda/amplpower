set N;
set G;
set L;
param OPF_TYPE symbolic;
param CONNECTIVITY symbolic;
param BASEMVA;
param MAXVOL;
param MINVOL;
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
param TAP {L};
param SHIFT {L};
param BR_STATUS {L};
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

var genp {g in G} >= PMIN[g], <= PMAX[g]:= PG0[g];
var genq {g in G} >= QMIN[g], <= QMAX[g]:= QG0[g];
var flowpf {l in L} >= -RATE_A[l], <=RATE_A[l] := PF0[l];
var flowpt {l in L} >= -RATE_A[l], <=RATE_A[l] := PT0[l];
var flowqf {l in L} >= -RATE_A[l], <=RATE_A[l] := QF0[l];
var flowqt {l in L} >= -RATE_A[l], <=RATE_A[l] := QT0[l];
var vol {n in N} >= VMIN[n], <= VMAX[n] := VOL0[n];
var volr {n in N} >= 0, <=VMAX[n] := VOLR0[n];
var voli {n in N} >= -VMAX[n], <=VMAX[n] := VOLI0[n];
var vol2 {n in N} >= VMIN[n]^2, <= VMAX[n]^2 := VOL0[n]^2;
var ang {n in N} >= AMIN[n], <= AMAX[n] := ANG0[n];
var cosft {l in L} >= 0, <= VMAX[F_BUS[l]]*VMAX[T_BUS[l]] := VOLR0[F_BUS[l]]*VOLR0[T_BUS[l]] + VOLI0[F_BUS[l]]*VOLI0[T_BUS[l]];
var sinft {l in L} >= -VMAX[F_BUS[l]]*VMAX[T_BUS[l]], <= VMAX[F_BUS[l]]*VMAX[T_BUS[l]] := VOLI0[F_BUS[l]]*VOLR0[T_BUS[l]] - VOLR0[F_BUS[l]]*VOLI0[T_BUS[l]];
# TODO: I think these bounds can be improved considering that sinft is a small number
var flowpf_aux {l in L} >= min(PFLODC[l],PFLOAC[l]), <= max(PFUPDC[l],PFUPAC[l]) := PF0[l];
var flowpt_aux {l in L} >= PTLOAC[l], <= PTUPAC[l] := PT0[l];
var flowqf_aux {l in L} >= QFLOAC[l], <= QFUPAC[l] := QF0[l];
var flowqt_aux {l in L} >= QTLOAC[l], <= QTUPAC[l] := QT0[l];
var u {n in N} >= 1, <= card(N);
var status {l in L} binary;
var statusf {l in L} binary;
var statust {l in L} binary;

########## OBJECTIVE FUNCTION ##########

minimize total_cost:
    sum {g in G} (COST_2[g] * (BASEMVA*genp[g])^2 + COST_1[g] * (BASEMVA*genp[g]) + COST_0[g]);

########## POWER BALANCE ##########

subject to active_power_balance_1 {n in N:OPF_TYPE='dc'}:
    sum {g in G} CG[g,n] * genp[g] - PD[n] = sum {l in L} (CF[l,n] * flowpf[l] + CT[l,n] * flowpt[l]);

subject to active_power_balance_2 {n in N:OPF_TYPE='acpolar'}:
    sum {g in G} CG[g,n] * genp[g] - PD[n] = GS[n]*vol[n]*vol[n] + sum {l in L} (CF[l,n] * flowpf[l] + CT[l,n] * flowpt[l]);

subject to reactive_power_balance_2 {n in N:OPF_TYPE='acpolar'}:
    sum {g in G} CG[g,n] * genq[g] - QD[n] = -BS[n]*vol[n]*vol[n] + sum {l in L} (CF[l,n] * flowqf[l] + CT[l,n] * flowqt[l]);

subject to active_power_balance_3_4 {n in N: OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    sum {g in G} CG[g,n] * genp[g] - PD[n] = GS[n]*vol2[n] + sum {l in L} (CF[l,n] * flowpf[l] + CT[l,n] * flowpt[l]);

subject to reactive_power_balance_3_4 {n in N: OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    sum {g in G} CG[g,n] * genq[g] - QD[n] = -BS[n]*vol2[n] + sum {l in L} (CF[l,n] * flowqf[l] + CT[l,n] * flowqt[l]);

########## POWER FLOW DEFINITIONS (BR_STATUS = 0) ##########

subject to active_flow_from_0 {l in L:BR_STATUS[l] == 0}:
    flowpf[l] = 0;

subject to active_flow_to_0 {l in L:BR_STATUS[l] == 0}:
    flowpt[l] = 0;

subject to reactive_flow_from_0 {l in L:BR_STATUS[l] == 0}:
    flowqf[l] = 0;

subject to reactive_flow_to_0 {l in L:BR_STATUS[l] == 0}:
    flowqt[l] = 0;

########## POWER FLOW DEFINITIONS (BR_STATUS = 1) ##########

subject to active_flow_from_1 {l in L:OPF_TYPE='dc' and BR_STATUS[l] == 1}:
    flowpf[l] = (1 / BR_X[l]) * (ang[F_BUS[l]] - ang[T_BUS[l]]);

subject to active_flow_to_1 {l in L:OPF_TYPE='dc' and BR_STATUS[l] == 1}:
    flowpt[l] = (1 / BR_X[l]) * (ang[T_BUS[l]] - ang[F_BUS[l]]);

subject to active_flow_from_2 {l in L:OPF_TYPE='acpolar' and BR_STATUS[l] == 1}:
    flowpf[l] = GFF[l]*vol[F_BUS[l]]*vol[F_BUS[l]] + vol[F_BUS[l]]*vol[T_BUS[l]]*(GFT[l]*cos(ang[F_BUS[l]]-ang[T_BUS[l]])+BFT[l]*sin(ang[F_BUS[l]]-ang[T_BUS[l]]));

subject to active_flow_to_2 {l in L:OPF_TYPE='acpolar' and BR_STATUS[l] == 1}:
    flowpt[l] = GTT[l]*vol[T_BUS[l]]*vol[T_BUS[l]] + vol[F_BUS[l]]*vol[T_BUS[l]]*(GTF[l]*cos(ang[T_BUS[l]]-ang[F_BUS[l]])+BTF[l]*sin(ang[T_BUS[l]]-ang[F_BUS[l]]));

subject to reactive_flow_from_2 {l in L:OPF_TYPE='acpolar' and BR_STATUS[l] == 1}:
    flowqf[l] = -BFF[l]*vol[F_BUS[l]]*vol[F_BUS[l]] - vol[F_BUS[l]]*vol[T_BUS[l]]*(BFT[l]*cos(ang[F_BUS[l]]-ang[T_BUS[l]])-GFT[l]*sin(ang[F_BUS[l]]-ang[T_BUS[l]]));

subject to reactive_flow_to_2 {l in L:OPF_TYPE='acpolar' and BR_STATUS[l] == 1}:
    flowqt[l] = -BTT[l]*vol[T_BUS[l]]*vol[T_BUS[l]] - vol[F_BUS[l]]*vol[T_BUS[l]]*(BTF[l]*cos(ang[T_BUS[l]]-ang[F_BUS[l]])-GTF[l]*sin(ang[T_BUS[l]]-ang[F_BUS[l]]));

subject to active_flow_from_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_STATUS[l] == 1}:
    flowpf[l] = GFF[l]*vol2[F_BUS[l]] + GFT[l]*cosft[l] + BFT[l]*sinft[l];

subject to active_flow_to_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_STATUS[l] == 1}:
    flowpt[l] = GTT[l]*vol2[T_BUS[l]] + GTF[l]*cosft[l] - BTF[l]*sinft[l];

subject to reactive_flow_from_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_STATUS[l] == 1}:
    flowqf[l] = -BFF[l]*vol2[F_BUS[l]] - BFT[l]*cosft[l] + GFT[l]*sinft[l];

subject to reactive_flow_to_3_4 {l in L:(OPF_TYPE='acrect' or OPF_TYPE = 'acjabr') and BR_STATUS[l] == 1}:
    flowqt[l] = -BTT[l]*vol2[T_BUS[l]] - BTF[l]*cosft[l] - GTF[l]*sinft[l];

########## POWER FLOW DEFINITIONS (BR_STATUS = 2) ##########

subject to active_flow_from_1_switch {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 2}:
    flowpf[l] = status[l] * (1 / BR_X[l]) * (ang[F_BUS[l]] - ang[T_BUS[l]]);

subject to active_flow_to_1_switch {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 2}:
    flowpt[l] = status[l] * (1 / BR_X[l]) * (ang[T_BUS[l]] - ang[F_BUS[l]]);

subject to active_flow_from_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_STATUS[l] == 2}:
    flowpf[l] = status[l] * (GFF[l]*vol[F_BUS[l]]*vol[F_BUS[l]] + vol[F_BUS[l]]*vol[T_BUS[l]]*(GFT[l]*cos(ang[F_BUS[l]]-ang[T_BUS[l]])+BFT[l]*sin(ang[F_BUS[l]]-ang[T_BUS[l]])));

subject to active_flow_to_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_STATUS[l] == 2}:
    flowpt[l] = status[l] * (GTT[l]*vol[T_BUS[l]]*vol[T_BUS[l]] + vol[F_BUS[l]]*vol[T_BUS[l]]*(GTF[l]*cos(ang[T_BUS[l]]-ang[F_BUS[l]])+BTF[l]*sin(ang[T_BUS[l]]-ang[F_BUS[l]])));

subject to reactive_flow_from_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_STATUS[l] == 2}:
    flowqf[l] = status[l] * (-BFF[l]*vol[F_BUS[l]]*vol[F_BUS[l]] - vol[F_BUS[l]]*vol[T_BUS[l]]*(BFT[l]*cos(ang[F_BUS[l]]-ang[T_BUS[l]])-GFT[l]*sin(ang[F_BUS[l]]-ang[T_BUS[l]])));

subject to reactive_flow_to_2_switch {l in L: OPF_TYPE == 'acpolar' and BR_STATUS[l] == 2}:
    flowqt[l] = status[l] * (-BTT[l]*vol[T_BUS[l]]*vol[T_BUS[l]] - vol[F_BUS[l]]*vol[T_BUS[l]]*(BTF[l]*cos(ang[T_BUS[l]]-ang[F_BUS[l]])-GTF[l]*sin(ang[T_BUS[l]]-ang[F_BUS[l]])));

subject to active_flow_from_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 2}:
    flowpf[l] = status[l] * (GFF[l]*vol2[F_BUS[l]] + GFT[l]*cosft[l] + BFT[l]*sinft[l]);

subject to active_flow_to_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 2}:
    flowpt[l] = status[l] * (GTT[l]*vol2[T_BUS[l]] + GTF[l]*cosft[l] - BTF[l]*sinft[l]);

subject to reactive_flow_from_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 2}:
    flowqf[l] = status[l] * (-BFF[l]*vol2[F_BUS[l]] - BFT[l]*cosft[l] + GFT[l]*sinft[l]);

subject to reactive_flow_to_3_4_switch {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 2}:
    flowqt[l] = status[l] * (-BTT[l]*vol2[T_BUS[l]] - BTF[l]*cosft[l] - GTF[l]*sinft[l]);

########## POWER FLOW DEFINITIONS (BR_STATUS = 3) ##########

subject to active_flow_from_1_bigm {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpf_aux[l] = (1 / BR_X[l]) * (ang[F_BUS[l]] - ang[T_BUS[l]]);

subject to active_flow_to_1_bigm {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpt_aux[l] = (1 / BR_X[l]) * (ang[T_BUS[l]] - ang[F_BUS[l]]);

subject to flowpf_aux_lower_1 {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    PFLODC[l] * (1 - status[l]) <= -flowpf[l] + flowpf_aux[l];

subject to flowpf_aux_upper_1 {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    -flowpf[l] + flowpf_aux[l] <= PFUPDC[l] * (1 - status[l]);

subject to flowpt_aux_lower_1 {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    PFLODC[l] * (1 - status[l]) <= -flowpt[l] + flowpt_aux[l];

subject to flowpt_aux_upper_1 {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    -flowpt[l] + flowpt_aux[l] <= PFUPDC[l] * (1 - status[l]);

subject to active_flow_from_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowpf_aux[l] = GFF[l] * vol2[F_BUS[l]] + GFT[l] * cosft[l] + BFT[l] * sinft[l];

subject to active_flow_to_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowpt_aux[l] = GTT[l] * vol2[T_BUS[l]] + GTF[l] * cosft[l] - BTF[l] * sinft[l];

subject to reactive_flow_from_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowqf_aux[l] = -BFF[l] * vol2[F_BUS[l]] - BFT[l] * cosft[l] + GFT[l] * sinft[l];

subject to reactive_flow_to_2_bigm {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowqt_aux[l] = -BTT[l] * vol2[T_BUS[l]] - BTF[l] * cosft[l] - GTF[l] * sinft[l];

subject to flowpf_aux_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    PFLOAC[l] * (1 - status[l]) <= -flowpf[l] + flowpf_aux[l];

subject to flowpf_aux_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    -flowpf[l] + flowpf_aux[l] <= PFUPAC[l] * (1 - status[l]);

subject to flowpt_aux_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    PTLOAC[l] * (1 - status[l]) <= -flowpt[l] + flowpt_aux[l];

subject to flowpt_aux_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    -flowpt[l] + flowpt_aux[l] <= PTUPAC[l] * (1 - status[l]);

subject to flowqf_aux_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    QFLOAC[l] * (1 - status[l]) <= -flowqf[l] + flowqf_aux[l];

subject to flowqf_aux_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    -flowqf[l] + flowqf_aux[l] <= QFUPAC[l] * (1 - status[l]);

subject to flowqt_aux_lower_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    QTLOAC[l] * (1 - status[l]) <= -flowqt[l] + flowqt_aux[l];

subject to flowqt_aux_upper_2 {l in L: (OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    -flowqt[l] + flowqt_aux[l] <= QTUPAC[l] * (1 - status[l]);

########## POWER FLOW LIMITS ##########

subject to flowf_limits_1 {l in L:OPF_TYPE='dc' and (BR_STATUS[l] == 1 or BR_STATUS[l] == 2)}:
    -RATE_A[l] <= flowpf[l] <= RATE_A[l];

subject to flowt_limits_1 {l in L:OPF_TYPE='dc' and (BR_STATUS[l] == 1 or BR_STATUS[l] == 2)}:
    -RATE_A[l] <= flowpt[l] <= RATE_A[l];

subject to flowf_limits_dc_lower {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpf[l] >= -RATE_A[l] * status[l];

subject to flowf_limits_dc_upper {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpf[l] <= RATE_A[l] * status[l];

subject to flowt_limits_dc_lower {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpt[l] >= -RATE_A[l] * status[l];

subject to flowt_limits_dc_upper {l in L: OPF_TYPE == 'dc' and BR_STATUS[l] == 3}:
    flowpt[l] <= RATE_A[l] * status[l];

subject to flow_limits_from_23_4 {l in L:(OPF_TYPE='acpolar' or OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr') and (BR_STATUS[l] == 1 or BR_STATUS[l] == 2)}:
    flowpf[l]^2 + flowqf[l]^2 <= RATE_A[l]^2;

subject to flow_limits_to_23_4 {l in L:(OPF_TYPE='acpolar' or OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr') and (BR_STATUS[l] == 1 or BR_STATUS[l] == 2)}:
    flowpt[l]^2 + flowqt[l]^2 <= RATE_A[l]^2;

subject to flow_limits_from_acrect {l in L: (OPF_TYPE == 'acpolar' or OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowpf[l]^2 + flowqf[l]^2 <= RATE_A[l]^2 * status[l];

subject to flow_limits_to_acrect {l in L: (OPF_TYPE == 'acpolar' or OPF_TYPE == 'acrect' or OPF_TYPE == 'acjabr') and BR_STATUS[l] == 3}:
    flowpt[l]^2 + flowqt[l]^2 <= RATE_A[l]^2 * status[l];

########## RECTANGULAR DEFINITIONS ##########

subject to eq_vol_squared {n in N:OPF_TYPE = 'acrect'}:
    vol2[n] == volr[n]*volr[n] + voli[n]*voli[n];

subject to eq_cosft {l in L:OPF_TYPE = 'acrect'}:
    cosft[l] == volr[F_BUS[l]]*volr[T_BUS[l]] + voli[F_BUS[l]]*voli[T_BUS[l]];

subject to eq_sinft {l in L:OPF_TYPE = 'acrect'}:
    sinft[l] == voli[F_BUS[l]]*volr[T_BUS[l]] - volr[F_BUS[l]]*voli[T_BUS[l]];

########## JABR RELAXATION ##########

subject to jabr_relaxation_ft {l in L:OPF_TYPE = 'acrect' or OPF_TYPE = 'acjabr'}:
    cosft[l]^2 + sinft[l]^2 <= vol2[F_BUS[l]]*vol2[T_BUS[l]];

# TODO: Add the relaxation by Muñoz where losses are positive?

########## SLACK BUS ##########

subject to eq_slack:
    ang[0] == 0;

subject to eq_slack_imag:
    voli[0] == 0;

########## CONNECTIVITY CONSTRAINTS ##########

subject to status_split {l in L: CONNECTIVITY = 'on' and (BR_STATUS[l] == 2 or BR_STATUS[l] == 3)}:
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

subject to fix_status_0 {l in L: BR_STATUS[l] == 0}:
    status[l] == 0;

subject to fix_status_1 {l in L: BR_STATUS[l] == 1}:
    status[l] == 1;

# TODO: Split AMPL models?
