set PLANTS;
set HOURS;
param Capacity{PLANTS};
param Cost{PLANTS};
param CO2_Emission{PLANTS};
param Demand{HOURS};
param RampRate{PLANTS};
var Gen{PLANTS, HOURS} >= 0;
# Objective: Minimize cost while considering emissions
minimize TotalCost:
    sum {p in PLANTS, t in HOURS} Cost[p] * Gen[p,t];
# Meet demand with 10% spinning reserve
subject to MeetDemand {t in HOURS}:
    sum {p in PLANTS} Gen[p,t] >= 1.1 * Demand[t];
# Enforce capacity limits
subject to CapacityLimit {p in PLANTS, t in HOURS}:
    Gen[p,t] <= Capacity[p];
# Ramp rate constraint (except for first hour)
subject to RampConstraint {p in PLANTS, t in HOURS: t > 1}:
    abs(Gen[p,t] - Gen[p,t-1]) <= RampRate[p];
# Solar and Wind availability constraints
subject to SolarConstraint {t in HOURS: t < 6 or t > 18}:
    Gen['Solar_1', t] = 0;
subject to WindConstraint {t in HOURS}:
    Gen['Wind_1', t] <= 300;  # Assuming wind variation
