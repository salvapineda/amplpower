# File test to check Ampl is working
import pytest
from amplpy import AMPL
import pandas as pd

def test_simple():
    ampl = AMPL()
    assert ampl.get_value('3.14') == pytest.approx(3.14)
    
# Data for example
def example_data():
    # Power plant characteristics and demand profile
    power_plants = df = pd.DataFrame(
        {
            "Capacity_MW": [500, 600, 400, 450, 1000, 200, 300],
            "Cost_per_MW": [50, 55, 40, 42, 30, 0, 0],  # Solar & Wind have no fuel cost
            "CO2_Emission_per_MW": [1.2, 1.1, 0.5, 0.4, 0, 0, 0],  # Carbon emissions per MW
            "Ramp_Rate": [100, 120, 80, 90, 50, 300, 400],  # Max MW change per hour
        },
        index=["Coal_1", "Coal_2", "Gas_1", "Gas_2", "Nuclear_1", "Solar_1", "Wind_1"],
    )
    demand = pd.DataFrame(
        {
            "Demand_MW": [
                800,
                750,
                700,
                680,
                660,
                700,
                850,
                1000,
                1200,
                1300,
                1400,
                1350,
                1250,
                1150,
                1100,
                1080,
                1050,
                1100,
                1150,
                1200,
                1300,
                1200,
                1000,
                900,
            ]
        },
        index=range(24),
    )

    return power_plants, demand

def example_solution():
    return {('Coal_1', 0): 0, ('Coal_1', 1): 0, ('Coal_1', 2): 0, ('Coal_1', 3): 0, ('Coal_1', 4): 0, ('Coal_1', 5): 0, ('Coal_1', 6): 0, ('Coal_1', 7): 0, ('Coal_1', 8): 0, ('Coal_1', 9): 0, ('Coal_1', 10): 0, ('Coal_1', 11): 0, ('Coal_1', 12): 0, ('Coal_1', 13): 0, ('Coal_1', 14): 0, ('Coal_1', 15): 0, ('Coal_1', 16): 0, ('Coal_1', 17): 0, ('Coal_1', 18): 0, ('Coal_1', 19): 35, ('Coal_1', 20): 0, ('Coal_1', 21): 0, ('Coal_1', 22): 0, ('Coal_1', 23): 0, \
    ('Coal_2', 0): 0, ('Coal_2', 1): 0, ('Coal_2', 2): 0, ('Coal_2', 3): 0, ('Coal_2', 4): 0, ('Coal_2', 5): 0, ('Coal_2', 6): 0, ('Coal_2', 7): 0, ('Coal_2', 8): 0, ('Coal_2', 9): 0, ('Coal_2', 10): 0, ('Coal_2', 11): 0, ('Coal_2', 12): 0, ('Coal_2', 13): 0, ('Coal_2', 14): 0, ('Coal_2', 15): 0, ('Coal_2', 16): 0, ('Coal_2', 17): 0, ('Coal_2', 18): 0, ('Coal_2', 19): 0, ('Coal_2', 20): 0, ('Coal_2', 21): 0, ('Coal_2', 22): 0, ('Coal_2', 23): 0, \
    ('Gas_1', 0): 0, ('Gas_1', 1): 5, ('Gas_1', 2): 0, ('Gas_1', 3): 0, ('Gas_1', 4): 0, ('Gas_1', 5): 0, ('Gas_1', 6): 0, ('Gas_1', 7): 30, ('Gas_1', 8): 110, ('Gas_1', 9): 190, ('Gas_1', 10): 225, ('Gas_1', 11): 145, ('Gas_1', 12): 65, ('Gas_1', 13): 5, ('Gas_1', 14): 0, ('Gas_1', 15): 0, ('Gas_1', 16): 0, ('Gas_1', 17): 0, ('Gas_1', 18): 5, ('Gas_1', 19): 85, ('Gas_1', 20): 165, ('Gas_1', 21): 120, ('Gas_1', 22): 40, ('Gas_1', 23): 0, \
    ('Gas_2', 0): 0, ('Gas_2', 1): 0, ('Gas_2', 2): 0, ('Gas_2', 3): 0, ('Gas_2', 4): 0, ('Gas_2', 5): 0, ('Gas_2', 6): 0, ('Gas_2', 7): 0, ('Gas_2', 8): 90, ('Gas_2', 9): 70, ('Gas_2', 10): 95, ('Gas_2', 11): 70, ('Gas_2', 12): 0, ('Gas_2', 13): 0, ('Gas_2', 14): 0, ('Gas_2', 15): 0, ('Gas_2', 16): 0, ('Gas_2', 17): 0, ('Gas_2', 18): 0, ('Gas_2', 19): 90, ('Gas_2', 20): 105, ('Gas_2', 21): 90, ('Gas_2', 22): 0, ('Gas_2', 23): 0, \
    ('Nuclear_1', 0): 580, ('Nuclear_1', 1): 520, ('Nuclear_1', 2): 470, ('Nuclear_1', 3): 448, ('Nuclear_1', 4): 426, ('Nuclear_1', 5): 470, ('Nuclear_1', 6): 520, ('Nuclear_1', 7): 570, ('Nuclear_1', 8): 620, ('Nuclear_1', 9): 670, ('Nuclear_1', 10): 720, ('Nuclear_1', 11): 770, ('Nuclear_1', 12): 810, ('Nuclear_1', 13): 760, ('Nuclear_1', 14): 710, ('Nuclear_1', 15): 688, ('Nuclear_1', 16): 660, ('Nuclear_1', 17): 710, ('Nuclear_1', 18): 760, ('Nuclear_1', 19): 810, ('Nuclear_1', 20): 860, ('Nuclear_1', 21): 810, ('Nuclear_1', 22): 760, ('Nuclear_1', 23): 710, \
     ('Solar_1', 0): 0, ('Solar_1', 1): 0, ('Solar_1', 2): 0, ('Solar_1', 3): 0, ('Solar_1', 4): 0, ('Solar_1', 5): 0, ('Solar_1', 6): 200, ('Solar_1', 7): 200, ('Solar_1', 8): 200, ('Solar_1', 9): 200, ('Solar_1', 10): 200, ('Solar_1', 11): 200, ('Solar_1', 12): 200, ('Solar_1', 13): 200, ('Solar_1', 14): 200, ('Solar_1', 15): 200, ('Solar_1', 16): 200, ('Solar_1', 17): 200, ('Solar_1', 18): 200, ('Solar_1', 19): 0, ('Solar_1', 20): 0, ('Solar_1', 21): 0, ('Solar_1', 22): 0, ('Solar_1', 23): 0, \
     ('Wind_1', 0): 300, ('Wind_1', 1): 300, ('Wind_1', 2): 300, ('Wind_1', 3): 300, ('Wind_1', 4): 300, ('Wind_1', 5): 300, ('Wind_1', 6): 300, ('Wind_1', 7): 300, ('Wind_1', 8): 300, ('Wind_1', 9): 300, ('Wind_1', 10): 300, ('Wind_1', 11): 300, ('Wind_1', 12): 300, ('Wind_1', 13): 300, ('Wind_1', 14): 300, ('Wind_1', 15): 300, ('Wind_1', 16): 300, ('Wind_1', 17): 300, ('Wind_1', 18): 300, ('Wind_1', 19): 300, ('Wind_1', 20): 300, ('Wind_1', 21): 300, ('Wind_1', 22): 300, ('Wind_1', 23): 300}


def test_power_plants_example():
    solver = 'gurobi'
    ampl = AMPL()
    ampl.reset()
    # Read model
    ampl.read('./tests/test_data/example.mod')
    # Get data
    power_plants, demand = example_data()
    # Load pandas data into AMPL
    ampl.set["PLANTS"] = power_plants.index
    ampl.set["HOURS"] = demand.index
    ampl.param["Capacity"] = power_plants["Capacity_MW"]
    ampl.param["Cost"] = power_plants["Cost_per_MW"]
    ampl.param["CO2_Emission"] = power_plants["CO2_Emission_per_MW"]
    ampl.param["Demand"] = demand["Demand_MW"]
    ampl.param["RampRate"] = power_plants["Ramp_Rate"]
        
    ampl.solve(solver=solver, mp_options="outlev=0")
      
    assert ampl.solve_result == 'solved'
    assert ampl.get_current_objective().value() == pytest.approx(549930)
        # Next code may not work since there are multiple solutions
        # Check solution:
        #solution = ampl.var["Gen"].to_dict()
        #manual_solution = example_solution()
        
        #for i in power_plants.index:
        #    for j in demand.index:
        #        assert solution[i,j] == pytest.approx(manual_solution[i,j])
        
