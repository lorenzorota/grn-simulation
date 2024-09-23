# Gene Regulatory Network Simulation

This is some code I wrote a long time ago for the course Modelling and Simulation at the University of Groningen, which may be of use to someone. We compare a boolean and kinetic model for simulating a GRN to see if they are equivalent in simulating causal relations among genes. The boolean model should be correct, but the kinetic model has a mistake somewhere, most likely in the ODE solver.

## Instructions

1. Create a virtualenv

2. Make sure to install all the dependencies
   - `pip install -r requirements.txt`

3. Test out either the kinetic or boolean model
   - `python main.py --mode [kinetic|boolean]`

4. Run the experiment
   - `python main.py --experiment`

## Contributors

- Lorenzo Rota
- Andrei Stoica