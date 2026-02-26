# Windkessel Cardiac Modeling Project 

The Windkessel Effect can be used to model the hemodynamics of the arterial system. The goal of this project is to use the two and three element Windkessel models to model the blood pressure and flow of one cardiac cycle and compare the two models. 

## Requirements 

This project used MATLAB and Signal Processing Toolbox Version 25.2 (R2025b)

## Mathematical Derivations 
### Two-Element Windkessel 

Using a circuit analogy and Kirchoff’s Law ($\sum I_{in} = \sum I_{out}$), the following equation states that cardiac output, $Q(t)$, is equal to flow through the resistive ($Q_R$ small arteries) and compliant ($Q_c$ large arteries) compartments:

$$
Q(t) = Q_R(t) + Q_C(t)
$$

Substituting Ohm’s Law for fluid flow ($Q(t) = \frac{P(t)}{R}$ where $R$ is fluid resistance) and the definition of vascular compliance ($Q(t) = C\frac{dP}{dt}$ where $C$ is capacitance):

$$
Q(t) = \frac{P(t)}{R} ​+ C\frac{dP}{dt}
$$

This equation is solved by estimating $\frac{dP}{dt}$ from the measured pressure and optimizing the $R$ and $C$ parameters to model flow. 

Blood pressure is described by rearranging the equation:

$$
P(t) = RQ(t) – CR\frac{dP}{dt}
$$

To solve with an ordinary differential equation solver, the equation can be written as:

$$
\frac{dP}{dt} = \frac{Q(t)}{C} – \frac{P(t)}{RC} 
$$

The equation is solved at $t$ and $P(t)$ to model blood pressure
### Three-Element Windkessel 

The Windkessel Model can also include an additional resistive element ($R$) representing the characteristic impedance of the proximal aorta. Using the circuit analogy (voltage in series where voltage is analogous to pressure) and Ohm’s Law for voltage drop ($V = IR$ where voltage is analogous to pressure, current is analogous to flow, and $R$ is resistance), pressure and flow can be modeled as follows:

$$
P(t) = RQ(t) + P_d(t)
$$

$$
Q(t) = \frac{P_d(t)}{R_d} ​+ C\frac{dP_d}{​dt}
$$

Because $P_d(t) = P(t) - Z_cQ(t)$:

$$
Q(t) = \frac{P(t)}{R_d} - \frac{RQ(t)}{R_d} ​+ C\frac{d(P(t) - RQ(t))}{dt} = \frac{P(t)}{R_d} - \frac{RQ(t)}{R_d} ​+ C\frac{dP}{dt} - CR\frac{dQ}{dt}
$$

Rearranging to solve pressure with ode45:

$$
\frac{dP}{dt} = \frac{RQ(t)}{R_dC} – \frac{P(t)}{R_dC} + R\frac{dQ}{dt} – \frac{Q(t)}{C}
$$

And rearranging to solve flow with ode45:

$$
\frac{dQ}{dt} = \frac{P(t)}{R_dCR} – \frac{Q(t)}{R_dC} - \frac{Q(t)}{CR} + \frac{1}{R}\frac{dP}{dt}
$$

