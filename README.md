# noise-suppression2025-scripts

This repository contains Matlab scripts for generating figures 3, 4, and 5 in the manuscript

Harry Dankowicz, Steven W Shaw, Oriel Shoshani (2025) *Improving frequency stability using slowly modulated adaptive feedback*, in review.

## Content
This repository contains:

- `Gen_Figures.m`: script to perform five simulations (numbered from 1 to 5) of three different cases:
	1) Standard Duffing Oscillator (stops after 100 time units).
	2) Oscillator with zero amplitude-to-frequency noise conversion (stops after 100 time units).
	3) Oscillator with theoretically zero phase diffusion (stops after 100 time units).
	4) Oscillator with theoretically zero phase diffusion (stops after 50 time units).
	3) Oscillator with theoretically zero phase diffusion (stops after 200 time units).
	
	Figure 5 of the manuscript can also be generated with this code by defining D_phi = sigma^2 and using different values of That and k.

## Usage

1. execute individual sections in the script `Gen_Figures.m`
