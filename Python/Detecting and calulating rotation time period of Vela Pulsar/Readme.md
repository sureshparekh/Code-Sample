# Analyzing the Time Period of Vela Pulsar

Welcome to the project on analyzing the time period of the Vela Pulsar. In this project, we have implemented our basic understanding of Pulsar Astronomy to calculate the Time Period of Vela Pulsar. Our choice of pulsar rests on the fact that it is the brightest object in the high energy gamma-ray sky. The simplistic data set consisting of only voltage signals makes our preliminary attempt as closely accurate as possible. The observations had been made at 326.5 MHz through a cylindrically paraboloid telescope at Ooty. A higher frequency creates a much lower delay in the arrival time of pulses and makes our calculations even more accurate.

## Introduction:

Pulsar timing is a method used in radio astronomy to study the properties of pulsars. Pulsars are highly magnetized, rotating neutron stars that emit beams of electromagnetic radiation out of their magnetic poles. As the pulsar rotates, these beams sweep across the sky, and when they intersect the Earth's line of sight, they appear as regular pulses of radiation. Vela Pulsar, also known as PSR B0833-45, is one of the most studied pulsars due to its stability and brightness in high-energy gamma-ray observations.

## About Ooty Radio Telescope:

The observations for this project were conducted using a cylindrically paraboloid telescope located at Ooty Radio Astronomy Observatory, India. The Ooty Radio Telescope is an important facility for radio astronomy research in India, operating at frequencies ranging from 10 MHz to 1.5 GHz.

## Process of Data Analysis:

### 1. Statistical Characteristics of the Signal:

We began our analysis by performing statistical evaluation of the raw voltage signals to verify some of its expected properties. We expect the signals to have a Gaussian distribution. To do this, we randomly and uniformly selected 100,000 voltage samples from both the north and the south arrays and plotted the histogram. As expected, the voltage signals demonstrate Gaussian distribution.

#### 1.1 Voltage Signal Characteristics:

Further analysis of voltage signal characteristics confirmed the expected Gaussian distribution.

#### 1.2 Power Signal Characteristics:

We studied the distribution of power signals, which are the square of the voltage signals, and confirmed an exponential distribution.

### 2. Properties of Signal in Time-Frequency Domain:

#### 2.1 Voltage Power Spectrum:

The Power Spectrum of a signal describes the power present in the signal as a function of frequency. We used the Fast Fourier Transform (FFT) to plot the power distribution as a function of frequency. The voltage power spectrum for both the northern and southern arrays were analyzed.

#### 2.2 Dynamic Spectrum:

The Dynamic Spectrum is a color-coded graph showing the relationship between Frequency (MHz) and Time (ms). It enables detection of pulsar signal indicators and dispersion effects in the interstellar medium.

### 3. Dispersion Measure and Frequency Delay:

Dispersion Measure (DM) is a parameter that shows up in observations as the broadening of an otherwise sharp pulse. It is calculated based on the pulse arrival time and frequency. The delay in the propagation of light due to electrostatic interaction between radio waves and charged particles in the Interstellar Medium is characterized by the Dispersion Measure.

## Additional Analysis:

### 4. Distance to Pulsar:

The distance to the pulsar (S) is given by the equation:

S = DM / ne

Where ne is the mean electron density between the pulsar and earth and is equal to 0.23 per cc. Hence, S = 294 pc.

### 5. Dedispersed Time series:

We eliminated the frequency-dependent time delays using the DM. The dedispersed time series showed significant peaks, indicating an increase in Signal to Noise Ratio (SNR).

### 6. Time Period of the pulsar:

Curve fitting was used to estimate the arrival times of individual pulses, and the time period of the pulsar was calculated as 89.3 ms. This value is in agreement with the currently accepted value.

## Conclusion:

Analyzing the time period of Vela Pulsar provides valuable insights into the properties of pulsars and the interstellar medium. By studying the statistical characteristics, power spectrum, and dispersion effects, we can deepen our understanding of pulsar astronomy and radio astronomy techniques.

This project showcases the application of basic pulsar astronomy principles and data analysis techniques to study the time period of Vela Pulsar, contributing to the broader field of astrophysics and radio astronomy.

An important thing to note is that all the quantitative figures are estimates with some amount of uncertainty. This can be due to uncertainties in other parameters such as the dispersion measure, low exposure time, small dataset etc. Longer observation time would significantly reduce uncertainties in the data. Since the source is a compact object, the visibility should remain constant as a function of time. Hence, the Fourier Transform of brightness distribution should not change with a change in baseline. Over the course of the project, we gained invaluable insights into the working of pulsars and how to decipher information from mere observations.
