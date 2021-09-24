# SAW Oscillator Phase Noise Calculator
Tool for calculating the phase/frequency noise in a closed-loop surface acoustic wave (SAW) oscillator where the SAW device can be a (two-port) delay line or a resonator

# Screenshots
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/01_start.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/02_saw_delay_line_block_diagram.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/03_saw_delay_line_components_phase_noise.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/04_saw_delay_line_oscillator_phase_noise.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/05_saw_delay_line_oscillator_frequency_noise.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/06_saw_resonator_block_diagram.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/07_saw_resonator_components_phase_noise.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/08_saw_resonator_oscillator_phase_noise.png)
![](https://github.com/yildi1337/SAWOscillatorPhaseNoiseCalculator/blob/main/screenshots/09_saw_resonator_oscillator_frequency_noise.png)

# Further information
The tool is mainly based on the [PySimpleGUI](https://pypi.org/project/PySimpleGUI/) and the [matplotlib](https://matplotlib.org/) libraries. For the underlying math, please refer to the following publications:

* P. Durdaut et al., Equivalence of Open-Loop and Closed-Loop Operation of SAW Resonators and Delay Lines, Sensors 19, 1, 185, 2019. https://doi.org/10.3390/s19010185
* P. Durdaut et al., Noise Analysis and Comparison of Phase- and Frequency-Detecting Readout Systems: Application to SAW Delay Line Magnetic Field Sensor, IEEE Sensors Journal, 19, 18, 8000-8008, 2019. https://doi.org/10.1109/JSEN.2019.2914965
* E. Rubiola, Phase noise and frequency stability in oscillators, Cambridge University Press, 2009. https://doi.org/10.1017/CBO9780511812798