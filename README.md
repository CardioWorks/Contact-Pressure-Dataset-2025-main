Welcome to the "Contact Pressure-PPG Dataset Research". This package offers the perfusion index and feature calculated from PPG for technical validation. It operates on MathWorks MATLAB R2020b.

"PI_ppg.m": main program, it calls the code in "functions" folder. First divide PPG into each cycle and calculate perfusion index and b/a feature on each cycle. Then, averaged at a trimmed ratio of 30%. Finally, calculated the trimmed averaged perfusion at the pressure of 30, 40, 50, 60, 70 and 80 mmHg.

"PI_data": contains some statistical data of perfusion index calculation.

Authors: Ziyi Wang, Yongbo Liang and Mohamed Elgendi
First uploade date: June 25, 2025