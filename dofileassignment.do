#Are death and visual impairment significant?#

. tab died vimp, col chi

#Death and visual impairment are significant#

#Are age, sex, BMI, education, occupation... potential confounders?

tab died agegrp, col chi   #yes
tab died sex, col chi  #no
tab died ethnic, col chi #yes
tab died religion, col chi  #yes
tab died occupation, col chi #yes
tab died education, col chi  #yes
tab died mfpermg, col chi  #continuous (load per mg)
tab died systolic, col chi #continuous (mmHg)
tab died diastolic, col chi #continuous (mmHg)
tab died pulse, col chi #continuous (pulse rate/minute)
tab died weight, col chi #continuous (kg)
tab died height, col chi #continuous (cms)
tab died map, col chi #continuous (mean arterial pressure)
tab died bmigrp, col chi #yes
tab died mfgrp, col chi #no
tab died mfpos, col chi #yes
tab died enter, col chi #continuous (dates of entry)
tab died exit, col chi #continuous (dates of exit)
