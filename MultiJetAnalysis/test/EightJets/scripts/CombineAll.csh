#! /bin/csh

set HISTO="MLP"

combine -M Asymptotic datacard_0400_133_040_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 400
combine -M Asymptotic datacard_0500_167_050_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 500
combine -M Asymptotic datacard_0600_200_060_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 600
combine -M Asymptotic datacard_0700_233_070_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 700
combine -M Asymptotic datacard_0800_267_080_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 800
combine -M Asymptotic datacard_0900_300_090_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 900
combine -M Asymptotic datacard_1000_333_100_$HISTO.txt --rMax 1   -n Coloron$HISTO -m 1000
combine -M Asymptotic datacard_1100_367_110_$HISTO.txt --rMax 2   -n Coloron$HISTO -m 1100
combine -M Asymptotic datacard_1200_400_120_$HISTO.txt --rMax 5   -n Coloron$HISTO -m 1200
combine -M Asymptotic datacard_1300_433_130_$HISTO.txt --rMax 10  -n Coloron$HISTO -m 1300
combine -M Asymptotic datacard_1400_467_140_$HISTO.txt --rMax 50  -n Coloron$HISTO -m 1400
combine -M Asymptotic datacard_1500_500_150_$HISTO.txt --rMax 100 -n Coloron$HISTO -m 1500
