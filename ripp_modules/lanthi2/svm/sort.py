# -*- coding: utf-8 -*-

import csv, operator
input_csv = "lanthi2_training_set.csv"
a = zip(*csv.reader(open(input_csv, "r")))
csv.writer(open("temp_out.csv", 'w')).writerows(a)


data = csv.reader(open("temp_out.csv"))
sortedlist = sorted(data, key=operator.itemgetter(0))


input_csv = "temp_out.csv"
a = zip(*sortedlist)
csv.writer(open("temp_out2.csv", 'w')).writerows(a)

