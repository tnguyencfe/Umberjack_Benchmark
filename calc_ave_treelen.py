"""
Calculate the average tree length fed to the indelbile partitions, put into csv
"""
import os
import glob
import csv


output_csv = "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/collate_treelen_partition.csv"

with open(output_csv, 'w') as fh_out:
    writer = csv.DictWriter(fh_out, fieldnames=["Name", "AveTreeLen"])
    writer.writeheader()
    for part_csv in glob.glob("/home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/*/fullpopn/partition.csv"):
        print part_csv
        with open(part_csv, 'rU') as fh_in:
            reader = csv.DictReader(fh_in)

            # TreeFile	TreeLen	Codons

            total_sites = 0.0
            weighted_total_treelen = 0.0
            for row in reader:
                weighted_total_treelen += (float(row["TreeLen"]) * float(row["Codons"]))
                total_sites +=  float(row["Codons"])

            weighted_ave_treelen = weighted_total_treelen / total_sites

            name = part_csv.replace("/home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/", "").replace("/fullpopn/partition.csv", "")
            outrow = dict(Name=name, AveTreeLen=weighted_ave_treelen)
            writer.writerow(outrow)








