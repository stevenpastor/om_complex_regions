import numpy as np
#from matplotlib import pyplot as plt
import csv
import argparse
from sklearn.cluster import KMeans


def open_matrix(in_matrix):
    """
    doc: todo
    """

    ids = []
    data = []
    with open(in_matrix, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            ids.append(line[0])
            data.append( [float(x) for x in line[1:]] )
        
        data = np.array(data)

    return [ids, data]

def kmeans_matrix(ids, data, clusters, threshold, outfile):
    """
    doc: todo
    """

    kmeans = KMeans(n_clusters=clusters)
    kmeans.fit(data)
    y_kmeans = kmeans.predict(data)

    # find proportion of group_one; no need to find group_two since group_two = 1 - group_one
    # NOTE: would need to change this if using k != 2:
    group_one = [ids[ind] for ind,i in enumerate(y_kmeans) if i == 0]
    group_two = [ids[ind] for ind,i in enumerate(y_kmeans) if i == 1]
    group_one_len = len(group_one) / (len(group_one) + len(group_two))

    # If does not meet threshold proportion, collapse into 1 group:
    threshold = float(threshold)
    if group_one_len >= threshold or group_one_len <= (1 - threshold):
        collapsed_group = [ids[ind] for ind,i in enumerate(y_kmeans) if i == 0] + [ids[ind] for ind,i in enumerate(y_kmeans) if i == 1]
        print("Collapsed into one group, since not met "+str(threshold)+" threshold:")
        print(collapsed_group)
        print("Wrote out these results to this file, one line as the collapsed group: "+outfile)
        myfile = open(outfile, 'w')
        for i in collapsed_group:
            myfile.write("%s\t" % i)
        myfile.write("\n")
        myfile.close()
    # Otherwise, keep as 2 separate groups:
    else:
        print("Group 1 Molecule IDs:")
        print(group_one)
        print("Group 2 Molecule IDs:")
        print(group_two)
        print("Wrote out these results to this file, 2 lines, one per group as shown: "+outfile)
        myfile = open(outfile, 'w')
        for i in group_one:
            myfile.write("%s\t" % i)
        myfile.write("\n")
        for i in group_two:
            myfile.write("%s\t" % i)
        myfile.write("\n")
        myfile.close()

    return

def main():
    """
    Main controller function
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, 
                        help="Matrix file input")
    parser.add_argument("-o", "--output", type=str, 
                        help="K Means groups output")
    parser.add_argument("-c", "--clusters", type=int,
                        nargs='?', default=2,
                        help="Number of clusters wanted, defaults to 2")
    parser.add_argument("-t", "--threshold", type=float,
                        nargs='?', default=0.8,
                        help="Group proportion threshold - if not met, collapses to one group")
    args = parser.parse_args()

    # open matrix data only:
    ids, data = open_matrix(args.input)

    # k means on the matrix, prints the <clusters> groups:
    kmeans_matrix(ids, data, args.clusters, args.threshold, args.output)


if __name__ == "__main__": main()

