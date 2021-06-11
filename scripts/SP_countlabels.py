import csv
def find_labels(xmap, contig_start_label, contig_end_label, contig_orientation):
  """
  this does not account for split mapped molecules 
  i.e., assumes one molecule ID per line in XMAP
  fix later
  """
  # swap if -:
  #print(contig_start_label, contig_end_label)
  if contig_orientation == "-":
    contig_start_label,contig_end_label = contig_end_label,contig_start_label
  #print(contig_start_label, contig_end_label)
  print("mol_id"+"\t"+"five_prime_labels"+"\t"+"three_prime_labels")
  # cvs module:
  with open(xmap, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    xmap_lines = [line for line in reader if "#" not in line[0]]
    for x in xmap_lines:
      # tuples, replace paren, remove trailing comma, comma split, and every other ele (contig labels):
      contig_labels_only = x[13].replace("(", "").replace(")", ",")[:-1].split(',')[::2]
      contig_labels_only_int = [int(i) for i in contig_labels_only] # make int from str!
      # check that molecule crosses into SDs:
      if contig_orientation == "-":
        smaller_five_prime = [i for i in contig_labels_only_int if i > contig_start_label] # 5' of SDs
        larger_five_prime = [i for i in contig_labels_only_int if i <= contig_start_label] # into SDs
        smaller_three_prime = [i for i in contig_labels_only_int if i >= contig_end_label] # into SDs
        larger_three_prime = [i for i in contig_labels_only_int if i < contig_end_label] # 3' of SDs
        # if crosses either side of SDs, count labels outside SDs:
        if len(smaller_five_prime) and len(larger_five_prime) and len(smaller_three_prime) and len(larger_three_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i > contig_start_label]), len([i for i in contig_labels_only_int if i < contig_end_label]))
        elif len(smaller_five_prime) and len(larger_five_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i > contig_start_label]), len([i for i in contig_labels_only_int if i < contig_end_label]))
        elif len(smaller_three_prime) and len(larger_three_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i > contig_start_label]), len([i for i in contig_labels_only_int if i < contig_end_label]))
      # check that molecule crosses into SDs (note the change in < or > below):
      elif contig_orientation == "+":
        smaller_five_prime = [i for i in contig_labels_only_int if i < contig_start_label] # 5' of SDs
        larger_five_prime = [i for i in contig_labels_only_int if i >= contig_start_label] # into SDs
        smaller_three_prime = [i for i in contig_labels_only_int if i <= contig_end_label] # into SDs
        larger_three_prime = [i for i in contig_labels_only_int if i > contig_end_label] # 3' of SDs
        # if crosses either side of SDs, count labels outside SDs:
        if len(smaller_five_prime) and len(larger_five_prime) and len(smaller_three_prime) and len(larger_three_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i < contig_start_label]), len([i for i in contig_labels_only_int if i > contig_end_label]))
        elif len(smaller_five_prime) and len(larger_five_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i < contig_start_label]), len([i for i in contig_labels_only_int if i > contig_end_label]))
        elif len(smaller_three_prime) and len(larger_three_prime) > 1:
          print(x[1], len([i for i in contig_labels_only_int if i < contig_start_label]), len([i for i in contig_labels_only_int if i > contig_end_label]))
  return
find_labels('results/11029B_initial_genome_check/11029B_fullContig211_molecules.xmap', 2214, 2253, "-")