=begin
1. Apply the lineage cutoff to filter low abundance variants
2. Obtain lineages from each well
3. Put all lineage sequences in one file and count number of identical lineags and their abundances
4. Make neighborjoining trees of lineages from all wells of one subject
=end

#open file in fasta format, turn into a sequence hash

def fasta_to_hash(infile)
  f=File.open(infile,"r")
  return_hash = {}
  name = ""
  while line = f.gets do
    if line =~ /^\>/
      name = line.chomp
      return_hash[name] = ""
    else
      return_hash[name] += line.chomp
    end
  end
  f.close
  return return_hash
end


class Sequence
  def initialize (name = ">sequence",dna_sequence ="")
    @name = name
    @dna_sequence = dna_sequence.upcase
    @aa_sequence = ""
    @aa_array = []
  end
  attr_accessor :name, :dna_sequence, :aa_sequence, :aa_array
  def rev_complement
    @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def rev_complement!
    @dna_sequence = @dna_sequence.reverse.upcase.tr('ATCG','TAGC')
  end
  def get_aa_sequence(initial_position = 0)
    @aa_sequence = ""
    require_sequence = @dna_sequence[initial_position..-1]
    base_array = []
    require_sequence.each_char {|base| base_array<<base}
    while (base_array.length>=3) do
      base_3= ""
      3.times{base_3 += base_array.shift}
      @aa_sequence<< amino_acid(base_3)
    end
  end

  def amino_acid (bases)
    case bases
    when /^TT[TCY]$/
      return "F"
    when /^TT[AGR]$/
      return "L"
    when /^CT.$/
      return "L"
    when /^AT[TCAHYWM]$/
      return "I"
    when "ATG"
      return "M"
    when /^GT.$/
      return "V"
    when /^TC.$/
      return "S"
    when /^CC.$/
      return "P"
    when /^AC.$/
      return "T"
    when /^GC.$/
      return "A"
    when /^TA[TCY]$/
      return "Y"
    when /^TA[AGR]$/ 
      return "*"
    when /^T[GR]A$/
      return "*"
    when /^CA[TCY]$/
      return "H"
    when /^CA[AGR]$/
      return "Q"
    when /^AA[TCY]$/
      return "N"
    when /^AA[AGR]$/
      return "K"
    when /^GA[TCY]$/
      return "D"
    when /^GA[AGR]$/
      return "E"
    when /^TG[TCY]$/
      return "C"
    when "TGG"
      return "W"
    when /^CG.$/
      return "R"
    when /^AG[TCY]$/
      return "S"
    when /^[AM]G[AGR]$/
      return "R"
    when /^GG.$/
      return "G"
    when /^[ATW][CGS][CTY]$/
      return "S"
    when /^[TCY]T[AGR]$/
      return "L"
    else
      return "#"
    end
  end
end

def count(array)
  hash = Hash.new(0)
  array.each do |element|
    hash[element] +=1
  end
  return hash
end


#pre-determined lineage abundance cutoff
lineage_cut_off = 0.025

#input is a directory of all sequences from all wells of one subject (one file per well)
indir = ARGV[0]
subject = File.basename(indir)
outdir = indir + "_out"
Dir.mkdir(outdir) unless File.directory?(outdir)

lineage_dir = outdir + "/lineage"
Dir.mkdir(lineage_dir) unless File.directory?(lineage_dir)
lineage_aa_dir = outdir + "/lineage_aa"
Dir.mkdir(lineage_aa_dir) unless File.directory?(lineage_aa_dir)
lineage_tree_dir = outdir + "/lineage_tree"
Dir.mkdir(lineage_tree_dir) unless File.directory?(lineage_tree_dir)

libs = []
all_sample_lineages = {}
all_sample_aa_lineages = {}
Dir.chdir(indir) {libs = Dir.glob("*")}
libs.each do |sample_well|
  puts sample_well
  path = indir + "/" + sample_well
  combined_seq = fasta_to_hash(path)
  seq_freq = count(combined_seq.values)
  sample_size = combined_seq.size
  lineage_out = File.open((lineage_dir + "/" + sample_well), "w")
  lineage_aa_out = File.open((lineage_aa_dir + "/" + sample_well), "w")
  seq_n = 1
  seq_freq.each do |seq,number|
    abund = number/sample_size.to_f
    if abund > lineage_cut_off
      seq_name = ">" + sample_well + "_" + seq_n.to_s + "_" + number.to_s + "_" + sample_size.to_s
      seq_n += 1
      all_sample_lineages[seq_name] = seq
      lineage_out.puts seq_name
      lineage_out.puts seq
      lineage_aa_out.puts seq_name
      s = Sequence.new(seq_name,seq)
      
      #amino acid sequences from nucleotide sequences
      s.get_aa_sequence(0)
      lineage_aa_out.puts s.aa_sequence
      all_sample_aa_lineages[seq_name] = s.aa_sequence
    end
  end
  lineage_out.close
  lineage_aa_out.close
end

#output lineage nucleotide sequences in fasta format for all wells
all_lineage = lineage_tree_dir + "/" + subject
all_lineage_out = File.open(all_lineage,"w")
all_sample_lineages.each do |k,v|
  all_lineage_out.puts k
  all_lineage_out.puts v
end
all_lineage_out.close

#output lineage amino acid sequences in fasta format for all wells. 
all_aa_lineage = lineage_tree_dir + "/" + subject + "_aa"
all_aa_lineage_out = File.open(all_aa_lineage,"w")
all_sample_aa_lineages.each do |k,v|
  all_aa_lineage_out.puts k
  all_aa_lineage_out.puts v
end
all_aa_lineage_out.close


#make lineage trees, require MUSCLE
puts `muscle -in #{all_lineage} -out #{all_lineage + ".fa"} -tree2 #{all_lineage + ".tre"} -cluster neighborjoining -quiet`
lineage_count = count(all_sample_lineages.values)

#output lineage abundance .csv file for IUPM calculation
n = 1
lineage_count_f = lineage_tree_dir + "/" + subject + ".csv"
lineage_count_out = File.open(lineage_count_f,"w")
lineage_count_out.puts "Lineage,Count"
lineage_count.values.each do |v|
  lineage_count_out.puts "L" + n.to_s + "," + v.to_s
  n += 1
end
lineage_count_out.close