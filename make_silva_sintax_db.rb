#!/usr/bin/env ruby

Signal.trap("PIPE", "EXIT")

require "abort_if"
require "fileutils"
require "parse_fasta"
require "optimist"
require "set"

include AbortIf
include AbortIf::Assert

good_ranks = [
  "superkingdom", # domain
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species",
]


def check_and_replace_no_break_space str
  bytes = str.bytes

  if (i1 = bytes.index(194)) &&
     (i2 = bytes.index(160)) &&
     (i1 = i2 - 1)

    bytes.map! do |byte|
      if byte == 194
        32
      elsif byte == 160
        nil
      else
        byte
      end
    end

    bytes.compact!

    bytes.pack("U*")
  else
    str
  end
end

# silva_tax_f = "/Users/moorer/projects/stec_cattle_microbiome_16s/external_data/silva_ltp_132_ssu/LTPs132_SSU.csv"
# names_f = "tax_dump/names.dmp"
# nodes_f = "tax_dump/nodes.dmp"
# fasta_f = "/Users/moorer/projects/stec_cattle_microbiome_16s/external_data/silva_ltp_132_ssu/LTPs132_SSU_compressed.fasta"
# fasta_tax_f = "TODO.fasta"
# lineages_with_na_f = "lineages_with_NA.txt"

opts = Optimist.options do
  banner <<-EOS

  Get silva taxonomy stuff.  I try and guess when things are wonky.

  Options:
  EOS

  opt(:silva_tax,
      "Input file",
      type: :string)
  opt(:names,
      "names.dmp",
      type: :string)
  opt(:nodes,
      "nodes.dmp",
      type: :string)
  opt(:silva_ltp,
      "silva ltp proj fasta",
      type: :string)
  opt(:outdir,
      "Output directory",
      type: :string,
      default: ".")
end

# Handle outdir arg
outdir = opts[:outdir]
FileUtils.mkdir_p outdir

#### sequence ids => taxa

silva_tax_f = opts[:silva_tax]
names_f = opts[:names]
nodes_f = opts[:nodes]
fasta_f = opts[:silva_ltp]

dirname = File.dirname fasta_f
extname = File.extname fasta_f
basename = File.basename fasta_f

fasta_tax_f = File.join outdir, "#{basename}.sintax.fasta"
lineages_with_na_f = File.join outdir, "#{basename}.sintax.lineages_with_missing_vals.txt"

################ read the nodes.dmp

AbortIf.logger.info { "Read nodes dump" }

TAXID2RANK = {}
TAXID2PARENT = {}
File.open(nodes_f, "rt").each_line do |line|
  tax_id, parent_tax_id, rank, *rest = line.downcase.sub("\t|\n", "").split "\t|\t"

  #abort_if TAXID2RANK.has_key?(tax_id), "'#{tax_id}' repeated in nodes_f '#{nodes_f}'"

  # rank = rank == "superkingdom" ? "domain" : rank

  # if good_ranks.include? rank

  TAXID2RANK[tax_id] = rank
  TAXID2PARENT[tax_id] = parent_tax_id
  # end
end

# set '1' => nil
TAXID2PARENT["1"] = nil

################ read the names.dmp

  #  Also sometimes things are
  # spelled differently in NCBI than in silva.  like thiobacillaceae
  # has an extra i in silva thiobacilliaceae.  Also lysobacteriales.


AbortIf.logger.info { "Reading names dump" }

NAME2TAXIDS = Hash.new do |ht, k|
  ht[k] = Set.new
end

TAXID2NAME = Hash.new do |ht, k|
  ht[k] = Set.new
end

File.open(names_f, "rt").each_line do |line|
  tax_id, name, unique_name, type = line.downcase.chomp("\t|\n").split("\t|\t")

  if TAXID2RANK.has_key? tax_id
    [name, unique_name].each do |name_id|
      unless name_id.empty?
        NAME2TAXIDS[name_id] << tax_id
        TAXID2NAME[tax_id] << {name: name, type: type} << {name: unique_name, type: type}
      end
    end
  end
end

def get_lineage taxid

  lineage = [taxid]
  parent = TAXID2PARENT[taxid]

  while parent
    lineage << parent
    parent = TAXID2PARENT[parent]
  end

  lineage.reverse
end

ALL_TAX_NAMES = Hash.new do |ht, rank|
  ht[rank] = Hash.new 0
end


SEQID2LINEAGE = {}
File.open(lineages_with_na_f, "w") do |f|

  File.open(silva_tax_f, "rt").each_line do |line|
    ary = line.downcase.chomp.split "\t"

    seq_id = ary[0]
    silva_lineage = ary[10].split(";").map do |name|
      name.
        sub("unclassified", "").
        sub("subgroup", "").
        # Silva and NCBI have different spellings for these.
        sub("thiobacilliaceae", "thiobacillaceae").
        sub("lysobacteriales", "lysobacterales").
        strip
    end

    # e.g., Erwinia gerundensis
    full_species = check_and_replace_no_break_space ary[4]
    # some might have 'subs. asoiten' or something like that on the end
    genus, species, *rest = full_species.split " "

    # just in case
    genus_lookup = genus.sub("unclassified", "").sub("subgroup", "").strip

    new_lineage = Hash[good_ranks.zip(Array.new(good_ranks.length, "NA"))]

    found_silva_names = Set.new

    lineage_lookups = []

    if NAME2TAXIDS.has_key? genus_lookup
      taxids = NAME2TAXIDS[genus_lookup]
      #assert !taxids.empty?

      taxids.each do |id|
        lineage = get_lineage id
        # this is the ncbi lineage lookup
        lineage_lookup = Hash[lineage.map {|id| [TAXID2RANK[id], TAXID2NAME[id]] }]
        lineage_lookups << lineage_lookup

        good_ranks.each do |rank|
          if lineage_lookup.has_key? rank
            names = lineage_lookup[rank]

            silva_lineage.each do |silva_name|
              # names is [{name: 'asoiten", type: 'arstoensten'}, ...]
              if names.map {|ht| ht[:name] }.include? silva_name
                if !new_lineage.has_key?(rank) || new_lineage[rank] == "NA"
                  new_lineage[rank] = silva_name
                  ALL_TAX_NAMES[rank][silva_name] += 1
                  found_silva_names << silva_name
                end
              end
            end
          end
        end
      end
    else
      AbortIf.logger.error { "Missing #{genus_lookup}" }
    end

    # Check if any of the silva names were not found
    missing_silva_names = silva_lineage.to_set - found_silva_names # TODO Is a chance we might miss some doubled names?
    missing_ranks = new_lineage.select { |rank, name| name == "NA" }

    # The numbers of these may be different as the silva lineage can
    # include things like subclass or other taxonomic levels that we
    # don't care about.
    if missing_ranks.has_key? "family"
      # check if any missing names end in aceae
      aceae_names = missing_silva_names.select { |name| name.end_with? "aceae" }
      if aceae_names.count == 1 # if there are more than one, it's ambiguous so don't guess
        new_lineage["family"] = aceae_names.first
        missing_silva_names.delete aceae_names.first
        missing_ranks.delete "family"

        AbortIf.logger.warn do
          "Gussing that rank 'family' in lineage " \
          "#{silva_lineage.inspect} for species " \
          "'#{genus} #{species}' is #{aceae_names.first}"
        end
      end
    end

    if missing_ranks.has_key? "order"
      ales_names = missing_silva_names.select { |name| name.end_with? "ales" }
      if ales_names.count == 1
        new_lineage["order"] = ales_names.first
        missing_silva_names.delete ales_names.first
        missing_ranks.delete "order"

        AbortIf.logger.warn do
          "Gussing that rank 'order' in lineage " \
          "#{silva_lineage.inspect} for species " \
          "'#{genus} #{species}' is #{ales_names.first}"
        end
      end
    end

    # Sometimes the silva lineage doesn't have each level we need but
    # the NCBI lineage does.

    missing_ranks.each do |rank, _|
      # all names should be "NA" by design still

      lineage_lookups.each do |lineage_lookup|
        if lineage_lookup.has_key?(rank)

          # get sci name
          this_rank_names = lineage_lookup[rank]
          this_sci_names = this_rank_names.select { |ht| ht[:type] == "scientific name" }

          if this_sci_names.count == 1 # didn't see any of these in the actual output
            AbortIf.logger.warn { "(one sci name) Using NCBI name '#{this_sci_names.first}' for rank '#{rank}' for species '#{genus} #{species}'.  It's silva lineage was '#{silva_lineage.inspect}'" }

            new_lineage[rank] = this_sci_names.first[:name]
          else
            # there are multiple valid sci names, let's check to see if
            # any of them have been seen before, and if they have, pick
            # the one that has been seen the most.
            all_names = this_sci_names.map do |ht|
              if ht[:name].empty? # TODO remove empty naems at the beginning.
                nil
              else

                [ht[:name], ALL_TAX_NAMES[rank][ht[:name]]]
              end
            end.compact.sort_by { |name, count| count }

            highest_count_name = all_names.last
            # returns [name, times_seen]

            # the highest count name still might never have been seen, but
            # that's okay too.
            new_lineage[rank] = highest_count_name.first

            AbortIf.logger.warn { "(multi sci names) Using NCBI name '#{highest_count_name.first}' (seen before #{highest_count_name.last} times. all_names #{all_names.inspect}) for rank '#{rank}' for species '#{genus} #{species}'.  It's silva lineage was '#{silva_lineage.inspect}'" }

          end
        else
          unless rank == "species"
            AbortIf.logger.warn { "Couldn't find '#{rank}' in the lineage lookup. (silva lineage was #{silva_lineage.inspect}" }
          end
        end
      end
    end

    new_lineage["species"] = species

    if new_lineage.values.any? { |val| val == "NA" }
      f.puts new_lineage.inspect
    end

    SEQID2LINEAGE[seq_id] = new_lineage
  end
end

def lineage_ht_to_s lineage_ht
  lineage_ht.map do |rank, name|
    name = name.gsub /\W+/, "_"

    if rank == "superkingdom"
      "d:#{name}"
    else
      "#{rank[0]}:#{name}"
    end
  end.join ","
end


File.open(fasta_tax_f, "w") do |f|

  ParseFasta::SeqFile.open(fasta_f).each_record do |rec|
    if SEQID2LINEAGE.has_key? rec.id.downcase
      lineage = SEQID2LINEAGE[rec.id.downcase]
      new_header = "#{rec.id};tax=#{lineage_ht_to_s(lineage)};"

      f.puts ">#{new_header}\n#{rec.seq}"
    else
      AbortIf.logger.error { "No tax info for #{rec.id}" }
    end
  end
end
